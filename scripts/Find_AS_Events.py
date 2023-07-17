#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2023.07.16

This is a script designed to enumerate all alternative splicing events that may be associated with a set of 
transcripts provided as input. Specifically, this script requires the following two files as input:
    * A GTF file containing genomic coordinates of input transcripts
    * A GTF file of reference gene annotations (e.g., from GENCODE)

Briefly, for each input transcript, our script will pull out the canonical transcript isoform of the corresponding
gene. We will next enumerate local differences in transcript structure between the input transcript and the 
canonical transcript isoform. Each local difference will be classified into one of the following alternative splicing
event types (relative to the canonical transcript isoform):
    * Exon skipping
    * Exon inclusion
    * Alternative 5' splice site
    * Alternative 3' splice site
    * Mutually exclusive exons
    * Intron retention
    * Exitron splicing
    * Alternative first exon
    * Alternative last exon
    * Complex splicing

Several edge cases (listed below) may arise when we compare the structure of a given transcript with the structure
of the canonical transcript isoform of the corresponding gene:
    * Case 1: The input transcript is intergenic (no canonical transcript is available to serve as a reference)
    * Case 2: The input transcript is the canonical transcript isoform of the corresponding gene
    * Case 3: The input transcript does not overlap at all with the canonical transcript isoform
    * Case 4: The input transcript only differs in transcript ends relative to the canonical transcript isoform

Our script will generate an output file containing the following pieces of information for each input transcript:
    * Chromosome
    * Strand
    * Gene ID
    * Transcript ID (input transcript)
    * Transcript ID (canonical transcript)
    * Comma-separated list of genomic coordinates describing alternative splicing event (input transcript)
    * Comma-separated list of genomic coordinates describing alternative splicing event (canonical transcript)
    * Alternative splicing event type
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import sys, argparse, re
from datetime import datetime
import networkx as nx
import pandas as pd
import numpy as np

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def outlog(myString):
    '''
    This is a function to print some output message string (myString) with the date and time
    '''
    print('[', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), '] ', myString, sep='', flush=True)

def checkOverlap(t1, t2):
    '''
    This is a function to check if two tuples t1 and t2 overlap
    '''
    # Make sure t1 and t2 are sorted
    t1, t2 = tuple(sorted(t1)), tuple(sorted(t2))
    return min(t1[1], t2[1]) >= max(t1[0], t2[0])

def bluntTxEnds(tx1, tx2):
    '''
    This is a function to "blunt" the ends of two transcripts (tx1 and tx2) if their terminal exons overlap
    '''
    # Check if the leftmost exons of tx1 and tx2 overlap
    if checkOverlap(tx1[0], tx2[0]):
        tx1 = [(max(tx1[0][0], tx2[0][0]), tx1[0][1])] + tx1[1:]
        tx2 = [(max(tx1[0][0], tx2[0][0]), tx2[0][1])] + tx2[1:]

    # Check if the rightmost exons of tx1 and tx2 overlap
    if checkOverlap(tx1[-1], tx2[-1]):
        tx1 = tx1[:-1] + [(tx1[-1][0], min(tx1[-1][1], tx2[-1][1]))]
        tx2 = tx2[:-1] + [(tx2[-1][0], min(tx1[-1][1], tx2[-1][1]))]
    
    return tx1, tx2

def reverseTupleList(myList):
    '''
    This is a function to reverse a list of tuples and their individual elements
    '''
    return [(tup[1], tup[0]) for tup in reversed(myList)]

def getIntronCoord(exonCoord):
    '''
    This is a function that returns a list of tuples representing intron coordinates
    relative to exonCoord
    '''
    return [(exonCoord[i][1], exonCoord[i+1][0]) for i in range(len(exonCoord)-1)]

def buildSpliceGraph(tx1, tx2, strand):
    '''
    This is a function to build a splice graph given exon coordinates from two transcripts
    '''
    # Initialize a directed graph
    spliceGraph = nx.DiGraph()

    # Reverse exon coordinates if strand is negative
    tx1, tx2 = (reverseTupleList(tx1), reverseTupleList(tx2)) if strand == '-' else (tx1, tx2)

    # Declare root and sink nodes; determine root and sink node values based on strand
    rootVal, sinkVal = (0, np.inf) if strand == '+' else (np.inf, 0)
    auxNodes = [('sink', {'attribute': 'sink', 'status': 1, 'value': sinkVal}), ('root', {'attribute': 'root', 'status': 1, 'value': rootVal})]
    spliceGraph.add_nodes_from(auxNodes)

    # Add nodes/edges from tx1 and tx2 to spliceGraph
    # Nodes unique to tx1 will have a status of 0; nodes unique to tx2 will have a status of 1
    for idx, tx in enumerate([tx1, tx2]):
        myNodes = [(str(exon[0]) + '_S', {'attribute': 'start', 'status': idx, 'value': exon[0]}) for exon in tx]
        myNodes += [(str(exon[1]) + '_E', {'attribute': 'end', 'status': idx, 'value': exon[1]}) for exon in tx]
        myEdges = [('root', str(tx[0][0]) + '_S'), (str(tx[-1][1]) + '_E', 'sink')]
        myEdges += [(str(edge[0]) + '_S', str(edge[1]) + '_E') for edge in tx]
        myEdges += [(str(edge[0]) + '_E', str(edge[1]) + '_S') for edge in getIntronCoord(tx)]
        spliceGraph.add_nodes_from(myNodes)
        spliceGraph.add_edges_from(myEdges)
            
    return spliceGraph

def findBubbles(spliceGraph):
    '''
    This is a function that will return all "bubbles" within a spliceGraph
    '''
    # Instantiate bubbles, currentNode, and bubbleStart
    bubbles, currentNode, bubbleStart = [], 'root', 'root'

    # Repeat the following algorithm until bubbleStart is set to 'sink'
    while not bubbleStart == 'sink':
        # Check if currentNode has an indegree of two
        if spliceGraph.in_degree(currentNode) == 2:
            # We have reached the end of a bubble; enumerate all simple paths between bubbleStart and currentNode
            bubbles.append(list(nx.all_simple_paths(spliceGraph, source=bubbleStart, target=currentNode)))

        if spliceGraph.out_degree(currentNode) == 2:
            # We have reached the beginning of a bubble; set bubbleStart to the currentNode and pick a child of
            # bubbleStart to be the currentNode
            bubbleStart = currentNode
            currentNode = list(spliceGraph.successors(currentNode))[0]
        elif spliceGraph.out_degree(currentNode) == 1:
            # Set the child of currentNode to currentNode
            currentNode = list(spliceGraph.successors(currentNode))[0]
        else:
            # Our currentNode is the sink node
            bubbleStart = currentNode
        
    return bubbles

def classifyBubble(bubble, attributes, status):
    '''
    This is a function that will classify bubbles as different alternative splicing event types
    '''
    # Figure out which path in bubble corresponds to tx1 and tx2 based on the status of the second
    # node in the longer path
    longPathIdx = pd.Series([len(path) for path in bubble]).idxmax()
    path1, path2 = (bubble[longPathIdx], bubble[longPathIdx^1]) if status[bubble[longPathIdx][1]] == 0 else (bubble[longPathIdx^1], bubble[longPathIdx])

    if min(len(path1), len(path2)) == 2 and max(len(path1), len(path2)) == 4:
        '''
        Check if bubble involves exon skipping or exon inclusion
            * bubble start node has attribute 'end' 
            * bubble end node has attribute 'start'
            * exon skipping: path2 has a length of 2 nodes
            * exon inclusion: path2 has a length of 4 nodes
        Check if bubble involves intron retention or exitron splicing
            * bubble start node has attribute 'start'
            * bubble end node has attribute 'end'
            * intron retention: path2 has a length of 2 nodes
            * exitron skipping: path2 has a length of 4 nodes
        '''
        if attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'start':
            event = 'exon_skipping' if len(path2) == 2 else 'exon_inclusion'
        elif attributes[path1[0]] == 'start' and attributes[path1[-1]] == 'end':
            event = 'intron_retention' if len(path2) == 2 else 'exitron_splicing'
        else:
            event = 'complex_splicing'
    
    elif len(path1) == 3 and len(path2) == 3:
        '''
        Check if bubble involves alternative 5' splice site
            * both bubble start and end nodes have attribute 'start'
        Check if bubble involves alternative 3' splice site
            * both bubble start and end nodes have attribute 'end'
        '''
        if attributes[path1[0]] == 'start' and attributes[path1[-1]] == 'start':
            event = 'alternative_5ss'
        elif attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'end':
            event = 'alternative_3ss'
        else:
            event = 'complex_splicing'
    
    elif len(path1) == 4 and len(path2) == 4:
        '''
        Check if bubble involves mutually exclusive exons:
            * bubble start node has attribute 'end'
            * bubble end node has attribute 'start'
            * internal exons are not overlapping
        Check if bubble involves an alternative first exon:
            * bubble start node is the root
            * bubble end node has attribute 'start'
            * first exons are not overlapping
        Check if bubble involves an alternative last exon:
            * bubble start node has attribute 'end'
            * bubble end node is the sink
            * last exons are not overlapping
        '''
        if attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'start' and not checkOverlap(tuple(path1[1:3]), tuple(path2[1:3])):
            event = 'mutually_exclusive_exons'
        elif attributes[path1[0]] == 'root' and attributes[path1[-1]] == 'start' and not checkOverlap(tuple(path1[1:3]), tuple(path2[1:3])):
            event = 'alternative_first_exon'
        elif attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'sink' and not checkOverlap(tuple(path1[1:3]), tuple(path2[1:3])):
            event = 'alternative_last_exon'
        else:
            event = 'complex_splicing'
    
    else:
       event = 'complex_splicing' 
            
    return event, path1, path2

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def getExonCoord(annoDF):
    '''
    This is a function that will return the genomic coordinates of exons for all transcripts represented
    in annoDF. Genomic coordinates of exons will be reported as a sorted list of tuples.
    '''
    # Pull out exon-level annotations from annoDF; then extricate transcript_id and gene_id
    exonDF = annoDF[annoDF[2] == 'exon'][[0,3,4,6,8]].reset_index(drop=True)
    exonDF['transcript_id'] = exonDF[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'transcript_id' in item][0])
    exonDF['gene_id'] = exonDF[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'gene_id' in item][0])

    # Group genomic coordinates of exons by transcript, chromosome, and strand
    exonDF = exonDF.groupby(['transcript_id', 'gene_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    exonDF['exonCoord'] = exonDF.apply(lambda x: [(x[3][i], x[4][i]) for i in range(len(x[3]))], axis=1)
    exonDF = exonDF[[0, 6, 'gene_id', 'transcript_id', 'exonCoord']]
    exonDF.columns = ['chr', 'strand', 'gene_id', 'transcript_id', 'exon_coord']

    return exonDF

def getCanonicalTx(annoDF):
    '''
    This is a function that will retrieve canonical transcripts of all genes represented in annoDF
    '''
    # Pull out transcript-level annotations from annoDF; then filter for annotations from canonical transcripts
    txDF = annoDF[annoDF[2] == 'transcript'][[8]].reset_index(drop=True)
    txDF = txDF[txDF[8].str.contains('Ensembl_canonical')].reset_index(drop=True)

    # Pull out gene_id and transcript_id information from txDF
    txDF['gene_id'] = txDF[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'gene_id' in item][0])
    txDF['transcript_id'] = txDF[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'transcript_id' in item][0])
    txDF = txDF.drop([8], axis=1)

    # Build dictionary mapping gene IDs to IDs of canonical transcript isoforms
    canonDict = dict(zip(txDF['gene_id'], txDF['transcript_id']))

    return canonDict

def findASEvents(myExonCoord, annoExonDict, canonDict):
    '''
    This is a function that will report all alternative splicing events associated with input transcripts
    '''
    output, colnames = [], ['chrom', 'strand', 'geneID', 'txID', 'canonical_txID', 'tx_path', 'canonical_path', 'event_type']

    # Iterate over rows of myExonCoord
    for idx, row in myExonCoord.iterrows():
        # Check that row['gene_id'] is defined
        if row['gene_id'] != 'NA':
            # Split row['gene_id'] by comma and iterate over each gene_id
            for gene_id in row['gene_id'].split(','):
                # Check if the input transcript is the canonical transcript
                if row['transcript_id'] != canonDict[gene_id]:
                    # Retrieve exon-level coordinates for canonical transcript
                    canonExonCoord = annoExonDict[canonDict[gene_id]]
                    # Check that the input transcript and canonical transcript have some degree of overlap
                    if checkOverlap((row['exon_coord'][0][0], row['exon_coord'][-1][1]), (canonExonCoord[0][0], canonExonCoord[-1][1])):
                        # Blunt the ends of the input transcript and canonical transcript
                        bluntTx, bluntCanon = bluntTxEnds(row['exon_coord'], canonExonCoord)

                        # Build a splice graph from bluntTx and bluntCanon
                        spliceGraph = buildSpliceGraph(bluntCanon, bluntTx, row['strand'])

                        # Generate dictionaries mapping nodes to attributes/statuses/values
                        attributes, status, valueDict = nx.get_node_attributes(spliceGraph, 'attribute'), nx.get_node_attributes(spliceGraph, 'status'), nx.get_node_attributes(spliceGraph, 'value')

                        # Enumerate all bubbles (alternative splicing events) in spliceGraph
                        bubbles = findBubbles(spliceGraph)

                        # Check if there are bubbles in spliceGraph
                        if len(bubbles) > 0:
                            classified = [classifyBubble(bubble, attributes, status) for bubble in bubbles]

                            # Iterate over classified bubbles
                            for bub in classified:
                                tx_path = ','.join([item.split('_')[0] for item in bub[2]])
                                canon_path = ','.join([item.split('_')[0] for item in bub[1]])
                                output.append([row['chr'], row['strand'], gene_id, row['transcript_id'], canonDict[gene_id], tx_path, canon_path, bub[0]])

                        else:
                            # The input transcript and canonical transcript only differ in transcript ends
                            output.append([row['chr'], row['strand'], gene_id, row['transcript_id'], canonDict[gene_id], 'NA', 'NA', 'transcript_end_difference'])
                    else:
                        # The input transcript and canonical transcript do not overlap
                        output.append([row['chr'], row['strand'], gene_id, row['transcript_id'], canonDict[gene_id], 'NA', 'NA', 'no_transcript_overlap'])
                else:
                    # The input transcript is identical to the canonical transcript
                    output.append([row['chr'], row['strand'], gene_id, row['transcript_id'], canonDict[gene_id], 'NA', 'NA', 'is_canonical'])
        else:
            # The input transcript is intergenic; can't really say what alternative splicing events are associated with it
            output.append([row['chr'], row['strand'], row['gene_id'], row['transcript_id'], 'NA', 'NA', 'NA', 'intergenic_transcript'])

        # Send update message for every 100 entries processed
        if (idx+1) % 100 == 0:
            outlog('Finished processing ' + str(idx+1) + ' out of ' + str(myExonCoord.shape[0]) + ' transcripts...')
        
    # Construct a dataframe from output
    outDF = pd.DataFrame(output, columns=colnames)

    return outDF

def main():
    message = 'Finds alternative splicing events associated with a set of input transcripts'
    parser = argparse.ArgumentParser(description=message)

    # Add arguments
    parser.add_argument('-i', metavar='/path/to/input/transcripts/GTF', required=True,
        help='path to GTF file containing genomic coordinates of input transcripts')
    parser.add_argument('-r', metavar='/path/to/reference/GTF', required=True,
        help='path to GTF file containing reference gene annotations')
    parser.add_argument('-o', metavar='/path/to/output/file', required=True,
        help='path to output file')

    # Parse command-line arguments
    args = parser.parse_args()
    infile, refAnno, outfile = args.i, args.a, args.r, args.o

    # Read in infile and refAnno and Pandas dataframes
    outlog('Reading in GTF files for input transcripts and reference gene annotations...')
    inDF = pd.read_csv(infile, sep='\t', header=None, comment='#')
    annoDF = pd.read_csv(refAnno, sep='\t', header=None, comment='#')

    # Retrieve IDs of canonical transcripts in annoDF
    outlog('Retrieving canonical transcript isoforms of genes provided as input...')
    canonDict = getCanonicalTx(annoDF)

    # Retrieve genomic coordinates of transcripts represented in inDF and annoDF
    outlog('Retrieving genomic coordinates of input transcripts and corresponding canonical transcripts...')
    myExonCoord, annoExonCoord = getExonCoord(inDF), getExonCoord(annoDF)

    # Generate a dictionary mapping transcript_id to exon_coord from annoExonCoord
    annoExonDict = dict(zip(annoExonCoord['transcript_id'], annoExonCoord['exon_coord']))

    # Identify AS events associated with transcripts in myExonCoord
    outlog('Identifying alternative splicing events associated with input transcripts...')
    outDF = findASEvents(myExonCoord, annoExonDict, canonDict)
    
    # Print outDF to outfile
    outDF.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    main()
