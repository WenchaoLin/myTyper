#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Date    : 2018-04-11
# Author  : Wenchao Lin (linwenchao@yeah.net)


import os
import sys
import argparse
import logging

from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from itertools import groupby
import operator
from core.sequence import reverseComplement
from tqdm import tqdm
from core import FastaReader,FastaWriter,FastaRecord
from core.matrix import DNAFULL
from core.align import aligner

FORMAT = '%(asctime)s [%(levelname)-8s] %(message)s'
logging.basicConfig(format=FORMAT, datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)





def read_barcode(fh):
    """read in the primer mapping file and
    example:

    #SampleID   ForwardBarcode  ReverseBarcode
    s5915   TCAGACGATGCGTCAT    GATCTCTACTATATGC
    s5916   TCAGACGATGCGTCAT    ACAGTCTATACTGCTG
    s5921   TCAGACGATGCGTCAT    ATGATGTGCTACATCT
    s5927   TCAGACGATGCGTCAT    CTGCGTGCTCTACGAC
    s5930   CTATACATGACTCTGC    GATCTCTACTATATGC
    s5931   CTATACATGACTCTGC    ACAGTCTATACTGCTG
   
    Args:
        fh (TYPE): file handle
    """

    barcode_dict = {}
    
    with open(fh) as f:
        rows = [ x.strip().split() for x in f if x.strip().startswith("#") == False]

    for row in rows:
        if barcode_dict.has_key(row[0]):
            log.error("Duplicated record %s found!" % row[0])
            exit(0)
        if row[1:] in barcode_dict.values():
            log.error("Duplicated record %s found!" % row[0])
            exit(0)
        barcode_dict[row[0]] = row[1:]

    return barcode_dict





def barcode_split(reads,bcs,mismatch=1,mode='slow'):
    """split reads into files by barcodes and primers
    
    Args:
        reads (TYPE): FastaReader reads object
        bcs (TYPE): barcodes object
        pms (TYPE): primers object
        mismatch (int, optional): not in use
        mode (str, optional): mode ='fast'; for match=0 fast barcode search.
    
    Returns:
        TYPE: Description
    """

    bcs_len = len(bcs.values()[0][0])
    check = int(bcs_len) * 2 - mismatch

    result = []

    reads_format = reads.sequence[:bcs_len] + '...' + reads.sequence[-bcs_len:]
    reads_barcode_forward = str(reads.sequence[:bcs_len])
    reads_barcode_reverse = reads.reverseComplement().sequence[:bcs_len]

    reads_revcom = reads.reverseComplement().sequence


    # name[0] is forward barcode name[1] is reverse barcode
    for name in bcs:
        # barcode完全匹配的快速搜索模式
        if mode == 'fast':
            if reads_barcode_forward == bcs[name][0] and reads_barcode_reverse == bcs[name][1]:
                result.append([reads.id,name,reads.sequence, bcs[name],'F',reads_format,bcs_len,bcs_len])
                continue
            elif reads_barcode_forward == bcs[name][1] and reads_barcode_reverse == bcs[name][0]:
                result.append([reads.id,name,reads_revcom, bcs[name],'R',reads_format,bcs_len,bcs_len])
                continue
        else:

            bc_alignmentsFF = pairwise2.align.localxx(reads_barcode_forward,bcs[name][0])
            bc_alignmentsFR = pairwise2.align.localxx(reads_barcode_reverse,bcs[name][1])
            bc_alignmentsRF = pairwise2.align.localxx(reads_barcode_forward,bcs[name][1])
            bc_alignmentsRR = pairwise2.align.localxx(reads_barcode_reverse,bcs[name][0])

            try:
                #找到有mistach个mismatch的barcode
                if int(bc_alignmentsFF[0][2]) + int(bc_alignmentsFR[0][2]) >= check:
                    # print( "%s : %s : %s : forward" % ( reads_format ,name, bcs[name]))
                    # print(format_alignment(*bc_alignmentsFF[0]))
                    # print(format_alignment(*bc_alignmentsFR[0]))
                    result.append([reads.id,name,reads.sequence, bcs[name],'F',reads_format,bc_alignmentsFF[0][2],bc_alignmentsFR[0][2]])
                    # result.append([reads.id,name])
                    continue
                elif int(bc_alignmentsRF[0][2]) + int(bc_alignmentsRR[0][2]) >= check:
                    # print( "%s : %s : %s : reverse" % (reads_format ,name, bcs[name]))
                    # print(format_alignment(*bc_alignmentsRF[0]))
                    # print(format_alignment(*bc_alignmentsRR[0]))
                    result.append([reads.id,name,reads_revcom, bcs[name],'R',reads_format,bc_alignmentsRF[0][2],bc_alignmentsRR[0][2]])
                    # result.append([reads.id,name])
                    continue
                else:
                    continue
            except:
                # log.error("barcode search Error, please check [%s] in your barcode file." % name)
                pass

    return result





def primer_split(reads,pms,mismatch=3):
    """split reads into files by barcodes and primers
    
    Args:
        reads (TYPE): FastaReader reads object
        bcs (TYPE): barcodes object
        mismatch (int, optional): not in use
        
        sample:
        ['m180326_085254_42208_c101249792550000001823299208051850_s1_p0/56/ccs', 's5927', 'TCAGACGATGCGTCATATGAAAGATGCAAAGCGCCTGAATTTTCTGACTCTTCCCATCAGACCCCCAAAGACACATGTGACCCACCACCCCATCTCTGACCATGAGGTCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGCGAGGACCAAACTCAGGACACTGAGCTTGTGGAGACCAGACCAGCAGGAGATAGAACCTTCCAGAAGTGGCAGCTGTGGTGGTGCCTTCTGGAGAAGAGCAGAGATACACATGCCATGTACAGCATGAGGGGCTGCCGAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGGGATGAGGGGTCATATCTGTCGTAGAGCACGCAG', ['TCAGACGATGCGTCAT', 'CTGCGTGCTCTACGAC'], 'F', 'TCAGACGATGCGTCAT...GTCGTAGAGCACGCAG', 16, 16]

    Returns:
        TYPE: list['readsID','barcodID','sequence',['barcodeF','barcodeR'],'Strand','sequence(head...tail)','numBarcodeMatchF','numBarcodeMatchR','bestPrimerID','bestLeftPrimerMismatch','bestRightPrimerMismatch','bestLeftPrimerScore','bestRightPrimerScore']
    """

    bcd_len = len(reads[3][0])

    hit_score = 0

    for x in pms:
        left_primer_reads = reads[2][bcd_len:bcd_len+len(pms[x][0])]
        right_primer_reads = reverseComplement(reads[2])[bcd_len:bcd_len+len(pms[x][1])]
        alignL = aligner(pms[x][0],left_primer_reads,method='global',matrix=DNAFULL, max_hits=1)
        alignR = aligner(pms[x][1],right_primer_reads,method='global',matrix=DNAFULL, max_hits=1)

        # ['count', 'end1', 'end2', 'index', 'n_gaps1', 'n_gaps2', 'n_mismatches', 'score', 'seq1', 'seq2', 'start1', 'start2']
        l_mismatches = alignL[0].n_mismatches
        r_mismatches = alignR[0].n_mismatches
        l_score = alignL[0].score
        r_score = alignR[0].score
        if l_score + r_score > hit_score:
            hit_score = l_score + r_score 
            hit_name = x
            hit_l_mismatches = l_mismatches
            hit_r_mismatches = r_mismatches
            hit_l_score = l_score
            hit_r_score = r_score

    reads += [hit_name,hit_l_mismatches,hit_r_mismatches,hit_l_score,hit_r_score]
    return reads




def runner(args,barcodes,primers):
    """Run the Typer 
    
    Args:
        args (TYPE): Description
        barcodes (TYPE): Description
        primers (TYPE): Description
    """
    obj  = FastaReader(args.fasta)

    barcode_result = []
    log.info("Start split reads with given barcodes")
    try:
        barcode_result = [barcode_split(record,barcodes,mismatch=args.mismatchB,mode=args.mode) for record in tqdm(list(obj)) if barcode_split(record,barcodes,mismatch=args.mismatchB,mode=args.mode)]
        barcode_result = [x[0] for x in barcode_result]
    except ValueError:
        log.error("Invalid FASTA file!")
        exit(0)


    log.info("%s reads with match given barcodes found" % len(barcode_result))

    log.info("Start split exons with given primers")
    primer_result = [primer_split(item,primers)for item in tqdm(barcode_result)]
    log.info("Start split exons by primer done!")

    primer_result = sorted(primer_result,key=lambda x:x[1])

    result_dic = {}
    for m,n in groupby(primer_result,key=lambda x:x[1]):
        result_dic[m] = list(n)

    for key in result_dic:
        log.debug("Write result of sample: %s" % key)
        output_sample(key,result_dic[key])





def output_sample(sampleID,content):

    folder_path = os.path.join(args.output,sampleID)

    try:
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
    except OSError:
        if not os.path.isdir(folder_path):
            raise

    record_this = {}
    try:
        for reads in content:
            try:
                record_this[reads[8]] += [FastaRecord(reads[0],reads[2])]
            except:
                record_this[reads[8]] = [FastaRecord(reads[0],reads[2])]
        
        for key in record_this:
            fasta_file = os.path.join(folder_path,key+'.fasta')
            write_records(FastaWriter,record_this[key],fasta_file)

        return True
    except:
        log.error("Error in write output file %s" % sampleID)
        exit(0)





def write_records(writer_class, records, file_name):
    n = 0
    with writer_class(file_name) as w:
        for record in records:
            w.writeRecord(record)
            n += 1

    log.debug("Completed writing %s records into %s" % (n,file_name))





def get_parser(parser=None):

    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('-f','--fasta', metavar='FASTA',
        help="Fasta of Amplicon Analysis or Reads of Insert output")
    add('-b', '--barcode', metavar='BARCODE', help="barcode mapping file" )
    add('--mismatchB', type=int,default=1,help="barcode mismatch")
    add('--mismatchM', type=int,default=3,help="primer mismatch")
    add('-p','--primer', metavar="PRIMER",help="primer mapping file")
    add('--mode', default='slow',help="[slow|fast] mode for barcode search")
    add('-o','--output', default='output',help="output prefix")
    add('--debug', action='store_true', help="Flag to enable Debug mode")

    args = parser.parse_args()

    return args




if __name__ == '__main__':

    args = get_parser()


    # read in barcode
    barcodes = read_barcode(args.barcode)
    log.info("%s barcodes load" % len(barcodes))

    # read in primers
    primers=read_barcode(args.primer)
    log.info("%s primers load" % len(primers))

    # start the main program
    runner(args,barcodes,primers)