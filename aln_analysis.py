#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Summary
"""



from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import Seq
from Bio import SeqRecord
import sys
import numpy as np
import argparse
import logging





def split_paired_barcode():

    with open(sys.argv[1]) as f:
        f.readline()
        rows = [x.strip().split() for x in f]
        d = dict((x[0],x[1:]) for x in rows)

    rs = SeqIO.parse(sys.argv[2],'fasta')

    out = {}

    for r in rs:
        barcode_forward = str(r.seq[:16])
        barcode_reverse = str(r.seq[-16:].reverse_complement())
        for k,v in d.items():
            if barcode_forward == v[0] and barcode_reverse == v[1]:
                item = (r.id + "_" + barcode_forward + ":" + barcode_reverse, r.seq[16:-16])
            elif barcode_forward == v[1] and barcode_reverse == v[0]:
                item = (r.id+"_" + barcode_forward + ":" + barcode_reverse+'_reverse',r.seq.reverse_complement()[16:-16])
            else:
                item = ()
            if item:
                if out.get(k):
                    out[k].append(item)
                else:
                    out[k] = [item]
                break

    duiying = open("raw-vs-new.txt",'w')
    newid = 0

    for k,v in out.items():
        with open(k+'.fa','w') as handle:
            for r in v:
                handle.write('>seq_%s\n%s\n' % (newid,r[1]))
                duiying.write('%s\tseq_%s\n' % (r[0],newid))
                newid +=1









def aln_to_matrix(summary_align):
    """Summary
    
    Args:
        summary_align (TYPE): Description
    
    Returns:
        TYPE: Description
    """
    consensus = summary_align.dumb_consensus()
    pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['N'])
    return pssm


def check_real_insert(thispos,threshold = .7):
   

    sum_insert = thispos['-']

    return sum_insert/sum(thispos.values())>=threshold




def main():
    """Summary
    """

    
    split_paired_barcode(sys.argv[1],sys.argv[2])

    
    aln = AlignIO.read(sys.argv[1],'fasta')
    summary_align = AlignInfo.SummaryInfo(aln)
    pssm = aln_to_matrix(summary_align)

    real_insert_pos = set()

    
    
    for record in list(aln):
        seq_changed = ''
        for i in range(aln.get_alignment_length()):
        # for i in np.arange(180,190):

            # the residue letter at the specified position.
            residue = pssm.get_residue(i)
            # the residue letter frequency at the specified position.
            res_freq = pssm[i]
            # print i+1,res_freq
            
            residue_this_pos = str(record.seq[i])
            # print("%s\t%s\t%s\t%s" % (record.id,i,residue,residue_this_pos))
            if check_real_insert(res_freq):
                real_insert_pos.add(i)
                next
            else:
               # print ("%s -> %s" % (residue,residue_this_pos))
                if residue != 'X':
                    seq_changed += residue
                else:
                    seq_changed += residue_this_pos


        # print
        record.seq = Seq.Seq(seq_changed)
        print record.seq
        print len(record.seq)

    # real_insert_pos 存储被认为是错误的indel位点位置
    # print real_insert_pos

            # seq.seq[0] = residue
            # print seq.seq




def get_parser(parser=None):

    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('-f','--fasta', metavar='INPUT',
        help="Fasta of Amplicon Analysis or Reads of Insert output")
    add('-b', '--barcode', metavar='BARCODE', 
        help="barcode mapping file" )
    add('--debug', action='store_true',
        help="Flag to enable Debug mode")
    args = parser.parse_args()

    return args





if __name__ == '__main__':

    args = get_parser()

    print args

    # main()

