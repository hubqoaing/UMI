#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os

import cPickle as pickle
import numpy as np
from optparse   import OptionParser

import htsaminfo


def prepare_optparser():
    usage ="""usage: %s [options] 
This file can read the result of HTseq result and then get the number of UMI.
With the count of UMI, the TPM value could be calculated for each gene in the
given gene list as well as the ercc list.

Using -h or --help for more information

Pre-processing:
samtools view -o /data/Analysis/huboqiang/b2_bar.sam                        \\
    /date/dongji/project/Huyuqiong/task10_2015_11_25/mismatch0/01.Tophat/b2_bar2/accepted_hits.bam &&\\
python /data/Analysis/huboqiang/software/anaconda/bin/htseq-count           \\
    -s no -f sam -a 10 -o /data/Analysis/huboqiang/b2_bar.umi.sort.gene.sam \\
    /data/Analysis/huboqiang/b2_bar.sam                                     \\
    /date/dongji/database/Database_RNA_v2/mm10/refGene.gtf                  \\
    >/data/Analysis/huboqiang/b2_bar.deseq.xls

Example:
    python %s -g genelist.xls -e list_gene_ercc.xls input.sam
   
    """ % (sys.argv[0],sys.argv[0])

    description = "TPM value for a given HTseq result file. "
    
    optparser = OptionParser(version="%s v0.2 20141130" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )
    optparser.add_option(
        "-g", "--geneList", default="list_gene.xls",
        help="\nGene list with given order. [default: %default]"
    )
    optparser.add_option(
        "-e", "--erccList", default="list_ercc.xls",
        help="\nERCC list with given order. [default: %default]"
    )

    optparser.add_option("-h","--help", action="help",
        help="\nShow this help message and exit."
    )
    return optparser


def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        input_sam = args[0]
        geneList = options.geneList
        erccList = options.erccList
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)

    m_tpmInfo = htsaminfo.HTSeqSamInfo(input_sam, geneList, erccList)
    m_tpmInfo.s1_SamToUMI()
    m_tpmInfo.s2_outUMIList()
    m_tpmInfo.s3_UMI2TPM()


if __name__ == '__main__':
    main()