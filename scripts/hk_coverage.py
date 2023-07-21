#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:46:20 2020

@author: jbyoung
"""

from argparse import ArgumentParser
import os
import sys
import pysam

def parse_command(command):
    """Sets up and parses input arguemnts"""
    parser = ArgumentParser(description='counts number of reads which overlap the housekeeping gene exons from the  \n'
                                        'fusionplex lung panel panel. only primary alignments are considered  \n')
    parser.add_argument('-b', '--bamfiles', nargs='+', required=True, help='list of bam files')
    parser.add_argument('-o', '--outfile', required=True, help='name of output tsv file')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return args

def count_hk_cov(bamfile, exons):
    """Parse through bam file, count only reads which are aligned to one of
    the HK exon regions, are mapped in proper pairs and not secondary or
    supplementary alignments. divide counts by two for final coverage"""
    samfile = pysam.AlignmentFile(bamfile, "rb")
    readcount = 0
    for exon in exons:
        for read in samfile.fetch(exon[0], exon[1], exon[2]):
            if read.is_secondary:
                pass
            elif read.is_supplementary:
                pass
            elif read.is_unmapped:
                pass
            elif read.mate_is_unmapped:
                pass
            elif not read.is_proper_pair:
                pass
            else:
                readcount += 1
                fragcount = round(readcount / 2)

    return (bamfile, readcount, fragcount)


def main():
    """main function"""
    input_command = parse_command(sys.argv[1:])
    bamfiles = input_command.bamfiles
    outfile_name = input_command.outfile

    exon_pos = [("chr19", 59063626, 59063805), ("chr19", 59063422, 59063552),
                ("chr19", 34890112, 34890240), ("chr19", 34890461, 34890536),
                ("chr3", 128525215, 128525433), ("chr3", 128526386, 128526514),
                ("chr9", 35059490, 35059798), ("chr9", 35059061, 35059216)]

    results = []
    for file in bamfiles:
        results.append(count_hk_cov(file, exon_pos))

    with open(outfile_name, 'w+') as out:
        out.write("Sample\tReadcount\tFragcount\n")
        for sample in results:
            out.write("%s\t%d\t%d\n" % (os.path.basename(sample[0]), sample[1], sample[2]))

if __name__ == "__main__":
    main()
