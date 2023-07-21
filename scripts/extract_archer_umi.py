################################################################################################
# extracts UMI barcodes from paired end reads derived from UMI tagged archer fusionplex panel
# 
# Expected UMI structure: 
# Forward read 5'- 8bp UMI - 13bp spacer - target sequence - 3'
# Reverse read 5'- target sequence - 3'
# The final UMI which is appened to both read headers consists of the 8bp UMI sequence from the 
# forward read and the first 10bp from the forward target sequence.
# 
# Requires python >= 3.0 and biopython.
# 
# usage: python extract_archer_umi.py -f1 read_1.fastq -f2 read_2.fastq -o output_prefix
#
################################################################################################

import gzip
import os, sys, re
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from argparse import ArgumentParser

def parse_command(command):
    parser = ArgumentParser(description='extracts UMI barcodes from paired end reads derived from UMI tagged archer \n'
                                        'fusionplex panel.The final UMI which is appened to both read headers consists  \n'
                                        'of the 8bp UMI sequence from the forward read and the first 10bp from the \n'
                                        'forward target sequence')
    parser.add_argument('-f1', '--fastq_1', required=True, help='fastq.gz R1 file')
    parser.add_argument('-f2', '--fastq_2', required=True, help='fastq.gz R2 file')
    parser.add_argument('-o', '--output', required=True, help='output prefix for umi tagged reads')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)

input_command = parse_command(sys.argv[1:])

input_file = input_command.fastq_1
input_file2 = input_command.fastq_2
output_filename = input_command.output + "_R1_marked.fastq.gz"
output_filename2 = input_command.output + "_R2_marked.fastq.gz"

output_handle = gzip.open(output_filename, "wt")
output_handle2 = gzip.open(output_filename2, "wt")
reads = FastqGeneralIterator(gzip.open(input_file, 'rt'))
reads2 = FastqGeneralIterator(gzip.open(input_file2, 'rt'))

for (f_id, f_seq, f_qual), (r_id, r_seq, r_qual) in zip(reads, reads2):
    umi = f_seq[:8]
    umi2 = str(f_seq[21:31])
    out_f_seq = str(f_seq[21:])
    out_f_qual = str(f_qual[21:])
    out_id = r_id.split()[0] + "_" + umi + umi2
    output_handle.write("@%s\n%s\n+\n%s\n" % (out_id, out_f_seq, out_f_qual))
    output_handle2.write("@%s\n%s\n+\n%s\n"% (out_id, r_seq, r_qual))
    
output_handle.close()
output_handle2.close()
