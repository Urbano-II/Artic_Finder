#!/usr/bin/env python3
# encoding: utf-8

# We import modules

import os
import sys
import pandas as pd
from Bio import SeqIO

# -------Functions-----------


def blaster(query, subject, identity, coverage, file):
    '''
    Does blast with a given subject and query and filters by given coverage,
    and identity cutoffs.
    Creates and overwrites file given.
    '''

    # Blast comand line for bash
    cmd = "blastp -query " + query + " -subject " + subject + " -evalue 0.00001" \
        ' -outfmt "6 sseqid pident qcovs sseq" -out /tmp/blast_to_filter.txt'
    os.system(cmd)

    # Read blast file with pandas and filter with identity and coverage
    to_filter = pd.read_csv('/tmp/blast_to_filter.txt',
                            header=None, index_col=False, sep='\t')
    filtered = to_filter.loc[(to_filter[1] >= identity)
                             & (to_filter[2] >= coverage)]

    # Erases previous blast file
    cmd2 = "rm " + file + " 2>/dev/null"
    os.system(cmd2)

    # Read query file and write contents onto new file (file)
    query_file = open(query, "r")
    query_seq = query_file.read().strip()
    query_file.close()
    file_file = open(file, "w")
    file_file.write(query_seq+"\n\n")
    file_file.close()

    # Write filtered contents onto file with a fasta format
    for index, row in filtered.iterrows():
        to_write = open(file, 'a')
        to_write.write(">" + row[0] + "\n")
        to_write.write(row[3] + "\n\n")
        to_write.close()


def look_in_gbk(file_in_blast, converted_gbk, file_out):
    '''
    Looks for blast results in converted gbk file, and creates a new file
    with the complete protein sequence.
    '''

    # Erases previous file
    cmd2 = "rm " + file_out + " 2>/dev/null"
    os.system(cmd2)

    # File with complete sequence form gbk converted
    complete_seq = open(file_out, 'a')

    # Count to read only once first blast record (query sequence)
    count = 0

    with open(converted_gbk, 'r') as complete_in:

        for record_complete in SeqIO.parse(complete_in, "fasta"):

            with open(file_in_blast, 'r') as blast_fasta:

                for record_blast in SeqIO.parse(blast_fasta, "fasta"):

                    # Read query from original blast and write contents onto
                    # new file (file)
                    if count == 0:
                        complete_seq.write(">" +
						                   str(record_blast.id) +
                                           "\n" +
										   str(record_blast.seq) +
										   "\n\n")
                        count = 1

                    # Copies complete sequence if ID is the same as blast hit
                    if str(record_blast.id) == str(record_complete.id):
                        complete_seq.write(">" +
						                   str(record_complete.id) +
                                           "\n" +
										   str(record_complete.seq) +
										   "\n\n")

    complete_seq.close()
