#!/usr/bin/env python3
# encoding: utf-8

# We import the modules
import os
import pathlib
import sys
from Bio import Seq
from Bio import SeqIO

# --------- Functions--------------


def directory_iterator(directory):
    '''
    Iterates over a given directory and converts multiple gbk with .gbff 
        termination into a fasta file using gbk_to_fasta.
    Fasta file generated has ">Locus_tag@Species_name" header.
    Creates a data/gbk_to_fasta directory, then saves the output fasta file 
        there.
    Returns converted file path.
    '''

    # Make new dirctories to put converted gbk files there
    new_dir = "data/gbk_to_fasta/"
    filename = "gbk_converted_to_fasta"
    os.makedirs(new_dir, exist_ok=True)

    # In order not to append to an existing file with the same name, we erase
    # the previous one.
    try:
        os.remove(new_dir + filename)
    except:
        pass

    # Iterate over .gbff files in given directory and convert them to fasta
    for entry in os.scandir(directory):

        if entry.path.endswith(".gbff") and entry.is_file():
            gbk_to_fasta(entry.path, new_dir + filename)

    return (new_dir+filename)


def gbk_to_fasta(input_gbk, file):
    '''
    Converts genbank to fasta. Input_file is given .gbff file to convert.
    Output is appended to "file" that is cretated if it didn't exist before
    '''

    with open(input_gbk, "r") as input_handle:

        # Save species name without blanks to use
        for record in SeqIO.parse(input_handle, "genbank"):
            species = record.annotations["organism"].replace(" ", "_")

            for feature in record.features:

                if feature.type == 'CDS':

                    try:

                        # If there is a protein sequence, save in fasta format
                        # (file)
                        if feature.qualifiers['translation'][0] != "":
                            file_out = open(file, "a")
                            file_out.write(
                                ">" +
                                feature.qualifiers['locus_tag'][0] +
                                "@" +
                                species +
                                "\n")
                            file_out.write(
                                feature.qualifiers['translation'][0] +
                                " \n\n")
                            file_out.close()

                    except:
                        pass
