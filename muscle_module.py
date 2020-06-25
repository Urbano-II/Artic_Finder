#!/usr/bin/env python3
# encoding: utf-8

# We import the module

import os


# --------Functions------------

def muscler(file_in, file_out):
    '''
    Function that uses muscle to align a fasta file (file_in), then makes a 
        Neigbour Joining tree (file_out) 
    Creates a temporary alignment file /tmp/muscle_out

    '''
    # Creates alignment
    cmd = "muscle -in " + file_in + " -out /tmp/muscle_out 2> /dev/null"
    os.system(cmd)

    # Creates tree
    cmd2 = "muscle -maketree -in /tmp/muscle_out -out " + \
        file_out + " -cluster neighborjoining 2> /dev/null"
    os.system(cmd2)
