#!/usr/bin/env python3
# encoding: utf-8

# We import modules
import os
import re
from Bio.ExPASy import Prosite, Prodoc
from Bio import Seq
from Bio import SeqIO

# ---Functions---------


def prositer(query_open, protein_domains, file_out):
    '''
    Function that finds protein domains in a query protein (from a given file
    query_open).
    Protein domains are extracted from a prosite database (protein_domains)
    (.dat file)
    Makes a file with information about the prosite domain and the blast domain
	(file_out)
    '''

    try:
        os.remove(file_out)
    except:
        pass

    with open(query_open, "r") as query_handle:

        # Open file with fasta format, we are going to look for the prosite
		# domains in this file
        for record_fasta in SeqIO.parse(query_handle, "fasta"):
            query_fasta = str(record_fasta.seq)  # Query sequence

            # Open prosite file and parse it for search
            handle = open(protein_domains, "r")
            records = Prosite.parse(handle)

            for record in records:
                # We traduce from prosite format to something readable for re
				# module
                protein_dom = record.pattern.replace(".", "").replace("-", "")
                protein_dom = protein_dom.replace("{", "[^").replace("}", "]")
                protein_dom = protein_dom.replace("(", "{").replace(")", "}")
                protein_dom = protein_dom.replace("x", ".")
                protein_dom = str(protein_dom)

                # If domain is found by search, information is saved
                if re.search(protein_dom, query_fasta) and protein_dom != "":
                    prot_info = open(file_out, 'a')
                    prot_info.write("Blast hit:"+record_fasta.id)
                    prot_info.write("\nname:" + record.name)
                    prot_info.write("\naccession:" + record.accession)
                    prot_info.write("\ndescription:" + record.description)
                    prot_info.write("\npattern:" + protein_dom)
                    prot_info.write("\nDomain:" +
                                    str(re.findall(protein_dom, query_fasta)[0]) +
                                    "\n\n")
                    prot_info.close()
