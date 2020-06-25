#!/usr/bin/env python3
#
# Manuel Fernández 2020

# We import the modules
import sys
import re
import os
from Bio import SeqIO

# We import the package modules
import gbk_to_fasta_module as gbk
import blast_module as bl
from muscle_module import muscler
from prosite_module import prositer

# ----- functions -------

# Function for error controlling


def help():

    message = (

        """
	#_Antartic Finder_#

	Usage: python main.py <genebank_directory> <query_fasta_file> 
		   <prosite_database.dat> [identity_cutoff] [coverage_cutoff] [-h] 

		-h: help
		- genebank_directory: a directory that has one or more genebank files
		  with .gbff termination to be used as subject
		- query_fasta_file: a single file that has one or more query in fasta
		  format
		- prosite_database.dat: a file that has protein domains in prosite
		  format 
		- identity_cutoff: a number, may be decimal, blast hits with lower
		  identity are filtered
		- coverage_cutoff: a number, may be decimal, blast hits with lower
		  coverage are filtered 

	Description: 

		This program does Blastp onto multiple genbank files. Hit results are
		then alligned with muscle.
		Alignment is used to create a Neighbour Joining tree with muscle.
		Also, it searches for domains provided form a prosite database within
		the blast hits.
		Makes "results" directory which has a directory for each query protein
		(same name) and contains the following:
		- A fasta file with the query sequence.
		- Blast hits in a fasta file with name "blast_query.fa".
		- Complete sequence looked onto converted gbk from blast hits with name
		  "blast_complete_query.fa"
		- A N-J tree in newick format with the name "tree_query.fa".
		- A file with all the found protein domains with the name 
		  "prosite_query.txt". Contains blast hit name, prosite name,
		  accession, description, pattern and blast hit domain.

		Also makes a data/gbk_to_fasta directory that has all genebank files
		fused into a file with fasta format (gbk_converted_to_fasta)
		Alignment is stored in the temporary directory /tmp/muscle_out
	
	Made with tears by Manuel Fernández



	"""
    )

    print(message)
    sys.exit()


def warning():

    message = (

        """
	Usage: python main.py <genebank_directory> <query_fasta_file>
	       <prosite_database.dat> [identity_cutoff] [coverage_cutoff] [-h] 

		-h: help
		- genebank_directory: a directory that has one or more genebank files
		  with .gbff termination to be used as subject
		- query_fasta_file: a single file that has one or more query in fasta
		  format
		- prosite_database.dat: a file that has protein domains in prosite
		  format 
		- identity_cutoff: a number, may be decimal, blast hits with lower
		  identity are filtered
		- coverage_cutoff: a number, may be decimal, blast hits with lower
		  coverage are filtered

	"""

    )

    print(message)
    sys.exit()


def is_fasta(file_to_check):
    """
    Controls that the query input is a fasta, if not, gives a warning message.
    """
    checking = open(file_to_check, "r")
    line = checking.read()
    checking.close()

    if not line.startswith(">"):

        print("\nError: provided " + file_to_check + " is not in fasta format!")
        warning()


def is_gbk(directory_to_check):
    """
    Controls that directory for genebank files is a directory. Only controls
    files with .gbff termination 
    Checks if each file has genebank format. 
    If not, gives a warning message.
    """
    check = 0  # Used later to check if directory has files in it

    # Check if the path given is a directory
    if os.path.isdir(directory_to_check) == False:

        print("\nError: " + directory_to_check + " is not a directory")
        warning()

    # Check if the files within the directory are genebank format (start with
    # LOCUS)
    # If a file does no end with .gbff, it is ignored and the program 
    # continues,but a warning message is printed
    for file_to_check in os.scandir(directory_to_check):

        if file_to_check.path.endswith(".gbff") and file_to_check.is_file():

            check = 1  # This indicates that directory has at least one file

            checking = open(file_to_check, "r")
            line = checking.read()
            checking.close()

            if not line.startswith("LOCUS"):

                print("\nError: provided " + file_to_check.path +
                      " is not in genebank format!")
                warning()
        else:

            print("\n" + file_to_check.path +
                  "does not have .gbff termination, so it was ignored")

    # Checks that directory has files in it
    if check == 0:
        print("No genebank files were found in given in " + directory_to_check)
        warning()


def is_prosite(file_to_check):
    """
    Controls that the query input has prosite format, if not, gives a warning
	message.
    """
    checking = open(file_to_check, "r")
    line = checking.read()
    checking.close()

    if not line.startswith("CC"):

        print("\nError:provided " + file_to_check +
              " is not in prosite format!")
        warning()


def input_control():
    """
    Controls all exceptions. 

    """

    arguments = sys.argv
    # We iterate over all the arguments
    for argument in arguments:

        if argument == "-h":
            help()

    # Check if correct number of arguments given
    if (len(arguments) < 4 or len(arguments) > 6):
        print("\nError: Icorrect number of arguments, must at least be 3, " +
              "not greater than 5")
        warning()

    # If there is no argument -h, we assign arguments to variable
    gbk_directory_file_path = arguments[1]
    query_proteins_file_path = arguments[2]
    prosite_database_file_path = arguments[3]

    # If cutoff values are given, we assign arguments to variable; if not,
    # default values are used
    try:
        identity_cutoff = arguments[4]
        default_id = False
    except:
        default_id = True
        identity_cutoff = 30

    try:
        coverage_cutoff = arguments[5]
        default_coverage = False
    except:
        default_coverage = True
        coverage_cutoff = 50

    # Check if identity and coverage cutoffs are numbers
    if (default_id == False) and (default_coverage == False):

        if ((identity_cutoff.replace(".", "", 1).isdigit() == False) 
            and (coverage_cutoff.replace(".", "", 1).isdigit() == False)):
            print("\nError: identity and coverage cutoffs are not a number")
            warning()

    if (default_id == False):

        if (identity_cutoff.replace(".", "", 1).isdigit() == False):
            print("\nError: identity cutoff is not a number")
            warning()

    if (default_coverage == False):

        if (coverage_cutoff.replace(".", "", 1).isdigit() == False):
            print("\nError: coverage cutoff is not a number")
            warning()

    # Check if the files are the correct type
    is_gbk(gbk_directory_file_path)  # Checks if genebank format
    is_fasta(query_proteins_file_path)  # Checks if fasta format
    is_prosite(prosite_database_file_path)  # Checks if prosite format

    return(gbk_directory_file_path, query_proteins_file_path,
           prosite_database_file_path, identity_cutoff, coverage_cutoff)


def completion(query_list, identity_cutoff, coverage_cutoff):
    '''
    Prints completion message
    '''
    print("Program completed!\n")

    print("Results directories generated:")
    for entry in query_list:
        print("\t" + entry)

    print("\nIdentity cutoff used: " + str(identity_cutoff))
    print("Coverage cutoff used: " + str(coverage_cutoff))


def waiting():
    '''
    Prints waiting message
    '''
    print("")
    # Waiting message
    print("\nPlease wait, this may take a while, patience is a virtue and " +
          "you are virtuous\n")


# ---Main Function

def main():

    # To control arguments
    arguments_from_input = input_control()

    # We retrieve arguments from input_control
    gbk_directory_file_path = arguments_from_input[0]
    query_proteins_file_path = arguments_from_input[1]
    prosite_database_file_path = arguments_from_input[2]
    identity_cutoff = float(arguments_from_input[3])
    coverage_cutoff = float(arguments_from_input[4])

    # Tells us to wait
    waiting()

    # We convert gbk to fasta format, and save destination path onto variable
    # gbk_converter
    # Bear in mind that it only converts files with .gbff termination, the 
    # rest are not even read
    gbk_coverted = gbk.directory_iterator(gbk_directory_file_path)
    print("Genebank conversion completed")
    print("Saved in: " + gbk_coverted + "\n")

    # Make a new directoy "results" to make subdirectories for each query
    # protein
    try:
        os.mkdir("results")
    except:
        pass

    # List for the names of the query directories created
    query_list = []

    # We iterate over each protein in the query file.
    # Do blast, alignment and domain search for each one.
    # Save the resuslts in a directory named after individual query protein.
    with open(query_proteins_file_path, "r") as query_prot:

        for record_query in SeqIO.parse(query_prot, "fasta"):
            query_protein_name = str(record_query.id)

            # Make direcory for each protein
            protein_directory_path = "results/" + query_protein_name
            query_list.append(protein_directory_path)  # Save dircetory name
            try:
                os.mkdir("results/" + query_protein_name)

            except:
                pass

            # Save protein query sequence in a file named query_file
            query_file = open(protein_directory_path + "/query_file.fa", "w")
            query_file.write(">" + query_protein_name +
                             "\n" + str(record_query.seq))
            query_file.close()

            # We do blast of protein query in genebank converted to fasta file
            # (gbk_converted)
            path_to_blast = protein_directory_path + "/query_file.fa"
            path_out_blast = protein_directory_path + \
                             "/blast_" + query_protein_name + ".fa"
            bl.blaster(path_to_blast, gbk_coverted, identity_cutoff,
                       coverage_cutoff, path_out_blast)

            # We look for blast hist in converted gbk
            path_out_looker = protein_directory_path + \
                "/blast_complete_" + query_protein_name + ".fa"
            bl.look_in_gbk(path_out_blast, gbk_coverted, path_out_looker)

            # We do the alignment and tree using muscle
            try:
                path_out_muscle = protein_directory_path + "/tree_" + \
                                  query_protein_name + ".nw"
                muscler(path_out_blast, path_out_muscle)

            except:
                print("No tree made for " + query_protein_name +
                      "Probably beacause there are no blast hits")

            # We do the domain search of blast hits and query using prosite
            # database
            path_out_prositer = protein_directory_path + "/prosite_" + \
                                query_protein_name + ".txt"
            prositer(path_out_looker, prosite_database_file_path,
                     path_out_prositer)

        # Program completion message
        completion(query_list, identity_cutoff, coverage_cutoff)


if __name__ == '__main__':
    main()
