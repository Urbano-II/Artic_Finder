# Artic_Finder

Requirements:
		Pyhton 3, Blast, Muscle, pandas module.

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
