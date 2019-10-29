# PPI-Network-Analysis
Protein-protein interaction network constructed with STRING database

To read the details of this Random-walk-with-restart(RWR) implementation and report regarding the analysis,<br/> please refer to Report.pdf

I added a Propagation algorithm as well, references to it can be found in the report as well.

Briefly, the algorithms generate the top 100 functional genes that are implicated in two diseases.

The rationale for doing so is to understand potential molecular links between the two diseases, allowing for better understanding of the two diseases and consequently, more efficient treatment of one/both diseases.

## Download data and setting up
Download the full links and aliases files from STRING database https://string-db.org/cgi/download.pl

The current latest version of the files are:

- protein.links.full.v11.0.txt.gz 
- protein.aliases.v11.0.txt.gz

Then extract out the files into text files

You will need at least python 3.5 and the following libraries:
- numpy
- csv
- argsparse

For the disease genes text, you will need two rows (row 1 for one disease and row 2 for the other)

Each row is a list of STRING identifiers for the genes of the particular disease, separated by commas

Refer to the sample_genes.txt where the first row is the PD related genes and second row is the T2D related genes

## Network Architecture

The network is stored in an adjacency matrix:
- proteins --> nodes
- links --> edges
- confidence scores --> weights of edges

## Running the algorithm
The algorithm has to be run on the command line using the run.py file.
The arguments are as follows:

**positional arguments**:
- links (Protein links text file from STRING)
- alias (Protein alias text file from STRING)
- type_of_analysis (Type of analysis: RWR or Propagation)
- disease_genes (Disease genes, see sample_genes.txt)
- param (Parameter for analysis; restart for RWR, alpha for Propagation)
- output_file (Output filename; csv file)

**optional arguments**:
- -h, --help (show this help message and exit)
  
An example:

python run.py protein.links.full.v11.0.txt protein.aliases.v11.0.txt RWR sample_genes.txt output.csv

This will run the RWR algorithm to find high ranking genes in PD and T2D and write the top 100 genes into output.csv

## Customization
To customize your own analysis, simply inherit the PPI_Network into your new class <br/> and add your own algorithms to the new class.

The arguments in the run.py has been set for RWR and Propagation only. Please edit accordingly.
