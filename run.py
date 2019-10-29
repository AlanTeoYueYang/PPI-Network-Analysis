import argparse
import RWR_analysis
import Propagation_analysis
import read_disease_genes
import csv

def main(links, aliases, type_of_analysis, disease_genes, param, output_file):
	network = None
	if type_of_analysis == 'RWR':
		network = RWR_analysis(links, aliases)
	elif type_of_analysis == 'Propagation':
		network = Propagation_analysis(links, aliases)
	else:
		print("Type of analysis not available")
	if network != None:
		genes = open(disease_genes, 'r').readline()
		genes1 = genes[0].split(',')
		genes2 = genes[1].split(',')
		output_data = network.intersection(genes1, genes2, param)
		output = open(output_file, 'w+')
		writer = csv.writer(output)
		writer.writerow(['Gene Name','Rank'])
		for row in output_data:
			writer.writerow(row)
		writer.close()
		print('Completed')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='PPI Network Analysis')
	parser.add_argument('protein_links', metavar='links', type=str, help='Protein links text file from STRING')
	parser.add_argument('protein_aliases', metavar='alias', type=str, help='Protein alias text file from STRING')
	parser.add_argument('type_of_analysis', type=str, help='Type of analysis: RWR or Propagation')
	parser.add_argument('disease_genes', type=str, help='Disease genes, see sample_genes.txt')
	parser.add_argument('param', type=int, help='Parameter for analysis; restart for RWR, alpha for Propagation')
	parser.add_argument('output_file', type=str, help='Output filename; csv')
    args = parser.parse_args()

    main(args.links, args.protein_aliases, args.type_of_analysis, args.disease_genes, args.param, args.output_file)