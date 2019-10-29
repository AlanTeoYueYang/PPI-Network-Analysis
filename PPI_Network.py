import numpy as np
import Protein

class PPI_Network:
	def __init__(self, protein_links, protein_aliases): 
		# download protein.links.full.v__.txt and protein.aliases.v__.txt from https://string-db.org/cgi/download.pl
		self.network = PPI_Network.create_network(self, protein_links, protein_aliases)
		self.matrix = PPI_Network.create_matrix(self)
		self.RWR_matrix = PPI_Network.creating_RWR_matrix(self)

	def protein_id(id_string):
		return id_string

	def create_network(self, protein_links, protein_aliases):
		print("Creating network")
		protein_network = {}
		protein_indices = []
		data = open(protein_links)
		names = open(protein_aliases)
		index = 0
		headers = True
		for i in data:
			if headers == False:
				temp = i.split(" ") # spilt data into list
				score = 1 
				for j in temp [3:15]: # setting channels to channel of interest; can be altered accordingly
					score *= (1-int(j)/1000) # calculate confidence scores of links
				score = 1 - score		
				protein_a_id = PPI_Network.protein_id(temp[0]) 
				protein_b_id = PPI_Network.protein_id(temp[1])
				if protein_a_id not in protein_network:
					protein_network[protein_a_id] = Protein(protein_a_id, index)
					protein_indices.append(protein_a_id)
					index+=1
				if protein_b_id not in protein_network:
					protein_network[protein_b_id] = Protein(protein_b_id, index)
					protein_indices.append(protein_b_id)
					index+=1
				if score >= 0.4:
					protein_network[protein_a_id].add_neighbour(protein_b_id,score)
					protein_network[protein_b_id].add_neighbour(protein_a_id,score)
			else:
				headers = False
		self.indices = protein_indices
		print("Adding names")
		headers = True
		for k in names:
			if headers == False:
				temp = k.split("\t")
				k_id = PPI_Network.protein_id(temp[0])
				curr_name = temp[1]
				if k_id in protein_network:
					# annotating names to proteins, criteria can be set accordingly
					if "HUMAN" in curr_name and type(protein_network[k_id].name) == int:
						protein_network[k_id].name = curr_name
					if temp[2].split(" ")[0] == "BLAST_UniProt_ID":
						protein_network[k_id].name = curr_name
			else:
				headers = False
		return np.array(protein_network)

	def create_matrix(self):
		print("Creating matrix")
		matrix = np.zeros(shape=(len(self.network),len(self.network)))
		for p_id in self.network:
			p = self.network[p_id]
			i = p.index  
			for neighbour_id in p.neighbours:
				neighbour = self.network[neighbour_id]
				score = p.neighbours[neighbour_id]
				j = neighbour.index
				matrix[i,j] = score
		return matrix
