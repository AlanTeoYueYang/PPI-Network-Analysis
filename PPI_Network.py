import numpy as np

class PPI_Network:
	def __init__(self, protein_links, protein_aliases): 
		# download protein.links.full.v__.txt and protein.aliases.v__.txt from https://string-db.org/cgi/download.pl
		self.matrix = np.zeros((1,1))
		self.indices = {}
		self.names = []
		self.add_network(protein_links, protein_aliases)

	def protein_id(id_string):
		return id_string

	def add_protein_to_matrix(self):
		col = np.zeros((len(self.matrix), 1))
		self.matrix = np.concatenate((self.matrix, col), 1)
		row = np.zeros((1, len(self.matrix[0])))
		self.matrix = np.concatenate((self.matrix, row), 0)

	def add_network(self, protein_links, protein_aliases):
		print("Creating network")
		data = open(protein_links, 'r')
		names = open(protein_aliases, 'r')
		index = 0
		headers = True
		for i in data:
			if headers == False:
				temp = i.split(" ") # spilt data into list
				score = 1 
				for j in temp [3:15]: # setting channels to channel of interest; can be altered accordingly
					score *= (1-int(j)/1000) # calculate confidence scores of links
				score = 1 - score		
				protein_a_id = temp[0] 
				protein_b_id = temp[1]
				if protein_a_id not in self.indices:
					self.indices[protein_a_id] = index
					protein_a_index = index
					if index > 0: self.add_protein_to_matrix()
					index+=1
					self.names.append([])
				if protein_b_id not in self.indices:
					self.indices[protein_b_id] = index
					protein_b_index = index
					self.add_protein_to_matrix()
					index += 1
					self.names.append([])
				if score >= 0.4:
					self.matrix[protein_a_index][protein_b_index] = score
			else:
				headers = False
		print("Adding names")
		headers = True
		for k in names:
			if headers == False:
				temp = k.split("\t")
				curr_name = temp[1]
				if temp[0] in self.indices:
					# annotating names to proteins, criteria can be set accordingly
					idx = self.indices[temp[0]]
					if "HUMAN" in curr_name:
						self.names[idx] = curr_name
					if temp[2].split(" ")[0] == "BLAST_UniProt_ID":
						self.names[idx] = curr_name
			else:
				headers = False
