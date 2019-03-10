import numpy as np
import math
import collections
import matplotlib.pyplot as plt
import networkx as nx

class Protein:
	def __init__(self, id, index):
		self.id = id
		self.index = index
		self.neighbours = {}
		self.name = id

	def add_neighbour(self, neighbour, score):
		self.neighbours[neighbour] = score

class PPI_Network:
	def __init__(self):
		self.network = PPI_Network.create_network(self)
		self.matrix = PPI_Network.create_matrix(self)
		self.RWR_matrix = PPI_Network.creating_RWR_matrix(self)

	def protein_id(id_string):
		try:
			return int(id_string[13:])
		except:
			return id_string

	def create_network(self):
		print("Creating network")
		protein_network = {}
		protein_indices = []
		data = open("C:/Users/Teo/Documents/NetBeansProjects/LKLNetworkBiology/stringdb_v10/9606.protein.links.full.v10.5.txt")
		names = open("C:/Users/Teo/Documents/NetBeansProjects/LKLNetworkBiology/stringdb_v10/9606.protein.aliases.v10.5.txt")
		index = 0
		headers = True
		for i in data:
			if headers == False:
				temp = i.split(" ")
				score = 1
				for j in temp[3:15]:
					score *= (1-int(j)/1000)
				score = 1 - score
				#score = int(temp[15][:len(temp[15])])/1000			
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
					if "HUMAN" in curr_name and type(protein_network[k_id].name) == int:
						protein_network[k_id].name = curr_name
					if temp[2].split(" ")[0] == "BLAST_UniProt_ID":
						protein_network[k_id].name = curr_name
			else:
				headers = False
		return protein_network

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

	def creating_RWR_matrix(self):
		print("Creating RWR matrix")
		temp_matrix = self.matrix.copy()
		D = list(map(lambda x: 1/len(self.network[x].neighbours) if len(self.network[x].neighbours) > 0 else 0, self.indices))
		temp_matrix *= np.array(D)
		return temp_matrix

	def creating_RWR_matrix_1(self):
		print("Creating RWR matrix")
		temp_matrix = self.matrix.copy()
		D = list(map(lambda x: len(test.network[x].neighbours), test.indices))
		for r in range(len(temp_matrix)):
			for c in range(len(temp_matrix)):
				norm = math.sqrt(D[r]*D[c])
				if norm != 0:
					temp_matrix[r][c] /= norm
		return temp_matrix

	def RWR(self, candidate_genes, restart):		
		print("Performing RWR")
		RWR_matrix = self.RWR_matrix
		candidate_genes_ids = list(map(lambda x: PPI_Network.protein_id(x),candidate_genes))
		p0 = np.zeros(len(self.network))
		start_prob = 1/len(candidate_genes)
		for gene_id in candidate_genes_ids:
			p0[self.network[gene_id].index] = start_prob
		L2_norm = 1
		p_next = p0.copy()
		iteration = 0
		while (L2_norm > 0.000001):
			p_prev = p_next.copy()
			p_next = (1-restart)*np.matmul(RWR_matrix, p_prev) + restart*(p0)
			L2_norm = sum((p_next-p_prev)**2)
			iteration += 1
			if (iteration > 1000000):
				break
		d = {}
		lst = list(p_next)
		for i in range(len(lst)):
			if lst[i] not in d:
				d[lst[i]] = []
			d[lst[i]].append(self.indices[i])
		sorted_d = collections.OrderedDict(sorted(d.items(), reverse = True))
		rank = []
		for m in sorted_d:
			for n in sorted_d[m]:
				if n not in candidate_genes_ids:
					rank.append([n,m])
		return rank

	def intersection_RWR(self,disease1_genes, disease2_genes, restart):
		RWR_matrix = self.RWR_matrix
		rank1 = PPI_Network.RWR(self,disease1_genes, restart, RWR_matrix)
		rank2 = PPI_Network.RWR(self,disease2_genes, restart, RWR_matrix)
		print("Finding intersection")
		d = {}
		for i in range(len(rank1)):
			d[rank1[i][0]] = i
		for j in range(len(rank2)):
			if rank2[j][0] in d:
				d[rank2[j][0]] += j
			else:
				d[rank2[j][0]] = j
		sorted_d = collections.OrderedDict(sorted(d.items(), key=lambda x:x[1]))
		intersectn = []
		for k in sorted_d:
			intersectn.append([test.network[k].name, sorted_d[k]])
		return intersectn
		

	def intersection_RWR_revised(self, disease1_genes, disease2_genes):
		
		combined_genes = disease1_genes + disease2_genes
		rank1 = PPI_Network.RWR(self,combined_genes, 0.75, RWR_matrix)
		rank2 = PPI_Network.RWR(self,disease2_genes, 0.75, RWR_matrix)
		max1 = max(list(map(lambda x:x[1], rank1)))
		rank1 = list(map(lambda x:[x[0],x[1]/max1], rank1))
		max2 = max(list(map(lambda x:x[0], rank2)))
		rank2 = list(map(lambda x:[x[0],x[1]/max2], rank2))
		print("Finding intersection")
		d = {}
		for i in rank1:
			d[i[0]] = i[1]
		for j in rank2:
			if j[0] in d and j[1] > d[j[0]]:
				d[j[0]] /= j[1]
		sorted_d = collections.OrderedDict(sorted(d.items(), key=lambda x:x[1], reverse=True))
		intersectn = []
		for k in sorted_d:
			intersectn.append([test.network[k].name, sorted_d[k]])
		return intersectn


	def creating_flow_matrix(self):
		print("Creating propogation matrix")
		temp_matrix = self.matrix.copy()
		D = temp_matrix.sum(axis = 1)
		for r in range(len(temp_matrix)):
			for c in range(len(temp_matrix)):
				norm = math.sqrt(D[r]*D[c])
				if norm != 0:
					temp_matrix[r][c] /= norm
		return temp_matrix	

	def propogation(self, candidate_genes, alpha, flow_matrix):
		print("Performing propogation")
		candidate_genes_ids = list(map(lambda x: PPI_Network.protein_id(x), candidate_genes))
		Y0 = np.zeros(len(self.network))
		for gene_id in candidate_genes_ids:
			Y0[self.network[gene_id].index] = 1
		L2_norm = 1
		iteration = 0
		Y_next = Y0.copy()
		while (L2_norm > 0.000001):
			Y_prev = Y_next.copy()
			Y_next = (alpha)*np.matmul(flow_matrix, Y_prev) + (1-alpha)*(Y0)
			L2_norm = sum((Y_next-Y_prev)**2)
			iteration += 1
			if (iteration > 1000000):
				break
		d = {}
		lst = list(Y_next)
		for i in range(len(lst)):
			if lst[i] not in d:
				d[lst[i]] = []
			d[lst[i]].append(self.indices[i])
		sorted_d = collections.OrderedDict(sorted(d.items(), reverse = True))
		rank = []
		for m in sorted_d:
			for n in sorted_d[m]:
				if n not in candidate_genes_ids:
					rank.append(n)
		return rank

	def intersection_propogation(self,disease1_genes, disease2_genes, alpha):
		flow_matrix = PPI_Network.creating_flow_matrix(self)
		rank1 = PPI_Network.propogation(self,disease1_genes, alpha, flow_matrix)
		rank2 = PPI_Network.propogation(self,disease2_genes, alpha, flow_matrix)
		print("Finding intersection")
		d = {}
		for i in range(len(rank1)):
			d[rank1[i]] = i
		for j in range(len(rank2)):
			if rank2[j] in d:
				d[rank2[j]] += j
		sorted_d = collections.OrderedDict(sorted(d.items(), key=lambda x: x[1]))
		intersectn = []
		for k in sorted_d:
			intersectn.append([test.network[k].name, sorted_d[k]])
		return intersectn

test = PPI_Network()
print("Network constructed")
PD_genes = [338345, 355865, 364204, 340278, 298910, 258080, 327214, 333142, 266087,299138]
T2D_genes = [304895, 385995, 264717, 417970, 415011, 419361, 418735, 427108, 276925, 365735, 155840, 377233, 257570, 257555]

ans = []

for gene in PD_genes:
	temp_genes = PD_genes.copy()
	temp_genes.remove(gene)
	lst = test.RWR(temp_genes,0.75)
	for i in range(len(lst)):
		if lst[i][0] == gene:
			ans.append([test.network[lst[i][0]].name, i, lst[i]])

for k in ans:
	print(k)

'''
combined_genes = T2D_genes + PD_genes
colors = ["orange"]*len(T2D_genes) + ["green"]*len(PD_genes)
nodes = combined_genes.copy()
edges = []
for i in combined_genes:
	if len(test.network[i].neighbours) < 3:
		print(i)
	for j in test.network[i].neighbours:
		edges.append((i,j))
		if i in T2D_genes:
			colors.append("red")
		if i in PD_genes:
			colors.append("blue")

G=nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)
nx.draw(G, node_size = 5, node_color=colors)
plt.show()
'''