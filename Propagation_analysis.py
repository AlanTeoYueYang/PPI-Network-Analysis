import PPI_Network

class Propagation_analysis(PPI_Network):
	def __init__(self, protein_links, protein_aliases):
		# download protein.links.full.v__.txt and protein.aliases.v__.txt from https://string-db.org/cgi/download.pl
		super(Propagation_analysis, self).__init__(protein_links, protein_aliases)

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

	def intersection(self, disease1_genes, disease2_genes, alpha):
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
