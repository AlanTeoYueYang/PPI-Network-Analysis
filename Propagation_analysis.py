from PPI_Network import PPI_Network

class Propagation_analysis(PPI_Network):
	def __init__(self, protein_links, protein_aliases):
		# download protein.links.full.v__.txt and protein.aliases.v__.txt from https://string-db.org/cgi/download.pl
		super().__init__(protein_links, protein_aliases)

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

	def propogation(self, candidate_genes_id, alpha, flow_matrix):
		print("Performing propogation")
		Y0 = np.zeros(len(self.network))
		for gene_id in candidate_genes_ids:
			Y0[self.indices[gene_id]] = 1
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
		for i, v in enumerate(Y_next):
			Y_next[i] = [self.names[i],v]
		rank = sorted(Y_next, key=lambda x: x[1], reverse=True)
		return rank

		return rank

	def intersection(self, disease1_genes, disease2_genes, alpha):
		flow_matrix = PPI_Network.creating_flow_matrix(self)
		ranks = {}
		print("Finding intersection")
		for dis_genes in diseases_genes:
			rank = self.propogation(dis_genes, alpha)
			for gene, rank_no in rank:
				if gene not in ranks:
					ranks[gene] = rank_no
				else:
					ranks[gene] += rank_no
		sorted_ranks = collections.OrderedDict(sorted(d.items(), key=lambda x:x[1]))
		return sorted_ranks
