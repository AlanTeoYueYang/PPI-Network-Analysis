from PPI_Network import PPI_Network 

class RWR_analysis(PPI_Network):
	def __init__(self, protein_links, protein_aliases):
		# download protein.links.full.v__.txt and protein.aliases.v__.txt from https://string-db.org/cgi/download.pl
		super().__init__(self, protein_links, protein_aliases)

	def RWR(self, candidate_genes_ids, restart):		
		print("Performing RWR")
		RWR_matrix = np.array(list(map(lambda x: x/len(x[x>0]), self.matrix))) # normalize scores
		p0 = np.zeros(len(self.matrix))
		start_prob = 1/len(candidate_genes_ids)
		for gene_id in candidate_genes_ids:
			p0[self.indices[gene_id]] = start_prob
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
		for i, v in enumerate(p_next):
			p_next[i] = [self.names[i],v]
		rank = sorted(p_next, key=lambda x: x[1], reverse=True)
		return rank

	def intersection(self, diseases_genes, restart):
		RWR_matrix = self.RWR_matrix
		ranks = {}
		print("Finding intersection")
		for dis_genes in diseases_genes:
			rank = self.RWR(dis_genes, restart)
			for gene, rank_no in rank:
				if gene not in ranks:
					ranks[gene] = rank_no
				else:
					ranks[gene] += rank_no
		sorted_ranks = collections.OrderedDict(sorted(d.items(), key=lambda x:x[1]))
		return sorted_ranks
