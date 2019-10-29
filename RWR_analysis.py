import PPI_Network

class RWR_analysis(PPI_Network):
	def __init__(self, protein_links, protein_aliases):
		# download protein.links.full.v__.txt and protein.aliases.v__.txt from https://string-db.org/cgi/download.pl
		super(RWR_analysis, self).__init__(protein_links, protein_aliases)

	def creating_RWR_matrix(self):
		print("Creating RWR matrix")
		temp_matrix = self.matrix.copy()
		# normalize scores
		D = list(map(lambda x: 1/len(self.network[x].neighbours) if len(self.network[x].neighbours) > 0 else 0, self.indices))
		temp_matrix *= np.array(D)
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

	def intersection(self,disease1_genes, disease2_genes, restart):
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