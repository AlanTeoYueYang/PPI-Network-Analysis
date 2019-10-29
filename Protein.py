class Protein:
	def __init__(self, id, index):
		self.id = id
		self.index = index
		self.neighbours = {}
		self.name = id

	def add_neighbour(self, neighbour, score):
		self.neighbours[neighbour] = score