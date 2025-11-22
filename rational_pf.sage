def to_exp_nozero(p):
	t = p.to_exp()
	t.reverse()
	while t and t[-1] == 0:
		t.pop()
	t.reverse()
	return t

def p_to_dw(p, h, v):
	if p:
		q = p.conjugate()
		horizontal = [h - q[0]] + to_exp_nozero(p)
		vertical = [v - p[0]] + to_exp_nozero(q)
		horizontal = horizontal[::-1]
		t = []
		for i in range(len(horizontal)):
			t += vertical[i]*[1] + horizontal[i]*[0]
		return t
	else:
		return [1]*v+[0]*h

class RationalPF:
	def __init__(self, diagram, h, v, labels):
		self.diagram = diagram
		# self.dyckword = p_to_dw(diagram, h, v)
		self.labels = labels
		self.horizontal = h
		self.vertical = v
		if diagram:
			self.full = [v - diagram[0]] + to_exp_nozero(diagram.conjugate())
		else:
			self.full = [v]

	def area(self):
		h = self.horizontal
		v = self.vertical
		staircase = Partition([floor((h-i)*v/h) for i in range(1,h)])
		return SkewPartition([staircase, self.diagram]).size()

	def rank(self):
		r = {}
		h = self.horizontal
		v = self.vertical
		g = gcd(h,v)
		p = self.diagram
		counter = 0
		if p:
			q = p.conjugate()
			horizontal = [h - q[0]] + to_exp_nozero(p)
			horizontal = horizontal[::-1]
			vertical = self.full
			for i in range(len(vertical)):
				for j in range(vertical[i]):
					try:
						r[self.labels[i][j]] = counter
					except:
						pass
					counter += h // g
				counter -= (v // g) * horizontal[i]
			return r
		else:
			return dict([(i,i) for i in range(len(self.labels[0]))])

	def labelling_permutation(self):
		return Permutation([i for j in self.labels for i in j])

	def dinv_pairs(self):
		r = self.rank()
		w = self.labelling_permutation().inverse()
		for i in w:
			for j in w:
				if i < j and w(i) < w(j) and r[i] == r[j]:
					yield (i,j)
				if i < j and w(i) > w(j) and r[i] == r[j] + 1:
					yield (i,j)

	def dinv(self):
		return len(list(self.dinv_pairs()))

	def trunc(self):
		t = []
		for i in range(len(self.full)):
			tt = []
			for j in range(self.full[i]):
				try:
					tt.append(self.labels[i][j])
				except:
					tt.append(' ')
			t.append(tt)
		return t

	# def skewtab(self):
	# 	t = []
	# 	p = self.diagram
	# 	if p:
	# 		for i in range(len(p)):
	# 			t.append([None]*p[i] + self.labels[i][::-1])
	# 		t.append(self.labels[-1][::-1])
	# 	else:
	# 		t.append(self.labels[0][::-1])
	# 	return SkewTableau(t).conjugate()

	# def pp(self):
	# 	self.skewtab().pp()

	def skewtab(self):
		t = []
		p = self.diagram
		tr = self.trunc()
		if p:
			for i in range(len(p)):
				t.append([None]*p[i] + tr[i][::-1])
			t.append(tr[-1][::-1])
		else:
			t.append(tr[0][::-1])
		return SkewTableau(t).conjugate()

	def pp(self):
		self.skewtab().pp()

def rational_pf(h,v):
	staircase = Partition([floor((h-i)*v/h) for i in range(1,h)])
	for n in range(staircase.size()+1):
		for x in Partitions(n, outer = staircase):
			if x:
				vertical = [v - x[0]] + x.conjugate().to_exp()
			else:
				vertical = [v]
			for osp in OrderedSetPartitions(v,vertical):
				t = [sorted(i) for i in osp]
				yield RationalPF(x, h, v, t)

def rpf(n,k):
	staircase = Partition([(k-i)*(n-k+1) for i in range(1,k)])
	yield RationalPF(Partition([]), k, k*(n-k+1), [list(range(1,k*(n-k+1)+1))])
	for m in range(staircase.size()+1):
		for x in Partitions(m, outer = staircase):
			if x:
				vertical = [k*(n-k+1) - x[0]] + to_exp_nozero(x.conjugate())
			else:
				vertical = [k*(n-k+1)]
			for y in Compositions(n, outer = vertical, inner = [i-n+k for i in vertical]):
				for osp in OrderedSetPartitions(n,y):
					# t = []
					# big = n+1
					# for i in range(len(y)):
					# 	diff = vertical[i] - y[i]
					# 	t.append(sorted(osp[i]) + list(range(big, big+diff)))
					# 	big += diff
					yield RationalPF(x, k, k*(n-k+1), [sorted(i) for i in osp])

