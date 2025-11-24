def to_exp_nozero(p):
	return [i for i in p.to_exp() if i != 0]

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
		self.diagram = Partition(diagram)
		# self.dyckword = p_to_dw(diagram, h, v)
		self.labels = labels
		self.horizontal = h
		self.vertical = v
		if diagram:
			self.fullv = [v - self.diagram[0]] + to_exp_nozero(self.diagram.conjugate())
			self.fullh = [h - self.diagram.conjugate()[0]] + to_exp_nozero(self.diagram)
			self.fullh.reverse()
		else:
			self.fullh = [h]
			self.fullv = [v]

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
		dr = {i:[] for i in range(v)}
		if p:
			for i in range(len(self.fullv)):
				for j in range(self.fullv[i]):
					try:
						r[self.labels[i][j]] = counter
						dr[counter].append(self.labels[i][j])
					except:
						#pass
						dr[counter].append(None)
					counter += h // g
				counter -= (v // g) * self.fullh[i]
			return r, dr
		else:
			return dict([(i,i) for i in range(len(self.labels[0]))]), dr

	def diagonal_reading(self):
		return self.rank()[1]

	def dr_set(self):
		dr = self.diagonal_reading()
		t = []
		for i in range(len(dr)):
			if dr[i]:
				fr = frozenset(dr[i])
				t.append(frozenset({fr,len(dr[i])-len(fr)}))
		return tuple(t)

	def labelling_permutation(self):
		return Permutation([i for j in self.labels for i in j])

	def dinv_pairs(self):
		r = self.rank()[0]
		w = self.labelling_permutation().inverse()
		for i in w:
			for j in w:
				if i < j and w(i) < w(j) and r[i] == r[j]:
					yield (i,j)
				if i < j and w(i) > w(j) and r[j] == r[i] + 1:
					yield (i,j)

	def dinv(self):
		return len(list(self.dinv_pairs()))

	def trunc(self):
		t = []
		for i in range(len(self.fullv)):
			tt = []
			for j in range(self.fullv[i]):
				try:
					tt.append(self.labels[i][j])
				except:
					tt.append(' ')
			t.append(tt)
		return t

	def skewtab(self):
		t = []
		p = self.diagram
		tr = self.trunc()
		if p:
			temp = self.vertical
			for i in range(len(p)):
				if p[i] < temp:
					t.append([None]*p[i]+tr[i][::-1])
					temp = p[i]
				else:
					t.append([None]*p[i])
			t.append(tr[-1][::-1])
		else:
			t.append(tr[0][::-1])
		return SkewTableau(t).conjugate()

	def pp(self):
		self.skewtab().pp()

	def latex(self):
		s = '\\begin{tikzpicture}[scale=0.5]\n'
		s += f'\\draw[dotted] (0,0) grid ({self.horizontal:d},{self.vertical:d});\n'
		s += '\\draw[thick] (0,0)'
		m = ''
		i,j = 0,0
		inc = 0.5
		for k in range(len(self.fullv)):
			h = j
			try:
				for l in self.labels[k]:
					m += f'\\node at ({i+inc:.1f},{h+inc:.1f}) {{{l:d}}};\n'
					h += 1
			except:
				pass
			j += self.fullv[k]
			s += f'--({i:d},{j:d})'
			i += self.fullh[k]
			s += f'--({i:d},{j:d})'
		s += ';\n'
		print(s + m + '\\end{tikzpicture}')

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

def weak_comp(n, k, outer, inner):
	for y in IntegerVectors(n, k, outer = outer):
		if all([y[i] > inner[i]-1 for i in range(k)]):
			yield Composition(y)

def rpf(n,k):
	K = k*(n-k+1)
	staircase = Partition([(k-i)*(n-k+1) for i in range(1,k)])
	for m in range(staircase.size()+1):
		for x in Partitions(m, outer = staircase):
			if x:
				vertical = [K - x[0]] + to_exp_nozero(x.conjugate())
				vertical = vertical + [0]*(k-len(vertical))
			else:
				vertical = [K] + [0]*(k-1)
			inside = [max(i-n+k, 0) for i in vertical]
			for y in weak_comp(n, k, vertical, inside):
				for osp in OrderedSetPartitions(n, y):
					yield RationalPF(x, k, K, [sorted(i) for i in osp])

def test(n,k,a):
	s = set()
	d = dict()
	R.<q> = QQ['q']
	for pf in rpf(n,k):
		if pf.area() == a:
			fs = pf.dr_set()
			s.add(fs)
			# d.setdefault(fs, [])
			# d[fs].append(pf)
			d.setdefault(fs, 0)
			d[fs] += q**(pf.dinv())
	return d
