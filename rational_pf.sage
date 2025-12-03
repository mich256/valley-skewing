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
		self.slope = v/h
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

	def pathdinv_boxes(self):
		La = self.diagram
		for (i,j) in La.cells():
			a = La.arm_length(i,j)
			l = La.leg_length(i,j)
			if l/(a+1) <= self.slope and a/(l+1) < 1/self.slope:
				yield (i,j)

	def pathdinv(self):
		return len(list(self.pathdinv_boxes()))

	def rank(self):
		r = {}
		h = self.horizontal
		v = self.vertical
		g = gcd(h,v)
		p = self.diagram
		counter = 0
		dr = {i:[] for i in range(v)}
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

	def diagonal_reading(self):
		return self.rank()[1]

	def lowest(self):
		return frozenset(self.diagonal_reading()[0])

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
		dinvs = []
		attacks = []
		r = self.rank()[0]
		w = self.labelling_permutation().inverse()
		for i in w:
			for j in w:
				if i < j and r[i] == r[j]:
					attacks.append((i,j))
					if w(i) < w(j):
						dinvs.append((i,j))
				if i < j and r[j] == r[i] + 1:
					attacks.append((i,j))
					if w(i) > w(j):
						dinvs.append((i,j))
		return dinvs, attacks

	def tdinv(self):
		return len(self.dinv_pairs()[0])

	def maxtdinv(self):
		return len(self.dinv_pairs()[1])

	def dinv(self):
		return self.pathdinv() + self.tdinv() - self.maxtdinv()

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

def fr_pp(tuple_of_frozensets):
	for i in tuple_of_frozensets:
		for j in i:
			if type(j) == int:
				if j > 0:
					tt = [None]*(j)
				else:
					tt = []
			else:
				t = list(j)
		yield tt + t

def latex_fs(list_of_fs):
	s = ''
	for i in list_of_fs:
		counter = 0
		for j in i:
			if j:
				s += str(j)
			else:
				counter += 1
		s += '\\emptyset'*counter+'|'
	return s[:-1]

def test(n,k,a):
	d = dict()
	R.<q> = QQ['q']
	print('\\begin{array}{|c|c|}\\hline')
	for pf in rpf(n,k):
		if pf.area() == a:
			fs = pf.dr_set()
			d.setdefault(fs, [])
			d[fs].append(pf)
			# d.setdefault(fs, 0)
			# d[fs] += q**(pf.tdinv())
	d2 = []
	for fs in d.keys():
		t = list(fr_pp(fs))
		cc = sum([q**(pf.tdinv()) for pf in d[fs]])
		d2.append((t, cc))
		print(latex_fs(t) + ' &' + latex(factor(cc)) + ' \\\\ \\hline')
	print('\\end{array}')
	return d, d2

def lowest(n,k,a):
	d = dict()
	R.<q> = QQ['q']
	for pf in rpf(n,k):
		if pf.area() == a:
			fs = pf.lowest()
			d.setdefault(fs, 0)
			d[fs] += q**(pf.tdinv())
	return d

load('osp-rational.sage')
def test_function(n):
	for k in range(1,n):
		for a in range(binomial(n,2)-binomial(k+1,2)):
			assert lowest_unm(n,k,a) == lowest(n,n-k,a)
