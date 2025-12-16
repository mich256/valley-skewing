from collections import defaultdict
import pyperclip

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

def latex_fs(fs):
	# list_of_fs = list(fr_pp(fs))
	s = ''
	for i in fs:
		counter = 0
		for j in i:
			if j:
				s += str(j)
			else:
				counter += 1
		s += '\\emptyset'*counter+'|'
	return s[:-1]

class RationalPF:
	def __init__(self, diagram, h, v, labels):
		self.diagram = Partition(diagram)
		# self.dyckword = p_to_dw(diagram, h, v)
		self.labels = labels
		self.horizontal = h
		self.vertical = v
		self.n = sum(len(i) for i in self.labels)
		self.slope = v/h
		if diagram:
			self.fullv = [v - self.diagram[0]] + to_exp_nozero(self.diagram.conjugate())
			self.fullh = [h - self.diagram.conjugate()[0]] + to_exp_nozero(self.diagram)
			self.fullh.reverse()
		else:
			self.fullh = [h]
			self.fullv = [v]

	def spf(self):
		load('osp-rational.sage')
		k = self.horizontal
		n = self.vertical//k -1 + k
		p = self.diagram
		if p:
			stack = [self.vertical-p[0]]
			index = 0
			label = [set(self.labels[index])]
			for i in range(k-1):
				if i < len(p):
					t = p[i] - p[i+1]
					stack.append(t)
					if t == 0:
						label.append(set())
					else:
						index += 1
						label.append(set(self.labels[index]))
				elif i == len(p)-1:
					stack.append(p[i])
					label.append(set(self.labels[-1]))
				else:
					stack.append(0)
					label.append(set())
			for i in range(k):
				stack[i] = n-k+1-(stack[i]-len(label[i]))
		else:
			stack = [self.vertical-n] + [n-k+1]*(k-1)
		return StackedPF(stack, label)

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
		ps = 0
		for i in range(len(self.fullv)):
			for j in range(self.fullv[i]):
				try:
					r[self.labels[ps][j]] = counter
					dr[counter].append(self.labels[ps][j])
				except:
					pass
					# dr[counter].append(None)
				counter += h // g
			counter -= (v // g) * self.fullh[i]
			ps += self.fullh[i]
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
				# t.append(frozenset({fr,len(dr[i])-len(fr)}))
				t.append(fr)
		return tuple(reversed(t))

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
		m += f'\\node at ({float(self.horizontal/2):.1f}, -0.5) {{${latex_fs(self.dr_set())}$}};\n'
		return s + m + '\\end{tikzpicture}'

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
	load('osp-rational.sage')
	for spf in stdstackpf(n,k):
		yield spf.rpf()

def lowest(n,k,a):
	d1 = defaultdict(list)
	d2 = defaultdict(int)
	R.<q> = QQ['q']
	for pf in rpf(n,k):
		if pf.area() == a:
			fs = pf.lowest()
			d1[fs].append(pf)
			d2[fs] += q**(pf.tdinv())
	return dict(d1),dict(d2)

def rpf_table(n,k,a):
	R.<q> = QQ['q']
	d = defaultdict(R)
	dl = defaultdict(set)
	s = '\\begin{table}[H]\n\\['
	s += '\\begin{array}{|c|c|c|}\\hline\n'
	for pf in rpf(n,k):
		if pf.area() == a:
			fs = pf.dr_set()
			d[fs] += q**(pf.tdinv())
			dl[pf.lowest()].add(fs)
	for low in dl.keys():
		switch = True
		for fs in dl[low]:
			if switch:
				s += ''.join(str(i) for i in sorted(low))
				switch = False
			s += ' &' + latex_fs(fs) 
			s += ' &' + latex(factor(d[fs])) + ' \\\\ \n'
		s += '\\hline\n'
	s += '\\end{array}\n\\]\n'
	s += f'\\caption{{$n={n},k={k},\\area={a}$.}}\n'
	s += '\\end{table}\n'
	return pyperclip.copy(s)

def test_function(n):
	load('mpf.sage')
	for k in range(1,n):
		for a in range(binomial(n,2)-binomial(k+1,2)):
			assert lowest_unm(n,k,a) == lowest(n,n-k,a)