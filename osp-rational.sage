def osp_to_permutation(osp):
	m = 0
	w = []
	for s in osp[::-1]:
		small = sorted([i for i in s if i > m])
		w = small + sorted(s-set(small)) + w
		m = w[0]
	return Permutation(w)

def osp_minimaj(osp):
	return osp_to_permutation(osp).major_index()

def perm_area(w):
	a = [len([j for j in range(i,len(w)) if j+1 in w.descents()]) for i in range(len(w))]
	a.reverse()
	return a

def osp_to_markings(osp):
	m = []
	for s in osp:
		m += [0] + [1]*(len(s)-1)
	return m

def perm_schedule(w):
	n = w.size()
	r = w.runs()
	r.append([0])
	v = w.inverse()
	s = [0]*n
	for j in range(len(r)-1):
		for i in r[j]:
			s[v(i)-1] = len([k for k in r[j] if k > i]) + len([k for k in r[j+1] if k < i])
	return s

def osp_dinv(osp):
	k = len(osp)
	temp = 0
	l = [sorted(i) for i in osp]
	for i in range(k):
		for j in range(i+1,k):
			oi = l[i]
			oj = l[j]
			for h in range(len(oi)):
				try:
					if oi[h] > oj[h]:
						temp += 1
					if oi[h] < oj[h+1]:
						temp += 1
				except:
					continue
	return temp

def osp_schedule(osp):
	w = osp_to_permutation(osp)
	sch = [0] * len(w)
	runs = w.runs()
	z = dict(zip(w,osp_to_markings(osp)))
	if len(runs) == 1:
		for c in runs[0]:
			if z[c] == 0:
				sch[c-1] = len([i for i in runs[0] if i > c and z[i] == 0]) + 1
			else:
				sch[c-1] = len([i for i in runs[0] if i < c and z[i] == 0])
		return sch
	for c in runs[0]:
		if z[c] == 0:
			sch[c-1] = len([i for i in runs[0] if i > c and z[i] == 0]) + len([i for i in runs[1] if i < c and z[i] == 0])
		else:
			sch[c-1] = len([i for i in runs[0] if i < c and z[i] == 0])
	for j in range(1,len(runs)-1):
		for c in runs[j]:
			if z[c] == 0:
				sch[c-1] = len([i for i in runs[j] if i > c and z[i] == 0]) + len([i for i in runs[j+1] if i < c and z[i] == 0])
			else:
				sch[c-1] = len([i for i in runs[j] if i < c and z[i] == 0]) + len([i for i in runs[j-1] if i > c and z[i] == 0])
	for c in runs[-1]:
		if z[c] == 0:
			sch[c-1] = len([i for i in runs[-1] if i > c and z[i] == 0]) + 1
		else:
			sch[c-1] = len([i for i in runs[-1] if i < c and z[i] == 0]) + len([i for i in runs[-2] if i > c and z[i] == 0])
	return sch

def minus_one_sch(osp):
	return sum([max(i-1,0) for i in osp_schedule(osp)])

def gale_compare(l1, l2):
	return all([l1[i] >= l2[i] for i in range(min(len(l1),len(l2)))])

def mosp_to_composition(m):
	return [len(j) for j in m]

def partial_sums(m):
	return [sum(m[:i+1]) for i in range(len(m))]

class StackedPF:
	def __init__(self, stack, label):
		self.stack = Composition(stack)
		self.label = label

	def __repr__(self):
		return str(self.stack) + str(self.label)

	def area(self):
		m = self.stack
		n = m.size()
		mm = partial_sums(mosp_to_composition(self.label)) + [n]*(len(self.stack) - len(self.label))
		return SkewPartition([mm[::-1], m.partial_sums()[::-1]]).column_lengths()

	def height(self):
		p1 = 0
		p2 = 0
		p3 = len(self.label[0])
		temp = []
		for j in range(len(self.label)):
			temp += list(range(p2 - p1, p3 - p1))
			p1 += self.stack[j]
			p2 += len(self.label[j])
			try:
				p3 += len(self.label[j+1])
			except:
				continue
		return temp

	def hdinv(self):
		pairs = []
		l = [j for i in self.label for j in sorted(i)]
		h = self.height()
		n = len(l)
		for i in range(n):
			for j in range(i+1,n):
				if (h[i] == h[j] and l[i] < l[j]) or (h[i] == h[j]+1 and l[i] > l[j]):
					pairs.append((i,j))
		return len(pairs)

	def wdinv(self):
		pairs = []
		l = [j for i in self.label for j in sorted(i)]
		a = self.area()
		diag = [0] + [i for i in self.stack.partial_sums()[:-1]]
		for i in diag:
			for j in range(i+1,len(l)):
				if (a[i] == a[j] and l[i] < l[j]) or (a[i] == a[j]+1 and l[i] > l[j]):
					pairs.append((i,j))
		return len(pairs) - len(l) + len(self.stack)

	def pp(self):
		m = self.stack
		osp = self.label
		mm = partial_sums(mosp_to_composition(self.label))
		l = [sorted(osp[0])]
		for i in range(len(mm)-1):
			l.append([' ']*mm[i]+sorted(osp[i+1]))
		l.reverse()
		SkewPartition([m.partial_sums()[::-1], m.partial_sums()[::-1][1:]]).pp()
		Tableau(l).pp()

	def rise(self):
		return ParkingFunction(labelling = [j for i in self.label for j in sorted(i)], area_sequence = self.height())

	def valley(self):
		return ParkingFunction(labelling = [j for i in self.label for j in sorted(i)], area_sequence = self.area())

def stdstackpf(n,k):
	for m in Compositions(n, min_length = k, max_length = k):
		for mm in Compositions(n, max_length = k):
			if gale_compare(mm.partial_sums(), m.partial_sums()):
				for osp in OrderedSetPartitions(n, mm):
					yield StackedPF(m, osp)

def stdstackpfn(n):
	for m in Compositions(n):
		k = len(m)
		for mm in Compositions(n, max_length = k):
			if gale_compare(mm.partial_sums(), m.partial_sums()):
				for osp in OrderedSetPartitions(n, mm):
					yield StackedPF(m, osp)

def maj_d(w,d):
	assert d < len(w)
	md = 0
	for i,j in w.inversions():
		if j-i < d:
			md += 1
		elif j-i == d:
			md += i
	return md

