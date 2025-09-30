def osp_to_permutation(osp):
	m = 0
	w = []
	for s in osp[::-1]:
		small = sorted([i for i in s if i > m])
		w = small + sorted(s-set(small)) + w
		m = w[0]
	return Permutation(w)

def osp_to_markings(osp):
	m = []
	for s in osp:
		m += [0] + [1]*(len(s)-1)
	return m

def osp_schedule(osp):
	w = osp_to_permutation(osp)
	sch = [0] * len(w)
	runs = w.runs()
	z = dict(zip(w,osp_to_markings(osp)))
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

def gale_compare(l1, l2):
	assert len(l1) == len(l2)
	return all([l1[i] >= l2[i] for i in range(len(l1))])

class StackedPF:
	def __init__(self, stack, label):
		self.stack = stack
		self.label = label

	def __repr__(self):
		return str(self.stack) + str(self.label)

	def area(self):
		m = self.stack
		mm = self.label.to_composition()
		return SkewPartition([mm.partial_sums()[::-1], m.partial_sums()[::-1]]).column_lengths()

	def height(self):
		m = self.stack
		p1 = self.stack.partial_sums()
		p2 = self.label.to_composition().partial_sums()
		temp = list(range(p1[0])) + [p1[0]] * (p2[0]-p1[0])
		for i in range(len(m)-1):
			temp += list(range(p2[i] - p1[i], m[i+1])) + [m[i+1]]*(p2[i+1]-p1[i+1])
		return temp

	def hdinv(self):
		pairs = []
		l = [j for i in self.label for j in sorted(i)]
		h = self.height()
		n = len(l)
		for i in range(n):
			for j in range(i,n):
				if (h[i] == h[j] and l[i] < l[j]) or (h[i] == h[j]+1 and l[i] > l[j]):
					pairs.append((i,j))
		return len(pairs)

	def wdinv(self):
		pairs = []
		l = [j for i in self.label for j in sorted(i)]
		a = self.area()
		diag = [0] + [i for i in self.stack.partial_sums()[:-1]]
		for i in diag:
			for j in range(i,len(l)):
				if (a[i] == a[j] and l[i] < l[j]) or (a[i] == a[j]+1 and l[i] > l[j]):
					pairs.append((i,j))
		return len(pairs) - len(l) + len(self.stack)

	def pp(self) -> None:
		m = self.stack
		osp = self.label
		mm = osp.to_composition().partial_sums()
		l = [sorted(osp[0])]
		for i in range(len(mm)-1):
			l.append([' ']*mm[i]+sorted(osp[i+1]))
		l.reverse()
		SkewPartition([m.partial_sums()[::-1], m.partial_sums()[::-1][1:]]).conjugate().pp()
		print('\n')
		Tableau(l).conjugate().pp()
		print('\n\n')

def stackpf(n,k):
	for m in Compositions(n, min_length = k, max_length = k):
		for mm in Compositions(n, min_length = k, max_length = k):
			if gale_compare(mm.partial_sums(), m.partial_sums()):
				for osp in OrderedSetPartitions(n, mm):
					yield StackedPF(m, osp)


