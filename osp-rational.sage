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

