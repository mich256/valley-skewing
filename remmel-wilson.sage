def op_to_perm(osp):
	return [j for i in osp for j in reversed(sorted(i))]

def op_to_star(osp):
	temp = []
	for i in osp:
		temp += [1]*(len(i)-1)
		temp.append(0)
	return temp

def op_major(osp):
	star = op_to_star(osp)
	w = op_to_perm(osp)
	temp = 0
	for i in range(osp.size()-1):
		if w[i] > w[i+1]:
			temp += i + 1 - sum(star[:i+1])
	return temp