def rises(dw):
	s = 1
	for i in range(1,len(dw)):
		s += dw[i]
		if dw[i] == 1 and dw[i-1] == 1:
			yield s

def valleys(pf):
	dw = pf.to_dyck_word()
	s = 1 + dw[1]
	w = pf.to_labelling_permutation()
	for i in range(2,len(dw)):
		s += dw[i]
		if dw[i] == 1 and dw[i-1] == 0:
			if w(s) > w(s-1):
				yield s
			if dw[i-2] == 0:
				yield s

class MPF:
	def __init__(self, pf, mark):
		self.pf = pf
		self.dw = pf.to_dyck_word()
		self.mark = mark
		self.marked_cars = {self.pf.to_labelling_permutation()(i) for i in self.mark}
		if len(mark) == 0:
			self.type = 'normal'
		elif next(iter(mark)) in set(rises(self.dw)):
			self.type = 'rise'
		elif next(iter(mark)) in set(valleys(pf)):
			self.type = 'valley'

	def area(self):
		if self.type == 'rise':
			return sum([self.pf.to_area_sequence()[i] for i in range(len(self.pf)) if i+1 not in self.mark])
		else:
			return self.dw.area()

	def dinv_code1(self):
		a = self.pf.to_area_sequence()
		n = len(self.pf)
		w = self.pf.to_labelling_permutation()
		for i in range(len(self.pf)):
			temp = 0
			for j in range(i):
				if self.type == 'valley':
					if j+1 not in self.mark:
						if a[j] == a[i] and w[i] > w[j]:
							temp += 1
						elif a[j] == a[i]+1 and w[i] < w[j]:
							temp += 1
				else:
					if a[j] == a[i] and w[i] > w[j]:
						temp += 1
					elif a[j] == a[i]+1 and w[i] < w[j]:
						temp += 1
			if self.type == 'valley' and i+1 in self.mark:
				yield temp-1
			else:
				yield temp

	def dinv_code2(self):
		a = self.pf.to_area_sequence()
		n = len(self.pf)
		w = self.pf.to_labelling_permutation()
		if self.type == 'valley':
			for i in range(len(self.pf)):
				temp = 0
				if i+1 not in self.mark:
					for j in range(i+1,n):
						if j+1 not in self.mark:
							if a[j] == a[i] and w[i] < w[j]:
								temp += 1
							elif a[j] == a[i]-1 and w[i] > w[j]:
								temp += 1
					yield temp
				elif i+1 in self.mark:
					for j in range(i):
						if j+1 not in self.mark:
							if a[j] == a[i] and w[i] > w[j]:
								temp += 1
							elif a[j] == a[i]+1 and w[i] < w[j]:
								temp += 1
					yield temp-1
			return
		else:
			for i in range(len(self.pf)):
				temp = 0
				for j in range(i+1,n):
					if a[j] == a[i] and w[i] < w[j]:
						temp += 1
					elif a[j] == a[i]-1 and w[i] > w[j]:
						temp += 1
				yield temp

	def dinv_pairs(self):
		if self.type == 'valley':
			a = self.pf.to_area_sequence()
			n = len(self.pf)
			w = self.pf.to_labelling_permutation()
			temp = []
			for i in range(len(self.pf)):
				if i+1 in self.mark:
					for j in range(i):
						if j+1 not in self.mark:
							if a[j] == a[i] and w[i] > w[j]:
								temp.append((j+1,i+1))
							elif a[j] == a[i] + 1 and w[i] < w[j]:
								temp.append((j+1,i+1))
				elif i+1 not in self.mark:
					for j in range(i):
						if j+1 not in self.mark:
							if a[j] == a[i] and w[i] > w[j]:
								temp.append((j+1,i+1))
							elif a[j] == a[i] + 1 and w[i] < w[j]:
								temp.append((j+1,i+1))
			return set(temp)
		else:
			return set([(i+1,j+1) for (i,j) in self.pf.dinversion_pairs()])

	def dinv_pairs_label(self):
		w = self.pf.to_labelling_permutation()
		return set([tuple(sorted((w(i),w(j)))) for (i,j) in self.dinv_pairs()])

	def dinv_code3(self):
		w = self.pf.to_labelling_permutation()
		v = w.inverse()
		sett = self.dinv_pairs()
		for i in range(1,len(self.pf)+1):
			temp = 0
			for j in range(1,i):
				if tuple(sorted((v(i),v(j)))) in sett:
					temp += 1
			if self.type == 'valley' and v(i) in self.mark:
				yield temp - 1
			else:
				yield temp

	def dinv(self):
		if self.type == 'valley':
			return sum(self.dinv_code2())
		else:
			return self.pf.dinv()

	def pp(self):
		self.pf.pretty_print()
		print(self.marked_cars)

	def __repr__(self):
		return self.pf, self.marked_cars

	def latex(self, aa = False, dd = False):
		w = self.pf.to_labelling_permutation()
		dw = self.dw
		n = dw.semilength()
		res = '\\begin{tikzpicture}[scale=0.5]\n'
		res += f'\\draw[dotted] (0,0) grid ({n:d},{n:d});\n'
		res += f'\\draw[thick,dotted] (0,0) -- ({n:d},{n:d});\n'
		res += '\\draw[thick] (0,0)'
		label = f'\\draw node at (0.5,0.5) {w[0]:d};\n'
		mark = ''
		stats = ''
		coord = [0,0]
		inc = 0.5
		for i in range(len(dw)):
			if dw[i] == 1:
				coord[1] += 1
				label += f'\\draw node at ({coord[0]+inc:.1f},{coord[1]-inc:.1f}) {w(coord[1]):d};\n'
				if coord[1] in self.mark:
					mark += f'\\draw node at ({coord[0]-inc:.1f},{coord[1]+inc:.1f}) {{*}};\n'
			if dw[i] == 0:
				coord[0] += 1
			res += '--({},{})'.format(coord)
		if aa:
			if self.type == 'rise':
				stats += f'\\draw node at (1,-.5) {{area-: {self.area():d}}};\n'
			else:
				stats += f'\\draw node at (1,-.5) {{area: {self.area():d}}};\n'
		if dd:
			if self.type == 'valley':
				stats += f'\\draw node at (1,-1.5) {{dinv-: {self.dinv():d}}};\n'
			else:
				stats += f'\\draw node at (1,-1.5) {{dinv: {self.dinv():d}}};\n'
		print(res + ';\n' + label + mark + stats + '\\end{tikzpicture}')

def riseMPF(n,k):
	for pf in ParkingFunctions(n):
		for marks in Subsets(set(rises(pf.to_dyck_word())),k):
			yield MPF(pf,marks)

def valleyMPF(n,k):
	for pf in ParkingFunctions(n):
		for marks in Subsets(set(valleys(pf)),k):
			yield MPF(pf,marks)


