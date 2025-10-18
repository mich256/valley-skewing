def rises(dw):
	s = 1
	for i in range(1,len(dw)):
		s += dw[i]
		if dw[i] == 1 and dw[i-1] == 1:
			yield s

def valleys(pf):
	dw = pf.to_dyck_word()
	s = 1 + dw[1]
	w = pf.cars_permutation()
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
		self.marked_cars = {self.pf.cars_permutation()(i) for i in self.mark}
		if len(mark) == 0:
			self.type = 'normal'
		elif next(iter(mark)) in set(rises(self.dw)):
			self.type = 'rise'
		elif next(iter(mark)) in set(valleys(pf)):
			self.type = 'valley'

	def area(self):
		if type == 'rise':
			return sum([self.dw.to_area_sequence()[i] for i in range(len(pf)) 
				if i not in self.mark])
		else:
			return self.dw.area()

	def dinv_code1(self):
		a = self.dw.to_area_sequence()
		n = len(self.pf)
		w = self.pf.cars_permutation()
		for i in range(len(self.pf)):
			temp = 0
			for j in range(i):
				if j+1 not in self.mark:
					if a[j] == a[i] and w[i] > w[j]:
						temp += 1
					elif a[j] == a[i]+1 and w[i] < w[j]:
						temp += 1
			if self.type == 'valley' and i+1 in self.mark:
				yield temp-1
			else:
				yield temp

	def dinv_code2(self):
		a = self.dw.to_area_sequence()
		n = len(self.pf)
		w = self.pf.cars_permutation()
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
					yield temp
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

def riseMPF(n,k):
	for pf in ParkingFunctions(n):
		for marks in Subsets(set(rises(pf.to_dyck_word())),k):
			yield MPF(pf,marks)

def valleyMPF(n,k):
	for pf in ParkingFunctions(n):
		for marks in Subsets(set(valleys(pf)),k):
			yield MPF(pf,marks)


