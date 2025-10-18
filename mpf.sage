
def rises(dw):
	s = 1
	for i in range(1,dw.semilength()):
		s += dw[i]
		if dw[i] == 1 and dw[i-1] == 1:
			yield s

def valleys(pf):
	s = 1
	dw = pf.to_dyck_word()
	w = pf.cars_permutation()
	for i in range(1,dw.semilength()):
		s += dw[i]
		if dw[i] == 1 and dw[i-1] == 0:
			if w[s] > w[s-1]:
				yield s
			try:
				if dw[i-2] == 0:
					yield s

class MPF:
	def __init__(self, pf, mark):
		self.pf = pf
		self.dw = 
		self.mark = mark
		if len(mark) == 0:
			self.type = 'unmarked'