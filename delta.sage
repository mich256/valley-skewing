R.<q,t> = QQ['q,t']
F = FractionField(R)
Sym = SymmetricFunctions(F)
e = Sym.e()
m = Sym.m()
Ht = Sym.macdonald().Ht()
QSym = QuasiSymmetricFunctions(F)

def bmu(m):
	m = Partition(m)
	return sum(q^i * t^j for i,j in m.cells())

def delta(f, g):
	return sum(f(bmu(i)-1)*j*Ht(i) for i,j in Ht(g))