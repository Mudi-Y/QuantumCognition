from ising_model import IsingModel

a = IsingModel(q=2, t=1, sigmaSq=2, delta=4)
lin, quad = a.isingModel

print(lin)
print(quad)

res = a.annealQPU()
print(res)
