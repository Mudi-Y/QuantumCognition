import pickle
import matplotlib.pyplot as plt
import pandas as pd

tmax = 15
drift = 4
diff = 2

with open(f'stepped_anneal_{tmax}_drift{drift}_diff{diff}.pickle', 'rb') as handle:
    b = pickle.load(handle)

#extract values into DF 
vals = []
for i in range(len(b)):
    vals.append(b[i][1])

vals = pd.DataFrame(vals)
vals = vals.div(vals.sum(axis=1), axis=0)


#plot values in 4x4 grid
fig = plt.figure(figsize=(40,30))
ticks = [1,2,3,4]

for i, (name, row) in enumerate(vals.iterrows()):
    ax = plt.subplot(4,4, i+1)
    ax.set_title(f"anneal time { 1 if row.name == 0 else row.name * tmax//15} {chr(956)}s", fontsize = 40)
    ax.set_xticks(ticks, fontsize = 1)
    ax.set_ylim(bottom=0, top=1)
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.bar(ticks, row)

fig.suptitle(f"QRW states with drift = {drift}, diffusion = {diff}, step = {tmax//15} {chr(956)}s", fontsize=100)
plt.savefig(f'/workspace/QuantumCognition/SteppedAnneal_{tmax}_drift{drift}_diff{diff}.png')