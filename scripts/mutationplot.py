#/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filenames = ['stats.dat']

for filename in filenames:
    with open(filename) as f:
        columnnames = f.readline().split()
        if columnnames[-1]=='\\':
            columnnames = columnnames[:-1] + f.readline().split()
            #columnnames += f.readline().split()
        data = f.readlines()
    dataf = [ item.split() for item in data ]
    for item in dataf:
        item.remove('|')
    df = pd.DataFrame(dataf, columns=columnnames).astype(float)

print df

####################################
## plot 2: diversity, nPool, nLib ##
####################################

fig = plt.figure()
ax = fig.add_subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

y_pos = df['gen']
patch_handles  = []
colors = ['red', 'goldenrod', 'green', 'turquoise', 'blue', 'magenta', 'grey']
colors = ['red', 'goldenrod', 'green', 'slategrey', 'blue', 'magenta', 'grey']
C = ['nAdd', 'nAddArRing', 'nBreak', 'nFlip', 'nNewRing', 'nRemove', 'nNoMutation']
for i in range(len(C))[:-1]:
    # add Fail
    C.insert(2*i+1, C[2*i]+'Fail')
    # add a dark color
    colors.insert(2*i, 'dark' + colors[2*i])
print "C:", C
print "colors:", colors
MuData = df[C]
for c, cfail in zip(C[::2], C[1::2]):
    #MuData[c] = MuData[c] - MuData[cfail]
    MuData.loc[:, c] = MuData[c] - MuData[cfail]
MuData = MuData.values.T

left = np.zeros(len(df))

for i, (d, c) in enumerate(zip(MuData, C)):
    patch_handles.append(
            ax.barh(y_pos, d, left=left, label=c, color=[colors[i]]))
    left += d

ax.set_yticks(y_pos)
ax.set_yticklabels(y_pos)
ax.set_xlabel('nMutations')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.legend(bbox_to_anchor=(1, 0.5), loc='center left')
plt.show()
