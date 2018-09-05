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
## plot 1: diversity, nPool, nLib ##
####################################
plt.figure(1)

plt.subplot(311)
plt.plot(df['gen'], df['diversity'])
plt.title('diversity plot')

plt.subplot(312)
plt.plot(df['gen'], df['nPool'])
plt.title('nPool')

plt.subplot(313)
plt.plot(df['gen'], df['nLib'])
plt.title('nLib')

plt.show()

