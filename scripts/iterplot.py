#/usr/bin/env python

from rdkit import Chem
import os

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA

from ACSESS import distance
# initiate distance
distance.metric = 'autocorr2d'
distance.dimRed = 40
n_components=3

distance.Init()

def set_nxy(n):
    xys = { '1':(1,1), '2':(2,1), '3':(1,3), '4':(2,2), '5':(3,2), 
            '6':(3,2), '7':(4,2), '8':(4,2), '9':(3,3), '10':(4,3)
          }
    return xys[str(n)]

def pcaplot3d(coords, gen):
    ax = fig.add_subplot(nx, ny, i+1, projection='3d')
    X_reduced = PCA(n_components=n_components).fit_transform(coords)
    ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2],
    #                  c=y,
                       cmap=plt.cm.Set1, edgecolor='k', s=40)

    ax.set_title("1st 3 PCA directions of gen: {}".format(str(gen)))
    ax.set_xlabel("1st eigenvector")
    ax.w_xaxis.set_ticklabels([])
    ax.set_ylabel("2nd eigenvector")
    ax.w_yaxis.set_ticklabels([])
    ax.set_zlabel("3rd eigenvector")
    ax.w_zaxis.set_ticklabels([])

def pcaplot2d(coords, gen):
    ax = fig.add_subplot(nx, ny, i+1)
    X_reduced = PCA(n_components=n_components).fit_transform(coords)
    ax.scatter(X_reduced[:, 0], X_reduced[:, 1],
    #                  c=y,
                       cmap=plt.cm.Set1)

    ax.set_title("1st 2 PCA directions of gen: {}".format(str(gen)))
    ax.set_xlabel("1st eigenvector")
    #ax.w_xaxis.set_ticklabels([])
    ax.set_ylabel("2nd eigenvector")
    #ax.w_yaxis.set_ticklabels([])

def pcaplot2d(coords, gen):
    ax = fig.add_subplot(nx, ny, i+1)



miniter = 5
maxiter = 10
iterfiles = []

# . get a list of it*.smi files
for i in range(miniter, maxiter+1):
    filename = "it{}.smi".format(i)
    if os.path.isfile(filename):
        iterfiles.append((i, filename))

n = len(iterfiles)
if n>10:
    raise NotImplementedError('too many files')
print "there are {} files".format(n)


fig = plt.figure()
nx, ny = set_nxy(n)
for i, (gen, iterf) in enumerate(iterfiles):
    print "we will use:", iterf

    # get the smiles from the file:
    supplier = Chem.SmilesMolSupplier(iterf, sanitize=True)
    lib = [mol for mol in supplier]
    print "{} SMILES strings were read".format(len(lib))


    # set coords
    passmols, coords = distance.HandleMolCoords(lib, _noDimRed=False)
    print coords[1], type(coords[1]), coords[1].shape

    # do a PCA plot
    if n_components==3:
        pcaplot3d(coords, gen)
    else:
        assert n_components==2
        pcaplot2d(coords, gen)

plt.tight_layout()
plt.show()
