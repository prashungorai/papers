from pylada.crystal import read,write,neighbors
import numpy as np
from glob import glob
import os

def myneighbors(strc=None, howmany=12, pos=None, tol=0.2):
    """
    A function that replaces the original pylada
    neighbors function that seems to be broken.
    It still uses the original one, just makes sure
    the output is correct. Also, it returns only the 
    first coordination shell.

    The tolerance is in fraction of the distance to 
    the closest atom and is 20% by default (test first)
    """

    nghs=[]

    # Loop over pylada neighbors (arbitrary small tolerance)
    for ng in neighbors(strc,howmany,pos,0.2):
        # If empty just add ng and the first one should have the shortest distance 
        if len(nghs)==0:
            nghs.append(ng)
            continue

        # Make sure the distance is within the tolerance
        elif nghs[0][-1]<=ng[-1]<=(1.+tol)*nghs[0][-1]:
            nghs.append(ng)

    return nghs

'''
This script assumes that it is the cations that are moving during ferroelectric switching since ionic radius for anionis are much larger
'''
# Anion list; make sure all the anions are added 
anions = ['N','O','S','Se']

# loop over polar structures of AlN-alloys
for path in sorted(glob('relaxed-polar/0.*')):
    print(path)
    # read structure
    strc = read.poscar(path + 'CONTCAR')
    scale = (strc.volume / 844.237)**(1/3)  # assume 72-atom cells; the number is for 332 supercell of GaN 
   #scale = 1.0
    bond1 = 1.966 * scale # basal planar ones 
    bond2 = 1.973 * scale # polar direction 
    name = path.replace('relaxed-polar/','data_bond_diff_scaling/')
    os.system('mkdir -p %s'%name)

    # find all the cation-to-anion vectors
    diff = []
    for atom in strc:
        if atom.type == 'Ga':  #only look at Al-N bonds
            nghs = myneighbors(strc,12,atom.pos,0.2)  # based on our knowledge for wurtzite-like structure
            if len(nghs)!= 4 : continue
            for ngh in nghs:
                if ngh[0].type in anions:
                    if np.dot(ngh[1],np.array([0,0,1])) > 0.9*ngh[2]: 
                        print('polar direction bond')
                        diff.append(float(ngh[2])-bond2)
                    else: 
                        print('basal direction bond')
                        diff.append(float(ngh[2])-bond1)
    print(len(diff))
    np.savetxt('%s/bond_decompose.dat'%(name),diff,header='%s'%name)
