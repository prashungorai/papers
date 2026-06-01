import os
from pylada.crystal import read,write,neighbors
from copy import deepcopy
import numpy as np
from glob import iglob
import pandas as pd
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
    for ng in neighbors(strc,howmany,pos,0.1):
        # If empty just add ng and the first one should have the shortest distance
        if len(nghs)==0:
            nghs.append(ng)
            continue

        # Make sure the distance is within the tolerance
        elif nghs[0][-1]<=ng[-1]<=(1.+tol)*nghs[0][-1]:
            nghs.append(ng)

    return nghs

# Anion list; make sure all the anions are added
anions = ['N','S','Se','O','I','C','Te']

# read structure

for name in sorted(iglob('GaN_w/band_gap_bowing/AlScN_HSE/static_HSE_222/0.222_*')):

    strc = read.poscar(name+'/CONTCAR')
    rows = []
    cation_no = 0
    for atom in strc:
        if atom.type not in anions:
            nghs = myneighbors(strc,4,atom.pos,0.25)
            cation_no+=1
            rows.append({
                'cation':               atom.type,
                'cation_no':            cation_no,
                'coordination_number':  len(nghs),
                'neighbour_species':    [ng[0].type for ng in nghs],
                'neighbour_distances':  [ng[-1] for ng in nghs],
            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(name, 'coordination_analysis.csv')
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

