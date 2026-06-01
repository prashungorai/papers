from pylada.crystal import read,neighbors,coordination_shells,gruber,supercell
from glob import iglob,glob
from copy import deepcopy
import spglib
#from vladan.format_spglib import *
import numpy as np
import math
import os
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

    # Loop ove pylada neighbors (arbitrary small tolerance)
    for ng in neighbors(strc,howmany,pos,0.2):
        # If empty just add ng and the first one should have the shortest distance
        if len(nghs)==0:
            nghs.append(ng)
            continue

        # Make sure the distance is within the tolerance
        elif nghs[0][-1]<=ng[-1]<=(1.+tol)*nghs[0][-1]:
            nghs.append(ng)

    return nghs

# Anion list; make sure all the anions are added
anions = ['N','O','S','Se']

if not os.path.exists('data_angle_diff'):
    os.makedirs('data_angle_diff')

angle1 = 108.20449  # angle associated with axial bonds
angle2 = 110.70770  # angle only associated with basal bonds 
for path in sorted(glob('MgHfAlN_NEB/39_images/MgHfAlN_39_U/end_images/POSCAR_*/POSCAR_p/POSCAR')):
    strc = read.poscar(path)
    name = path.replace('MgHfAlN_NEB/39_images/MgHfAlN_39_U/end_images/','').replace('/POSCAR_p/POSCAR','') 
    diffangle = []
    for i in range(len(strc)):
        if strc[i].type !='Al': continue # only look at the distortion to N-Al-N angle
        nghs = myneighbors(strc=strc, howmany=12, pos=strc[i], tol=0.15)  # same tolerance as bond length
        if len(nghs)!=4 : continue
        for j in range(len(nghs)):
            for k in range(len(nghs)):
                if j>=k: continue
                v1 = nghs[j][1]
                v2 = nghs[k][1]
                value = np.round(np.dot(v1,v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)),8)
                angle = math.degrees(np.arccos(value)) 
                v1norm = np.linalg.norm(v1)
                v2norm = np.linalg.norm(v2)
                v1proj = np.absolute(np.dot(v1,np.array([0,0,1])))
                v2proj = np.absolute(np.dot(v2,np.array([0,0,1])))
                if v1proj < 0.9*v1norm and v2proj < 0.9*np.linalg.norm(v2):
                    diffangle.append(angle-angle2)
                else: 
                    diffangle.append(angle-angle1)
    
    print(name,len(diffangle)) # check if the number is expected
    np.savetxt('data_angle_diff/%s.dat'%(name),diffangle,header='%s'%name)
