import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis import contacts, distances

import sys
psf = sys.argv[1]
dcd = sys.argv[2]
output=sys.argv[3]
#psf="../fus214/fus_SOP_sc_frags1.psf"
#dcd="../fus214/all_production.dcd"
u=mda.Universe(psf,dcd)
core_region = u.select_atoms("resname CA and bynum 0:428")
print(core_region)



coordinates=[]
skip=1000
for ts in u.trajectory[0::skip]:
    coordinates.append(core_region.positions)
#print(coordinates)
print(len(coordinates))
print(type(coordinates))
print (np.ndim(coordinates))
print (np.shape(coordinates))


def compute_distances(frame,i,j):
    dist=np.linalg.norm(coordinates[frame][i]-coordinates[frame][j])
    return dist
print (coordinates[0][1],coordinates[0][2])
compute_distances(0,1,2)
#print(coordinates[0],coordinates[1])

natoms=len(core_region)
z=[0]*natoms

x=[[[0]*natoms]*natoms]*(1+int(len(u.trajectory)/skip))
count=[0]*natoms
print(len(x))
ts_id=0
z1=[0]*natoms
for ts in u.trajectory[0::skip]:
    ts_id +=1
    for i in range(natoms):
        for j in range(i,natoms):
            ts_int=int(ts.frame/skip)
            
            x[ts_int][i][j]=compute_distances(ts_int,i,j)
            z[abs(i-j)] = z[abs(i-j)]+x[ts_int][i][j]
            count[abs(i-j)] =count[abs(i-j)]+1
            
        if (count[abs(i-j)] !=0):
            z1[abs(i-j)]=z[abs(i-j)]/count[i-j]
#            print(count[abs(i-j)])
#        z[abs(i-j)]=z[abs(i-j)]/count[abs(i-j)]
#            print(x[ts_int][i][j],ts_int,i,j,count[i-j],z[abs(i-j)]/count[i-j],z1)

for i in range(natoms):
    for j in range(natoms):
        z1[abs(i-j)]=z[abs(i-j)]/count[i-j]

z1=np.array(z1)
np.savetxt('polymer_scaling'+output+'.dat',z1)        
import matplotlib.pyplot as plt
print(z1)
rijvsij=[]
#for i in range(natoms):
#    for j in range(natoms):
#        rijvsij[]=(z[abs(i-j)],abs(i-j))
plt.plot((z1),'o')
plt.savefig('out_polymer_scaling.pdf')
plt.show()
