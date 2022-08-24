import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import figure, axes, pie, title, show, cm
import MDAnalysis as mda
from MDAnalysis.analysis import contacts, distances
from scipy import cluster, linalg, spatial

import MDAnalysis as mda

psf="/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus163/fus163_SOP_sc_frags1.psf"
dcd="/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus163/all_fus163.dcd"
#psf="/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus_SOP_sc_frags1.psf"
#dcd="/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/all_production.dcd"

u=mda.Universe(psf,dcd)

chi_row=[]
ref=open("/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus108/core1/all_core1","r")
#list=ref.readline().split()
for line in ref:
    line = line.strip()
    columns = line.split()        
    chi_row.append(float(columns[0]))
    
print (len(chi_row)) 



nstar=[]
#for ts in u.trajectory[0::skip]:
for x in range(0,2000010):
#    ts.frame = int(ts.frame/skip)
    if chi_row[x]>0.215:
        nstar.append(x)



skip=1
#print(nstar)
new_list_cluster = [(i) * skip for i in nstar]

c_all = u.select_atoms("(bynum 1:326)")
#protein = u.select_atoms("c_all")
with mda.Writer("nstar_fus163.pdb", c_all.n_atoms) as W:
    for ts in u.trajectory[new_list_cluster]:
        W.write(c_all)
