import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import figure, axes, pie, title, show, cm
import MDAnalysis as mda

chi_row=[]
ref=open("/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/core1/all_core1","r")
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
print(nstar)
print(len(nstar))


#u=mda.Universe('/Users/u0819028/mnt/lakers/project/histatin_5/gly_0.5_alasidechain/neutral_histidine/trial1/production/HISTATIN5.psf','/Users/u0819028/mnt/lakers/project/histatin_5/gly_0.5_alasidechain/neutral_histidine/trial1/production/histatin5_trial1.dcd')
u=mda.Universe('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus_SOP_sc_frags1.psf','/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/all_production.dcd')
c_all = u.select_atoms("(bynum 1:428)")

print(len(u.trajectory))
#nstar=[1,3]
k=-1
skip=1
#with open('filename.txt', 'w') as f:
import sys
sys. stdout = open("test.txt", "w")
new_list_nstar = [(i) * skip for i in nstar]
#print(new_list_nstar)
for nframe in new_list_nstar:
    coordinates=[]
    k=k+1
#    print(nframe)
    for ts in u.trajectory[nframe]:
        coordinates.append(c_all.positions)
#        print(len(coordinates[0]))
        def compute_distance(i,j):    
            dist=np.linalg.norm(coordinates[0][i-1]-coordinates[0][j-1])
            return dist
#        compute_distance(2,12)
    


    def read_psf():

        psf=open('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus_SOP_sc_frags1.psf')
        for _ in range(5):
            next(psf)
        k=0
        charge=[]
        for line in psf:        
            while k<428:
                k=k+1
    #            print(line)
                idx,a,b,c,d,q,m,e = ((psf.readline().strip()).split())
    #            print(idx)
                charge.append(tuple([idx, q]))
        return(charge)
    #    print(charge) 
    read_psf()


    def read_pair_debye():
        debye_list=[]
        deb=open('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus214.debye')
        for line in deb:
    #        print(line)
            columns= ((line.strip()).split())
            i=columns[1]
            j=columns[2]
    #        j= ((deb.readline().strip()).split())[2]
            debye_list.append([i, j])
        return (debye_list)
    #    print(debye_list)
    read_pair_debye()


    def read_pair_fene():
        fene_list=[]
        deb=open('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus214.bondparams')
        for line in deb:
    #        print(line)
            columns= ((line.strip()).split())
            i=columns[1]
            j=columns[2]
            k=columns[3]
            l=columns[4]
    #        j= ((deb.readline().strip()).split())[2]
            fene_list.append([i, j, k, l])
        return (fene_list)
    #    print(debye_list)
    read_pair_fene()

    def read_bond_fene():
        fene_bond_list=[]
        deb=open('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus214.data')
        for _ in range(879):
            next(deb)
        for line in deb:
    #        print(line)
            columns= ((line.strip()).split())
    #        print(columns)
            i=columns[0]
            j=columns[2]
            k=columns[3]
    #        j= ((deb.readline().strip()).split())[2]
            fene_bond_list.append([i, j, k])
        return (fene_bond_list)
    #    print(debye_list)
    x=read_bond_fene()
#    print(x[0])


    def read_pair_local():
        pair_local_list=[]
        deb=open('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus214.local')
        for line in deb:
    #        print(line)
            columns= ((line.strip()).split())
            i=columns[1]
            j=columns[2]
            k=columns[4]
            l=columns[5]
    #        j= ((deb.readline().strip()).split())[2]
            pair_local_list.append([i, j, k, l])
        return (pair_local_list)
    #    print(debye_list)
    read_pair_local()

    def read_pair_nonlocal():
        pair_nonlocal_list=[]
        deb=open('/Users/u0819028/Desktop/Coarse_grained_IDP/fus_peptide/trajectories/fus214/fus214.nonlocal')
        for line in deb:
    #        print(line)
            columns= ((line.strip()).split())
            i=columns[1]
            j=columns[2]
            k=columns[4]
            l=columns[5]
    #        j= ((deb.readline().strip()).split())[2]
            pair_nonlocal_list.append([i, j, k, l])
        return pair_nonlocal_list
    #    print(debye_list)
    read_pair_nonlocal()

    import math
    q=[0]* len(c_all)
    E=0
    def debye_force(i,j):
        C=(1.6*10**-19)**2/((4*math.pi*8.85*10**-12)*78)
    #    print(i,j)
        q[i]=float(read_psf()[i-1][1])
        q[j]=float(read_psf()[j-1][1])
        kappa=(1/(9.6133*10**-10))

        dist=compute_distance(i,j)*(10**-10)
    #    print(i,j,dist,q[i],q[j],C,coordinates[i-1],coordinates[j-1])
    #    E=0
    #    if q[i] != 0 and q[j] != 0:
        E =(C*q[i]*q[j]/dist)*np.exp(-kappa*dist)
    #    print(np.exp(-kappa*dist))    
        #conversion from Joule to kcal/mol
        return E *(6.022*10**23)*0.000239006

    #debye_force(2,10)




    import math
    E=0
    bond_list=read_bond_fene()
    pair_list=read_pair_fene()
    def fene_force(i):
        k=20
        dist=compute_distance(int(bond_list[i][1]),int(bond_list[i][2]))
        ro=2
        rosq=ro**2
        rref=pair_list[i][3]
    #    print(k,dist,rosq,rref)
        diffr=(float(dist) - float(rref))**2
    #    print(i,j,dist,q[i],q[j],C,coordinates[i-1],coordinates[j-1])
    #    E=0
    #    if q[i] != 0 and q[j] != 0:
        E=-0.5*k*(rosq)*math.log(1-(diffr/rosq))
    #    E =(C*q[i]*q[j]/dist)*np.exp(-kappa*dist)
    #    print(np.exp(-kappa*dist))    
        #conversion from Joule to kcal/mol
    #    return E *(6.022*10**23)*0.000239006
        return E

    #fene_force(1)



    import math
    E=0
    pair_list_nonlocal=read_pair_nonlocal()
    def non_local_force(i):
        eps=float(pair_list_nonlocal[i][2])
        dist=float(compute_distance(int(pair_list_nonlocal[i][0]),int(pair_list_nonlocal[i][1])))
        sigma=float(pair_list_nonlocal[i][3])
    #    print(eps,dist,sigma)
    #    rref=pair_list[i][3]
    #    print(k,dist,rosq,rref)
    #    diffr=(float(dist) - float(rref))**2
    #    print(i,j,dist,q[i],q[j],C,coordinates[i-1],coordinates[j-1])
    #    E=0
    #    if q[i] != 0 and q[j] != 0:
        if dist==0:
            E=0
        else:
            E=eps*((sigma/dist)**12-2*(sigma/dist)**6)
    #    E =(C*q[i]*q[j]/dist)*np.exp(-kappa*dist)
    #    print(np.exp(-kappa*dist))
        #conversion from Joule to kcal/mol
    #    return E *(6.022*10**23)*0.000239006
        return E

    #local_force(251)

    import math
    E=0
    pair_list_local=read_pair_local()
    def local_force(i):
        eps=float(pair_list_local[i][2])
        dist=float(compute_distance(int(pair_list_local[i][0]),int(pair_list_local[i][1])))
        sigma=float(pair_list_local[i][3])
    #    print(eps,dist,sigma)
    #    rref=pair_list[i][3]
    #    print(k,dist,rosq,rref)
    #    diffr=(float(dist) - float(rref))**2
    #    print(i,j,dist,q[i],q[j],C,coordinates[i-1],coordinates[j-1])
    #    E=0
    #    if q[i] != 0 and q[j] != 0:
        if dist==0:
            E=0
        else:
            E=eps*((sigma/dist)**6)
    #    E =(C*q[i]*q[j]/dist)*np.exp(-kappa*dist)
    #    print(np.exp(-kappa*dist))
        #conversion from Joule to kcal/mol
    #    return E *(6.022*10**23)*0.000239006
        return E

    #local_force(4)

    fene=0
    elocal=0
    enonlocal=0

    index_debye=(read_pair_debye())
    debtotal=0
    for i in range(len(read_pair_debye())):
        debtotal += debye_force(int(index_debye[i][0]),int(index_debye[i][1]))
#        print(int(index_debye[i][0]),int(index_debye[i][1]),coordinates)
#    print(debtotal)


    for i in range(427):
        fene +=fene_force(i)
#    print(fene)

    for i in range(2342):
        elocal +=local_force(i)
#    print(elocal)

    for i in range(89464):
        enonlocal +=non_local_force(i)
#    print(enonlocal)
#    print(elocal+enonlocal+debtotal+fene,elocal,enonlcoal,debtotal,fene)
    print(elocal+enonlocal+debtotal+fene,elocal,enonlocal,debtotal,fene)

    
#print((coordinates[0][0]),(coordinates[1][0]))
