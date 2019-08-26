import numpy as np
import sys,os,h5py
import Pk_library as PKL
import units_library as UL


rho_crit = UL.units().rho_crit #critical density of the Universe at z=0
##################################### INPUT #############################################
prefixes  = ['TNG100', 'TNG100DM', 'TNG300', 'TNG300DM']
BoxSizes  = [75.0,     75.0,       205.0,    205.0] #Mpc/h
redshifts = [0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0]
axis      = 0 #axis along which RSD are placed (no matter in our case)
MAS       = 'CIC'
threads   = 48 #openmp threads
#########################################################################################

# do a loop over the different simulations
for prefix,BoxSize in zip(prefixes,BoxSizes):

    # do a loop over the different redshifts
    for z in redshifts[::-1]:

        # find the name of the output file
        fout = 'Pk_%s_z=%.1f.txt'%(prefix,z)
        if os.path.exists(fout):  continue

        # read the file containing the mass per voxel
        fin = '%s_z=%.1f.hdf5'%(prefix,z)
        f = h5py.File(fin, 'r')
        delta = f['density'][:]
        f.close()

        # compute the value of Omega_m
        Omega_m = np.sum(delta, dtype=np.float64)/BoxSize**3/rho_crit
        print 'Omega_m = %.4f'%Omega_m

        # compute the overdensity field
        delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0
        print np.min(delta), np.mean(delta, dtype=np.float64), np.max(delta)

        # compute Pk and save results to file
        Pk = PKL.Pk(delta, BoxSize, axis, MAS, threads)
        np.savetxt(fout, np.transpose([Pk.k3D, Pk.Pk[:,0], Pk.Nmodes3D]))
