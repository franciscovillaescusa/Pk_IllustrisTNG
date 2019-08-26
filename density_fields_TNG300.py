import numpy as np
import MAS_library as MASL
import sys,os,h5py,time
import units_library as UL

U = UL.units();  rho_crit = U.rho_crit #h^2 Msun/Mpc^3
################################ INPUT ########################################
dims = 2048
MAS  = 'CIC'

root = '/simons/scratch/sgenel/Illustris_IllustrisTNG_public_data_release/L205n2500TNG/output/'
prefix_out = 'TNG300'
snapnums = np.array([21, 25, 33, 40, 50, 59, 72, 84, 99])

do_CDM   = True
do_gas   = True
do_stars = True
do_BH    = True
##############################################################################


# do a loop over the different redshifts
for snapnum in snapnums:

    # find the snapshot name
    snapshot = '%s/snapdir_%03d/snap_%03d.0.hdf5'%(root,snapnum,snapnum)
    if not(os.path.exists(snapshot)):  continue 

    # read header
    f        = h5py.File(snapshot, 'r')
    redshift = f['Header'].attrs[u'Redshift']
    BoxSize  = f['Header'].attrs[u'BoxSize']/1e3  #Mpc/h
    filenum  = f['Header'].attrs[u'NumFilesPerSnapshot']
    Omega_m  = f['Header'].attrs[u'Omega0']
    Omega_L  = f['Header'].attrs[u'OmegaLambda']
    h        = f['Header'].attrs[u'HubbleParam']
    Masses   = f['Header'].attrs[u'MassTable']*1e10  #Msun/h
    Nall     = f['Header'].attrs[u'NumPart_Total']
    f.close()

    # determine whether it is an N-body or hydrodynamic simulation
    if Nall[0]==0:  Nbody = True
    else:           Nbody = False

    # find output file name
    print 'Working with %s'%snapshot
    f_out = '%s_z=%.1f.hdf5'%(prefix_out,redshift)
    if os.path.exists(f_out):  continue
        
    # define the array hosting the density field
    delta = np.zeros((dims,dims,dims), dtype=np.float32)

    # do a loop over all subfiles in a given snapshot
    M_total, start = 0.0, time.time()
    for j in xrange(filenum):

        snapshot = '%s/snapdir_%03d/snap_%03d.%d.hdf5'%(root,snapnum,snapnum,j)
        f = h5py.File(snapshot, 'r')

        ### CDM ###
        if do_CDM:
            pos  = (f['PartType1/Coordinates'][:]/1e3).astype(np.float32)        
            mass = np.ones(pos.shape[0], dtype=np.float32)*Masses[1] #Msun/h
            MASL.MA(pos, delta, BoxSize, MAS, mass)  #CDM
            M_total += np.sum(mass, dtype=np.float64)

        ### GAS ###
        if do_gas and not(Nbody):
            pos  = (f['PartType0/Coordinates'][:]/1e3).astype(np.float32)
            mass = f['PartType0/Masses'][:]*1e10  #Msun/h
            MASL.MA(pos, delta, BoxSize, MAS, mass)  #gas
            M_total += np.sum(mass, dtype=np.float64)

        ### Stars ###
        if do_stars and not(Nbody):
            pos  = (f['PartType4/Coordinates'][:]/1e3).astype(np.float32)        
            mass = f['PartType4/Masses'][:]*1e10  #Msun/h
            MASL.MA(pos, delta, BoxSize, MAS, mass)  #stars
            M_total += np.sum(mass, dtype=np.float64)

        ### Black-holes ###
        if do_BH and not(Nbody):
            pos  = (f['PartType5/Coordinates'][:]/1e3).astype(np.float32)        
            mass = f['PartType5/Masses'][:]*1e10  #Msun/h
            MASL.MA(pos, delta, BoxSize, MAS, mass)  #black-holes
            M_total += np.sum(mass, dtype=np.float64)

        f.close()

        print '%03d -----> Omega_m = %.4f  : %6.0f s'\
            %(j, M_total/(BoxSize**3*rho_crit), time.time()-start)

    f = h5py.File(f_out,'w')
    f.create_dataset('density',  data=delta)
    f.close()

