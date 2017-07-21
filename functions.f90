module functions

    Implicit none
     
contains

subroutine read_data_CF3(path_to_datafile)
    use arrays
    Implicit none
    Integer*4 :: arrays_dimension,p
    Integer :: stat
    character(len=*) :: path_to_datafile

    open(11,file=path_to_datafile)

    read(11,*)

    arrays_dimension = 0

    Do 

        read(11,*,iostat=stat)

        If (stat .ne. 0) then

            exit

        Else

            arrays_dimension = arrays_dimension + 1 

        End If

    End Do

    close(11)

    allocate (redshift(1:arrays_dimension),luminosity_distance(1:arrays_dimension),galactic_longitude(1:arrays_dimension),&
    galactic_latitude(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    read(11,*)

    Do p=1,arrays_dimension

       read(11,*) redshift(p),luminosity_distance(p),galactic_longitude(p),galactic_latitude(p)

    End Do

    close(11)

end subroutine read_data_CF3

subroutine write_python_script_angular_distribution_galaxies()
    
    use fiducial
    Implicit none
!!$    Real*8:: 
!!$    Real*8:: 
!!$    Real*8,dimension(nbins):: 
!!$    logical :: 
!!$    character*16 :: 
!!$    character*16 :: 
!!$    Character(len=10) :: 
        
!!$    fmt = '(es16.10)' 
!!$    write(string_omega_b,fmt) param_omega_b
!!$    write(string_omega_cdm,fmt) param_omega_cdm
!!$    write(string_n_s,fmt) param_n_s
!!$    write(string_A_s,fmt) param_A_s
!!$    write(string_H0,fmt) param_H0
!!$    write(string_m_ncdm,fmt) param_m_ncdm
!!$    write(string_MG_beta2,fmt) param_MG_beta2
!!$    write(string_nc_bias_b0,fmt) param_nc_bias_b0

    open(10, file='./scripts/angular_distribution_galaxies_CF3.py')

    write(10,*) 'import numpy as np'
    write(10,*) 'import math'
    write(10,*) 'import matplotlib as mpl'
    write(10,*) "mpl.use('Agg')"
    write(10,*) 'import matplotlib.pyplot as py'
    write(10,*) 'from matplotlib.colors import LogNorm'
    write(10,*) 'import healpy as hp'

!!$marker = ['bo','rs','go','bv','r^','g<','b>','rD','g+','bx','r*','bo']
!!$
!!$exefile = open('execution_information.txt','w')
!!$
!!$# WE LOAD DATA FROM COSMICFLOWS-3 
!!$
!!$table1 = np.dtype([('The Catalogue of Principal Galaxies (PGC) Number',np.str_,7),('Luminosity distance (weighted average if more than one source)',np.float64),('Number of distance sources',np.int32),('Luminosity distance modulus (weighted average if more than one source)',np.float64),('1 sigma uncertainty in distance modulus (NOTE! NOT fractional error in distance as given in CF2)',np.float64),('Distance method: Cepheid Period-Luminosity relation',np.str_,1),('Distance method: Tip of the Red Giant Branch from collaboration HST program (Jacobs et al. 2009)',np.str_,1),('Distance method: Tip of the Red Giant Branch (literature)',np.str_,1),('Distance method: miscellaneous (RR Lyrae, Horizontal Branch, Eclipsing Binary)',np.str_,1),('Distance method: Surface Brightness Fluctuation',np.str_,1),('Distance method: type Ia supernovae',np.str_,1),('Distance method: optical (I band) spiral luminosity-rotation correlation (Tully-Fisher)',np.str_,1),('Distance method: Spitzer [3.6] band spiral luminosity-rotation correlation (Tully-Fisher)',np.str_,1),('Distance method: E/S0 Fundamental Plane - sources ENEAR, EFAR, SMAC',np.str_,1),('Luminosity distance modulus carried over from Cosmicflows-2.1',np.float64),('1 sigma uncertainty in distance modulus carried over from Cosmicflows-2.1',np.float64),('Supernova ID',np.str_,7),('Number of separate analyses of this supernova',np.int32),('Luminosity distance modulus of supernova averaged over all contributions',np.float64),('Luminosity distance modulus determined from luminosity-rotation correlation using Spitzer [3.6] photometry',np.float64),('1 sigma uncertainty in distance modulus determined from luminosity-rotation correlation using Spitzer [3.6] photometry',np.float64),('Luminosity distance modulus from 6dFGS (Springob et al. 2014) with CF3 zero point',np.float64),('1 sigma uncertainty in distance modulus from 6dFGS (Springob et al. 2014) with CF3 zero point',np.float64),('6dFGS morphological type code',np.str_,3),('Right Ascension (J2000)',np.str_,8),('Declination (J2000)',np.str_,9),('Galactic longitude',np.float64),('Galactic latitude',np.float64),('Supergalactic longitude',np.float64),('Supergalactic latitude',np.float64),('Morphological type; RC3 numeric code',np.float64),('Reddening at B band from Schlafly & Finkbeiner (2011,ApJ,737,103) dust maps',np.float64),('Total B magnitude from LEDA',np.float64),('2MASS Ks magnitude, extinction corrected from Huchra et al. 2012 or else Lavaux-Hudson 2011',np.float64),('Heliocentric velocity',np.int32),('Velocity in Galactic standard of rest (circular velocity at Sun of 239 km/s; total velocity 251 km/s toward l=90,b=0)',np.int32),('Velocity in Local Sheet standard of rest (Tully et al. 2008)',np.int32),('Velocity in CMB standard of rest (Fixsen et al. 1996)',np.int32),('Velocity in CMB standard of rest adjusted in accordance with a cosmological model with Omega_matter=0.27 and Omega_Lambda=0.73',np.int32),('Common name',np.str_,11),('2MASS nest group ID',np.int32),('Number of galaxies with distance measures in group',np.int32),('Luminosity distance modulus to group; weighted average of all contributions',np.float64),('1 sigma error in group distance modulus',np.float64),('Luminosity distance to group',np.float64),('Abell cluster ID; ASxxx are from southern supplement list',np.str_,5),('Alternate name for group or cluster',np.str_,9),('Number of galaxies with positions and velocities in group',np.int32),('PGC ID of dominant galaxy in group',np.int32),('Galactic latitude of group',np.float64),('Galactic longitude of group',np.float64),('Supergalactic longitude of group',np.float64),('Supergalactic latitude of group',np.float64),('Log summed K luminosity of group, adjusted by correction factor for lost light; assumed distance is Vmgp/75',np.float64),('Luminosity selection function correction factor',np.str_,5),('Projected velocity dispersion anticipated by corrected luminosity',np.str_,4),('Projected second turnaround radius anticipated by corrected intrinsic luminosity; assumed distance is Vmgp/75',np.float64),('Group heliocentric velocity',np.int32),('Group velocity in Galactic standard of rest (circular velocity at Sun of 239 km/s; total velocity 251 km/s toward l=90,b=0)',np.int32),('Group velocity in Local Sheet standard of rest (Tully et al. 2008)',np.int32),('Group velocity in CMB standard of rest (Fixsen et al. 1996)',np.int32),('Group velocity in CMB standard of rest adjusted in accordance with a cosmological model with Omega_matter=0.27 and Omega_Lambda=0.73',np.int32),('Group r.m.s. velocity dispersion',np.str_,4),('Group mass (x10^12) from virial theorem with bi-weight dispersion and radius parameters; assumed distance is Vmgp/75',np.float64),('Group mass (x10^12) based on corrected intrinsic luminosity and M/L prescription; assumed distance is Vmgp/75',np.float64),('Crook et al. (2007) low density group ID',np.int32),('Crook et al. (2007) high density group ID',np.int32),('Group ID from 2MASS++ catalog of Lavaux & Hudson (2011)',np.int32),('Makarov & Karachentsev (2011, MNRAS 412, 2498) group ID',np.int32),('Internal group ID',np.int32)])
!!$
!!$pgc,Dist,Nd,DM,eDM,C,T,L,M,S,N,H,I,F,DM2,eD2,SNIa,Ns,DMsn,DMsp,eDsp,DM6d,eD6d,Mt,RAJ,DeJ,Glon,Glat,SGL,SGB,Ty,Asf,Btot,Ks,Vhel,Vgsr,Vls,Vcmb,Vmod,Name,Nest,Ndgp,DMgp,eDgp,Dgp,Abell,GroupName,NV,PGC1,Glongp,Glatgp,SGLgp,SGBgp,lgLgp,cf,sigp,R2t,Vhgp,Vggp,Vlsgp,Vcgp,Vmgp,Vrms,bwMass12,L_Mass12,LDC,HDC,twoMplusplus,MKgp,Icnt = np.loadtxt('CF-3.txt',unpack=True,dtype=table1,skiprows=5,delimiter=',',usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69]) 
!!$
!!$Nside = 1 
!!$
!!$lmax = 3*Nside - 1
!!$
!!$map = np.arange(hp.nside2npix(Nside))
!!$
!!$map_nc = np.arange(hp.nside2npix(Nside),dtype=float)
!!$
!!$shot_noise = np.arange(hp.nside2npix(Nside),dtype=float)
!!$
!!$map[:] = 0
!!$
!!$map_nc[:] = 0.
!!$
!!$shot_noise[:] = 0.
!!$
!!$z_min = 0.01 
!!$
!!$z_max = 0.05 
!!$
!!$z_step = 0.005
!!$
!!$speed_of_light = 299792.458 # km/s
!!$
!!$exefile.write("THE NUMBER OF DISTANCE MEASUREMENTS IN THE CF3 CATALOGUE IS: "+str(len(Dist))+"\n")
!!$
!!$data = np.zeros((len(Dist),7))
!!$
!!$for index in range(len(Dist)):
!!$    data[index,0] = Vcmb[index]/speed_of_light # REDSHIFT 
!!$    data[index,1] = Dist[index]                # LUMINOSITY DISTANCE
!!$    data[index,2] = Glon[index]                # GALACTIC LONGITUDE
!!$    data[index,3] = Glat[index]                # GALACTIC LATITUDE
!!$    data[index,4] = RAJ[index]                 # RIGHT ASCENSION
!!$    data[index,5] = DeJ[index]                 # DECLINATION
!!$    data[index,6] = Vmod[index]/speed_of_light # CORRECTED REDSHIFT
!!$
!!$exefile.write("LOWEST UNCORRECTED REDSHIFT IN THE CATALOGUE IS: "+str(np.min(data[:,0]))+"\n")
!!$
!!$exefile.write("HIGHEST UNCORRECTED REDSHIFT IN THE CATALOGUE IS: "+str(np.max(data[:,0]))+"\n")
!!$
!!$z_max_temp = z_min + z_step
!!$
!!$z_min_temp = z_min
!!$
!!$counter = 0
!!$
!!$while (z_max_temp <= z_max):
!!$
!!$    hp.graticule(dpar=10.,dmer=10.,coord='G')
!!$
!!$    for index in range(len(Dist)):
!!$
!!$        if ( (data[index,0] > z_min_temp) & (data[index,0] <= z_max_temp) ):
!!$
!!$            hp.projplot(data[index,2],data[index,3],marker[counter],lonlat=True,coord='G')
!!$    
!!$    exefile.write("NUMBER OF GALAXIES IN REDSHIFT BIN Z_MIN "+str(z_min_temp)+" Z_MAX "+str(z_max_temp)+" IS: "+str(np.sum(map))+"\n")
!!$    
!!$    z_max_temp = z_max_temp + z_step
!!$
!!$    counter = counter + 1
!!$
!!$
!!$py.savefig('galaxy_distribution_redshift_bin_z_min_'+str(z_min_temp)+'_z_max'+str(z_max_temp)+'_CF3.pdf')

    write(10,*) 'exit()'

    
!!$    ! Background parameters and anisotropic stress
!!$                                                                                                        
!!$    write(10,'(a6, es16.10)') 'A_s = ', param_A_s  
!!$
!!$    write(10,'(a6, es16.10)') 'n_s = ', param_n_s
!!$ 
!!$    write(10,'(a5, es16.10)') 'H0 = ', param_H0
!!$
!!$    write(10,'(a10, es16.10)') 'omega_b = ', param_omega_b
!!$
!!$    write(10,'(a12, es16.10)') 'omega_cdm = ', param_omega_cdm
!!$
!!$    write(10,'(a11, es16.10)') 'tau_reio = ', param_tau_reio
!!$
!!$    write(10,'(a11, es16.10)') 'MG_beta2 = ', param_MG_beta2
!!$
!!$    write(10,'(a13, es16.10)') 'nc_bias_b0 = ', param_nc_bias_b0
!!$
!!$    ! Parameters for massive neutrinos                                                                                            
!!$
!!$    write(10,'(a7, f5.3)') 'N_ur = ', real(param_N_ur)
!!$
!!$    write(10,'(a9, f5.3)') 'N_ncdm = ', real(param_N_ncdm)
!!$
!!$    write(10,'(a11, f5.3)') 'deg_ncdm = ', real(param_deg_ncdm)
!!$
!!$    write(10,'(a9, es16.10)') 'm_ncdm = ', param_m_ncdm
!!$
!!$    ! Number counts in the output                                                                                            
!!$
!!$    write(10,'(a12)') 'output = nCl'
!!$    
!!$    !write(10,'(a20)') 'non linear = halofit'
!!$ 
!!$    write(10,'(a32)') 'dNdz_selection = analytic_euclid'
!!$
!!$    write(10,'(a32)') 'dNdz_evolution = analytic_euclid'
!!$
!!$    write(10,'(a20)') 'selection = gaussian'
!!$
!!$    write(10,'(a17, 4(f10.8, a1),f10.8)') 'selection_mean = ', z_bin_centers(1),',', z_bin_centers(2),',', z_bin_centers(3),',',&
!!$    z_bin_centers(4),',',z_bin_centers(5)!,',',z_bin_centers(6),',',z_bin_centers(7),',',z_bin_centers(8),',',&
!!$!    z_bin_centers(9),',',z_bin_centers(10)
!!$
!!$    write(10,'(a18, 4(f10.8, a1),f10.8)') 'selection_width = ', z_bin_widths(1),',',z_bin_widths(2),',',z_bin_widths(3),',',&
!!$    z_bin_widths(4),',',z_bin_widths(5)!,',',z_bin_widths(6),',',z_bin_widths(7),',',z_bin_widths(8),',',&
!!$ !   z_bin_widths(9),',',z_bin_widths(10)
!!$
!!$    write(10,'(a17, 4(f10.8, a1),f10.8)') 'selection_bias = ', z_bin_bias(1),',',z_bin_bias(2),',',z_bin_bias(3),',',&
!!$    z_bin_bias(4),',',z_bin_bias(5)!,',',z_bin_bias(6),',',z_bin_bias(7),',',z_bin_bias(8),',',&
!!$  !  z_bin_bias(9),',',z_bin_bias(10)
!!$
!!$    write(10,'(a31, 4(f10.8, a1),f10.8)') 'selection_magnification_bias = ', s_z_mag_bias(1),',',s_z_mag_bias(2),',',&
!!$         s_z_mag_bias(3),',',s_z_mag_bias(4),',',s_z_mag_bias(5)!,',',z_bin_bias(6),',',z_bin_bias(7),',',z_bin_bias(8),',',&
!!$!    z_bin_bias(9),',',z_bin_bias(10)
!!$
!!$    write(10,'(a15,i2)') 'non_diagonal = ',nbins-1
!!$
!!$    write(10,'(a13)') 'headers = yes'
!!$
!!$    write(10,'(a17)') 'bessel file = yes'
!!$
!!$    write(10,'(a12,i4)') 'l_max_lss = ', lmax_class
!!$
!!$    write(10,'(a8,i1)') 'l_min = ', lmin
!!$
!!$!    write(10,'(a27,f2.0)') 'selection_magnitude_bias = ', 0.
!!$
!!$    write(10,'(a14)') 'format = class'
!!$
!!$    write(10,'(a17)') 'gauge = newtonian'
!!$
!!$    ! PRECISION PARAMETERS
!!$
!!$    write(10,'(a40, f6.0)') 'l_switch_limber_for_cl_density_over_z = ', real(l_switch_limber_for_cl_density_over_z)
!!$
!!$    write(10,'(a28, f5.2)') 'selection_sampling_bessel = ', real(bessel)
!!$
!!$    write(10,'(a12, f5.1)') 'q_linstep = ', real(q)
!!$
!!$    write(10,'(a24, f5.2)') 'k_max_tau0_over_l_max = ', real(kmaxtau0)

    close(10)

end subroutine write_python_script_angular_distribution_galaxies

subroutine write_mpi_file(name_python_file)
    Implicit none
    character(len=*) :: name_python_file

    open(12,file='./scripts/'//trim(name_python_file)//'.mpi')

    write(12,'(a11)') '#!/bin/bash'
    write(12,'(a12)') '#SBATCH -N 1'
    write(12,'(a23)') '#SBATCH -J ANG-DIST-CF3'
    write(12,'(a12)') '#SBATCH -n 1'
    write(12,'(a21)') '#SBATCH -t 2-00:00:00'
    write(12,'(a21)') '#SBATCH -o job.%J.out'
    write(12,*) 
    write(12,*)'prun python ./scripts/'//trim(name_python_file)//'.py'

    close(12)

end subroutine write_mpi_file


end module functions
