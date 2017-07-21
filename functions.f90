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

    open(10, file='./scripts/angular_distribution_galaxies_CF3.py')

    write(10,'(a)') 'import numpy as np'
    write(10,'(a)') 'import math'
    write(10,'(a)') 'import matplotlib as mpl'
    write(10,'(a)') "mpl.use('Agg')"
    write(10,'(a)') 'import matplotlib.pyplot as py'
    write(10,'(a)') 'from matplotlib.colors import LogNorm'
    write(10,'(a)') 'import healpy as hp'
    write(10,'(a)') "marker = ['bo','rs','go','bv','r^','g<','b>','rD','g+','bx','r*','bo']"

    write(10,'(a)') "table1 = np.dtype([('The Catalogue of Principal Galaxies (PGC) Number',np.str_,7),('Luminosity distance &
         &(weighted average if more than one source)',np.float64),('Number of distance sources',np.int32),('Luminosity &
         &distance modulus (weighted average if more than one source)',np.float64),('1 sigma uncertainty in distance modulus &
         &(NOTE! NOT fractional error in distance as given in CF2)',np.float64),('Distance method: Cepheid Period-Luminosity &
         &relation',np.str_,1),('Distance method: Tip of the Red Giant Branch from collaboration HST program (Jacobs et al. &
         &2009)',np.str_,1),('Distance method: Tip of the Red Giant Branch (literature)',np.str_,1),('Distance method: &
         &miscellaneous (RR Lyrae, Horizontal Branch, Eclipsing Binary)',np.str_,1),('Distance method: Surface Brightness &
         &Fluctuation',np.str_,1),('Distance method: type Ia supernovae',np.str_,1),('Distance method: optical (I band) spiral &
         &luminosity-rotation correlation (Tully-Fisher)',np.str_,1),('Distance method: Spitzer [3.6] band spiral &
         &luminosity-rotation correlation (Tully-Fisher)',np.str_,1),('Distance method: E/S0 Fundamental Plane - sources ENEAR, &
         &EFAR, SMAC',np.str_,1),('Luminosity distance modulus carried over from Cosmicflows-2.1',np.float64),('1 sigma &
         &uncertainty in distance modulus carried over from Cosmicflows-2.1',np.float64),('Supernova ID',np.str_,7),('Number of &
         &separate analyses of this supernova',np.int32),('Luminosity distance modulus of supernova averaged over all &
         &contributions',np.float64),('Luminosity distance modulus determined from luminosity-rotation correlation using Spitzer &
         &[3.6] photometry',np.float64),('1 sigma uncertainty in distance modulus determined from luminosity-rotation correlation &
         &using Spitzer [3.6] photometry',np.float64),('Luminosity distance modulus from 6dFGS (Springob et al. 2014) with CF3 &
         &zero point',np.float64),('1 sigma uncertainty in distance modulus from 6dFGS (Springob et al. 2014) with CF3 zero &
         &point',np.float64),('6dFGS morphological type code',np.str_,3),('Right Ascension (J2000)',np.str_,8),('Declination &
         &(J2000)',&
         &np.str_,9),('Galactic longitude',np.float64),('Galactic latitude',np.float64),('Supergalactic longitude',np.float64),&
         &('Supergalactic latitude',np.float64),('Morphological type; RC3 numeric code',np.float64),('Reddening at B band from &
         &Schlafly & Finkbeiner (2011,ApJ,737,103) dust maps',np.float64),('Total B magnitude from LEDA',np.float64),('2MASS Ks &
         &magnitude, extinction corrected from Huchra et al. 2012 or else Lavaux-Hudson 2011',np.float64),('Heliocentric &
         &velocity',&
         &np.int32),('Velocity in Galactic standard of rest (circular velocity at Sun of 239 km/s; total velocity 251 km/s &
         &toward &
         &l=90,b=0)',np.int32),('Velocity in Local Sheet standard of rest (Tully et al. 2008)',np.int32),('Velocity in CMB &
         &standard &
         &of rest (Fixsen et al. 1996)',np.int32),('Velocity in CMB standard of rest adjusted in accordance with a cosmological &
         &model &
         &with Omega_matter=0.27 and Omega_Lambda=0.73',np.int32),('Common name',np.str_,11),('2MASS nest group ID',np.int32),&
         &('Number of galaxies with distance measures in group',np.int32),('Luminosity distance modulus to group; weighted &
         &average &
         &of all contributions',np.float64),('1 sigma error in group distance modulus',np.float64),('Luminosity distance to &
         &group',&
         &np.float64),('Abell cluster ID; ASxxx are from southern supplement list',np.str_,5),('Alternate name for group or &
         &cluster',&
         &np.str_,9),('Number of galaxies with positions and velocities in group',np.int32),('PGC ID of dominant galaxy in &
         &group',&
         &np.int32),('Galactic latitude of group',np.float64),('Galactic longitude of group',np.float64),('Supergalactic &
         &longitude of &
         &group',np.float64),('Supergalactic latitude of group',np.float64),('Log summed K luminosity of group, adjusted &
         &by correction &
         &factor for lost light; assumed distance is Vmgp/75',np.float64),('Luminosity selection function correction &
         &factor',np.str_,5),&
         &('Projected velocity dispersion anticipated by corrected luminosity',np.str_,4),('Projected second turnaround radius &
         &anticipated &
         &by corrected intrinsic luminosity; assumed distance is Vmgp/75',np.float64),('Group heliocentric velocity',np.int32),&
         &('Group velocity in Galactic standard of rest (circular velocity at Sun of 239 km/s; total velocity 251 km/s &
         &toward l=90,b=0)'&
         &,np.int32),('Group velocity in Local Sheet standard of rest (Tully et al. 2008)',np.int32),('Group velocity &
         &in CMB standard of &
         &rest (Fixsen et al. 1996)',np.int32),('Group velocity in CMB standard of rest adjusted in accordance with &
         &a cosmological model &
         &with Omega_matter=0.27 and Omega_Lambda=0.73',np.int32),('Group r.m.s. velocity dispersion',np.str_,4),&
         &('Group mass (x10^12) &
         &from virial theorem with bi-weight dispersion and radius parameters; assumed distance is Vmgp/75',np.float64),&
         &('Group mass &
         &(x10^12) based on corrected intrinsic luminosity and M/L prescription; assumed distance is Vmgp/75',np.float64),&
         &('Crook et al. &
         &(2007) low density group ID',np.int32),('Crook et al. (2007) high density group ID',np.int32),('Group ID &
         &from 2MASS++ catalog &
         &of Lavaux & Hudson (2011)',np.int32),('Makarov & Karachentsev (2011, MNRAS 412, 2498) group ID',np.int32),&
         &('Internal group ID',np.int32)])"

    write(10,'(a)') "pgc,Dist,Nd,DM,eDM,C,T,L,M,S,N,H,I,F,DM2,eD2,SNIa,Ns,DMsn,DMsp,eDsp,DM6d,eD6d,Mt,RAJ,DeJ,Glon,&
         &Glat,SGL,SGB,Ty,Asf,&
         &Btot,Ks,Vhel,Vgsr,Vls,Vcmb,Vmod,Name,Nest,Ndgp,DMgp,eDgp,Dgp,Abell,GroupName,NV,PGC1,Glongp,Glatgp,SGLgp,&
         &SGBgp,lgLgp,cf,sigp,&
         &R2t,Vhgp,Vggp,Vlsgp,Vcgp,Vmgp,Vrms,bwMass12,L_Mass12,LDC,HDC,twoMplusplus,MKgp,Icnt = np.loadtxt('../data/CF-3.txt',&
         &unpack=True,dtype=&
         &table1,skiprows=5,delimiter=',',usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,&
         &25,26,27,28,29,30,&
         &31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,&
         &67,68,69])"

    write(10,'(a,i4)') 'Nside = ', nsmax

    write(10,'(a)') 'lmax = 3*Nside - 1'

    write(10,'(a,es16.10)') 'z_min = ', redshift_min 

    write(10,'(a,es16.10)') 'z_max = ', redshift_max 

    write(10,'(a,es16.10)') 'z_step = ', redshift_step

    write(10,'(a)') 'speed_of_light = 299792.458 # km/s'

    write(10,'(a)') 'data = np.zeros((len(Dist),7))'

    write(10,'(a)') 'for index in range(len(Dist)):'
    write(10,'(a)') '    data[index,0] = Vcmb[index]/speed_of_light # REDSHIFT' 
    write(10,'(a)') '    data[index,1] = Dist[index]                # LUMINOSITY DISTANCE'
    write(10,'(a)') '    data[index,2] = Glon[index]                # GALACTIC LONGITUDE'
    write(10,'(a)') '    data[index,3] = Glat[index]                # GALACTIC LATITUDE'
    write(10,'(a)') '    data[index,4] = RAJ[index]                 # RIGHT ASCENSION'
    write(10,'(a)') '    data[index,5] = DeJ[index]                 # DECLINATION'
    write(10,'(a)') '    data[index,6] = Vmod[index]/speed_of_light # CORRECTED REDSHIFT'

    write(10,'(a)') 'z_max_temp = z_min + z_step'

    write(10,'(a)') 'z_min_temp = z_min'

    write(10,'(a)') 'counter = 0'

    write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

    write(10,'(a)') 'while (z_max_temp <= z_max):'

    write(10,'(a)') '    for index in range(len(Dist)):'

    write(10,'(a)') '        if ( (data[index,0] > z_min_temp) & (data[index,0] <= z_max_temp) ):'

    write(10,'(a)') "            hp.projplot(data[index,2],data[index,3],marker[counter],lonlat=True,coord='G')"
    
    write(10,'(a)') '    z_max_temp = z_max_temp + z_step'

    write(10,'(a)') '    counter = counter + 1'

    write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_whole_CF3_dataset.pdf')"

    write(10,'(a)') 'py.close()'

    write(10,'(a)') 'z_max_temp = z_min + z_step'

    write(10,'(a)') 'z_min_temp = z_min'

    write(10,'(a)') 'counter = 0'

    write(10,'(a)') 'while (z_max_temp <= z_max):'

    write(10,'(a)') "    hp.graticule(dpar=10.,dmer=10.,coord='G')"

    write(10,'(a)') '    for index in range(len(Dist)):'

    write(10,'(a)') '        if ( (data[index,0] > z_min_temp) & (data[index,0] <= z_max_temp) ):'

    write(10,'(a)') "            hp.projplot(data[index,2],data[index,3],marker[counter],lonlat=True,coord='G')"
    
    write(10,'(a)') '    z_max_temp = z_max_temp + z_step'

    write(10,'(a)') '    counter = counter + 1'

    write(10,'(a)') "    py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
         &str(z_min_temp)+'_z_max'+str(z_max_temp)+'_CF3.pdf')"

    write(10,'(a)') '    py.close()'

    write(10,'(a)') 'z_max_temp = z_min + z_step'

    write(10,'(a)') 'z_min_temp = z_min'

    write(10,'(a)') 'counter = 0'

    write(10,'(a)') 'while (z_max_temp <= z_max):'

    write(10,'(a)') "    hp.graticule(dpar=10.,dmer=10.,coord='G')"

    write(10,'(a)') '    for index in range(len(Dist)):'

    write(10,'(a)') '        if ( (data[index,0] > z_min_temp) & (data[index,0] <= z_max_temp) ):'

    write(10,'(a)') "            hp.projplot(data[index,2],data[index,3],marker[counter],lonlat=True,coord='G')"

    write(10,'(a)') '    z_min_temp = z_max_temp'
    
    write(10,'(a)') '    z_max_temp = z_max_temp + z_step'

    write(10,'(a)') '    counter = counter + 1'

    write(10,'(a)') "    py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
         &str(z_min_temp)+'_z_max'+str(z_max_temp)+'_CF3.pdf')"

    write(10,'(a)') '    py.close()'

    write(10,'(a)') 'exit()'

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
    write(12,'(a)') 'prun python ./scripts/'//trim(name_python_file)//'.py'

    close(12)

end subroutine write_mpi_file


end module functions
