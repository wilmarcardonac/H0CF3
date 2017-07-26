module functions

    Implicit none
     
contains

subroutine read_data_CF3(path_to_datafile)

  use arrays
  use fiducial

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

  number_galaxies_in_CF3 = arrays_dimension

  allocate (redshift(1:arrays_dimension),luminosity_distance(1:arrays_dimension),galactic_longitude(1:arrays_dimension),&
       galactic_latitude(1:arrays_dimension),stat=status1)

  open(11,file=path_to_datafile)

  read(11,*)

  Do p=1,arrays_dimension

     read(11,*) redshift(p),luminosity_distance(p),galactic_longitude(p),galactic_latitude(p)

  End Do

  close(11)

end subroutine read_data_CF3

subroutine write_python_script_angular_distribution_galaxies(option)
    
    use fiducial
    Implicit none

    Integer*4 :: option

    If (option .eq. 1) then

       open(10, file='./scripts/ang_dist_gal_whole_CF3_dataset.py')

    Else If (option .eq. 2) then

       open(10, file='./scripts/ang_dist_gal_redshift_bins_CF3_dataset.py')
       
    Else If (option .eq. 3) then

       open(10, file='./scripts/ang_dist_gal_redshift_bins_cumulative_CF3_dataset.py')

    Else

       write(UNIT_EXE_FILE,*) 'SUBROUTINE TO WRITE PYTHON SCRIPT ONLY ACCEPTS THREE OPTIONS: 1, 2, OR 3.'&
             &'DIFFERENT OPTIONS WAS GIVEN'
       
       stop

    End If

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

    If (option .eq. 1) then

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[0],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+z_step) & (data[index,0] <= z_min+2.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[1],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+2.*z_step) & (data[index,0] <= z_min+3.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[2],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+3.*z_step) & (data[index,0] <= z_min+4.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[3],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+4.*z_step) & (data[index,0] <= z_min+5.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[4],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+5.*z_step) & (data[index,0] <= z_min+6.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[5],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+6.*z_step) & (data[index,0] <= z_min+7.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[6],lonlat=True,coord='G')"

       write(10,'(a)') '    elif ( (data[index,0] > z_min+7.*z_step) & (data[index,0] <= z_min+8.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[7],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_whole_CF3_dataset.pdf')"

       write(10,'(a)') 'py.close()'

    Else If (option .eq. 2) then

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[0],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+z_step) & (data[index,0] <= z_min+2.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[1],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+z_step)+'_z_max'+str(z_min+2.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+2.*z_step) & (data[index,0] <= z_min+3.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[2],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+2.*z_step)+'_z_max'+str(z_min+3.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+3.*z_step) & (data[index,0] <= z_min+4.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[3],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+3.*z_step)+'_z_max'+str(z_min+4.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+4.*z_step) & (data[index,0] <= z_min+5.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[4],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+4.*z_step)+'_z_max'+str(z_min+5.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+5.*z_step) & (data[index,0] <= z_min+6.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[5],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+5.*z_step)+'_z_max'+str(z_min+6.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+6.*z_step) & (data[index,0] <= z_min+7.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[6],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+6.*z_step)+'_z_max'+str(z_min+7.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min+7.*z_step) & (data[index,0] <= z_min+8.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[7],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min+7.*z_step)+'_z_max'+str(z_min+8.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

    Else If (option .eq. 3) then

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[0],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+2.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[1],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+2.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+3.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[2],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+3.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+4.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[3],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+4.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+5.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[4],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+5.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+6.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[5],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+6.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+7.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[6],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+7.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

       write(10,'(a)') "hp.graticule(dpar=10.,dmer=10.,coord='G')"

       write(10,'(a)') 'for index in range(len(Dist)):'

       write(10,'(a)') '    if ( (data[index,0] > z_min) & (data[index,0] <= z_min+8.*z_step) ):'

       write(10,'(a)') "        hp.projplot(data[index,2],data[index,3],marker[7],lonlat=True,coord='G')"

       write(10,'(a)') "py.savefig('../figures/angular_distribution_galaxies_redshift_bin_z_min_'+&
            &str(z_min)+'_z_max'+str(z_min+8.*z_step)+'_CF3.pdf')"

       write(10,'(a)') 'py.close()'

    End If

    write(10,'(a)') 'exit()'

    close(10)

end subroutine write_python_script_angular_distribution_galaxies

subroutine write_mpi_file(name_python_file,option)
    
  Implicit none

  character(len=*) :: name_python_file

  Integer*4 :: option

  open(12,file='./scripts/'//trim(name_python_file)//'.mpi')

  write(12,'(a11)') '#!/bin/bash'
  write(12,'(a12)') '#SBATCH -N 1'
  If (option .eq. 1) then
  write(12,'(a23)') '#SBATCH -J 1-ANG-DIST-CF3'
  Else If (option .eq. 2) then
  write(12,'(a23)') '#SBATCH -J 2-ANG-DIST-CF3'
  Else If (option .eq. 3) then
  write(12,'(a23)') '#SBATCH -J 3-ANG-DIST-CF3'
  End If
  write(12,'(a12)') '#SBATCH -n 1'
  write(12,'(a21)') '#SBATCH -t 2-12:00:00'
  write(12,'(a21)') '#SBATCH -o job.%J.out'
  write(12,*) 
  write(12,'(a)') 'prun python '//trim(name_python_file)//'.py'

  close(12)

end subroutine write_mpi_file

subroutine compute_number_counts_map(zmin,zmax,jackknife,data_index_excluded,dip_amplitude,dip_latitude,dip_longitude)

  use healpix_types
  use udgrade_nr, only: udgrade_ring, udgrade_nest
  use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring,remove_dipole, ang2pix_ring, vec2ang
  use fitstools, only: getsize_fits, input_map, output_map
  use head_fits
  use fiducial
  use arrays

  Implicit none

  Character(len=80),dimension(1:60) :: header
  Character(len=4) :: Ns
  Character(len=5) :: xmin,xmax

  Real*8 :: zmin,zmax,mean_nc,dip_amplitude,dip_latitude,dip_longitude
  Real(kind=DP),dimension(0:DEGREE_REMOVE_DIPOLE*DEGREE_REMOVE_DIPOLE-1) :: multipoles ! SAVES MONOPOLE AND DIPOLE OF CMB MAP 
  Real(kind=DP),dimension(1:2) :: zbounds ! BOUNDS TO COMPUTE DIPOLE AND MONOPOLE
  Real(kind=DP) :: theta,phi ! COLATITUDE AND LONGITUDE

  Integer*8 :: data_index, index_i,data_index_excluded
  Integer(kind=I8B) :: ipring ! NUMBERS PIXELS IN NUMBER COUNTS MAP

  Logical :: jackknife
  Logical :: exist

  zbounds(:) = 0.d0

  If (.not.jackknife) then

     If (data_index_excluded .lt. 0) then

        continue 

     Else

        write(UNIT_EXE_FILE,*) 'WHEN NOT DOING JACKKNIFE ANALYSIS SET "data_index_excluded" TO A NEGATIVE INTEGER'&
             &'SUBROUTINE "compute_number_counts_map" '

        stop

     End If

  End If

  write(xmin,'(F5.2)') zmin
  write(xmax,'(F5.2)') zmax
  write(Ns,'(I4.4)') nsmax

  allocate (map(0:npixC-1,1:1), map_nc(0:npixC-1,1:1), shot_noise(0:npixC-1,1:1),& 
         mask(0:npixC-1,1:1), stat = status1)

  call write_minimal_header(header, 'MAP', nside = nsmax, ordering = ORDERING_NC_MAPS, coordsys = SYS_COORD) ! HEADER OF NUMBER COUNTS MAPS

  map(0:npixC-1,1) = 0.d0 ! INITIALIZATION OF NUMBER COUNTS MAP
  map_nc(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF NORMALISED NUMBER COUNTS MAP
  shot_noise(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF SHOT NOISE MASK
  mask(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF NUMBER COUNTS MASK

  Do data_index=1,number_galaxies_in_CF3

     If ( (redshift(data_index) .gt. zmin) .and. (redshift(data_index) .le. zmax) ) then

        If (jackknife) then

           If (data_index .eq. data_index_excluded) then

              continue 

           Else

              If ( galactic_latitude(data_index) .lt. 0.d0 ) then

                 call ang2pix_ring(nsmax,abs(galactic_latitude(data_index))*DEG2RAD+HALFPI,&
                      galactic_longitude(data_index)*DEG2RAD,ipring)

                 map(ipring,1) = map(ipring,1) + 1

              Else

                 call ang2pix_ring(nsmax,-abs(galactic_latitude(data_index))*DEG2RAD+HALFPI,&
                      galactic_longitude(data_index)*DEG2RAD,ipring)

                 map(ipring,1) = map(ipring,1) + 1

              End If

           End If

        Else

           If ( galactic_latitude(data_index) .lt. 0.d0 ) then

              call ang2pix_ring(nsmax,abs(galactic_latitude(data_index))*DEG2RAD+HALFPI,&
                   galactic_longitude(data_index)*DEG2RAD,ipring)

              map(ipring,1) = map(ipring,1) + 1

           Else

              call ang2pix_ring(nsmax,-abs(galactic_latitude(data_index))*DEG2RAD+HALFPI,&
                   galactic_longitude(data_index)*DEG2RAD,ipring)

              map(ipring,1) = map(ipring,1) + 1

           End If

        End If

     Else

        continue

     End If

  End Do

  mean_nc = sum(map)/npixC

  Do index_i=0,npixC-1

     If ( map(index_i,1) .ne. 0.d0 ) then

        mask(index_i,1) = 1.d0

        shot_noise(index_i,1) = 1.d0/map(index_i,1)

        map_nc(index_i,1) = ( map(index_i,1) - mean_nc )/mean_nc

     End If

  End Do

  If (jackknife) then

     continue

  Else

     inquire(file = './output/map_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//'_Nside_'//trim(Ns)//'.fits',exist=exist)

     If (exist) then

        call system("rm ./output/map_zmin_"//trim(xmin)//"_zmax_"//trim(xmax)//"_Nside_"//trim(Ns)//".fits")

     Else

        continue

     End If

     call output_map(map,header,'./output/map_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//'_Nside_'//trim(Ns)//'.fits')

     inquire(file = './output/mask_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//'_Nside_'//trim(Ns)//'.fits',exist=exist)

     If (exist) then

        call system("rm ./output/mask_zmin_"//trim(xmin)//"_zmax_"//trim(xmax)//"_Nside_"//trim(Ns)//".fits")

     Else

        continue

     End If

     call output_map(mask,header,'./output/mask_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//'_Nside_'//trim(Ns)//'.fits')

     inquire(file = './output/shot_noise_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//&
          '_Nside_'//trim(Ns)//'.fits',exist=exist)

     If (exist) then

        call system("rm ./output/shot_noise_zmin_"//trim(xmin)//"_zmax_"//trim(xmax)//&
          "_Nside_"//trim(Ns)//".fits")

     Else

        continue

     End If

     call output_map(shot_noise,header,'./output/shot_noise_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//&
          '_Nside_'//trim(Ns)//'.fits')

     inquire(file = './output/number_counts_map_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//&
          '_Nside_'//trim(Ns)//'.fits',exist=exist)

     If (exist) then

        call system("rm ./output/number_counts_map_zmin_"//trim(xmin)//"_zmax_"//trim(xmax)//&
          "_Nside_"//trim(Ns)//".fits")

     Else

        continue

     End If

     call output_map(map_nc,header,'./output/number_counts_map_zmin_'//trim(xmin)//'_zmax_'//trim(xmax)//&
          '_Nside_'//trim(Ns)//'.fits')

  End If

  call remove_dipole(nsmax,map_nc(0:npixC-1,1),RING_ORDERING,DEGREE_REMOVE_DIPOLE,multipoles,zbounds,HPX_DBADVAL)

  dip_amplitude = sqrt(multipoles(1)**2 + multipoles(2)**2 + multipoles(3)**3)

  call vec2ang(multipoles(1:3),theta,phi)

  dip_longitude = phi

  dip_latitude = HALFPI - theta

  deallocate (map, shot_noise, mask, map_nc, stat = status1)

end subroutine compute_number_counts_map

subroutine jackknife_analysis(zmin,zmax)

!  use healpix_types
!  use udgrade_nr, only: udgrade_ring, udgrade_nest
!  use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring,remove_dipole, ang2pix_ring, vec2ang
!  use fitstools, only: getsize_fits, input_map, output_map
!  use head_fits
  use fiducial
  use arrays

  Implicit none

  Character(len=4) :: Ns
  Character(len=5) :: xmin,xmax

  Real*8 :: zmin,zmax

  Integer*8 :: counter_data_points, data_index, index_jackknife_analysis, dimension_jackknife_analysis

  write(xmin,'(F5.2)') zmin
  write(xmax,'(F5.2)') zmax
  write(Ns,'(I4.4)') nsmax

  counter_data_points = 0

  Do data_index=1,number_galaxies_in_CF3

     If ( (redshift(data_index) .gt. zmin) .and. (redshift(data_index) .le. zmax) ) then

        counter_data_points = counter_data_points + 1

     Else

        continue

     End If

  End Do

  allocate (jackknife_data_indices(1:counter_data_points),jackknife_dipole_amplitude(1:counter_data_points),&
       jackknife_galactic_longitude(1:counter_data_points),jackknife_galactic_latitude(1:counter_data_points),stat=status1)

  dimension_jackknife_analysis = counter_data_points

  counter_data_points = 1

  Do data_index=1,number_galaxies_in_CF3

     If ( (redshift(data_index) .gt. zmin) .and. (redshift(data_index) .le. zmax) ) then

        jackknife_data_indices(counter_data_points) = data_index

        counter_data_points = counter_data_points + 1

     Else

        continue

     End If

  End Do

  Do index_jackknife_analysis=1,dimension_jackknife_analysis

     call compute_number_counts_map(zmin,zmax,.true.,jackknife_data_indices(index_jackknife_analysis),&
          jackknife_dipole_amplitude(index_jackknife_analysis),jackknife_galactic_latitude(index_jackknife_analysis),&
          jackknife_galactic_longitude(index_jackknife_analysis))

  End Do

  open(UNIT_JACKKNIFE_FILE,file= PATH_TO_JACKKNIFE_ANALYSIS_OUTPUT//trim('jackknife_analysis_zmin')//'_'//trim(xmin)//&
       '_zmax_'//trim(xmax)//'_Nside_'//trim(Ns)//'.txt')

  Do index_jackknife_analysis=1,dimension_jackknife_analysis

     write(UNIT_JACKKNIFE_FILE,89) jackknife_dipole_amplitude(index_jackknife_analysis),&
          jackknife_galactic_latitude(index_jackknife_analysis),jackknife_galactic_longitude(index_jackknife_analysis)

89   Format(3E20.10)

  End Do

  close(UNIT_JACKKNIFE_FILE)

  deallocate(jackknife_data_indices,jackknife_dipole_amplitude,jackknife_galactic_longitude,&
       jackknife_galactic_latitude,stat=status1)

end subroutine jackknife_analysis

end module functions
