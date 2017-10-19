import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as py
from matplotlib.colors import LogNorm
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord

marker = ['bo','rs','go','bv','r^','g<','b>','rD','g+','bx','r*','bo']

exefile = open('execution_information.txt','w')

# WE LOAD DATA FROM COSMICFLOWS-3 

table1 = np.dtype([('The Catalogue of Principal Galaxies (PGC) Number',np.str_,7),('Luminosity distance (weighted average if more than one source)',np.float64),('Number of distance sources',np.int32),('Luminosity distance modulus (weighted average if more than one source)',np.float64),('1 sigma uncertainty in distance modulus (NOTE! NOT fractional error in distance as given in CF2)',np.float64),('Distance method: Cepheid Period-Luminosity relation',np.str_,1),('Distance method: Tip of the Red Giant Branch from collaboration HST program (Jacobs et al. 2009)',np.str_,1),('Distance method: Tip of the Red Giant Branch (literature)',np.str_,1),('Distance method: miscellaneous (RR Lyrae, Horizontal Branch, Eclipsing Binary)',np.str_,1),('Distance method: Surface Brightness Fluctuation',np.str_,1),('Distance method: type Ia supernovae',np.str_,1),('Distance method: optical (I band) spiral luminosity-rotation correlation (Tully-Fisher)',np.str_,1),('Distance method: Spitzer [3.6] band spiral luminosity-rotation correlation (Tully-Fisher)',np.str_,1),('Distance method: E/S0 Fundamental Plane - sources ENEAR, EFAR, SMAC',np.str_,1),('Luminosity distance modulus carried over from Cosmicflows-2.1',np.float64),('1 sigma uncertainty in distance modulus carried over from Cosmicflows-2.1',np.float64),('Supernova ID',np.str_,7),('Number of separate analyses of this supernova',np.int32),('Luminosity distance modulus of supernova averaged over all contributions',np.float64),('Luminosity distance modulus determined from luminosity-rotation correlation using Spitzer [3.6] photometry',np.float64),('1 sigma uncertainty in distance modulus determined from luminosity-rotation correlation using Spitzer [3.6] photometry',np.float64),('Luminosity distance modulus from 6dFGS (Springob et al. 2014) with CF3 zero point',np.float64),('1 sigma uncertainty in distance modulus from 6dFGS (Springob et al. 2014) with CF3 zero point',np.float64),('6dFGS morphological type code',np.str_,3),('Right Ascension (J2000)',np.str_,8),('Declination (J2000)',np.str_,9),('Galactic longitude',np.float64),('Galactic latitude',np.float64),('Supergalactic longitude',np.float64),('Supergalactic latitude',np.float64),('Morphological type; RC3 numeric code',np.float64),('Reddening at B band from Schlafly & Finkbeiner (2011,ApJ,737,103) dust maps',np.float64),('Total B magnitude from LEDA',np.float64),('2MASS Ks magnitude, extinction corrected from Huchra et al. 2012 or else Lavaux-Hudson 2011',np.float64),('Heliocentric velocity',np.int32),('Velocity in Galactic standard of rest (circular velocity at Sun of 239 km/s; total velocity 251 km/s toward l=90,b=0)',np.int32),('Velocity in Local Sheet standard of rest (Tully et al. 2008)',np.int32),('Velocity in CMB standard of rest (Fixsen et al. 1996)',np.int32),('Velocity in CMB standard of rest adjusted in accordance with a cosmological model with Omega_matter=0.27 and Omega_Lambda=0.73',np.int32),('Common name',np.str_,11),('2MASS nest group ID',np.int32),('Number of galaxies with distance measures in group',np.int32),('Luminosity distance modulus to group; weighted average of all contributions',np.float64),('1 sigma error in group distance modulus',np.float64),('Luminosity distance to group',np.float64),('Abell cluster ID; ASxxx are from southern supplement list',np.str_,5),('Alternate name for group or cluster',np.str_,9),('Number of galaxies with positions and velocities in group',np.int32),('PGC ID of dominant galaxy in group',np.int32),('Galactic latitude of group',np.float64),('Galactic longitude of group',np.float64),('Supergalactic longitude of group',np.float64),('Supergalactic latitude of group',np.float64),('Log summed K luminosity of group, adjusted by correction factor for lost light; assumed distance is Vmgp/75',np.float64),('Luminosity selection function correction factor',np.str_,5),('Projected velocity dispersion anticipated by corrected luminosity',np.str_,4),('Projected second turnaround radius anticipated by corrected intrinsic luminosity; assumed distance is Vmgp/75',np.float64),('Group heliocentric velocity',np.int32),('Group velocity in Galactic standard of rest (circular velocity at Sun of 239 km/s; total velocity 251 km/s toward l=90,b=0)',np.int32),('Group velocity in Local Sheet standard of rest (Tully et al. 2008)',np.int32),('Group velocity in CMB standard of rest (Fixsen et al. 1996)',np.int32),('Group velocity in CMB standard of rest adjusted in accordance with a cosmological model with Omega_matter=0.27 and Omega_Lambda=0.73',np.int32),('Group r.m.s. velocity dispersion',np.str_,4),('Group mass (x10^12) from virial theorem with bi-weight dispersion and radius parameters; assumed distance is Vmgp/75',np.float64),('Group mass (x10^12) based on corrected intrinsic luminosity and M/L prescription; assumed distance is Vmgp/75',np.float64),('Crook et al. (2007) low density group ID',np.int32),('Crook et al. (2007) high density group ID',np.int32),('Group ID from 2MASS++ catalog of Lavaux & Hudson (2011)',np.int32),('Makarov & Karachentsev (2011, MNRAS 412, 2498) group ID',np.int32),('Internal group ID',np.int32)])

pgc,Dist,Nd,DM,eDM,C,T,L,M,S,N,H,I,F,DM2,eD2,SNIa,Ns,DMsn,DMsp,eDsp,DM6d,eD6d,Mt,RAJ,DeJ,Glon,Glat,SGL,SGB,Ty,Asf,Btot,Ks,Vhel,Vgsr,Vls,Vcmb,Vmod,Name,Nest,Ndgp,DMgp,eDgp,Dgp,Abell,GroupName,NV,PGC1,Glongp,Glatgp,SGLgp,SGBgp,lgLgp,cf,sigp,R2t,Vhgp,Vggp,Vlsgp,Vcgp,Vmgp,Vrms,bwMass12,L_Mass12,LDC,HDC,twoMplusplus,MKgp,Icnt = np.loadtxt('CF-3.txt',unpack=True,dtype=table1,skiprows=5,delimiter=',',usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69]) 

# LOADING JLA DATA

table2 = np.dtype([('name',np.str_,9), ('zcmb',np.float64), ('zhel',np.float64), ('dz',np.float64), ('mb',np.float64), ('dmb',np.float64), ('x1',np.float64), ('dx1',np.float64), ('color',np.float64), ('dcolor',np.float64), ('3rdvar',np.float64), ('d3rdvar',np.float64), ('tmax',np.float64), ('dtmax',np.float64), ('cov_m_s',np.float64), ('cov_m_c',np.float64), ('cov_s_c',np.float64), ('set',np.float64), ('ra',np.float64), ('dec',np.float64), ('biascor',np.float64)])

name, zcmb, zhel, dz, mb, dmb, x1, dx1, color, dcolor, threerdvar, d3rdvar, tmax, dtmax, cov_m_s, cov_m_c, cov_s_c, set, ra, dec, biascor = np.loadtxt('jla_likelihood_v6/data/jla_lcparams.txt', unpack=True,dtype=table2,delimiter=' ') 

# WE SET A FEW CONSTANT VALUES FOR THE BINNING AND EXECUTION

CF3 = False

JLA = True 

z_min = 0.01 

z_max = 0.05 

z_step = 0.005

speed_of_light = 299792.458 # km/s

# WE PRINT OUT THE NUMBER OF DISTANCES IN THE CATALOGUE 

exefile.write("THE NUMBER OF DISTANCE MEASUREMENTS IN THE CF3 CATALOGUE IS: "+str(len(Dist))+"\n")

# WE DEFINE AN ARRAY CONTAINING REDSHIFT (INCLUDING NO CORRECTION), LUMINOSITY DISTANCE, GALACTIC LONGITUDE, GALACTIC LATITUDE, RIGHT ASCENSION, DECLINATION, CORRECTED REDSHIFT 

data = np.zeros((len(Dist),7))

data_jla = np.zeros((len(name),3))

cf3data = open('cosmic_flows_data.txt','w')

jladata = open('jla_data.txt','w')

jladata.write("# REDSHIFT    GALACTIC LONGITUDE     GALACTIC LATITUDE\n")
 
cf3data.write("# REDSHIFT    LUMINOSITY DISTANCE    GALACTIC LONGITUDE    GALACTIC LATITUDE\n")

for index in range(len(Dist)):
    data[index,0] = Vcmb[index]/speed_of_light # REDSHIFT 
    data[index,1] = Dist[index]                # LUMINOSITY DISTANCE
    data[index,2] = Glon[index]                # GALACTIC LONGITUDE
    data[index,3] = Glat[index]                # GALACTIC LATITUDE
    data[index,4] = RAJ[index]                 # RIGHT ASCENSION
    data[index,5] = DeJ[index]                 # DECLINATION
    data[index,6] = Vmod[index]/speed_of_light # CORRECTED REDSHIFT
    string_aux_0 = '%13.10E' % data[index,0]
    string_aux_1 = '%13.10E' % data[index,1]
    string_aux_2 = '%13.10E' % data[index,2]
    string_aux_3 = '%13.10E' % data[index,3]
    cf3data.write(string_aux_0+" "+string_aux_1+" "+string_aux_2+" "+string_aux_3+"\n")


for index in range(len(name)):
    data_jla[index,0] = zcmb[index]                   # REDSHIFT 
    c_icrs = SkyCoord(ra=ra[index]*u.degree, dec=dec[index]*u.degree, frame='icrs')
    data_jla[index,1] = c_icrs.galactic.l.value       # GALACTIC LONGITUDE
    data_jla[index,2] = c_icrs.galactic.b.value       # GALACTIC LATITUDE                
    string_aux_0 = '%13.10E' % data_jla[index,0]
    string_aux_1 = '%13.10E' % data_jla[index,1]
    string_aux_2 = '%13.10E' % data_jla[index,2]
    jladata.write(string_aux_0+" "+string_aux_1+" "+string_aux_2+"\n")

jladata.close()

cf3data.close()

exit()

exefile.write("LOWEST UNCORRECTED REDSHIFT IN THE CATALOGUE IS: "+str(np.min(data[:,0]))+"\n")

#print 'LOWEST CORRECTED REDSHIFT IN THE CATALOGUE IS: ', np.min(data[:,6])

exefile.write("HIGHEST UNCORRECTED REDSHIFT IN THE CATALOGUE IS: "+str(np.max(data[:,0]))+"\n")

#print 'HIGHEST CORRECTED REDSHIFT IN THE CATALOGUE IS: ', np.max(data[:,6])

# WE BIN DATA ACCORDING TO UNCORRECTED REDSHIFT REDSHIFT  

z_max_temp = z_min + z_step

z_min_temp = z_min

counter = 0

while (z_max_temp <= z_max):

#    datafile = open('redshift_bin_z_min_'+str(z_min_temp)+'_z_max'+str(z_max_temp)+'_CF3_dataset.txt','w')

    hp.graticule(dpar=10.,dmer=10.,coord='G')

#    datafile.write("#REDSHIFT    LUMINOSITY DISTANCE    GALACTIC LONGITUDE    GALACTIC LATITUDE    RIGHT ASCENSION    DECLINATION    CORRECTED REDSHIFT\n")

    for index in range(len(Dist)):

        if ( (data[index,0] > z_min_temp) & (data[index,0] <= z_max_temp) ):

#            datafile.write(str(data[index,0])+" "+str(data[index,1])+" "+str(data[index,2])+" "+str(data[index,3])+" "+str(data[index,4])+" "+str(data[index,5])+" "+str(data[index,6])+"\n" )

#            for index_i in range(len(map)):
                
#                if ( hp.ang2pix(Nside,data[index,2],data[index,3],lonlat=True) == index_i ):

#                    map[index_i] = map[index_i] + 1

            hp.projplot(data[index,2],data[index,3],marker[counter],lonlat=True,coord='G')

    # COMPUTE MEAN OF NUMBER OF GALAXIES PER PIXEL
#    mean_nc = np.mean(map)

    # COMPUTE SHOT NOISE IN CURRENT REDSHIFT BIN
#    for index_j in range(len(map)):

#        shot_noise[index_j] = 1./float(map[index_j])
    
    exefile.write("NUMBER OF GALAXIES IN REDSHIFT BIN Z_MIN "+str(z_min_temp)+" Z_MAX "+str(z_max_temp)+" IS: "+str(np.sum(map))+"\n")
    
    # COMPUTE NUMBER COUNTS MAP (N_i - mean_N)/mean_N WHERE N_i IS THE NUMBER OF GALAXIES IN PIXEL i 
#    for index_i in range(len(map)):

#        map_nc[index_i] = (float(map[index_i]) - mean_nc)/mean_nc

#    print map,mean_nc,(map[1] - mean_nc)/mean_nc,np.max(map_nc),np.min(map_nc),map_nc[0],map_nc[3]

#    exefile.write("MONOPOLE AND DIPOLE FOR CURRENT NUMBER COUNTS MAP IS "+str(hp.fit_dipole(map_nc))+"\n")

#    hp.mollview(map_nc,coord='G',title='Dipole amplitude is '+str(np.sqrt(hp.fit_dipole(map_nc)[1][0]**2 + hp.fit_dipole(map_nc)[1][1]**2 + hp.fit_dipole(map_nc)[1][2]**2)))

#    hp.projplot(hp.vec2ang(hp.fit_dipole(map_nc)[1],lonlat=True)[0][0],hp.vec2ang(hp.fit_dipole(map_nc)[1],lonlat=True)[1][0],marker[counter+1],lonlat=True,coord='G')

#    py.savefig('galaxy_distribution_redshift_bin_z_min_'+str(z_min_temp)+'_z_max_'+str(z_max_temp)+'_CF3.pdf')

#    hp.mollview(shot_noise,coord='G',title='Dipole amplitude is '+str(np.sqrt(hp.fit_dipole(shot_noise)[1][0]**2 + hp.fit_dipole(shot_noise)[1][1]**2 + hp.fit_dipole(shot_noise)[1][2]**2)))

#    hp.projplot(hp.vec2ang(hp.fit_dipole(shot_noise)[1],lonlat=True)[0][0],hp.vec2ang(hp.fit_dipole(shot_noise)[1],lonlat=True)[1][0],marker[counter+1],lonlat=True,coord='G')

#    py.savefig('shot_noise_galaxy_distribution_redshift_bin_z_min_'+str(z_min_temp)+'_z_max_'+str(z_max_temp)+'_CF3.pdf')

#    py.close()

#    cl = hp.anafast(map_nc,lmax=lmax)

#    nl = hp.anafast(shot_noise,lmax=lmax)

#    ell = np.arange(len(cl))

#    py.figure()

#    py.plot(ell,cl/nl,label='Signal to Noise')

#    py.plot(ell,cl,label='Signal')

#    py.plot(ell,nl,label='Noise')

#    py.xlabel(r'$\ell$')

#    py.yscale('log')

#    py.legend()

#    py.savefig('angular_power_spectra_galaxy_distribution_redshift_bin_z_min_'+str(z_min_temp)+'_z_max_'+str(z_max_temp)+'_CF3.pdf')

#    py.close()

#    exit()

#    z_min_temp = z_max_temp

    z_max_temp = z_max_temp + z_step

    counter = counter + 1

#    datafile.close()

#    map[:] = 0

#    map_nc[:] = 0.

#    shot_noise[:] = 0.

#hp.graticule(dpar=10.,dmer=10.,coord='G')

#for index in range(len(Dist)):

#    if ( (data[index,0] > 0.01) & (data[index,0] <= 0.011) ):

#        hp.projplot(data[index,2],data[index,3],marker[1],lonlat=True,coord='G')

py.savefig('galaxy_distribution_redshift_bin_z_min_'+str(z_min_temp)+'_z_max'+str(z_max_temp)+'_CF3.pdf')

exit() 

# WE COMPUTE THE NUMBER OF GALAXIES WITH TFR MEASUREMENTS IN CF2 

counter_galaxies = 0 

for index in range(len(VLS)):
    if (TFR[index] == 'H'):

        counter_galaxies = counter_galaxies + 1

print 'THERE ARE ',counter_galaxies, ' GALAXIES WITH TFR MEASUREMENTS IN CF-2.1 DATASET'

# WE DEFINE A NEW 2D ARRAY CONTAINING ONLY INFORMATION ABOUT GALAXIES WITH TFR MEASUREMENTS

data = np.zeros((counter_galaxies,2)) # 2D ARRAY CONTAINING VELOCITIES AND HUBBLE PARAMETER

counter = 0

for index in range(len(VLS)):

    if (TFR[index] == 'H'):

        data[counter,0] = VLS[index] #VLS VELOCITIES FOR GALAXIES WITH TFR DISTANCES 

        data[counter,1] = VLS[index]/D[index] # HUBBLE PARAMETER FOR THOSE GALAXIES
                                       
        counter = counter + 1

# WE MAKE A SCATTER PLOT HUBBLE PARAMETER VS LOCAL SHEET VELOCITIES

ax = py.gca()

ax.scatter(data[:,0],data[:,1],color='k',label='TFR in CF2.1')

ax.set_yscale('log')

py.ylim(4.,4.e2)

py.xlim(0.,3.e4)

# WE MAKE A PLOT OF BINNED H0 VALUES

low = np.min(data[:,0]) 

upper = np.max(data[:,0]) 

temp_low = low # TEMPORARY LOWER EDGE

step = 1.e3

temp_up = temp_low + step # TEMPORARY UPPER EDGE

means = np.zeros(int((upper-low)/step))

means_v = np.zeros(int((upper-low)/step))

std = np.zeros(int((upper-low)/step))

counter = 0

while (temp_up <= upper):

    current_array = data[np.where((data[:,0] <= temp_up) & (data[:,0] > temp_low))[0],1]

    current_array_v = data[np.where((data[:,0] <= temp_up) & (data[:,0] > temp_low))[0],0]
    
    means[counter] = np.mean(current_array)

    means_v[counter] = np.mean(current_array_v)

    std[counter] = np.std(current_array)/np.sqrt(float(len(current_array)))

    counter = counter + 1

    temp_low = temp_up

    temp_up = temp_up + step

H0_mean = np.mean(means)

print 'MEAN AND STANDARD DEVIATION FOR THE HUBBLE VALUES : ', np.mean(means), np.std(means)

weighted_mean = 0.

weight_norm = 0.

for index in range(len(means)):

    weight_norm = 1/std[index]**2 + weight_norm

for index in range(len(means)):
    
    weighted_mean = means[index]/std[index]**2/weight_norm + weighted_mean 
    
print 'WEIGHTED MEAN IS : ', weighted_mean

#py.hist(V,alpha=0.2,label='distances')

#py.savefig('histogram.pdf')

#py.close()

#fig = py.figure()


py.errorbar(means_v,means,xerr=None,yerr=std,ecolor='r',marker=marker[1],mfc='r',fmt='',ls='None',label='Binned H0')

#py.hist(V,alpha=0.2,label='distances')

Pe = np.linspace(4.e3,3.e4,num=50.)

be = np.ones(len(Pe))

be2 = np.ones(len(Pe))

be3 = np.ones(len(Pe))

for index in range(len(Pe)):
    be[index] = 74.6
    be2[index] = H0_mean
    be3[index] = weighted_mean

py.plot(Pe,be,markersize='small',color='b',label='CF2 fit: 74.6')

py.plot(Pe,be2,markersize='small',color='yellow',label='Mean of binned H0: '+str(H0_mean))

py.plot(Pe,be,markersize='small',color='magenta',label='Weighted mean of binned H0: '+str(weighted_mean))

py.title('Cosmicflows--2.1 including '+str(counter_galaxies)+' galaxies TFR')

py.ylabel(r'log Hubble Parameter [km/s/Mpc]')

py.xlabel(r'Velocity$_{LS}$ [km/s]')

py.legend(numpoints=1,loc=4)

py.savefig('H0-Fig10.pdf')

#py.close()


exit()

# We create files to write data 

dataA = open('dataA.txt','w')

dataB = open('dataB.txt','w')

dataAB = open('dataAB.txt','w')

dataC = open('dataC.txt','w')

dataABC = open('dataABC.txt','w')

# An array with header

#data.write("#Name    Period    H    Sigma_m    V    I\n")

# counter = 0

for i in range(len(n)):
    for j in range(len(N)):

        if ( (N[j] == n[i]) and ( math.log10(p[j])<=1.0 )):
            dataA.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( (N[j] == n[i]) and ((math.log10(p[j])>1.0) and (math.log10(p[j])< 1.8))):
            dataB.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( (N[j] == n[i]) and ( math.log10(p[j])>= 1.8 ) ):
            dataC.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( (N[j] == n[i]) and ( math.log10(p[j])< 1.8 ) ):
            dataAB.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( N[j] == n[i] ):
            dataABC.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )


dataA.close()
           
dataB.close()

dataAB.close()

dataC.close()

dataABC.close()
 
exit()
