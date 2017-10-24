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
 
cf3data.write("# REDSHIFT    LUMINOSITY DISTANCE    GALACTIC LONGITUDE    GALACTIC LATITUDE    CORRECTED REDSHIFT\n")

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
    string_aux_4 = '%13.10E' % data[index,6]
    cf3data.write(string_aux_0+" "+string_aux_1+" "+string_aux_2+" "+string_aux_3+" "+string_aux_4+"\n")


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

