import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py 

marker = ['bo','rs','go','bv','r^','g<','b>','rD','g+','bx','r*','bo']

z_mean, A, A_mean, AL68, AR68, AL95, AR95, LA, LA_mean, LAL68, LAR68, LAL95, LAR95, LO, LO_mean, LOL68, LOR68, LOL95, LOR95, A_noise, LA_noise, LO_noise = np.loadtxt('./output/confidence_limits.txt',unpack=True)

# DIPOLE AMPLITUDE

fig = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    py.errorbar(z_mean[index],A[index],yerr=[[A_mean[index]-AL68[index]],[AR68[index] - A_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Dipole Amplitude',elinewidth=3,ms=10)
    py.plot(z_mean[index],A_noise[index],color='blue',marker='*')

#py.xlim(65.,83.)
#py.ylim(-2.,83.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel('Dipole Amplitude', fontsize=18)
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)
py.savefig("./output/Dipole_amplitude.pdf",format="pdf", dpi=1000)

# DIPOLE LONGITUDE

fig2 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Dipole Longitude',elinewidth=3,ms=10)
    py.plot(z_mean[index],LA_noise[index],color='blue',marker='*')
#py.xlim(65.,83.)
#py.ylim(-2.,83.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel('Dipole Longitude', fontsize=18)
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)
py.savefig("./output/Dipole_longitude.pdf",format="pdf", dpi=1000)

# DIPOLE LATITUDE

fig3 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Dipole Latitude',elinewidth=3,ms=10)
    py.plot(z_mean[index],LO_noise[index],color='blue',marker='*')
#py.xlim(65.,83.)
#py.ylim(-2.,83.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel('Dipole Latitude', fontsize=18)
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)
py.savefig("./output/Dipole_latitude.pdf",format="pdf", dpi=1000)

exit()
