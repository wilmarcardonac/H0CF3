import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py 

marker = ['bo','rs','go','bv','r^','g<','b>','rD','g+','bx','r*','bo']

# NSIDE = 1

z_mean, A, A_mean, AL68, AR68, AL95, AR95, LA, LA_mean, LAL68, LAR68, LAL95, LAR95, LO, LO_mean, LOL68, LOR68, LOL95, LOR95 = np.loadtxt('./output/nside-1/confidence_limits.txt',unpack=True)

# NSIDE = 2

z_mean_2, A_2, A_mean_2, AL68_2, AR68_2, AL95_2, AR95_2, LA_2, LA_mean_2, LAL68_2, LAR68_2, LAL95_2, LAR95_2, LO_2, LO_mean_2, LOL68_2, LOR68_2, LOL95_2, LOR95_2 = np.loadtxt('./output/nside-2/confidence_limits.txt',unpack=True)

# NSIDE = 4

z_mean_4, A_4, A_mean_4, AL68_4, AR68_4, AL95_4, AR95_4, LA_4, LA_mean_4, LAL68_4, LAR68_4, LAL95_4, LAR95_4, LO_4, LO_mean_4, LOL68_4, LOR68_4, LOL95_4, LOR95_4 = np.loadtxt('./output/nside-4/confidence_limits.txt',unpack=True)

# NSIDE = 8 

z_mean_8, A_8, A_mean_8, AL68_8, AR68_8, AL95_8, AR95_8, LA_8, LA_mean_8, LAL68_8, LAR68_8, LAL95_8, LAR95_8, LO_8, LO_mean_8, LOL68_8, LOR68_8, LOL95_8, LOR95_8 = np.loadtxt('./output/nside-8/confidence_limits.txt',unpack=True)

# NSIDE = 16

z_mean_16, A_16, A_mean_16, AL68_16, AR68_16, AL95_16, AR95_16, LA_16, LA_mean_16, LAL68_16, LAR68_16, LAL95_16, LAR95_16, LO_16, LO_mean_16, LOL68_16, LOR68_16, LOL95_16, LOR95_16 = np.loadtxt('./output/nside-16/confidence_limits.txt',unpack=True)

# NSIDE = 32

z_mean_32, A_32, A_mean_32, AL68_32, AR68_32, AL95_32, AR95_32, LA_32, LA_mean_32, LAL68_32, LAR68_32, LAL95_32, LAR95_32, LO_32, LO_mean_32, LOL68_32, LOR68_32, LOL95_32, LOR95_32 = np.loadtxt('./output/nside-32/confidence_limits.txt',unpack=True)

# DIPOLE AMPLITUDE

fig = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],A[index],yerr=[[A_mean[index]-AL68[index]],[AR68[index] - A_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10)

        py.errorbar(z_mean_2[index],A_2[index],yerr=[[A_mean_2[index]-AL68_2[index]],[AR68_2[index] - A_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10)

        py.errorbar(z_mean_4[index],A_4[index],yerr=[[A_mean_4[index]-AL68_4[index]],[AR68_4[index] - A_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10)

        py.errorbar(z_mean_8[index],A_8[index],yerr=[[A_mean_8[index]-AL68_8[index]],[AR68_8[index] - A_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10)

        py.errorbar(z_mean_16[index],A_16[index],yerr=[[A_mean_16[index]-AL68_16[index]],[AR68_16[index] - A_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10)

        py.errorbar(z_mean_32[index],A_32[index],yerr=[[A_mean_32[index]-AL68_32[index]],[AR68_32[index] - A_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10)

    else:

        py.errorbar(z_mean[index],A[index],yerr=[[A_mean[index]-AL68[index]],[AR68[index] - A_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_2[index],A_2[index],yerr=[[A_mean_2[index]-AL68_2[index]],[AR68_2[index] - A_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_4[index],A_4[index],yerr=[[A_mean_4[index]-AL68_4[index]],[AR68_4[index] - A_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10)
        
        py.errorbar(z_mean_8[index],A_8[index],yerr=[[A_mean_8[index]-AL68_8[index]],[AR68_8[index] - A_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_16[index],A_16[index],yerr=[[A_mean_16[index]-AL68_16[index]],[AR68_16[index] - A_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_32[index],A_32[index],yerr=[[A_mean_32[index]-AL68_32[index]],[AR68_32[index] - A_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10)

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

py.legend(loc=0,numpoints=1,ncol=2)

py.savefig("./output/Dipole_amplitude_nside.pdf",format="pdf", dpi=1000)

# DIPOLE LATITUDE

fig2 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10)

        py.errorbar(z_mean_2[index],LA_2[index],yerr=[[LA_mean_2[index]-LAL68_2[index]],[LAR68_2[index] - LA_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10)

        py.errorbar(z_mean_4[index],LA_4[index],yerr=[[LA_mean_4[index]-LAL68_4[index]],[LAR68_4[index] - LA_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10)

        py.errorbar(z_mean_8[index],LA_8[index],yerr=[[LA_mean_8[index]-LAL68_8[index]],[LAR68_8[index] - LA_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10)

        py.errorbar(z_mean_16[index],LA_16[index],yerr=[[LA_mean_16[index]-LAL68_16[index]],[LAR68_16[index] - LA_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10)

        py.errorbar(z_mean_32[index],LA_32[index],yerr=[[LA_mean_32[index]-LAL68_32[index]],[LAR68_32[index] - LA_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10)

    else:

        py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_2[index],LA_2[index],yerr=[[LA_mean_2[index]-LAL68_2[index]],[LAR68_2[index] - LA_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_4[index],LA_4[index],yerr=[[LA_mean_4[index]-LAL68_4[index]],[LAR68_4[index] - LA_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_8[index],LA_8[index],yerr=[[LA_mean_8[index]-LAL68_8[index]],[LAR68_8[index] - LA_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_16[index],LA_16[index],yerr=[[LA_mean_16[index]-LAL68_16[index]],[LAR68_16[index] - LA_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_32[index],LA_32[index],yerr=[[LA_mean_32[index]-LAL68_32[index]],[LAR68_32[index] - LA_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10)


#py.xlim(65.,83.)
#py.ylim(-2.,83.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel(r'Dipole Latitude $[\degree]$', fontsize=18)
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)
py.legend(loc=0,numpoints=1,ncol=2)
py.savefig("./output/Dipole_latitude_nside.pdf",format="pdf", dpi=1000)

# DIPOLE LONGITUDE

fig3 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10)

        py.errorbar(z_mean_2[index],LO_2[index],yerr=[[LO_mean_2[index]-LOL68_2[index]],[LOR68_2[index] - LO_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10)

        py.errorbar(z_mean_4[index],LO_4[index],yerr=[[LO_mean_4[index]-LOL68_4[index]],[LOR68_4[index] - LO_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10)

        py.errorbar(z_mean_8[index],LO_8[index],yerr=[[LO_mean_8[index]-LOL68_8[index]],[LOR68_8[index] - LO_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10)

        py.errorbar(z_mean_16[index],LO_16[index],yerr=[[LO_mean_16[index]-LOL68_16[index]],[LOR68_16[index] - LO_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10)

        py.errorbar(z_mean_32[index],LO_32[index],yerr=[[LO_mean_32[index]-LOL68_32[index]],[LOR68_32[index] - LO_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10)

    else:

        py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_2[index],LO_2[index],yerr=[[LO_mean_2[index]-LOL68_2[index]],[LOR68_2[index] - LO_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_4[index],LO_4[index],yerr=[[LO_mean_4[index]-LOL68_4[index]],[LOR68_4[index] - LO_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_8[index],LO_8[index],yerr=[[LO_mean_8[index]-LOL68_8[index]],[LOR68_8[index] - LO_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_16[index],LO_16[index],yerr=[[LO_mean_16[index]-LOL68_16[index]],[LOR68_16[index] - LO_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10)

        py.errorbar(z_mean_32[index],LO_32[index],yerr=[[LO_mean_32[index]-LOL68_32[index]],[LOR68_32[index] - LO_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10)

#py.xlim(65.,83.)
#py.ylim(-2.,83.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel(r'Dipole Longitude $[\degree]$', fontsize=18)
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)
py.legend(loc=0,numpoints=1,ncol=2)
py.savefig("./output/Dipole_longitude_nside.pdf",format="pdf", dpi=1000)

exit()
