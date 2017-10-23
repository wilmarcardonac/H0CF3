import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py 

marker = ['o','s','*','o','r^','g<','b>','rD','g+','bx','r*','bo']

# NSIDE = 1

z_mean, A, A_mean, AL68, AR68, AL95, AR95, LA, LA_mean, LAL68, LAR68, LAL95, LAR95, LO, LO_mean, LOL68, LOR68, LOL95, LOR95, A_noise, LA_noise, LO_noise = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-1/confidence_limits.txt',unpack=True)

z_mean_jla, A_jla, A_mean_jla, AL68_jla, AR68_jla, AL95_jla, AR95_jla, LA_jla, LA_mean_jla, LAL68_jla, LAR68_jla, LAL95_jla, LAR95_jla, LO_jla, LO_mean_jla, LOL68_jla, LOR68_jla, LOL95_jla, LOR95_jla, A_noise_jla, LA_noise_jla, LO_noise_jla = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-1/confidence_limits_jla.txt',unpack=True)

# NSIDE = 2

z_mean_2, A_2, A_mean_2, AL68_2, AR68_2, AL95_2, AR95_2, LA_2, LA_mean_2, LAL68_2, LAR68_2, LAL95_2, LAR95_2, LO_2, LO_mean_2, LOL68_2, LOR68_2, LOL95_2, LOR95_2, A_noise_2, LA_noise_2, LO_noise_2 = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-2/confidence_limits.txt',unpack=True)

z_mean_2_jla, A_2_jla, A_mean_2_jla, AL68_2_jla, AR68_2_jla, AL95_2_jla, AR95_2_jla, LA_2_jla, LA_mean_2_jla, LAL68_2_jla, LAR68_2_jla, LAL95_2_jla, LAR95_2_jla, LO_2_jla, LO_mean_2_jla, LOL68_2_jla, LOR68_2_jla, LOL95_2_jla, LOR95_2_jla, A_noise_2_jla, LA_noise_2_jla, LO_noise_2_jla = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-2/confidence_limits_jla.txt',unpack=True)

# NSIDE = 4

z_mean_4, A_4, A_mean_4, AL68_4, AR68_4, AL95_4, AR95_4, LA_4, LA_mean_4, LAL68_4, LAR68_4, LAL95_4, LAR95_4, LO_4, LO_mean_4, LOL68_4, LOR68_4, LOL95_4, LOR95_4, A_noise_4, LA_noise_4, LO_noise_4 = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-4/confidence_limits.txt',unpack=True)

z_mean_4_jla, A_4_jla, A_mean_4_jla, AL68_4_jla, AR68_4_jla, AL95_4_jla, AR95_4_jla, LA_4_jla, LA_mean_4_jla, LAL68_4_jla, LAR68_4_jla, LAL95_4_jla, LAR95_4_jla, LO_4_jla, LO_mean_4_jla, LOL68_4_jla, LOR68_4_jla, LOL95_4_jla, LOR95_4_jla, A_noise_4_jla, LA_noise_4_jla, LO_noise_4_jla = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-4/confidence_limits_jla.txt',unpack=True)

# NSIDE = 8 

z_mean_8, A_8, A_mean_8, AL68_8, AR68_8, AL95_8, AR95_8, LA_8, LA_mean_8, LAL68_8, LAR68_8, LAL95_8, LAR95_8, LO_8, LO_mean_8, LOL68_8, LOR68_8, LOL95_8, LOR95_8, A_noise_8, LA_noise_8, LO_noise_8 = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-8/confidence_limits.txt',unpack=True)

z_mean_8_jla, A_8_jla, A_mean_8_jla, AL68_8_jla, AR68_8_jla, AL95_8_jla, AR95_8_jla, LA_8_jla, LA_mean_8_jla, LAL68_8_jla, LAR68_8_jla, LAL95_8_jla, LAR95_8_jla, LO_8_jla, LO_mean_8_jla, LOL68_8_jla, LOR68_8_jla, LOL95_8_jla, LOR95_8_jla, A_noise_8_jla, LA_noise_8_jla, LO_noise_8_jla = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-8/confidence_limits_jla.txt',unpack=True)

# NSIDE = 16

z_mean_16, A_16, A_mean_16, AL68_16, AR68_16, AL95_16, AR95_16, LA_16, LA_mean_16, LAL68_16, LAR68_16, LAL95_16, LAR95_16, LO_16, LO_mean_16, LOL68_16, LOR68_16, LOL95_16, LOR95_16, A_noise_16, LA_noise_16, LO_noise_16 = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-16/confidence_limits.txt',unpack=True)

z_mean_16_jla, A_16_jla, A_mean_16_jla, AL68_16_jla, AR68_16_jla, AL95_16_jla, AR95_16_jla, LA_16_jla, LA_mean_16_jla, LAL68_16_jla, LAR68_16_jla, LAL95_16_jla, LAR95_16_jla, LO_16_jla, LO_mean_16_jla, LOL68_16_jla, LOR68_16_jla, LOL95_16_jla, LOR95_16_jla, A_noise_16_jla, LA_noise_16_jla, LO_noise_16_jla = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-16/confidence_limits_jla.txt',unpack=True)

# NSIDE = 32

z_mean_32, A_32, A_mean_32, AL68_32, AR68_32, AL95_32, AR95_32, LA_32, LA_mean_32, LAL68_32, LAR68_32, LAL95_32, LAR95_32, LO_32, LO_mean_32, LOL68_32, LOR68_32, LOL95_32, LOR95_32, A_noise_32, LA_noise_32, LO_noise_32 = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-32/confidence_limits.txt',unpack=True)

z_mean_32_jla, A_32_jla, A_mean_32_jla, AL68_32_jla, AR68_32_jla, AL95_32_jla, AR95_32_jla, LA_32_jla, LA_mean_32_jla, LAL68_32_jla, LAR68_32_jla, LAL95_32_jla, LAR95_32_jla, LO_32_jla, LO_mean_32_jla, LOL68_32_jla, LOR68_32_jla, LOL95_32_jla, LOR95_32_jla, A_noise_32_jla, LA_noise_32_jla, LO_noise_32_jla = np.loadtxt('./output/zmin_0.01_zmax_0.05/nside-32/confidence_limits_jla.txt',unpack=True)

# DIPOLE AMPLITUDE

fig = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],A[index],yerr=[[A_mean[index]-AL68[index]],[AR68[index] - A_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10,marker=marker[0],markersize=0.6)

        py.errorbar(z_mean_jla[index],A_jla[index],yerr=[[A_mean_jla[index]-AL68_jla[index]],[AR68_jla[index] - A_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='JLA Nside = 1',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],A_2[index],yerr=[[A_mean_2[index]-AL68_2[index]],[AR68_2[index] - A_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10,marker=marker[0],markersize=0.5)

        py.errorbar(z_mean_2_jla[index],A_2_jla[index],yerr=[[A_mean_2_jla[index]-AL68_2_jla[index]],[AR68_2_jla[index] - A_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='JLA Nside = 2',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],A_4[index],yerr=[[A_mean_4[index]-AL68_4[index]],[AR68_4[index] - A_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10,marker=marker[0],markersize=0.4)

        py.errorbar(z_mean_4_jla[index],A_4_jla[index],yerr=[[A_mean_4_jla[index]-AL68_4_jla[index]],[AR68_4_jla[index] - A_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='JLA Nside = 4',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],A_8[index],yerr=[[A_mean_8[index]-AL68_8[index]],[AR68_8[index] - A_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10,marker=marker[0],markersize=0.3)

        py.errorbar(z_mean_8_jla[index],A_8_jla[index],yerr=[[A_mean_8_jla[index]-AL68_8_jla[index]],[AR68_8_jla[index] - A_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='JLA Nside = 8',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],A_16[index],yerr=[[A_mean_16[index]-AL68_16[index]],[AR68_16[index] - A_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10,marker=marker[0],markersize=0.2)

        py.errorbar(z_mean_16_jla[index],A_16_jla[index],yerr=[[A_mean_16_jla[index]-AL68_16_jla[index]],[AR68_16_jla[index] - A_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='JLA Nside = 16',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],A_32[index],yerr=[[A_mean_32[index]-AL68_32[index]],[AR68_32[index] - A_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10,marker=marker[0],markersize=0.1)

        py.errorbar(z_mean_32_jla[index],A_32_jla[index],yerr=[[A_mean_32_jla[index]-AL68_32_jla[index]],[AR68_32_jla[index] - A_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='JLA Nside = 32',elinewidth=2,ms=5,marker=marker[1],markersize=1)

    else:

        py.errorbar(z_mean[index],A[index],yerr=[[A_mean[index]-AL68[index]],[AR68[index] - A_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=0.6)

        py.errorbar(z_mean_jla[index],A_jla[index],yerr=[[A_mean_jla[index]-AL68_jla[index]],[AR68_jla[index] - A_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],A_2[index],yerr=[[A_mean_2[index]-AL68_2[index]],[AR68_2[index] - A_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=0.5)

        py.errorbar(z_mean_2_jla[index],A_2_jla[index],yerr=[[A_mean_2_jla[index]-AL68_2_jla[index]],[AR68_2_jla[index] - A_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],A_4[index],yerr=[[A_mean_4[index]-AL68_4[index]],[AR68_4[index] - A_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=0.4)

        py.errorbar(z_mean_4_jla[index],A_4_jla[index],yerr=[[A_mean_4_jla[index]-AL68_4_jla[index]],[AR68_4_jla[index] - A_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)
        
        py.errorbar(z_mean_8[index],A_8[index],yerr=[[A_mean_8[index]-AL68_8[index]],[AR68_8[index] - A_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=0.3)

        py.errorbar(z_mean_8_jla[index],A_8_jla[index],yerr=[[A_mean_8_jla[index]-AL68_8_jla[index]],[AR68_8_jla[index] - A_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],A_16[index],yerr=[[A_mean_16[index]-AL68_16[index]],[AR68_16[index] - A_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=0.2)

        py.errorbar(z_mean_16_jla[index],A_16_jla[index],yerr=[[A_mean_16_jla[index]-AL68_16_jla[index]],[AR68_16_jla[index] - A_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],A_32[index],yerr=[[A_mean_32[index]-AL68_32[index]],[AR68_32[index] - A_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=0.1)

        py.errorbar(z_mean_32_jla[index],A_32_jla[index],yerr=[[A_mean_32_jla[index]-AL68_32_jla[index]],[AR68_32_jla[index] - A_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

#py.xlim(65.,83.)
#py.ylim(-1.,8.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel('Dipole Amplitude', fontsize=18)
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)

py.legend(loc=0,numpoints=1,ncol=2,fontsize=10)

py.savefig("./output/Dipole_amplitude_nside.pdf",format="pdf", dpi=1000)

# DIPOLE LATITUDE

fig2 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LA_jla[index],yerr=[[LA_mean_jla[index]-LAL68_jla[index]],[LAR68_jla[index] - LA_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='JLA Nside = 1',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LA_2[index],yerr=[[LA_mean_2[index]-LAL68_2[index]],[LAR68_2[index] - LA_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LA_2_jla[index],yerr=[[LA_mean_2_jla[index]-LAL68_2_jla[index]],[LAR68_2_jla[index] - LA_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='JLA Nside = 2',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LA_4[index],yerr=[[LA_mean_4[index]-LAL68_4[index]],[LAR68_4[index] - LA_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LA_4_jla[index],yerr=[[LA_mean_4_jla[index]-LAL68_4_jla[index]],[LAR68_4_jla[index] - LA_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='JLA Nside = 4',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LA_8[index],yerr=[[LA_mean_8[index]-LAL68_8[index]],[LAR68_8[index] - LA_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LA_8_jla[index],yerr=[[LA_mean_8_jla[index]-LAL68_8_jla[index]],[LAR68_8_jla[index] - LA_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='JLA Nside = 8',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LA_16[index],yerr=[[LA_mean_16[index]-LAL68_16[index]],[LAR68_16[index] - LA_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LA_16_jla[index],yerr=[[LA_mean_16_jla[index]-LAL68_16_jla[index]],[LAR68_16_jla[index] - LA_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='JLA Nside = 16',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LA_32[index],yerr=[[LA_mean_32[index]-LAL68_32[index]],[LAR68_32[index] - LA_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LA_32_jla[index],yerr=[[LA_mean_32_jla[index]-LAL68_32_jla[index]],[LAR68_32_jla[index] - LA_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='JLA Nside = 32',elinewidth=2,ms=5,marker=marker[1],markersize=1)

    else:

        py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LA_jla[index],yerr=[[LA_mean_jla[index]-LAL68_jla[index]],[LAR68_jla[index] - LA_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LA_2[index],yerr=[[LA_mean_2[index]-LAL68_2[index]],[LAR68_2[index] - LA_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LA_2_jla[index],yerr=[[LA_mean_2_jla[index]-LAL68_2_jla[index]],[LAR68_2_jla[index] - LA_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LA_4[index],yerr=[[LA_mean_4[index]-LAL68_4[index]],[LAR68_4[index] - LA_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LA_4_jla[index],yerr=[[LA_mean_4_jla[index]-LAL68_4_jla[index]],[LAR68_4_jla[index] - LA_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LA_8[index],yerr=[[LA_mean_8[index]-LAL68_8[index]],[LAR68_8[index] - LA_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LA_8_jla[index],yerr=[[LA_mean_8_jla[index]-LAL68_8_jla[index]],[LAR68_8_jla[index] - LA_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LA_16[index],yerr=[[LA_mean_16[index]-LAL68_16[index]],[LAR68_16[index] - LA_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LA_16_jla[index],yerr=[[LA_mean_16_jla[index]-LAL68_16_jla[index]],[LAR68_16_jla[index] - LA_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LA_32[index],yerr=[[LA_mean_32[index]-LAL68_32[index]],[LAR68_32[index] - LA_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LA_32_jla[index],yerr=[[LA_mean_32_jla[index]-LAL68_32_jla[index]],[LAR68_32_jla[index] - LA_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)


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
py.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
py.savefig("./output/Dipole_latitude_nside.pdf",format="pdf", dpi=1000)

# DIPOLE LONGITUDE

fig3 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LO_jla[index],yerr=[[LO_mean_jla[index]-LOL68_jla[index]],[LOR68_jla[index] - LO_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='JLA Nside = 1',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LO_2[index],yerr=[[LO_mean_2[index]-LOL68_2[index]],[LOR68_2[index] - LO_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LO_2_jla[index],yerr=[[LO_mean_2_jla[index]-LOL68_2_jla[index]],[LOR68_2_jla[index] - LO_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='JLA Nside = 2',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LO_4[index],yerr=[[LO_mean_4[index]-LOL68_4[index]],[LOR68_4[index] - LO_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LO_4_jla[index],yerr=[[LO_mean_4_jla[index]-LOL68_4_jla[index]],[LOR68_4_jla[index] - LO_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='JLA Nside = 4',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LO_8[index],yerr=[[LO_mean_8[index]-LOL68_8[index]],[LOR68_8[index] - LO_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LO_8_jla[index],yerr=[[LO_mean_8_jla[index]-LOL68_8_jla[index]],[LOR68_8_jla[index] - LO_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='JLA Nside = 8',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LO_16[index],yerr=[[LO_mean_16[index]-LOL68_16[index]],[LOR68_16[index] - LO_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LO_16_jla[index],yerr=[[LO_mean_16_jla[index]-LOL68_16_jla[index]],[LOR68_16_jla[index] - LO_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='JLA Nside = 16',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LO_32[index],yerr=[[LO_mean_32[index]-LOL68_32[index]],[LOR68_32[index] - LO_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LO_32_jla[index],yerr=[[LO_mean_32_jla[index]-LOL68_32_jla[index]],[LOR68_32_jla[index] - LO_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='JLA Nside = 32',elinewidth=2,ms=5,marker=marker[1],markersize=1)

    else:

        py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LO_jla[index],yerr=[[LO_mean_jla[index]-LOL68_jla[index]],[LOR68_jla[index] - LO_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LO_2[index],yerr=[[LO_mean_2[index]-LOL68_2[index]],[LOR68_2[index] - LO_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LO_2_jla[index],yerr=[[LO_mean_2_jla[index]-LOL68_2_jla[index]],[LOR68_2_jla[index] - LO_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LO_4[index],yerr=[[LO_mean_4[index]-LOL68_4[index]],[LOR68_4[index] - LO_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LO_4_jla[index],yerr=[[LO_mean_4_jla[index]-LOL68_4_jla[index]],[LOR68_4_jla[index] - LO_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LO_8[index],yerr=[[LO_mean_8[index]-LOL68_8[index]],[LOR68_8[index] - LO_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LO_8_jla[index],yerr=[[LO_mean_8_jla[index]-LOL68_8_jla[index]],[LOR68_8_jla[index] - LO_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LO_16[index],yerr=[[LO_mean_16[index]-LOL68_16[index]],[LOR68_16[index] - LO_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LO_16_jla[index],yerr=[[LO_mean_16_jla[index]-LOL68_16_jla[index]],[LOR68_16_jla[index] - LO_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LO_32[index],yerr=[[LO_mean_32[index]-LOL68_32[index]],[LOR68_32[index] - LO_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LO_32_jla[index],yerr=[[LO_mean_32_jla[index]-LOL68_32_jla[index]],[LOR68_32_jla[index] - LO_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

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
py.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
py.savefig("./output/Dipole_longitude_nside.pdf",format="pdf", dpi=1000)

# PLOTS SIGNAL TO NOISE 
# DIPOLE AMPLITUDE

fig = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.plot(z_mean[index],A[index]/A_noise[index],ls='None',label='Nside = 1',marker=marker[2],color='blue')

        py.plot(z_mean_jla[index],A_jla[index]/A_noise_jla[index],color='blue',ls='None',label='JLA Nside = 1',marker=marker[3])

        py.plot(z_mean_2[index],A_2[index]/A_noise_2[index],color='red',ls='None',label='Nside = 2',marker=marker[2])

        py.plot(z_mean_2_jla[index],A_2_jla[index]/A_noise_2_jla[index],color='red',ls='None',label='JLA Nside = 2',marker=marker[3])

        py.plot(z_mean_4[index],A_4[index]/A_noise_4[index],color='magenta',ls='None',label='Nside = 4',marker=marker[2])

        py.plot(z_mean_4_jla[index],A_4_jla[index]/A_noise_4_jla[index],color='magenta',ls='None',label='JLA Nside = 4',marker=marker[3])

        py.plot(z_mean_8[index],A_8[index]/A_noise_8[index],color='yellow',ls='None',label='Nside = 8',marker=marker[2])

        py.plot(z_mean_8_jla[index],A_8_jla[index]/A_noise_8_jla[index],color='yellow',ls='None',label='JLA Nside = 8',marker=marker[3])

        py.plot(z_mean_16[index],A_16[index]/A_noise_16[index],color='green',ls='None',label='Nside = 16',marker=marker[2])

        py.plot(z_mean_16_jla[index],A_16_jla[index]/A_noise_16_jla[index],color='green',ls='None',label='JLA Nside = 16',marker=marker[3])

        py.plot(z_mean_32[index],A_32[index]/A_noise_32[index],color='black',ls='None',label='Nside = 32',marker=marker[2])

        py.plot(z_mean_32_jla[index],A_32_jla[index]/A_noise_32_jla[index],color='black',ls='None',label='JLA Nside = 32',marker=marker[3])

    else:

        py.plot(z_mean[index],A[index]/A_noise[index],ls='None',marker=marker[2],color='blue')

        py.plot(z_mean_jla[index],A_jla[index]/A_noise_jla[index],color='blue',ls='None',marker=marker[3])

        py.plot(z_mean_2[index],A_2[index]/A_noise_2[index],color='red',ls='None',marker=marker[2])

        py.plot(z_mean_2_jla[index],A_2_jla[index]/A_noise_2_jla[index],color='red',ls='None',marker=marker[3])

        py.plot(z_mean_4[index],A_4[index]/A_noise_4[index],color='magenta',ls='None',marker=marker[2])

        py.plot(z_mean_4_jla[index],A_4_jla[index]/A_noise_4_jla[index],color='magenta',ls='None',marker=marker[3])

        py.plot(z_mean_8[index],A_8[index]/A_noise_8[index],color='yellow',ls='None',marker=marker[2])

        py.plot(z_mean_8_jla[index],A_8_jla[index]/A_noise_8_jla[index],color='yellow',ls='None',marker=marker[3])

        py.plot(z_mean_16[index],A_16[index]/A_noise_16[index],color='green',ls='None',marker=marker[2])

        py.plot(z_mean_16_jla[index],A_16_jla[index]/A_noise_16_jla[index],color='green',ls='None',marker=marker[3])

        py.plot(z_mean_32[index],A_32[index]/A_noise_32[index],color='black',ls='None',marker=marker[2])

        py.plot(z_mean_32_jla[index],A_32_jla[index]/A_noise_32_jla[index],color='black',ls='None',marker=marker[3])


#py.xlim(65.,83.)
#py.ylim(1.e-5,8.)
py.xlabel(r'$\bar{z}$', fontsize=18)
py.tick_params(axis='both', which='major', labelsize=15)
py.ylabel('Dipole Amplitude', fontsize=18)
py.yscale('log')
#minor_ticks = np.arange(65, 85, 1)
#ax = fig.add_subplot(111)
#ax.tick_params(axis='both', which='minor', labelsize=12)
#ax.set_xticks(minor_ticks, minor = True)
#py.gca().axes.get_yaxis().set_visible(False)
#py.tight_layout(pad=0, h_pad=0, w_pad=0)

py.legend(loc=0,numpoints=1,ncol=2,fontsize=10)

py.savefig("./output/Dipole_amplitude_nside_signal_to_noise.pdf",format="pdf", dpi=1000)
exit()
# DIPOLE LATITUDE

fig2 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LA_jla[index],yerr=[[LA_mean_jla[index]-LAL68_jla[index]],[LAR68_jla[index] - LA_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='JLA Nside = 1',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LA_2[index],yerr=[[LA_mean_2[index]-LAL68_2[index]],[LAR68_2[index] - LA_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LA_2_jla[index],yerr=[[LA_mean_2_jla[index]-LAL68_2_jla[index]],[LAR68_2_jla[index] - LA_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='JLA Nside = 2',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LA_4[index],yerr=[[LA_mean_4[index]-LAL68_4[index]],[LAR68_4[index] - LA_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LA_4_jla[index],yerr=[[LA_mean_4_jla[index]-LAL68_4_jla[index]],[LAR68_4_jla[index] - LA_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='JLA Nside = 4',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LA_8[index],yerr=[[LA_mean_8[index]-LAL68_8[index]],[LAR68_8[index] - LA_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LA_8_jla[index],yerr=[[LA_mean_8_jla[index]-LAL68_8_jla[index]],[LAR68_8_jla[index] - LA_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='JLA Nside = 8',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LA_16[index],yerr=[[LA_mean_16[index]-LAL68_16[index]],[LAR68_16[index] - LA_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LA_16_jla[index],yerr=[[LA_mean_16_jla[index]-LAL68_16_jla[index]],[LAR68_16_jla[index] - LA_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='JLA Nside = 16',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LA_32[index],yerr=[[LA_mean_32[index]-LAL68_32[index]],[LAR68_32[index] - LA_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LA_32_jla[index],yerr=[[LA_mean_32_jla[index]-LAL68_32_jla[index]],[LAR68_32_jla[index] - LA_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='JLA Nside = 32',elinewidth=2,ms=5,marker=marker[1],markersize=1)

    else:

        py.errorbar(z_mean[index],LA[index],yerr=[[LA_mean[index]-LAL68[index]],[LAR68[index] - LA_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LA_jla[index],yerr=[[LA_mean_jla[index]-LAL68_jla[index]],[LAR68_jla[index] - LA_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LA_2[index],yerr=[[LA_mean_2[index]-LAL68_2[index]],[LAR68_2[index] - LA_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LA_2_jla[index],yerr=[[LA_mean_2_jla[index]-LAL68_2_jla[index]],[LAR68_2_jla[index] - LA_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LA_4[index],yerr=[[LA_mean_4[index]-LAL68_4[index]],[LAR68_4[index] - LA_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LA_4_jla[index],yerr=[[LA_mean_4_jla[index]-LAL68_4_jla[index]],[LAR68_4_jla[index] - LA_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LA_8[index],yerr=[[LA_mean_8[index]-LAL68_8[index]],[LAR68_8[index] - LA_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LA_8_jla[index],yerr=[[LA_mean_8_jla[index]-LAL68_8_jla[index]],[LAR68_8_jla[index] - LA_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LA_16[index],yerr=[[LA_mean_16[index]-LAL68_16[index]],[LAR68_16[index] - LA_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LA_16_jla[index],yerr=[[LA_mean_16_jla[index]-LAL68_16_jla[index]],[LAR68_16_jla[index] - LA_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LA_32[index],yerr=[[LA_mean_32[index]-LAL68_32[index]],[LAR68_32[index] - LA_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LA_32_jla[index],yerr=[[LA_mean_32_jla[index]-LAL68_32_jla[index]],[LAR68_32_jla[index] - LA_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)


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
py.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
py.savefig("./output/Dipole_latitude_nside.pdf",format="pdf", dpi=1000)

# DIPOLE LONGITUDE

fig3 = py.figure()
yerrV = np.zeros(2)

for index in range(len(z_mean)):

    if (index == 0):

        py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='Nside = 1',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LO_jla[index],yerr=[[LO_mean_jla[index]-LOL68_jla[index]],[LOR68_jla[index] - LO_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',label='JLA Nside = 1',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LO_2[index],yerr=[[LO_mean_2[index]-LOL68_2[index]],[LOR68_2[index] - LO_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='Nside = 2',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LO_2_jla[index],yerr=[[LO_mean_2_jla[index]-LOL68_2_jla[index]],[LOR68_2_jla[index] - LO_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',label='JLA Nside = 2',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LO_4[index],yerr=[[LO_mean_4[index]-LOL68_4[index]],[LOR68_4[index] - LO_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='Nside = 4',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LO_4_jla[index],yerr=[[LO_mean_4_jla[index]-LOL68_4_jla[index]],[LOR68_4_jla[index] - LO_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',label='JLA Nside = 4',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LO_8[index],yerr=[[LO_mean_8[index]-LOL68_8[index]],[LOR68_8[index] - LO_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='Nside = 8',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LO_8_jla[index],yerr=[[LO_mean_8_jla[index]-LOL68_8_jla[index]],[LOR68_8_jla[index] - LO_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',label='JLA Nside = 8',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LO_16[index],yerr=[[LO_mean_16[index]-LOL68_16[index]],[LOR68_16[index] - LO_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='Nside = 16',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LO_16_jla[index],yerr=[[LO_mean_16_jla[index]-LOL68_16_jla[index]],[LOR68_16_jla[index] - LO_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',label='JLA Nside = 16',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LO_32[index],yerr=[[LO_mean_32[index]-LOL68_32[index]],[LOR68_32[index] - LO_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='Nside = 32',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LO_32_jla[index],yerr=[[LO_mean_32_jla[index]-LOL68_32_jla[index]],[LOR68_32_jla[index] - LO_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',label='JLA Nside = 32',elinewidth=2,ms=5,marker=marker[1],markersize=1)

    else:

        py.errorbar(z_mean[index],LO[index],yerr=[[LO_mean[index]-LOL68[index]],[LOR68[index] - LO_mean[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_jla[index],LO_jla[index],yerr=[[LO_mean_jla[index]-LOL68_jla[index]],[LOR68_jla[index] - LO_mean_jla[index]]],xerr=None,ecolor='blue',mfc='blue',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_2[index],LO_2[index],yerr=[[LO_mean_2[index]-LOL68_2[index]],[LOR68_2[index] - LO_mean_2[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_2_jla[index],LO_2_jla[index],yerr=[[LO_mean_2_jla[index]-LOL68_2_jla[index]],[LOR68_2_jla[index] - LO_mean_2_jla[index]]],xerr=None,ecolor='red',mfc='red',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_4[index],LO_4[index],yerr=[[LO_mean_4[index]-LOL68_4[index]],[LOR68_4[index] - LO_mean_4[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_4_jla[index],LO_4_jla[index],yerr=[[LO_mean_4_jla[index]-LOL68_4_jla[index]],[LOR68_4_jla[index] - LO_mean_4_jla[index]]],xerr=None,ecolor='magenta',mfc='magenta',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_8[index],LO_8[index],yerr=[[LO_mean_8[index]-LOL68_8[index]],[LOR68_8[index] - LO_mean_8[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_8_jla[index],LO_8_jla[index],yerr=[[LO_mean_8_jla[index]-LOL68_8_jla[index]],[LOR68_8_jla[index] - LO_mean_8_jla[index]]],xerr=None,ecolor='yellow',mfc='yellow',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_16[index],LO_16[index],yerr=[[LO_mean_16[index]-LOL68_16[index]],[LOR68_16[index] - LO_mean_16[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_16_jla[index],LO_16_jla[index],yerr=[[LO_mean_16_jla[index]-LOL68_16_jla[index]],[LOR68_16_jla[index] - LO_mean_16_jla[index]]],xerr=None,ecolor='green',mfc='green',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

        py.errorbar(z_mean_32[index],LO_32[index],yerr=[[LO_mean_32[index]-LOL68_32[index]],[LOR68_32[index] - LO_mean_32[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=3,ms=10,marker=marker[0],markersize=1)

        py.errorbar(z_mean_32_jla[index],LO_32_jla[index],yerr=[[LO_mean_32_jla[index]-LOL68_32_jla[index]],[LOR68_32_jla[index] - LO_mean_32_jla[index]]],xerr=None,ecolor='black',mfc='black',fmt='',ls='None',elinewidth=2,ms=5,marker=marker[1],markersize=1)

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
py.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
py.savefig("./output/Dipole_longitude_nside.pdf",format="pdf", dpi=1000)


exit()
