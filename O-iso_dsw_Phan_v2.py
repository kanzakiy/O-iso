# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

outdir = '../oiso_output/'

dsw_d18_max = np.loadtxt(outdir + 'd18sw_dyn_max.txt')
dsw_d18_min = np.loadtxt(outdir + 'd18sw_dyn_min.txt')
dsw_d17_min = np.loadtxt(outdir + 'd17sw_dyn_min.txt')
dsw_d17_max = np.loadtxt(outdir + 'd17sw_dyn_max.txt')
dsw_capd17_min = np.loadtxt(outdir + 'capd17sw_dyn_min.txt')
dsw_capd17_max = np.loadtxt(outdir + 'capd17sw_dyn_max.txt')

dsw_d18_ss = np.loadtxt(outdir + 'd18sw_ss.txt')
dsw_d17_ss = np.loadtxt(outdir + 'd17sw_ss.txt')
dsw_capd17_ss = np.loadtxt(outdir + 'capd17sw_ss.txt')

figsize = (4,6)
fig = plt.figure(figsize=figsize)

nx =1 
ny =3 


axes = [plt.subplot2grid((ny,nx), (i,0)) for i in range(ny)]

grey = '0.5'

axes[0].fill_between(dsw_d18_ss[:,0],np.amin(dsw_d18_ss[:,1:],axis=1),np.amax(dsw_d18_ss[:,1:],axis=1)
    ,hatch ='xxxxxx',linestyle = 'solid',edgecolor = grey,facecolor='None',label = 'Steady state')
axes[0].fill_between(dsw_d18_max[:,0],np.amin(dsw_d18_min[:,1:],axis=1),np.amax(dsw_d18_max[:,1:],axis=1)
    ,hatch ='//////',linestyle = 'solid',edgecolor = 'k',facecolor='None',label = 'Dynamic')

axes[2].fill_between(dsw_capd17_ss[:,0],np.amin(dsw_capd17_ss[:,1:],axis=1),np.amax(dsw_capd17_ss[:,1:],axis=1)
    ,hatch ='xxxxxx',linestyle = 'solid',edgecolor = grey,facecolor='None')
axes[2].fill_between(dsw_capd17_min[:,0],np.amin(dsw_capd17_max[:,1:],axis=1),np.amax(dsw_capd17_min[:,1:],axis=1)
    ,hatch ='//////',linestyle = 'solid',edgecolor = 'k',facecolor='None')

axes[1].fill_between(dsw_d17_ss[:,0],np.amin(dsw_d17_ss[:,1:],axis=1),np.amax(dsw_d17_ss[:,1:],axis=1)
    ,hatch ='xxxxxx',linestyle = 'solid',edgecolor = grey,facecolor='None')
axes[1].fill_between(dsw_d17_min[:,0],np.amin(dsw_d17_max[:,1:],axis=1),np.amax(dsw_d17_min[:,1:],axis=1)
    ,hatch ='//////',linestyle = 'solid',edgecolor = 'k',facecolor='None')
    
axes[0].legend(facecolor = 'None', edgecolor = 'None',loc = 'lower right')

axes[2].set_ylim(-0.2,0.3)

axes[0].invert_xaxis()
axes[1].invert_xaxis()
axes[2].invert_xaxis()

axes[2].set_xlabel('Age (Ma)')

axes[0].set_xticklabels([])
axes[1].set_xticklabels([])

axes[0].set_ylabel(r'$\mathregular{\delta^{18}O}$'+ ' ('+u'$‰$'+')')
axes[2].set_ylabel(r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')')
axes[1].set_ylabel(r"$\mathregular{\delta^{17}O}$"+ ' ('+u'$‰$'+')')

fig.align_labels(axes)
fig.subplots_adjust(left=0.21,right= 0.97 ,bottom=0.1,top = 0.97,wspace=0.2,hspace=0.1)
outdir = 'C:/Users/YK/Desktop/Isotope/'

plt.savefig(outdir+'O_dsw_Phan.svg', transparent=True)
plt.savefig(outdir+'O_dsw_Phan.pdf', transparent=True)

plt.show()
plt.clf()
plt.close()