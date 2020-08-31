# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

outdir = '../oiso_output/'

iso = 'd18'
# iso = 'd17'
# iso = 'capd17'

dr_cw = np.loadtxt(outdir + iso+'r_cw_dyn.txt')
dr_ha1 = np.loadtxt(outdir + iso+'rp_ha_1_dyn.txt')
dr_ha2 = np.loadtxt(outdir + iso+'rp_ha_2_dyn.txt')
dr_ha3 = np.loadtxt(outdir + iso+'rp_ha_3_dyn.txt')

figsize = (8,5)
fig = plt.figure(figsize=figsize)

nx =2 
ny =2 


axes = [[plt.subplot2grid((ny,nx), (i,j)) for j in range(nx)] for i in range(ny)]

res = [
    [dr_cw,dr_ha1]
    ,[dr_ha2,dr_ha3]
    ]

titles = [
    ['Shale','Basalt']
    ,['Dike','Gabbro']
    ]

if iso == 'd18': ylabel = r"$\mathregular{\delta^{18}O}$"+ ' ('+u'$‰$'+')'
if iso == 'd17': ylabel = r"$\mathregular{\delta^{17}O}$"+ ' ('+u'$‰$'+')'
if iso == 'capd17': ylabel = r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')'

for i in range(ny):
    for j in range(nx):
        axes[i][j].fill_between(res[i][j][:,0],np.amin(res[i][j][:,1:3],axis=1),np.amax(res[i][j][:,1:3],axis=1)
            ,hatch ='//////',linestyle = 'solid',edgecolor = 'k',facecolor='None')
        axes[i][j].invert_xaxis()
        if i==ny-1:axes[i][j].set_xlabel('Age (Ma)')
        else:axes[i][j].set_xticklabels([])
        
        if j==0:axes[i][j].set_ylabel(ylabel)
        # else:axes[i][j].set_yticklabels([])
        
        axes[i][j].text(0.5,1.02,titles[i][j],horizontalalignment='center',transform=axes[i][j].transAxes)
        

fig.align_labels(axes)
fig.subplots_adjust(left=0.12,right= 0.97 ,bottom=0.1,top = 0.95,wspace=0.23,hspace=0.15)
outdir = 'C:/Users/YK/Desktop/Isotope/'

plt.savefig(outdir+'O_'+iso+'r_Phan.svg', transparent=True)
plt.savefig(outdir+'O_'+iso+'r_Phan.pdf', transparent=True)

plt.show()
plt.clf()
plt.close()