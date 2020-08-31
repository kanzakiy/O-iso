# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def d2dp(d):
    dp = 1e3*np.log(d*1e-3+1.)
    return dp

outdir = '../oiso_output/'

iso = 'd18'

d18r_cw = np.loadtxt(outdir + iso+'r_cw_dyn.txt')
d18r_ha1 = np.loadtxt(outdir + iso+'rp_ha_1_dyn.txt')
d18r_ha2 = np.loadtxt(outdir + iso+'rp_ha_2_dyn.txt')
d18r_ha3 = np.loadtxt(outdir + iso+'rp_ha_3_dyn.txt')

iso = 'capd17'
capd17r_cw = np.loadtxt(outdir + iso+'r_cw_dyn.txt')
capd17r_ha1 = np.loadtxt(outdir + iso+'rp_ha_1_dyn.txt')
capd17r_ha2 = np.loadtxt(outdir + iso+'rp_ha_2_dyn.txt')
capd17r_ha3 = np.loadtxt(outdir + iso+'rp_ha_3_dyn.txt')

c = d18r_cw[:,0]
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['black','red'])
mymap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mymap)
cmap.set_array([])

figsize = (8,5)
fig = plt.figure(figsize=figsize)

nx =2 
ny =2 


axes = [[plt.subplot2grid((ny,nx), (i,j)) for j in range(nx)] for i in range(ny)]

res18 = [
    [d18r_cw,d18r_ha1]
    ,[d18r_ha2,d18r_ha3]
    ]
    
res17 = [
    [capd17r_cw,capd17r_ha1]
    ,[capd17r_ha2,capd17r_ha3]
    ]

titles = [
    ['Shale','Basalt']
    ,['Dike','Gabbro']
    ]

for i in range(ny):
    for j in range(nx):
        for k in range(c.shape[0]):
            axes[i][j].plot(d2dp(res18[i][j][k,1]),res17[i][j][k,1],'o',c=cmap.to_rgba(c[k]))
            axes[i][j].plot(d2dp(res18[i][j][k,2]),res17[i][j][k,2],'o',c=cmap.to_rgba(c[k]))
            if i==ny-1:axes[i][j].set_xlabel(r"$\mathregular{\delta'^{18}O}$"+ ' ('+u'$‰$'+')')
            else:axes[i][j].set_xticklabels([])
            
            if j==0:axes[i][j].set_ylabel(r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')')
            
        axes[i][j].text(0.5,1.02,titles[i][j],horizontalalignment='center',transform=axes[i][j].transAxes)
        

# cbaxes = fig.add_axes([0.15, 0.1, 0.75, 0.02]) 
cbaxes = fig.add_axes([0.9, 0.15, 0.01, 0.75]) 
cticks = [500,400,300,200,100,0]
cbar = fig.colorbar(cmap, cax = cbaxes, ticks = cticks, orientation = 'vertical')
cbar.set_label('Age (Ma)', rotation=90)

fig.align_labels(axes)
fig.subplots_adjust(left=0.12,right= 0.88 ,bottom=0.1,top = 0.95,wspace=0.23,hspace=0.15)
outdir = 'C:/Users/YK/Desktop/Isotope/'

plt.savefig(outdir+'O_d18-capD17r_Phan.svg', transparent=True)
plt.savefig(outdir+'O_d18-capD17r_Phan.pdf', transparent=True)

plt.show()
plt.clf()
plt.close()