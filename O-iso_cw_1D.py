# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

outdir = '../oiso_output/'

sld_d18 = np.loadtxt(outdir + 'o_iso_cw_prof_sld_d18.txt')
pw_d18 = np.loadtxt(outdir + 'o_iso_cw_prof_pw_d18.txt')
sld_d17 = np.loadtxt(outdir + 'o_iso_cw_prof_sld_d17.txt')
pw_d17 = np.loadtxt(outdir + 'o_iso_cw_prof_pw_d17.txt')
sld_D17 = np.loadtxt(outdir + 'o_iso_cw_prof_sld_capD17.txt')
pw_D17 = np.loadtxt(outdir + 'o_iso_cw_prof_pw_capD17.txt')

# sld_D17[:,2:] = sld_D17[:,2:] - pw_D17[:,2:]  # to give the difference between sld and water phase D17O

n_lines = 11
c = np.linspace(-20, 0, n_lines)
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['black','red'])
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mymap)
cmap.set_array([])


figsize = (6.5,3)
fig = plt.figure(figsize=figsize)

nx =2 
ny =1 


axes = [plt.subplot2grid((ny,nx), (0,i)) for i in range(nx)]
res = [sld_d18, sld_D17]
xlabels = [  
    r'$\mathregular{\delta^{18}O}$'+ ' ('+u'$‰$'+')'
    ,r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')'
    ]

for j in range(nx):
    for i in range(n_lines):
        cbar_ax = axes[j].plot(res[j][1:,2+i],res[j][1:,1],c=cmap.to_rgba(c[i]))
        axes[j].set_ylim(0,1)
        axes[j].set_xlabel(xlabels[j])
        axes[j].invert_yaxis()
        if j==1: axes[j].set_yticklabels([])
        else: axes[j].set_ylabel(r'$\mathregular{\it{z}}$'+'/'+r'$\mathregular{\it{z}_{\rm{tot}}}$')
cbaxes = fig.add_axes([0.87, 0.18, 0.02, 0.75]) 
# cticks = [-20, -16, -12, -8, -4, 0]
cbar = fig.colorbar(cmap, cax = cbaxes, ticks = c)
cbar.set_label('Rainwater '+r'$\mathregular{\delta^{18}}$'+'O ('+u'$‰$'+')', rotation=90)
fig.subplots_adjust(left=0.12,right= 0.83 ,bottom=0.17,top = 0.95,wspace=0.15,hspace=0.1)

outdir = 'C:/Users/YK/Desktop/Isotope/'

plt.savefig(outdir+'O_cw-1d.svg', transparent=True)
plt.savefig(outdir+'O_cw-1d.pdf', transparent=True)

plt.show()
plt.clf()
plt.close()