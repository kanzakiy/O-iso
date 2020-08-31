# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

outdir = '../oiso_output/'

sld_1_d18 = np.loadtxt(outdir + 'o_iso_ha_d_sld_1_prof_d18.txt')
pw_1_d18 = np.loadtxt(outdir + 'o_iso_ha_d_pw_1_prof_d18.txt')
sld_2_d18 = np.loadtxt(outdir + 'o_iso_ha_d_sld_2_prof_d18.txt')
pw_2_d18 = np.loadtxt(outdir + 'o_iso_ha_d_pw_2_prof_d18.txt')
sld_3_d18 = np.loadtxt(outdir + 'o_iso_ha_d_sld_3_prof_d18.txt')
pw_3_d18 = np.loadtxt(outdir + 'o_iso_ha_d_pw_3_prof_d18.txt')

sld_1_d17 = np.loadtxt(outdir + 'o_iso_ha_d_sld_1_prof_d17.txt')
pw_1_d17 = np.loadtxt(outdir + 'o_iso_ha_d_pw_1_prof_d17.txt')
sld_2_d17 = np.loadtxt(outdir + 'o_iso_ha_d_sld_2_prof_d17.txt')
pw_2_d17 = np.loadtxt(outdir + 'o_iso_ha_d_pw_2_prof_d17.txt')
sld_3_d17 = np.loadtxt(outdir + 'o_iso_ha_d_sld_3_prof_d17.txt')
pw_3_d17 = np.loadtxt(outdir + 'o_iso_ha_d_pw_3_prof_d17.txt')

sld_1_D17 = np.loadtxt(outdir + 'o_iso_ha_d_sld_1_prof_capd17.txt')
pw_1_D17 = np.loadtxt(outdir + 'o_iso_ha_d_pw_1_prof_capd17.txt')
sld_2_D17 = np.loadtxt(outdir + 'o_iso_ha_d_sld_2_prof_capd17.txt')
pw_2_D17 = np.loadtxt(outdir + 'o_iso_ha_d_pw_2_prof_capd17.txt')
sld_3_D17 = np.loadtxt(outdir + 'o_iso_ha_d_sld_3_prof_capd17.txt')
pw_3_D17 = np.loadtxt(outdir + 'o_iso_ha_d_pw_3_prof_capd17.txt')

n_lines = 11
c = np.linspace(-20, 0, n_lines)
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['black','red'])
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mymap)
cmap.set_array([])


figsize = (7,6)
fig = plt.figure(figsize=figsize)

nx =2 
ny =3 


axes = [[plt.subplot2grid((ny,nx), (i,j)) for i in range(ny)] for j in range(nx)]
res = [
    [
    [sld_1_d18, pw_1_d18]
    ,[sld_2_d18, pw_2_d18]
    ,[sld_3_d18, pw_3_d18]
    ]
    ,[
    [sld_1_D17, pw_1_D17]
    ,[sld_2_D17, pw_2_D17]
    ,[sld_3_D17, pw_3_D17]
    ]
    ]
ylabels = [
    [  
    'Basalt '+r'$\mathregular{\delta^{18}O}$'+ ' ('+u'$‰$'+')'
    ,'Dike '+r'$\mathregular{\delta^{18}O}$'+ ' ('+u'$‰$'+')'
    ,'Gabbro '+r'$\mathregular{\delta^{18}O}$'+ ' ('+u'$‰$'+')'
    ]
    ,[  
    'Basalt '+r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')'
    ,'Dike '+r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')'
    ,'Gabbro '+r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')'
    ]
    ]
    
ylabels = [
    'Basalt'
    ,'Dike'
    ,'Gabbro'
    ]
titles = [
    r'$\mathregular{\delta^{18}O}$'+ ' ('+u'$‰$'+')'
    ,r"$\mathregular{\Delta'^{17}O}$"+ ' ('+u'$‰$'+')'
    ]

for j in range(nx):
    for i in range(ny):
        for k in range(n_lines):
            axes[j][i].plot(res[j][i][0][0:,1],res[j][i][0][0:,2+k],c=cmap.to_rgba(c[k]))
            axes[j][i].plot(res[j][i][1][0:,1],res[j][i][1][0:,2+k],'--',c=cmap.to_rgba(c[k]))
            axes[j][i].set_xlim(0,1)
            if j==0: axes[j][i].set_ylabel(ylabels[i])
            # else: axes[j][i].set_yticklabels([])
            if i==0 and k==0: 
                axes[j][i].text(0.5,1.1,titles[j],horizontalalignment='center',transform=axes[j][i].transAxes)
            if i!=ny-1: axes[j][i].set_xticklabels([])
            else: axes[j][i].set_xlabel(r'$\mathregular{\it{x}}$'+'/'+r'$\mathregular{\it{x}_{\rm{lim}}}$')
cbaxes = fig.add_axes([0.15, 0.1, 0.75, 0.02]) 
# cticks = [-20, -16, -12, -8, -4, 0]
cbar = fig.colorbar(cmap, cax = cbaxes, ticks = c, orientation = 'horizontal')
cbar.set_label('Seawater '+r'$\mathregular{\delta^{18}}$'+'O ('+u'$‰$'+')', rotation=0)
fig.align_labels(axes)
fig.subplots_adjust(left=0.1,right= 0.95 ,bottom=0.25,top = 0.92,wspace=0.2,hspace=0.1)
outdir = 'C:/Users/YK/Desktop/Isotope/'

plt.savefig(outdir+'O_ha-2d.svg', transparent=True)
plt.savefig(outdir+'O_ha-2d.pdf', transparent=True)

plt.show()
plt.clf()
plt.close()