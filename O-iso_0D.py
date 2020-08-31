# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

def rock_R(psi,kappa,alpha,Rsw,Rmc,beta_kin):
    Rr = (alpha*Rsw+alpha*Rmc/psi+Rmc/(kappa**beta_kin))/(1.+alpha/psi+1./(kappa**beta_kin))
    return Rr
    
def porewater_R(psi,kappa,alpha,Rsw,Rmc,beta_kin):
    Rp = (Rsw+Rmc/psi+Rsw/(kappa**beta_kin))/(1.+alpha/psi+1./(kappa**beta_kin))
    return Rp
    
def distance_eq(Rr, Rp, alpha):
    omega = Rr/alpha/Rp
    return omega

def buff_calc(psi, kappa, alpha, Rstd, phi,w,L,h,rhos,ms,beta_kin):
    buff = Rstd*1e-3*alpha*(1.-phi)*w*L*h*rhos*ms/(1.+alpha/psi+1./(kappa**beta_kin))
    return buff 

def d2r(d, rstd):
    r = (d/1e3 + 1.)*rstd
    return r

def r2d(r, rstd):
    d = (r/rstd - 1.)*1e3
    return d

def d_prime2r(d_prime, rstd):
    r = rstd*np.exp(d/1e3)
    return r

def r2d_prime(r, rstd):
    d_prime = 1e3*np.log(r/rstd)
    return d_prime

def d_prime2d(d_prime, rstd):
    r = d_prime2r(d_prime, rstd)
    d = r2d(r, rstd)
    return d

def d2d_prime(d, rstd):
    r = d2r(d_prime, rstd)
    d_prime = r2d_prime(r, rstd)
    return d_prime
    
# x --> psi
# y --> kappa

beta_eq = (1/16.-1/17.)/(1/16.-1/18.)  # Theoretial one by Young et al. (2002) 
# beta_eq = 0.51 
fact_kin = 1. # 1 to infinity assuming 0 < beta_kin < beta_eq
beta_kin = beta_eq/fact_kin

s_ref =  0.5305 #  Pack and Herwartz (2014)
c_ref =  0. #  Pack and Herwartz (2014)

dd = 0.2

dx, dy = dd, dd

x = np.arange(-5,5+dx/2,dx)
y = np.arange(-5,5+dy/2,dy)

X, Y = np.meshgrid(x, y)

Ts,  dsw, alpha, phi, kref,      E, tau,  h,   w,  L,  rhos, rhof, ms,  mf  = [
273, 0,   2,    0.1,  10**-8.5, 50, 1e6, 2e3,3e-2, 1e8, 3e3, 1e3, 31.3, 55.6
]   
alpha = np.exp(alpha/1e3)
alpha_17 = alpha**beta_eq

Rstd = 2.0052e-3
Rstd_17 = 3.799e-4 

Rsw = d2r(dsw,Rstd)
Rmc = d2r(5.7,Rstd)

Rsw_17 = d2r(dsw,Rstd_17)
Rmc_17 = d2r(2.86,Rstd_17)  # Pack and Herwartz (2014)
# Rmc_17 = d2r(3.,Rstd_17)  # Pack and Herwartz (2014)

Rr = rock_R(10.**X,10**Y,alpha,Rsw, Rmc,1.)
Rp = porewater_R(10.**X,10**Y,alpha,Rsw, Rmc,1.)
omega = distance_eq(Rr,Rp,alpha)
buff = buff_calc(10.**X, 10.**Y, alpha, Rstd, phi,w,L,h,rhos,ms,1.)

kin_fact = beta_kin
kin_fact = 1.  # this assume kinetics of 17O/16O exchange is not different from that of 18O/16O
Rr_17 = rock_R(10.**X,10**Y,alpha_17,Rsw_17, Rmc_17,kin_fact)
Rp_17 = porewater_R(10.**X,10**Y,alpha_17,Rsw_17, Rmc_17,kin_fact)
omega_17 = distance_eq(Rr_17,Rp_17,alpha_17)
buff_17 = buff_calc(10.**X, 10.**Y, alpha_17, Rstd_17, phi,w,L,h,rhos,ms,kin_fact)

dr = r2d(Rr,Rstd) 
dp = r2d(Rp,Rstd) 

dr_prime = r2d_prime(Rr,Rstd) 
dp_prime = r2d_prime(Rp,Rstd) 

dr_17 = r2d(Rr_17,Rstd_17) 
dp_17 = r2d(Rp_17,Rstd_17) 

dr_17_prime = r2d_prime(Rr_17,Rstd_17) 
dp_17_prime = r2d_prime(Rp_17,Rstd_17) 

capD17r = dr_17_prime - (s_ref*dr_prime + c_ref)
capD17p = dp_17_prime - (s_ref*dp_prime + c_ref)

cmap = 'gnuplot'

obj_plot = '18'
obj_plot = '17'
    
figsize = (7,5)
fig = plt.figure(figsize=figsize)

nx =2 
ny =2 

axes = [[plt.subplot2grid((ny,nx), (j,i)) for i in range(nx)] for j in range(ny)]

if obj_plot =='18':
    labels = [
        'Porewater '+r'$\mathregular{\delta ^{18}}$'+'O'+ ' ('+u'$‰$'+')'
        ,'Solid rock '+r'$\mathregular{\delta ^{18}}$'+'O'+ ' ('+u'$‰$'+')'
        ,r'$\mathregular{\Omega}$'
        ,'Buffering capacity\n(10'+r'$\mathregular{^{9}}$'+' mol yr'+r'$\mathregular{^{-1}}$'+' '+u'$‰$'+r'$\mathregular{^{-1}}$'+')'
        ]
    
    Zs = [
        dp
        ,dr
        ,omega
        ,buff/1e9
        ]
elif obj_plot =='17':
    labels = [
        'Porewater '+r"$\mathregular{\delta' ^{17}}$"+'O'+ ' ('+u'$‰$'+')'
        ,'Solid rock '+r"$\mathregular{\delta' ^{17}}$"+'O'+ ' ('+u'$‰$'+')'
        ,'Porewater '+r"$\mathregular{\Delta' ^{17}}$"+'O'+ ' ('+u'$‰$'+')'
        ,'Solid rock '+r"$\mathregular{\Delta' ^{17}}$"+'O'+ ' ('+u'$‰$'+')'
        ]
    
    Zs = [
        dp_17_prime
        ,dr_17_prime
        ,capD17p
        ,capD17r
        ]
    
for k in range(nx*ny):
    i = k%nx
    j = (k-i)/nx
    cf = axes[j][i].pcolormesh(X, Y, Zs[k], cmap=cmap)
    pp=fig.colorbar (cf,ax= axes[j][i]) 
    pp.set_label(labels[k])
    
    cs = axes[j][i].contour(X, Y, Zs[k], colors = 'w')
    axes[j][i].clabel(cs, inline=1, fontsize=10)

    axes[j][i].set_xlabel(r'$\mathregular{\log \it{\Psi}}$')
    axes[j][i].set_ylabel(r'$\mathregular{\log \kappa}$')
    
    axes[j][i].set_xticks([-4,-2,0,2,4])
    axes[j][i].set_yticks([-4,-2,0,2,4])

fig.tight_layout()

outdir = 'C:/Users/YK/Desktop/Isotope/'

plt.savefig(outdir+obj_plot+'O-0d.svg', transparent=True)
plt.savefig(outdir+obj_plot+'O-0d.pdf', transparent=True)
    
plt.show()
plt.clf()
plt.close()