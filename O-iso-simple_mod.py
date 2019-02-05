import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.misc import derivative
from scipy.optimize import fsolve
import time

start=time.clock()

def clear(infileA,infileB):
    a=np.zeros(infileA.shape[0])
    b=np.zeros(infileB.shape[0])
    cnt=0
    cnt2=1
    for n in range(infileA.shape[0]):
        if n==0:
            a[cnt]=infileA[n]
            b[cnt]=infileB[n]
            cnt+=1
            continue
        if infileA[n]==infileA[n-1]:
            a[cnt-1]+=infileA[n]
            b[cnt-1]+=infileB[n]
            cnt2+=1
            continue
        if infileA[n]!=infileA[n-1]:
            if n==1:
                a[cnt]=infileA[n]
                b[cnt]=infileB[n]
                cnt+=1
                continue
            if infileA[n-1]==infileA[n-2]:
                a[cnt-1]=a[cnt-1]/cnt2
                b[cnt-1]=b[cnt-1]/cnt2
                cnt2=1
                a[cnt]=infileA[n]
                b[cnt]=infileB[n]
                cnt+=1
                continue
            a[cnt]=infileA[n]
            b[cnt]=infileB[n]
            cnt+=1
            continue
    a=a[:cnt]
    b=b[:cnt]
    return a,b

def FRdel(sca):
    return 1e3*(sca/((1.0-sca)*0.0020052)-1.0)

def FRr(sca):
    return sca/(1.0-sca)

def Rfr(sca):
    return sca/(1.0+sca)

outfile = 'C:/Users/YK/Desktop/Isotope/wallmann_calc_full.txt'
infiledir = 'C:/Users/YK/Desktop/Earth history/'

fl=open(infiledir + 'fsp.txt')
A=np.loadtxt(fl)
fl.close()
A=np.transpose(A)
x,y=clear(-A[0],A[1])

fsp=interpolate.interp1d(-x,y,kind='linear')

fl=open(infiledir+'18O.txt')
A=np.loadtxt(fl)
fl.close()
A=np.transpose(A)
v,w=clear(-A[0],A[2])

obs18O=interpolate.interp1d(-v,w,kind='linear')

tmin=max(x.min(),v.min())
tmax=min(x.max(),v.max())

N=int(1e3)
t=np.linspace(tmin,tmax,N) #time [Myr]
delt=(tmax-tmin)/(N-1) #time discrete[Myr]


Fw=fsp(-t)*7e18 #[mol Myr^-1]
Falt=fsp(-t)*6e18 #[mol Myr^-1]
Fsp=fsp(-t)*20e18 #[mol Myr^-1]

rw=1.0
rup=1.0
rdeep=0.5
Fnew=rw*Fw #[mol Myr^-1]
Freup=rup*Falt #[mol Myr^-1]
Fredeep=rdeep*Fsp #[mol Myr^-1]

Fm=fsp(-t)*3e18 #[mol Myr^-1]

F18pw=fsp(-t)*2e14
# 18O loss during porewater formation and recycling [mol Myr^-1]
F18ps=2.017e-3*7.8e19*fsp(-t) # formation of sedimentary rock 18O [mol Myr^-1]
F18pup=2.0126e-3*1.2e20*fsp(-t)
# formation of upper crust 18O at spreading zones [mol Myr^-1]
F18pdeep=2.0126e-3*1.6e21*fsp(-t)
# formation of deep crust 18O at spreading zones [mol Myr^-1]
F18m=2.015e-3*Fm
# H218O release through mantle degassing [mol Myr^-1]

F18w=np.zeros(N)
F18alt=np.zeros(N)
F18sp=np.zeros(N)
F18wex=np.zeros(N)
F18upex=np.zeros(N)
F18deepex=np.zeros(N)
F18new=np.zeros(N)
F18reup=np.zeros(N)
F18redeep=np.zeros(N)
F18ms=np.zeros(N)
F18sup=np.zeros(N)
F18sdeep=np.zeros(N)

FRf=np.zeros(N)
FRw=np.zeros(N)
FRup=np.zeros(N)
FRdeep=np.zeros(N)

Rf=np.zeros(N)
Rw=np.zeros(N)
Rup=np.zeros(N)
Rdeep=np.zeros(N)

Sw=np.zeros(N)
Sup=np.zeros(N)
Sdeep=np.zeros(N)

d18ss=np.zeros(N)
d18sw=np.zeros(N)
d18w=np.zeros(N)
d18up=np.zeros(N)
d18deep=np.zeros(N)


fspv=np.zeros(N)

fspv = fsp(-t[:])

M=np.zeros(N) # 
M[0]=1.55e24/18.0 # 

alw=1.02
alup=1.015
aldeep=1.0
FRsi=2.017e-3
FRoc=2.0126e-3
FRm=2.015e-3

kw=1.2e17
kup=1.0e17
kdeep=4.9e17

Rf[0]=0.0020052*(1.0-8.4e-3)
FRf[0]=Rfr(Rf[0])

afunc=lambda x:(F18ps[0]+alw*Rf[0]*Fw[0]/(1.0+alw*Rf[0])+F18pw[0]
                -fsp(-t[0])*kw*(2.017e-3/((1.0-2.017e-3)*Rf[0]*alw)-1.0)
                -x*Fnew[0]
                -x*fsp(-t[0])*7.8e19-x*(1.0-rw)*Fw[0])

bfunc=lambda x:(F18pup[0]+alup*Rf[0]*Falt[0]/(1.0+alup*Rf[0])
                -fsp(-t[0])*kup*(2.0126e-3/((1.0-2.0126e-3)*Rf[0]*alup)-1.0)
                -x*Freup[0]-x*fsp(-t[0])*1.2e20-x*(1.0-rup)*Falt[0])

cfunc=lambda x:(F18pdeep[0]+aldeep*Rf[0]*Fsp[0]/(1.0+aldeep*Rf[0])
                -fsp(-t[0])*kdeep*(2.0126e-3/((1.0-1.0126e-3)*aldeep*Rf[0])-1.0)
                -x*Fredeep[0]-x*fsp(-t[0])*1.6e21-x*(1.0-rdeep)*Fsp[0])

asol=fsolve(afunc,0.5)
bsol=fsolve(bfunc,0.5)
csol=fsolve(cfunc,0.5)

# force initial 18O in sediments, up crust and deep crust
asol[0]=Rfr(0.0020052*(1.0+15.0e-3))
bsol[0]=Rfr(0.0020052*(1.0+8.0e-3))
csol[0]=Rfr(0.0020052*(1.0+5.0e-3))

print FRdel(asol[0]),FRdel(bsol[0]),FRdel(csol[0])

FRw[0]=asol[0]
FRup[0]=bsol[0]
FRdeep[0]=csol[0]

Rw[0]=FRr(FRw[0])
Rup[0]=FRr(FRup[0])
Rdeep[0]=FRr(FRdeep[0])

M18f=np.zeros(N)
M18w=np.zeros(N)
M18up=np.zeros(N)
M18deep=np.zeros(N)

Mw=7e22
Mup=1e22
Mdeep=1.6e23

M18f[0]=FRf[0]*M[0]
M18w[0]=FRw[0]*Mw
M18up[0]=FRup[0]*Mup
M18deep[0]=FRdeep[0]*Mdeep

F18w[0]=alw*Rf[0]*Fw[0]/(1.0+alw*Rf[0])
F18alt[0]=alup*Rf[0]*Falt[0]/(1.0+alup*Rf[0])
F18sp[0]=aldeep*Rf[0]*Fsp[0]/(1.0+aldeep*Rf[0])

F18wex[0]=fsp(-t[0])*kw*(FRr(FRsi)/(alw*Rf[0])-1.0)
F18upex[0]=fsp(-t[0])*kup*(FRr(FRoc)/(alup*Rf[0])-1.0)
F18deepex[0]=fsp(-t[0])*kdeep*(FRr(FRoc)/(aldeep*Rf[0])-1.0)

F18new[0]=FRw[0]*Fnew[0]
F18reup[0]=FRup[0]*Freup[0]
F18redeep[0]=FRdeep[0]*Fredeep[0]

F18ms[0]=FRw[0]*fsp(-t[0])*7.8e19+FRw[0]*(1.0-rw)*Fw[0]
F18sup[0]=FRup[0]*fsp(-t[0])*1.2e20+FRup[0]*(1.0-rup)*Falt[0]
F18sdeep[0]=FRdeep[0]*fsp(-t[0])*1.6e21+FRdeep[0]*(1.0-rdeep)*Fsp[0]

temp=np.zeros(N)
temp[0]=(16.9-4.38*(obs18O(-t[0])-FRdel(FRf[0]))
         +0.10*(obs18O(-t[0])-FRdel(FRf[0]))**2.0)

for n in range(1,N):
#    print n
    M[n]=(M[n-1]+delt*
          (Fm[n-1]+Fnew[n-1]+Freup[n-1]+Fredeep[n-1]
           -Fw[n-1]-Falt[n-1]-Fsp[n-1]))

    M18f[n]=(M18f[n-1]+delt*
          (F18new[n-1]+F18reup[n-1]+F18redeep[n-1]+F18m[n-1]
           +F18wex[n-1]+F18upex[n-1]+F18deepex[n-1]
           -F18w[n-1]-F18alt[n-1]-F18sp[n-1]-F18pw[n-1]))

    M18w[n]=(M18w[n-1]+delt*
             (F18ps[n-1]+F18w[n-1]+F18pw[n-1]-F18wex[n-1]
              -F18new[n-1]-F18ms[n-1]))

    M18up[n]=(M18up[n-1]+delt*
              (F18pup[n-1]+F18alt[n-1]-F18upex[n-1]
               -F18reup[n-1]-F18sup[n-1]))

    M18deep[n]=(M18deep[n-1]+delt*
                (F18pdeep[n-1]+F18sp[n-1]-F18deepex[n-1]
                 -F18redeep[n-1]-F18sdeep[n-1]))

    FRf[n]=M18f[n]/M[n]
    FRw[n]=M18w[n]/Mw
    FRup[n]=M18up[n]/Mup
    FRdeep[n]=M18deep[n]/Mdeep

    Rf[n]=FRf[n]/(1.0-FRf[n])
    Rw[n]=FRw[n]/(1.0-FRw[n])
    Rup[n]=FRup[n]/(1.0-FRup[n])
    Rdeep[n]=FRdeep[n]/(1.0-FRdeep[n])

    F18w[n]=alw*Rf[n]*Fw[n]/(1.0+alw*Rf[n])
    F18alt[n]=alup*Rf[n]*Falt[n]/(1.0+alup*Rf[n])
    F18sp[n]=aldeep*Rf[n]*Fsp[n]/(1.0+aldeep*Rf[n])

    F18wex[n]=fsp(-t[n])*kw*((FRsi/(1.0-FRsi))/(alw*Rf[n])-1.0)
    F18upex[n]=fsp(-t[n])*kup*((FRoc/(1.0-FRoc))/(alup*Rf[n])-1.0)
    F18deepex[n]=fsp(-t[n])*kdeep*((FRoc/(1.0-FRoc))/(aldeep*Rf[n])-1.0)

    F18new[n]=FRw[n]*Fnew[n]
    F18reup[n]=FRup[n]*Freup[n]
    F18redeep[n]=FRdeep[n]*Fredeep[n]

    F18ms[n]=FRw[n]*fsp(-t[n])*7.8e19+FRw[n]*(1.0-rw)*Fw[n]
    F18sup[n]=FRup[n]*fsp(-t[n])*1.2e20+FRup[n]*(1.0-rup)*Falt[n]
    F18sdeep[n]=FRdeep[n]*fsp(-t[n])*1.6e21+FRdeep[n]*(1.0-rdeep)*Fsp[n]

    temp[n]=(16.9-4.38*(obs18O(-t[n])-(1e3*(FRf[n]/((1.0-FRf[n])*0.0020052)-1.0)))
         +0.10*(obs18O(-t[n])-(1e3*(FRf[n]/((1.0-FRf[n])*0.0020052)-1.0)))**2.0)

d18sw[:]=FRdel(FRf[:])
d18w[:]=FRdel(FRw[:])
d18up[:]=FRdel(FRup[:])
d18deep[:]=FRdel(FRdeep[:])
d18ss[:]=(rw*Fw[:]*d18w[:]+rup*Falt[:]*d18up[:]+rup*Fsp[:]*d18deep[:]
          +1e3*((1.0-alw)*Fw[:]+(1.0-alup)*Falt[:]+(1.0-aldeep)*Fsp[:])
          +FRdel(FRm)*Fm[:]
          +(FRdel(FRsi)-1e3*np.log(alw))*kw/2e-3*fsp(-t[:])
          +(FRdel(FRoc)-1e3*np.log(alup))*kup/2e-3*fsp(-t[:])
          +(FRdel(FRoc)-1e3*np.log(aldeep))*kdeep/2e-3*fsp(-t[:])
          -1e3*F18pw[:]/2e-3)/(
              rw*Fw[:]+rup*Falt[:]+rdeep*Fsp[:]
              -((1.0-alw)*Fw[:]+(1.0-alup)*Falt[:]+(1.0-aldeep)*Fsp[:])
              +kw/2e-3*fsp(-t[:])+kup/2e-3*fsp(-t[:])+kdeep/2e-3*fsp(-t[:]))


end=time.clock()
print 'duration (min) =',(end-start)/60.0
fo=open(outfile, 'w')
for n in range(N):
    print>>fo,t[n],d18sw[n],d18ss[n]
fo.close()

print np.amax(M/fspv),np.amin(M/fspv)

plt.subplot(321)
plt.plot(t,M*18.0)
##plt.show()
##plt.clf()

plt.subplot(322)
plt.plot(t,d18sw)
plt.plot(t,d18ss)
plt.plot(t,obs18O(-t))
#;plt.show();plt.clf()

plt.subplot(323)
plt.plot(t,F18wex);plt.plot(t,F18upex);plt.plot(t,F18deepex)
#plt.show();plt.clf()

plt.subplot(324)
plt.plot(t,FRdel(FRw));plt.plot(t,FRdel(FRup))
plt.plot(t,FRdel(FRdeep))#;plt.show();plt.clf()

plt.subplot(325)
plt.plot(t,temp)
#plt.show();plt.clf()

##plt.subplot(326)
##plt.plot(t,(FRsi/(1.0-FRsi))/(alw*Rf));
##plt.plot(t,(FRoc/(1.0-FRoc))/(alup*Rf));
##plt.plot(t,(FRoc/(1.0-FRoc))/(aldeep*Rf));


##plt.subplot(326) # sed
##plt.plot(t,F18ps,label='F18ps');
##plt.plot(t,F18w,label='F18w');
##plt.plot(t,F18pw,label='F18pw');
##plt.plot(t,F18wex,label='F18wex');
##plt.plot(t,F18new,label='F18new');
##plt.plot(t,F18ms,label='F18ms');
##plt.legend();

##plt.subplot(326) # up
##plt.plot(t,F18pup,label='F18pup');
##plt.plot(t,F18alt,label='F18alt');
##plt.plot(t,F18upex,label='F18upex');
##plt.plot(t,F18reup,label='F18reup');
##plt.plot(t,F18sup,label='F18sup');
##plt.legend();

plt.subplot(326) # deep
plt.plot(t,F18pdeep,label='F18pdeep');
plt.plot(t,F18sp,label='F18sp');
plt.plot(t,F18deepex,label='F18deepex');
plt.plot(t,F18redeep,label='F18redeep');
plt.plot(t,F18sdeep,label='F18sdeep');
plt.legend();

plt.tight_layout()
plt.subplots_adjust
plt.show()
plt.clf()

plt.plot(t,F18new,':',label='F18new')
plt.plot(t,F18reup,':',label='F18reup')
plt.plot(t,F18redeep,':',label='F18redeep')
plt.plot(t,F18m,':',label='F18m')
plt.plot(t,F18wex,label='F18wex')
plt.plot(t,F18upex,label='F18upex')
plt.plot(t,F18deepex,label='F18deepex')
plt.plot(t,F18w,label='F18w')
plt.plot(t,F18alt,label='F18alt')
plt.plot(t,F18sp,label='F18sp')
plt.plot(t,F18pw,label='F18pw')
plt.legend()
plt.show()
plt.clf()

plt.plot(t,F18redeep+F18m+F18deepex-F18sp,label='sum-f,rep')
plt.plot(t,F18new+F18reup+F18redeep+F18m
           +F18wex+F18upex+F18deepex
           -F18w-F18alt-F18sp-F18pw,label='sum-f,tot')
plt.plot(t,F18ps+F18w+F18pw-F18wex
              -F18new-F18ms,label='sum-w,tot')
plt.plot(t,F18pup+F18alt-F18upex
               -F18reup-F18sup,label='sum-up,tot')
plt.plot(t,F18pdeep+F18sp-F18deepex
                 -F18redeep-F18sdeep,label='sum-deep,tot')
plt.legend()
plt.show()
plt.clf()
