module o_iso_ha_mod
use o_iso_ref_mod
implicit none 
! call o_iso_ha_sense_dsw_prof
! call o_iso_ha_sense_dsw_Phan
! endprogram o_iso_ha_main
contains
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************

subroutine o_iso_ha_sense_dsw_prof

implicit none

integer, parameter :: nz = 200
integer, parameter :: nps1 = 11  ! dsw

real(kind=8) :: L1 = 0.5d0 ! km
real(kind=8) :: L2 = 1.5d0     
real(kind=8) :: L3 = 4.0d0  
real(kind=8) :: fsp = 3.0d0 ! km^2 yr^-1
real(kind=8) :: T1 = 5.0d0 ! deg C
real(kind=8) :: T2 = 200.0d0
real(kind=8) :: T3 = 350.0d0  
real(kind=8) :: kref = 10.0d0**(-8.50d0) ! yr^-1
real(kind=8) :: E1 = 50.0d0 ! kj mol^-1
real(kind=8) :: E2 = 50.0d0
real(kind=8) :: E3 = 50.0d0
real(kind=8) :: tres = 1.0d6! yr
real(kind=8) :: alfa1 = exp(27.0d0/1d3)
real(kind=8) :: alfa2 = exp(7.0d0/1d3)
real(kind=8) :: alfa3 = exp(2.5d0/1d3)
real(kind=8) :: f1 = 1.0d16  ! kg yr^-1
real(kind=8) :: f2 = 1.0d13
real(kind=8) :: f3 = 1.0d13
real(kind=8) :: p1 = 0.25d0
real(kind=8) :: p2 = 0.05d0
real(kind=8) :: p3 = 0.02d0
real(kind=8) :: a = -2.0d0
real(kind=8) :: beta = 0.876002291d0 ! anddesite (Zhao and Zheng, 2003)
real(kind=8) :: lambda = 0.5305d0 
real(kind=8) :: gamma = 0d0

real(kind=8) k1, k2, k3
real(kind=8) L
real(kind=8) r1av,r2av,r3av
real(kind=8) rw1av,rw2av,rw3av
real(kind=8) jocnss
real(kind=8) avx,avy,sumx2,sumy2,sumxy,slp,itcpt,crr
real(kind=8),dimension(2) :: dmc, rsmow, dsw,fsw,fmc
real(kind=8),dimension(2,nps1) :: dswlist,jocnsslist
real(kind=8),dimension(nz) :: z
real(kind=8),dimension(2,nps1,nz) :: r1,r2,r3,rw1,rw2,rw3
real(kind=8),dimension(3) :: T_in,F_in,poro_in,L_in,beta_eq,beta_kin
real(kind=8),dimension(2,3) :: alpha_in,k_in

integer iz,is1,iiso,iiso2

real(kind=8) k_arrhenius,alpha_zhao03,r2f,d2r,r2d,f2r,t2beta_eq,mwl_Luz,dp2d,r2dp
character (50) beta_eq_model
character (2) isochr
!---------------------------

rsmow(1) = rsmow_18
rsmow(2) = rsmow_17 
dmc(1) = d18_mc
dmc(2) = d17_mc

k1 = k_arrhenius(kref,E1,T1) 
k2 = k_arrhenius(kref,E2,T2) 
k3 = k_arrhenius(kref,E3,T3) 

alfa1 = alpha_zhao03(T1,beta)
alfa2 = alpha_zhao03(T2,beta)
alfa3 = alpha_zhao03(T3,beta)

T_in = (/T1,T2,T3/)
k_in(1,:) = (/k1,k2,k3/)
poro_in = (/p1,p2,p3/)
alpha_in(1,:) = (/alfa1,alfa2,alfa3/)
F_in = (/f1,f2,f3/)
L_in = (/L1,L2,L3/)

beta_eq_model = 'sharp16'
beta_eq(:) = t2beta_eq(T_in(:)+273d0, beta_eq_model)
beta_kin(:) = 1d0
if (switch_KIE) beta_kin(:) = beta_eq(:)

k_in(2,:) = k_in(1,:)**beta_kin(:)
alpha_in(2,:) = alpha_in(1,:)**beta_eq(:)

do is1 = 1,nps1
    dsw(1) = -20.0d0 + 20.0d0*(is1-1.0d0)/(nps1-1.0d0)
    dsw(2) = dp2d(mwl_Luz(dsw(1),rsmow(1)),rsmow(2))! changing 18O to 17O 
    dswlist(:,is1) = dsw(:)
    fsw(1) = r2f(d2r(dsw(1),rsmow(1)),d2r(dsw(2),rsmow(2)))
    fsw(2) = r2f(d2r(dsw(2),rsmow(2)),d2r(dsw(1),rsmow(1)))
    fmc(1) = r2f(d2r(dmc(1),rsmow(1)),d2r(dmc(2),rsmow(2)))
    fmc(2) = r2f(d2r(dmc(2),rsmow(2)),d2r(dmc(1),rsmow(1)))

    do iiso= 1, 2        
        print *, is1, dsw(iiso)
        
        call o_iso_ha( &
            & fsw(iiso),nz,T_in,L_in,fsp,F_in,poro_in,tres,fmc(iiso),alpha_in(iiso,:),k_in(iiso,:),a  &! input 
            & ,r1(iiso,is1,:),r2(iiso,is1,:),r3(iiso,is1,:),rw1(iiso,is1,:),rw2(iiso,is1,:),rw3(iiso,is1,:)  &! output
            & ,r1av,r2av,r3av,rw1av,rw2av,rw3av,jocnss,L,z              &! output
            & ) 
        
        jocnsslist(iiso,is1) = jocnss
    
    enddo 

end do

do iiso = 1,2
    if (iiso ==1) then 
        isochr ='18'
        iiso2 = 2
    elseif (iiso ==2) then 
        isochr ='17'
        iiso2 = 1
    endif 

    open(unit=15,file = '../oiso_output/o_iso_ha_d_sld_1_prof_d'//isochr//'.txt',status = 'replace')
    open(unit=16,file = '../oiso_output/o_iso_ha_d_sld_2_prof_d'//isochr//'.txt',status = 'replace')
    open(unit=17,file = '../oiso_output/o_iso_ha_d_sld_3_prof_d'//isochr//'.txt',status = 'replace')
    open(unit=18,file = '../oiso_output/o_iso_ha_d_pw_1_prof_d'//isochr//'.txt',status = 'replace')
    open(unit=19,file = '../oiso_output/o_iso_ha_d_pw_2_prof_d'//isochr//'.txt',status = 'replace')
    open(unit=20,file = '../oiso_output/o_iso_ha_d_pw_3_prof_d'//isochr//'.txt',status = 'replace')
    do iz = 1, nz
        write(15,*) z(iz),z(iz)/L, (r2d(f2r(r1(iiso,is1,iz),r1(iiso2,is1,iz)),rsmow(iiso)),is1=1,nps1)
        write(16,*) z(iz),z(iz)/L, (r2d(f2r(r2(iiso,is1,iz),r2(iiso2,is1,iz)),rsmow(iiso)),is1=1,nps1)
        write(17,*) z(iz),z(iz)/L, (r2d(f2r(r3(iiso,is1,iz),r3(iiso2,is1,iz)),rsmow(iiso)),is1=1,nps1)
        write(18,*) z(iz),z(iz)/L, (r2d(f2r(rw1(iiso,is1,iz),rw1(iiso2,is1,iz)),rsmow(iiso)),is1=1,nps1)
        write(19,*) z(iz),z(iz)/L, (r2d(f2r(rw2(iiso,is1,iz),rw2(iiso2,is1,iz)),rsmow(iiso)),is1=1,nps1)
        write(20,*) z(iz),z(iz)/L, (r2d(f2r(rw3(iiso,is1,iz),rw3(iiso2,is1,iz)),rsmow(iiso)),is1=1,nps1)
    enddo
    close(15);close(16);close(17);close(18);close(19);close(20)


    avx = sum(dswlist(iiso,:))/nps1
    avy = sum(jocnsslist(iiso,:))/nps1
    sumx2 = sum((dswlist(iiso,:)-avx)**2.0d0)
    sumy2 = sum((jocnsslist(iiso,:)-avy)**2.0d0)
    sumxy = sum((dswlist(iiso,:)-avx)*(jocnsslist(iiso,:)-avy))
    slp = sumxy/sumx2
    itcpt = avy - slp*avx
    crr = sumxy**2.0d0/sumx2/sumy2
    
    print *, iiso,slp,itcpt,-itcpt/slp,crr
    
    if (crr < 0.999d0) then 
        print *, 'error', crr
        pause
    end if
enddo

isochr ='17'
iiso = 2; iiso2 = 1
open(unit=15,file = '../oiso_output/o_iso_ha_d_sld_1_prof_capd'//isochr//'.txt',status = 'replace')
open(unit=16,file = '../oiso_output/o_iso_ha_d_sld_2_prof_capd'//isochr//'.txt',status = 'replace')
open(unit=17,file = '../oiso_output/o_iso_ha_d_sld_3_prof_capd'//isochr//'.txt',status = 'replace')
open(unit=18,file = '../oiso_output/o_iso_ha_d_pw_1_prof_capd'//isochr//'.txt',status = 'replace')
open(unit=19,file = '../oiso_output/o_iso_ha_d_pw_2_prof_capd'//isochr//'.txt',status = 'replace')
open(unit=20,file = '../oiso_output/o_iso_ha_d_pw_3_prof_capd'//isochr//'.txt',status = 'replace')
do iz = 1, nz
    write(15,*) z(iz),z(iz)/L, (r2dp(f2r(r1(iiso,is1,iz),r1(iiso2,is1,iz)),rsmow(iiso))  &
        & -(lambda*r2dp(f2r(r1(iiso2,is1,iz),r1(iiso,is1,iz)),rsmow(iiso2)) + gamma) ,is1=1,nps1)
    write(16,*) z(iz),z(iz)/L, (r2dp(f2r(r2(iiso,is1,iz),r2(iiso2,is1,iz)),rsmow(iiso))&
        & -(lambda*r2dp(f2r(r2(iiso2,is1,iz),r2(iiso,is1,iz)),rsmow(iiso2)) + gamma) ,is1=1,nps1)
    write(17,*) z(iz),z(iz)/L, (r2dp(f2r(r3(iiso,is1,iz),r3(iiso2,is1,iz)),rsmow(iiso))&
        & -(lambda*r2dp(f2r(r3(iiso2,is1,iz),r3(iiso,is1,iz)),rsmow(iiso2)) + gamma) ,is1=1,nps1)
    write(18,*) z(iz),z(iz)/L, (r2dp(f2r(rw1(iiso,is1,iz),rw1(iiso2,is1,iz)),rsmow(iiso))&
        & -(lambda*r2dp(f2r(rw1(iiso2,is1,iz),rw1(iiso,is1,iz)),rsmow(iiso2)) + gamma) ,is1=1,nps1)
    write(19,*) z(iz),z(iz)/L, (r2dp(f2r(rw2(iiso,is1,iz),rw2(iiso2,is1,iz)),rsmow(iiso))&
        & -(lambda*r2dp(f2r(rw2(iiso2,is1,iz),rw2(iiso,is1,iz)),rsmow(iiso2)) + gamma) ,is1=1,nps1)
    write(20,*) z(iz),z(iz)/L, (r2dp(f2r(rw3(iiso,is1,iz),rw3(iiso2,is1,iz)),rsmow(iiso))&
        & -(lambda*r2dp(f2r(rw3(iiso2,is1,iz),rw3(iiso,is1,iz)),rsmow(iiso2)) + gamma) ,is1=1,nps1)
enddo
close(15);close(16);close(17);close(18);close(19);close(20)

endsubroutine o_iso_ha_sense_dsw_prof

!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************

subroutine o_iso_ha_sense_dsw_Phan

implicit none

integer, parameter :: nz = 200
integer, parameter :: nps1 = 11  ! o18wi
integer, parameter :: nry = 58, nrx = 25

real(kind=8) :: L1 = 0.5d0 ! km
real(kind=8) :: L2 = 1.5d0     
real(kind=8) :: L3 = 4.0d0  
real(kind=8) :: fsp = 3.0d0 ! km^2 yr^-1
real(kind=8) :: T1 = 5.0d0 ! deg C
real(kind=8) :: T2 = 200.0d0
real(kind=8) :: T3 = 350.0d0  
real(kind=8) :: kref = 10.0d0**(-8.50d0) ! yr^-1
real(kind=8) :: E1 = 50.0d0 ! kj mol^-1
real(kind=8) :: E2 = 50.0d0
real(kind=8) :: E3 = 50.0d0
real(kind=8) :: tres = 1.0d6! yr
real(kind=8) :: alfa1 = exp(27.0d0/1d3)
real(kind=8) :: alfa2 = exp(7.0d0/1d3)
real(kind=8) :: alfa3 = exp(2.5d0/1d3)
real(kind=8) :: f1 = 1.0d16  ! kg yr^-1
real(kind=8) :: f2 = 1.0d13
real(kind=8) :: f3 = 1.0d13
real(kind=8) :: p1 = 0.25d0
real(kind=8) :: p2 = 0.05d0
real(kind=8) :: p3 = 0.02d0
real(kind=8) :: a = -2.0d0
real(kind=8) :: beta = 0.876002291d0 ! anddesite (Zhao and Zheng, 2003)

real(kind=8) k1, k2, k3
real(kind=8) L
real(kind=8) r1av,r2av,r3av
real(kind=8) rw1av,rw2av,rw3av
real(kind=8) jocnss
real(kind=8),dimension(nz) :: r1,r2,r3,rw1,rw2,rw3,z
real(kind=8),dimension(nrx,nry) :: mry
real(kind=8),dimension(nry) :: fl,fa,fsr,rco2,temp,tp
real(kind=8),dimension(3) :: T_in,F_in,poro_in,L_in,beta_eq,beta_kin
real(kind=8),dimension(2,3) :: alpha_in,k_in
real(kind=8),dimension(2) :: fsw,fmc,dsw,rsmow,dmc
real(kind=8),dimension(2) :: avx,avy,sumx2,sumy2,sumxy,slp,itcpt,crr
real(kind=8),dimension(2,nps1) :: dswlist,jocnsslist

integer iry,is1,iiso

real(kind=8) k_arrhenius,alpha_zhao03,d2r,r2f,mwl_Luz,dp2d,t2beta_eq
character (50) beta_eq_model
!---------------------------

open(30, file ='./Royer2014.txt',status = 'old')
read(30,*) mry
close(30)

fl(:) = mry(6,:)   !  carbonate area (not used)
fa(:) = mry(7,:)    ! continental are (not used)
fsr(:) = mry(16,:) ! spreading rate
rco2(:) = mry(18,:) ! co2
temp(:) = mry(19,:) ! surface temperature 
tp(:) = mry(1,:)  ! age

rsmow(1) = rsmow_18
rsmow(2) = rsmow_17 
dmc(1) = d18_mc
dmc(2) = d17_mc

open (57,file='../oiso_output/GS-oca.txt',status='replace')

do iry = 1, nry

    kref = 10.0d0**(-8.50d0)!*fsr(iry)!*rco2(iry)**0.50d0
    fsp = 3.0d0*fsr(iry)
    tres = 1d6/fsr(iry)

    ! T1 = 5.0d0 + 3.0d0*(temp(iry) -temp(nry))
    T1 = 5.0d0 + 1.9d0*(temp(iry) -temp(nry))

    k1 = k_arrhenius(kref,E1,T1) 
    k2 = k_arrhenius(kref,E2,T2) 
    k3 = k_arrhenius(kref,E3,T3) 

    f1=1d16*(fsp/3.0d0)!**0.50d0
    f2=1d13*(fsp/3.0d0)!**0.50d0
    f3=1d13*(fsp/3.0d0)!**0.50d0

    alfa1 = alpha_zhao03(T1,beta)
    alfa2 = alpha_zhao03(T2,beta)
    alfa3 = alpha_zhao03(T3,beta)

    T_in = (/T1,T2,T3/)
    poro_in = (/p1,p2,p3/)
    F_in = (/f1,f2,f3/)
    L_in = (/L1,L2,L3/)
    
    k_in(1,:) = (/k1,k2,k3/)
    alpha_in(1,:) = (/alfa1,alfa2,alfa3/)
    
    beta_eq_model = 'sharp16'
    beta_eq(:) = t2beta_eq(T_in(:)+273d0, beta_eq_model)
    beta_kin(:) = 1d0
    if (switch_KIE) beta_kin(:) = beta_eq(:)
    
    k_in(2,:) = k_in(1,:)**beta_kin(:)
    alpha_in(2,:) = alpha_in(1,:)**beta_eq(:)

    do is1 = nps1,1,-1
        dsw(1) = -1.0d0*(1.0*is1-1.0d0)/(nps1 - 1.0d0)*20.0d0
        dsw(2) = dp2d(mwl_Luz(dsw(1),rsmow(1)),rsmow(2))! changing 18O to 17O 
        dswlist(:,is1) = dsw(:)
        fsw(1) = r2f(d2r(dsw(1),rsmow(1)),d2r(dsw(2),rsmow(2)))
        fsw(2) = r2f(d2r(dsw(2),rsmow(2)),d2r(dsw(1),rsmow(1)))
        fmc(1) = r2f(d2r(dmc(1),rsmow(1)),d2r(dmc(2),rsmow(2)))
        fmc(2) = r2f(d2r(dmc(2),rsmow(2)),d2r(dmc(1),rsmow(1)))
        

        do iiso= 1, 2        
            print *, is1, dsw(iiso)
        
            call o_iso_ha( &
                & fsw(iiso),nz,T_in,L_in,fsp,F_in,poro_in,tres,fmc(iiso),alpha_in(iiso,:),k_in(iiso,:),a  &! input 
                & ,r1,r2,r3,rw1,rw2,rw3   &! output
                & ,r1av,r2av,r3av,rw1av,rw2av,rw3av,jocnss,L,z              &! output
                & ) 
            
            jocnsslist(iiso,is1) = jocnss
            print*,jocnss
        enddo 

    end do
    
    do iiso = 1,2
        avx(iiso) = sum(dswlist(iiso,:))/nps1
        avy(iiso) = sum(jocnsslist(iiso,:))/nps1
        sumx2(iiso) = sum((dswlist(iiso,:)-avx(iiso))**2.0d0)
        sumy2(iiso) = sum((jocnsslist(iiso,:)-avy(iiso))**2.0d0)
        sumxy(iiso) = sum((dswlist(iiso,:)-avx(iiso))*(jocnsslist(iiso,:)-avy(iiso)))
        slp(iiso) = sumxy(iiso)/sumx2(iiso)
        itcpt(iiso) = avy(iiso) - slp(iiso)*avx(iiso)
        crr(iiso) = sumxy(iiso)**2.0d0/sumx2(iiso)/sumy2(iiso)

        if (crr(iiso) < 0.999d0) then 
            print *, 'error', crr(iiso), iiso
            pause
        end if
    enddo

    write(57,*) slp(1),itcpt(1),slp(2),itcpt(2)

end do

close(57)

endsubroutine o_iso_ha_sense_dsw_Phan

!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************

subroutine o_iso_ha_drdp_Phan

implicit none

integer, parameter :: nz = 200
integer, parameter :: nps1 = 2  ! o18wi
integer, parameter :: nry = 58, nrx = 25

real(kind=8) :: L1 = 0.5d0 ! km
real(kind=8) :: L2 = 1.5d0     
real(kind=8) :: L3 = 4.0d0  
real(kind=8) :: fsp = 3.0d0 ! km^2 yr^-1
real(kind=8) :: T1 = 5.0d0 ! deg C
real(kind=8) :: T2 = 200.0d0
real(kind=8) :: T3 = 350.0d0  
real(kind=8) :: kref = 10.0d0**(-8.50d0) ! yr^-1
real(kind=8) :: E1 = 50.0d0 ! kj mol^-1
real(kind=8) :: E2 = 50.0d0
real(kind=8) :: E3 = 50.0d0
real(kind=8) :: tres = 1.0d6! yr
real(kind=8) :: alfa1 = exp(27.0d0/1d3)
real(kind=8) :: alfa2 = exp(7.0d0/1d3)
real(kind=8) :: alfa3 = exp(2.5d0/1d3)
real(kind=8) :: f1 = 1.0d16  ! kg yr^-1
real(kind=8) :: f2 = 1.0d13
real(kind=8) :: f3 = 1.0d13
real(kind=8) :: p1 = 0.25d0
real(kind=8) :: p2 = 0.05d0
real(kind=8) :: p3 = 0.02d0
real(kind=8) :: a = -2.0d0
real(kind=8) :: beta = 0.876002291d0 ! anddesite (Zhao and Zheng, 2003)
real(kind=8) :: lambda = 0.5305d0
real(kind=8) :: gamma = 0d0

real(kind=8) k1, k2, k3
real(kind=8) L
real(kind=8) jocnss
real(kind=8),dimension(nz) :: z
real(kind=8),dimension(2,nps1,nz) :: r1,r2,r3,rw1,rw2,rw3
real(kind=8),dimension(nrx,nry) :: mry
real(kind=8),dimension(nry) :: fl,fa,fsr,rco2,temp,tp
real(kind=8),dimension(3) :: T_in,F_in,poro_in,L_in,beta_eq,beta_kin
real(kind=8),dimension(2,3) :: alpha_in,k_in
real(kind=8),dimension(2) :: fsw,fmc,dsw,rsmow,dmc
real(kind=8),dimension(2,nps1) :: r1av,r2av,r3av
real(kind=8),dimension(2,nps1) :: rw1av,rw2av,rw3av

integer iry,is1,iiso

real(kind=8) k_arrhenius,alpha_zhao03,d2r,r2f,mwl_Luz,dp2d,t2beta_eq,f2r,r2d,r2dp
character (50) beta_eq_model
real(kind=8) d18o_max(4,nry),d18o_min(4,nry) 
!---------------------------

open(30, file ='./Royer2014.txt',status = 'old')
read(30,*) mry
close(30)

open(30, file = '../oiso_output/d18sw_dyn_max.txt',status = 'old')
read(30,*) d18o_max
close(30)

open(30, file = '../oiso_output/d18sw_dyn_min.txt',status = 'old')
read(30,*) d18o_min
close(30)

fl(:) = mry(6,:)   !  carbonate area (not used)
fa(:) = mry(7,:)    ! continental are (not used)
fsr(:) = mry(16,:) ! spreading rate
rco2(:) = mry(18,:) ! co2
temp(:) = mry(19,:) ! surface temperature 
tp(:) = mry(1,:)  ! age

rsmow(1) = rsmow_18
rsmow(2) = rsmow_17 
dmc(1) = d18_mc
dmc(2) = d17_mc

open (57,file='../oiso_output/d18rp_ha_1_dyn.txt',status='replace')
open (58,file='../oiso_output/d17rp_ha_1_dyn.txt',status='replace')
open (59,file='../oiso_output/capd17rp_ha_1_dyn.txt',status='replace')

open (60,file='../oiso_output/d18rp_ha_2_dyn.txt',status='replace')
open (61,file='../oiso_output/d17rp_ha_2_dyn.txt',status='replace')
open (62,file='../oiso_output/capd17rp_ha_2_dyn.txt',status='replace')

open (63,file='../oiso_output/d18rp_ha_3_dyn.txt',status='replace')
open (64,file='../oiso_output/d17rp_ha_3_dyn.txt',status='replace')
open (65,file='../oiso_output/capd17rp_ha_3_dyn.txt',status='replace')

do iry = 1, nry

    kref = 10.0d0**(-8.50d0)!*fsr(iry)!*rco2(iry)**0.50d0
    fsp = 3.0d0*fsr(iry)
    tres = 1d6/fsr(iry)

    ! T1 = 5.0d0 + 3.0d0*(temp(iry) -temp(nry))
    T1 = 5.0d0 + 1.9d0*(temp(iry) -temp(nry))

    k1 = k_arrhenius(kref,E1,T1) 
    k2 = k_arrhenius(kref,E2,T2) 
    k3 = k_arrhenius(kref,E3,T3) 

    f1=1d16*(fsp/3.0d0)!**0.50d0
    f2=1d13*(fsp/3.0d0)!**0.50d0
    f3=1d13*(fsp/3.0d0)!**0.50d0

    alfa1 = alpha_zhao03(T1,beta)
    alfa2 = alpha_zhao03(T2,beta)
    alfa3 = alpha_zhao03(T3,beta)

    T_in = (/T1,T2,T3/)
    poro_in = (/p1,p2,p3/)
    F_in = (/f1,f2,f3/)
    L_in = (/L1,L2,L3/)
    
    k_in(1,:) = (/k1,k2,k3/)
    alpha_in(1,:) = (/alfa1,alfa2,alfa3/)
    
    beta_eq_model = 'sharp16'
    beta_eq(:) = t2beta_eq(T_in(:)+273d0, beta_eq_model)
    beta_kin(:) = 1d0
    if (switch_KIE) beta_kin(:) = beta_eq(:)
    
    k_in(2,:) = k_in(1,:)**beta_kin(:)
    alpha_in(2,:) = alpha_in(1,:)**beta_eq(:)

    do is1 = nps1,1,-1
        if (is1 == 1) dsw(1) = maxval(d18o_max(2:,iry))
        if (is1 == nps1) dsw(1) = minval(d18o_min(2:,iry)) 
        dsw(2) = dp2d(mwl_Luz(dsw(1),rsmow(1)),rsmow(2))! changing 18O to 17O 
        fsw(1) = r2f(d2r(dsw(1),rsmow(1)),d2r(dsw(2),rsmow(2)))
        fsw(2) = r2f(d2r(dsw(2),rsmow(2)),d2r(dsw(1),rsmow(1)))
        fmc(1) = r2f(d2r(dmc(1),rsmow(1)),d2r(dmc(2),rsmow(2)))
        fmc(2) = r2f(d2r(dmc(2),rsmow(2)),d2r(dmc(1),rsmow(1)))
        

        do iiso= 1, 2        
            print *, is1, dsw(iiso)
        
            call o_iso_ha( &
                & fsw(iiso),nz,T_in,L_in,fsp,F_in,poro_in,tres,fmc(iiso),alpha_in(iiso,:),k_in(iiso,:),a  &! input 
                & ,r1(iiso,is1,:),r2(iiso,is1,:),r3(iiso,is1,:),rw1(iiso,is1,:),rw2(iiso,is1,:),rw3(iiso,is1,:)   &! output
                & ,r1av(iiso,is1),r2av(iiso,is1),r3av(iiso,is1),rw1av(iiso,is1),rw2av(iiso,is1),rw3av(iiso,is1),jocnss,L,z &! output
                & ) 
        enddo 

    end do

    write(57,*) tp(iry), r2d(f2r(r1av(1,1),r1av(2,1)),rsmow(1)), r2d(f2r(r1av(1,2),r1av(2,2)),rsmow(1))  &
        & ,r2d(f2r(rw1av(1,1),rw1av(2,1)),rsmow(1)), r2d(f2r(rw1av(1,2),rw1av(2,2)),rsmow(1))
    write(58,*) tp(iry), r2d(f2r(r1av(2,1),r1av(1,1)),rsmow(2)), r2d(f2r(r1av(2,2),r1av(1,2)),rsmow(2))  &
        & ,r2d(f2r(rw1av(2,1),rw1av(1,1)),rsmow(2)), r2d(f2r(rw1av(2,2),rw1av(1,2)),rsmow(2))
    write(59,*) tp(iry), r2dp(f2r(r1av(2,1),r1av(1,1)),rsmow(2))-(lambda*r2dp(f2r(r1av(1,1),r1av(2,1)),rsmow(1))+gamma) &  
        & ,r2dp(f2r(r1av(2,2),r1av(1,2)),rsmow(2))-(lambda*r2dp(f2r(r1av(1,2),r1av(2,2)),rsmow(1))+gamma)  &
        & ,r2dp(f2r(rw1av(2,1),rw1av(1,1)),rsmow(2))-(lambda*r2dp(f2r(rw1av(1,1),rw1av(2,1)),rsmow(1))+gamma)  &
        & ,r2dp(f2r(rw1av(2,2),rw1av(1,2)),rsmow(2))-(lambda*r2dp(f2r(rw1av(1,2),rw1av(2,2)),rsmow(1))+gamma)  

    write(60,*) tp(iry), r2d(f2r(r2av(1,1),r2av(2,1)),rsmow(1)), r2d(f2r(r2av(1,2),r2av(2,2)),rsmow(1))  &
        & ,r2d(f2r(rw2av(1,1),rw2av(2,1)),rsmow(1)), r2d(f2r(rw2av(1,2),rw2av(2,2)),rsmow(1))
    write(61,*) tp(iry), r2d(f2r(r2av(2,1),r2av(1,1)),rsmow(2)), r2d(f2r(r2av(2,2),r2av(1,2)),rsmow(2))  &
        & ,r2d(f2r(rw2av(2,1),rw2av(1,1)),rsmow(2)), r2d(f2r(rw2av(2,2),rw2av(1,2)),rsmow(2))
    write(62,*) tp(iry), r2dp(f2r(r2av(2,1),r2av(1,1)),rsmow(2))-(lambda*r2dp(f2r(r2av(1,1),r2av(2,1)),rsmow(1))+gamma) &  
        & ,r2dp(f2r(r2av(2,2),r2av(1,2)),rsmow(2))-(lambda*r2dp(f2r(r2av(1,2),r2av(2,2)),rsmow(1))+gamma)  &
        & ,r2dp(f2r(rw2av(2,1),rw2av(1,1)),rsmow(2))-(lambda*r2dp(f2r(rw2av(1,1),rw2av(2,1)),rsmow(1))+gamma)  &
        & ,r2dp(f2r(rw2av(2,2),rw2av(1,2)),rsmow(2))-(lambda*r2dp(f2r(rw2av(1,2),rw2av(2,2)),rsmow(1))+gamma)  

    write(63,*) tp(iry), r2d(f2r(r3av(1,1),r3av(2,1)),rsmow(1)), r2d(f2r(r3av(1,2),r3av(2,2)),rsmow(1))  &
        & ,r2d(f2r(rw3av(1,1),rw3av(2,1)),rsmow(1)), r2d(f2r(rw3av(1,2),rw3av(2,2)),rsmow(1))
    write(64,*) tp(iry), r2d(f2r(r3av(2,1),r3av(1,1)),rsmow(2)), r2d(f2r(r3av(2,2),r3av(1,2)),rsmow(2))  &
        & ,r2d(f2r(rw3av(2,1),rw3av(1,1)),rsmow(2)), r2d(f2r(rw3av(2,2),rw3av(1,2)),rsmow(2))
    write(65,*) tp(iry), r2dp(f2r(r3av(2,1),r3av(1,1)),rsmow(2))-(lambda*r2dp(f2r(r3av(1,1),r3av(2,1)),rsmow(1))+gamma) &  
        & ,r2dp(f2r(r3av(2,2),r3av(1,2)),rsmow(2))-(lambda*r2dp(f2r(r3av(1,2),r3av(2,2)),rsmow(1))+gamma)  &
        & ,r2dp(f2r(rw3av(2,1),rw3av(1,1)),rsmow(2))-(lambda*r2dp(f2r(rw3av(1,1),rw3av(2,1)),rsmow(1))+gamma)  &
        & ,r2dp(f2r(rw3av(2,2),rw3av(1,2)),rsmow(2))-(lambda*r2dp(f2r(rw3av(1,2),rw3av(2,2)),rsmow(1))+gamma)  

end do

close(57)
close(58)
close(59)
close(60)
close(61)
close(62)
close(63)
close(64)
close(65)

endsubroutine o_iso_ha_drdp_Phan

!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************

subroutine o_iso_ha( &
    & o18wi_in,nz,T_in,L_in,fsp,F_in,poro_in,tres,o18ri,alpha_in,k_in,a  &! input 
    & ,r1,r2,r3,rw1,rw2,rw3,r1av,r2av,r3av,rw1av,rw2av,rw3av,jocnss,L,z            &! output
    & ) 

! 3 box oceanic crust alteration model 
! converted from old .F file (oxygen18-1d-mix2_ss_v4_20_inF_FLDC-disp-dense-fp_v2-GS)
! to facilitate application to triple oxygen isotopes (16O-17O-18O)

implicit none

real(kind=8),intent(in) :: o18wi_in,fsp,tres,o18ri,a  ! fsp km^2 yr^-1
real(kind=8),dimension(3), intent(in) :: T_in,L_in,F_in,poro_in,alpha_in,k_in
integer,intent(in) :: nz

real(kind=8),dimension(nz),intent(out) :: r1,r2,r3,rw1,rw2,rw3
real(kind=8),intent(out) :: r1av,r2av,r3av,rw1av,rw2av,rw3av,jocnss,L,z(nz)

real(kind=8) :: L1 = 0.5d0 ! km
real(kind=8) :: L2 = 1.5d0     
real(kind=8) :: L3 = 4.0d0  
real(kind=8) :: rl = 1d5 ! km  ridge length
real(kind=8) :: T1 = 5.0d0 ! deg C
real(kind=8) :: T2 = 200.0d0
real(kind=8) :: T3 = 350.0d0  
real(kind=8) :: E1 = 50.0d0 ! kj mol^-1
real(kind=8) :: E2 = 50.0d0
real(kind=8) :: E3 = 50.0d0
real(kind=8) :: ms = 0.5d0/16.0d0*1.0d3 ! mol kg^-1
real(kind=8) :: mw = 0.89d0/16.0d0*1.0d3
real(kind=8) :: dense = 2.7d12 ! kg km^-3
real(kind=8) :: alfa1 = exp(27.0d0/1d3)
real(kind=8) :: alfa2 = exp(7.0d0/1d3)
real(kind=8) :: alfa3 = exp(2.5d0/1d3)
real(kind=8) :: f1 = 1.0d16  ! kg yr^-1
real(kind=8) :: f2 = 1.0d13
real(kind=8) :: f3 = 1.0d13
real(kind=8) :: p1 = 0.25d0
real(kind=8) :: p2 = 0.05d0
real(kind=8) :: p3 = 0.02d0
real(kind=8) :: densef = 1.0d12 ! kg km^-3
real(kind=8) :: denses = 3.0d12 ! kg km^-3

real(kind=8) :: tol = 1d-8

real(kind=8) :: swad = 1.0d0    ! 1.0 if porewater is carried with rock by spreading
real(kind=8) :: swsp = 1.0d0    ! 1.0 if water flows single path (lateral mixing)
real(kind=8) :: swmtl = 0.0d0   !  1.0 if initial porewater has mantle isotopic comp.or that of seawater in eq. with rock
                                !  0.0 if seawater comp.
real(kind=8) :: swbndsw = 1.0d0  ! 1.0  if initial iso. comp. is that of seawater in single flow path     
real(kind=8) :: swdif = 1.0d0  !  switch for diffuse flow (vertical mixing)      

real(kind=8) dz
real(kind=8) rwi
real(kind=8) rmi
real(kind=8) rp1i,rp2i,rp3i
real(kind=8) w, q
real(kind=8) k1, k2, k3
real(kind=8) rri 
real(kind=8) r1x(nz),r2x(nz),r3x(nz)
real(kind=8) rw1x(nz),rw2x(nz),rw3x(nz)
real(kind=8) j1t,j2t,j3t
real(kind=8) j1ts,j2ts,j3ts
real(kind=8) jocn, jocn2, jocn3
real(kind=8) errit
real(kind=8) qz1(nz),qz2(nz),qz3(nz)
real(kind=8) b1,b2,b3
real(kind=8) f1c,f2c,f3c
real(kind=8) dsp1(nz),dsp2(nz),dsp3(nz)
real(kind=8) dense1, dense2, dense3

integer row
integer iz

integer ipiv(3*2*(nz))
integer info

external DGESV

real(kind=8) amx1(3*2*(nz),3*2*(nz)), ymx1(3*2*(nz))
real(kind=8) ymx2(3*2*(nz))

!---------------------------------------------------------------------------

T1 = T_in(1); T2 = T_in(2); T3 = T_in(3)
L1 = L_in(1); L2 = L_in(2); L3 = L_in(3)
f1 = F_in(1); f2 = F_in(2); f3 = F_in(3)
p1 = poro_in(1); p2 = poro_in(2); p3 = poro_in(3)
alfa1 = alpha_in(1); alfa2 = alpha_in(2); alfa3 = alpha_in(3)
k1 = k_in(1); k2 = k_in(2); k3 = k_in(3)


rri = o18ri 
rwi = o18wi_in 

rp1i = rwi
rp2i = rwi*alfa2
rp3i = rwi*alfa3

w = fsp/rl  ! spreading rate / ridge length 

L = w*tres  ! km

b1 = log10(f1*a*log(10.0d0)/rl/L/(10.0d0**a-1.0d0))
b2 = log10(f2*a*log(10.0d0)/rl/L/(10.0d0**a-1.0d0))
b3 = log10(f3*a*log(10.0d0)/rl/L/(10.0d0**a-1.0d0))

do iz = 1,nz
    z(iz) = (iz)*L/(nz)
end do

qz1 = 0.0d0
qz2 = 0.0d0
qz3 = 0.0d0

do iz = 1, nz
    qz1(iz) = 10.0d0**(a*z(iz)/L+b1)
    qz2(iz) = 10.0d0**(a*z(iz)/L+b2)
    qz3(iz) = 10.0d0**(a*z(iz)/L+b3)
end do      

dz = z(2) - z(1)

f1c = dz*(0.50d0*(qz1(1)+qz1(nz))+sum(qz1(2:nz-1)))*rl
f2c = dz*(0.50d0*(qz2(1)+qz2(nz))+sum(qz2(2:nz-1)))*rl
f3c = dz*(0.50d0*(qz3(1)+qz3(nz))+sum(qz3(2:nz-1)))*rl

! if(abs((f1-f1c)/f1) > 1.0d-2) then
    ! write(*,*) log10(f1),log10(f1c),abs((f1-f1c)/f1) *1d2
! end if

! if(abs((f2-f2c)/f2) > 1.0d-2) then
    ! write(*,*) log10(f2),log10(f2c),abs((f2-f2c)/f2) *1d2
! end if

! if(abs((f3-f3c)/f3) > 1.0d-2) then
    ! write(*,*) log10(f3),log10(f3c),abs((f3-f3c)/f3) *1d2 
! end if

do iz = 1, nz
    dsp1(iz) = 1d-3*0.15d0*(1d3*(-0.00088d0*(z(iz)/w*1d-6)**2.0d0+0.322d0*(z(iz)/w*1d-6)+2.0d0))**0.61d0    &! Logngitudial macrodispersitivity (km)
        & *qz1(iz)/p1/L1/densef/dz/dz*(-0.00088d0*(z(iz)/w*1d-6)**2.0d0+0.322d0*(z(iz)/w*1d-6)+2.0d0)
    dsp2(iz) = 1d-3*0.15d0*(1d3*(-0.00088d0*(z(iz)/w*1d-6)**2.0d0+0.322d0*(z(iz)/w*1d-6)+2.0d0))**0.61d0    &! Logngitudial macrodispersitivity (km)
        & *qz2(iz)/p2/L2/densef/dz/dz*(-0.00088d0*(z(iz)/w*1d-6)**2.0d0+0.322d0*(z(iz)/w*1d-6)+2.0d0)
    dsp3(iz) = 1d-3*0.15d0*(1d3*(-0.00088d0*(z(iz)/w*1d-6)**2.0d0+0.322d0*(z(iz)/w*1d-6)+2.0d0))**0.61d0    &! Logngitudial macrodispersitivity (km)
        & *qz3(iz)/p3/L3/densef/dz/dz*(-0.00088d0*(z(iz)/w*1d-6)**2.0d0+0.322d0*(z(iz)/w*1d-6)+2.0d0)
end do

j1ts = 0.0d0
j1t = 0.0d0
j2ts = 0.0d0
j2t = 0.0d0
j3ts = 0.0d0
j3t = 0.0d0
jocn = 0.0d0
jocn2 = 0.0d0
jocn3 = 0.0d0

r1x = rri
rw1x = rwi
r2x = rri
rw2x = rwi
r3x = rri
rw3x = rwi

dense1 = densef*p1 + denses*(1.0d0-p1)
dense2 = densef*p2 + denses*(1.0d0-p2)
dense3 = densef*p3 + denses*(1.0d0-p3)

errit = 1d4

do while (errit > tol)

    amx1 = 0.0d0
    ymx1 = 0.0d0
    ymx2 = 0.0d0

    do iz = 1, nz

        row = 1 + 6*(iz-1)

        if (iz == 1) then
            !-- box 1 ---
            amx1(row,row) =  -mw*k1*((1.0d0-rw1x(iz))-alfa1*(-1.0d0)*rw1x(iz))/(1.0d0-p1)*dense1/denses &
                & - w/dz/1.0d0
            ymx1(row) = &
                & -w*(r1x(iz)-rri)/dz/1.0d0  &
                & -mw*k1*((1.0d0-rw1x(iz))*r1x(iz)-alfa1*(1.0d0-r1x(iz))*rw1x(iz))/(1.0d0-p1)*dense1/denses
            amx1(row, row+1) = -mw*k1*((-1.0d0)*r1x(iz)-alfa1*(1.0d0-r1x(iz)))/(1.0d0-p1)*dense1/denses

            ymx1(row + 1) = &
                & +qz1(iz)*(rwi-rw1x(iz))/L1/densef/p1*swdif  &
                & +qz2(iz)*(rw2x(iz)-rw1x(iz))/L1/densef/p1*swdif &
                & + dsp1(iz)*(rwi-rw1x(iz))*swsp  &
                & + dsp1(iz+1)*(rw1x(iz+1)-rw1x(iz))*swsp &
                & -w*(rw1x(iz)-rwi)/dz/1.0d0*swad *(1.0d0-swmtl)  &
                & -w*(rw1x(iz)-rp1i)/dz/1.0d0*swad *swmtl &
                & +ms*k1*((1.0d0-rw1x(iz))*r1x(iz)-alfa1*(1.0d0-r1x(iz))*rw1x(iz))/(p1)*dense1/densef

            amx1(row+1,row + 1) =  &
                & +qz1(iz)*(-1.0d0)/L1/densef/p1*swdif  &
                & +qz2(iz)*(-1.0d0)/L1/densef/p1*swdif &
                & + dsp1(iz)*(-1.0d0)*swsp &
                & + dsp1(iz+1)*(-1.0d0)*swsp &
                & -w*(1.0d0)/dz/1.0d0*swad *(1.0d0-swmtl)  &
                & -w*(1.0d0)/dz/1.0d0*swad *swmtl &
                & +ms*k1*((-1.0d0)*r1x(iz)-alfa1*(1.0d0-r1x(iz)))/(p1)*dense1/densef
            amx1(row+1,row +7) = + dsp1(iz+1)*(1.0d0)*swsp
            amx1(row+1,row +3) = +qz2(iz)*(1.0d0)/L1/densef/p1*swdif
            amx1(row + 1, row) = +ms*k1*((1.0d0-rw1x(iz))-alfa1*(-1.0d0)*rw1x(iz))/(p1)*dense1/densef

            !-- box 2 ---

            amx1(row+2,row+2) =  &
                & -mw*k2*((1.0d0-rw2x(iz))-alfa2*rw2x(iz)*(-1.0d0))/(1.0d0-p2)*dense2/denses  &
                & - w/dz/1.0d0
            ymx1(row+2) =  &
                & -w*(r2x(iz)-rri)/dz/1.0d0  &
                & -mw*k2*((1.0d0-rw2x(iz))*r2x(iz)-alfa2*rw2x(iz)*(1.0d0-r2x(iz)))/(1.0d0-p2)*dense2/denses
            amx1(row+2, row+3) = &
                & -mw*k2*((-1.0d0)*r2x(iz)-alfa2*(1.0d0-r2x(iz)))/(1.0d0-p2)*dense2/denses

            ymx1(row + 3) =  &
                & +qz2(iz)*(rw1x(iz)-rw2x(iz))/L2/densef/p2*swdif  &
                & +qz3(iz)*(rw3x(iz)-rw2x(iz))/L2/densef/p2*swdif  &
                & + dsp2(iz)*(rwi-rw2x(iz))*swsp  &
                & + dsp2(iz+1)*(rw2x(iz+1)-rw2x(iz))*swsp  &
                & -w*(rw2x(iz)-rwi)/dz/1.0d0*swad *(1.0d0-swmtl)  &
                & -w*(rw2x(iz)-rp2i)/dz/1.0d0*swad *swmtl  &
                & +ms*k2*((1.0d0-rw2x(iz))*r2x(iz)-alfa2*(1.0d0-r2x(iz))*rw2x(iz))/(p2)*dense2/densef
            amx1(row+3,row + 3) =  &
                & +qz2(iz)*(-1.0d0)/L2/densef/p2*swdif  &
                & +qz3(iz)*(-1.0d0)/L2/densef/p2*swdif  &
                & + dsp2(iz)*(-1.0d0)*swsp  &
                & + dsp2(iz+1)*(-1.0d0)*swsp  &
                & -w*(1.0d0)/dz/1.0d0*swad *(1.0d0-swmtl)  &
                & -w*(1.0d0)/dz/1.0d0*swad *swmtl  &
                & +ms*k2*((-1.0d0)*r2x(iz)-alfa2*(1.0d0-r2x(iz)))/(p2)*dense2/densef
            amx1(row+3,row +9) = + dsp2(iz+1)*(1.0d0)*swsp
            amx1(row+3,row +5) = +qz3(iz)*(1.0d0)/L2/densef/p2*swdif
            amx1(row+3,row +1) = +qz2(iz)*(1.0d0)/L2/densef/p2*swdif
            amx1(row + 3, row + 2) = +ms*k2*((1.0d0-rw2x(iz))-alfa2*(-1.0d0)*rw2x(iz))/(p2)*dense2/densef

            !-- box 3 ---

            amx1(row+4,row+4) = &
                & -w*(1.0d0)/dz/1.0d0 &
                & -mw*k3*((1.0d0-rw3x(iz))-alfa3*rw3x(iz)*(-1.0d0))/(1.0d0-p3)*dense3/denses
            ymx1(row+4) = &
                & -w*(r3x(iz)-rri)/dz/1.0d0  &
                & -mw*k3*((1.0d0-rw3x(iz))*r3x(iz)-alfa3*rw3x(iz)*(1.0d0-r3x(iz)))/(1.0d0-p3)*dense3/denses
            amx1(row+4, row+5) = -mw*k3*((-1.0d0)*r3x(iz)-alfa3*(1.0d0-r3x(iz)))/(1.0d0-p3)*dense3/denses

            ymx1(row + 5) = &
                & +qz3(iz)*(rw2x(iz)-rw3x(iz))/L3/densef/p3*swdif  &
                & + dsp3(iz)*(rwi-rw3x(iz))*swsp &
                & + dsp3(iz+1)*(rw3x(iz+1)-rw3x(iz))*swsp  &
                & -w*(rw3x(iz)-rwi)/dz/1.0d0*swad *(1.0d0-swmtl) &
                & -w*(rw3x(iz)-rp3i)/dz/1.0d0*swad *swmtl  &
                & +ms*k3*((1.0d0-rw3x(iz))*r3x(iz)-alfa3*(1.0d0-r3x(iz))*rw3x(iz))/(p3)*dense3/densef
            amx1(row+5,row + 5) = &
                & +qz3(iz)*(-1.0d0)/L3/densef/p3*swdif  &
                & + dsp3(iz)*(-1.0d0)*swsp  &
                & + dsp3(iz+1)*(-1.0d0)*swsp  &
                & -w*(1.0d0)/dz/1.0d0*swad *(1.0d0-swmtl)  &
                & -w*(1.0d0)/dz/1.0d0*swad *swmtl  &
                & +ms*k3*((-1.0d0)*r3x(iz)-alfa3*(1.0d0-r3x(iz)))/(p3)*dense3/densef
            amx1(row+5,row +11) = + dsp3(iz+1)*(1.0d0)*swsp
            amx1(row+5,row +3) = +qz3(iz)*(1.0d0)/L3/densef/p3*swdif
            amx1(row + 5, row + 4) = +ms*k3*((1.0d0-rw3x(iz))-alfa3*(-1.0d0)*rw3x(iz))/(p3)*dense3/densef

        else if (iz == nz) then
            !-- box 1 ---
            amx1(row,row) =  -mw*k1*((1.0d0-rw1x(iz))-alfa1*(-1.0d0)*rw1x(iz))/(1.0d0-p1)*dense1/denses  &
                & - w/dz/1.0d0
            amx1(row,row-6) = -w*(-1.0d0)/dz/1.0d0
            ymx1(row) = &
                & -w*(r1x(iz)-r1x(iz-1))/dz/1.0d0  &
                & -mw*k1*((1.0d0-rw1x(iz))*r1x(iz)-alfa1*(1.0d0-r1x(iz))*rw1x(iz))/(1.0d0-p1)*dense1/denses
            amx1(row, row+1) = -mw*k1*((-1.0d0)*r1x(iz)-alfa1*(1.0d0-r1x(iz)))/(1.0d0-p1)*dense1/denses

            ymx1(row + 1) = &
                & +qz1(iz)*(rwi-rw1x(iz))/L1/densef/p1*swdif  &
                & +qz2(iz)*(rw2x(iz)-rw1x(iz))/L1/densef/p1*swdif &
                & + dsp1(iz)*(rw1x(iz-1)-rw1x(iz))*swsp &
                & -w*(rw1x(iz)-rw1x(iz-1))/dz/1.0d0*swad  &
                & +ms*k1*((1.0d0-rw1x(iz))*r1x(iz)-alfa1*(1.0d0-r1x(iz))*rw1x(iz))/(p1)*dense1/densef

            amx1(row+1,row + 1) =  &
                & +qz1(iz)*(-1.0d0)/L1/densef/p1*swdif  &
                & +qz2(iz)*(-1.0d0)/L1/densef/p1*swdif  &
                & + dsp1(iz)*(-1.0d0)*swsp  &
                & -w*(1.0d0)/dz/1.0d0*swad   &
                & +ms*k1*((-1.0d0)*r1x(iz)-alfa1*(1.0d0-r1x(iz)))/(p1)*dense1/densef
            amx1(row+1,row - 5) =  &
                & + dsp1(iz)*(1.0d0)*swsp &
                & -w*(-1.0d0)/dz/1.0d0*swad 
            amx1(row+1,row +3) = +qz2(iz)*(1.0d0)/L1/densef/p1*swdif
            amx1(row + 1, row) = +ms*k1*((1.0d0-rw1x(iz))-alfa1*(-1.0d0)*rw1x(iz))/(p1)*dense1/densef

            !-- box 2 ---

            amx1(row+2,row+2) = -mw*k2*((1.0d0-rw2x(iz))-alfa2*rw2x(iz)*(-1.0d0))/(1.0d0-p2)*dense2/denses  &
                & - w/dz/1.0d0
            amx1(row+2,row-4) = -w*(-1.0d0)/dz/1.0d0
            ymx1(row+2) =  -w*(r2x(iz)-r2x(iz-1))/dz/1.0d0  &
                & -mw*k2*((1.0d0-rw2x(iz))*r2x(iz)-alfa2*rw2x(iz)*(1.0d0-r2x(iz)))/(1.0d0-p2)*dense2/denses
            amx1(row+2, row+3) = -mw*k2*((-1.0d0)*r2x(iz)-alfa2*(1.0d0-r2x(iz)))/(1.0d0-p2)*dense2/denses

            ymx1(row + 3) = &
                & +qz2(iz)*(rw1x(iz)-rw2x(iz))/L2/densef/p2*swdif  &
                & +qz3(iz)*(rw3x(iz)-rw2x(iz))/L2/densef/p2*swdif  &
                & + dsp2(iz)*(rw2x(iz-1)-rw2x(iz))*swsp  &
                & -w*(rw2x(iz)-rw2x(iz-1))/dz/1.0d0*swad   &
                & +ms*k2*((1.0d0-rw2x(iz))*r2x(iz)-alfa2*(1.0d0-r2x(iz))*rw2x(iz))/(p2)*dense2/densef
            amx1(row+3,row + 3) = &
                & +qz2(iz)*(-1.0d0)/L2/densef/p2*swdif  &
                & +qz3(iz)*(-1.0d0)/L2/densef/p2*swdif  &
                & + dsp2(iz)*(-1.0d0)*swsp   &
                & -w*(1.0d0)/dz/1.0d0*swad   &
                & +ms*k2*((-1.0d0)*r2x(iz)-alfa2*(1.0d0-r2x(iz)))/(p2)*dense2/densef
            amx1(row+3,row - 3) =  &
                & + dsp2(iz)*(1.0d0)*swsp  &
                & -w*(-1.0d0)/dz/1.0d0*swad 
            amx1(row+3,row +5) = +qz3(iz)*(1.0d0)/L2/densef/p2*swdif
            amx1(row+3,row +1) = +qz2(iz)*(1.0d0)/L2/densef/p2*swdif
            amx1(row + 3, row + 2) = +ms*k2*((1.0d0-rw2x(iz))-alfa2*(-1.0d0)*rw2x(iz))/(p2)*dense2/densef

            !-- box 3 ---

            amx1(row+4,row+4) = &
                & -w*(1.0d0)/dz/1.0d0  &
                & -mw*k3*((1.0d0-rw3x(iz))-alfa3*rw3x(iz)*(-1.0d0))/(1.0d0-p3)*dense3/denses
            amx1(row+4,row-2) = -w*(-1.0d0)/dz/1.0d0
            ymx1(row+4) =  &
                & -w*(r3x(iz)-r3x(iz-1))/dz/1.0d0  &
                & -mw*k3*((1.0d0-rw3x(iz))*r3x(iz)-alfa3*rw3x(iz)*(1.0d0-r3x(iz)))/(1.0d0-p3)*dense3/denses
            amx1(row+4, row+5) = -mw*k3*((-1.0d0)*r3x(iz)-alfa3*(1.0d0-r3x(iz)))/(1.0d0-p3)*dense3/denses

            ymx1(row + 5) = &
                & +qz3(iz)*(rw2x(iz)-rw3x(iz))/L3/densef/p3*swdif  &
                & + dsp3(iz)*(rw3x(iz-1)-rw3x(iz))*swsp  &
                & -w*(rw3x(iz)-rw3x(iz-1))/dz/1.0d0*swad   &
                & +ms*k3*((1.0d0-rw3x(iz))*r3x(iz)-alfa3*(1.0d0-r3x(iz))*rw3x(iz))/(p3)*dense3/densef
            amx1(row+5,row + 5) =   &
                & +qz3(iz)*(-1.0d0)/L3/densef/p3*swdif  &
                & + dsp3(iz)*(-1.0d0)*swsp  &
                & -w*(1.0d0)/dz/1.0d0*swad   &
                & +ms*k3*((-1.0d0)*r3x(iz)-alfa3*(1.0d0-r3x(iz)))/(p3)*dense3/densef
            amx1(row+5,row - 1) = &
            & + dsp3(iz)*(1.0d0)*swsp  &
            & -w*(-1.0d0)/dz/1.0d0*swad 
            amx1(row+5,row +3) = +qz3(iz)*(1.0d0)/L3/densef/p3*swdif
            amx1(row + 5, row + 4) = +ms*k3*((1.0d0-rw3x(iz))-alfa3*(-1.0d0)*rw3x(iz))/(p3)*dense3/densef

        else 
            !-- box 1 ---
            amx1(row,row) =  -mw*k1*((1.0d0-rw1x(iz))-alfa1*(-1.0d0)*rw1x(iz))/(1.0d0-p1)*dense1/denses &
                & - w/dz/1.0d0
            amx1(row,row-6) = -w*(-1.0d0)/dz/1.0d0
            ymx1(row) = &
                & -w*(r1x(iz)-r1x(iz-1))/dz/1.0d0  &
                & -mw*k1*((1.0d0-rw1x(iz))*r1x(iz)-alfa1*(1.0d0-r1x(iz))*rw1x(iz))/(1.0d0-p1)*dense1/denses
            amx1(row, row+1) = -mw*k1*((-1.0d0)*r1x(iz)-alfa1*(1.0d0-r1x(iz)))/(1.0d0-p1)*dense1/denses

            ymx1(row + 1) = &
                & +qz1(iz)*(rwi-rw1x(iz))/L1/densef/p1*swdif &
                & +qz2(iz)*(rw2x(iz)-rw1x(iz))/L1/densef/p1*swdif &
                & + dsp1(iz)*(rw1x(iz-1)-rw1x(iz))*swsp &
                & + dsp1(iz+1)*(rw1x(iz+1)-rw1x(iz))*swsp &
                & -w*(rw1x(iz)-rw1x(iz-1))/dz/1.0d0*swad  &
                & +ms*k1*((1.0d0-rw1x(iz))*r1x(iz)-alfa1*(1.0d0-r1x(iz))*rw1x(iz))/(p1)*dense1/densef

            amx1(row+1,row + 1) =  &
                & +qz1(iz)*(-1.0d0)/L1/densef/p1*swdif  &
                & +qz2(iz)*(-1.0d0)/L1/densef/p1*swdif  &
                & + dsp1(iz)*(-1.0d0)*swsp   &
                & + dsp1(iz+1)*(-1.0d0)*swsp   &
                & -w*(1.0d0)/dz/1.0d0*swad  &
                & +ms*k1*((-1.0d0)*r1x(iz)-alfa1*(1.0d0-r1x(iz)))/(p1)*dense1/densef
            amx1(row+1,row - 5) = + dsp1(iz)*(1.0d0)*swsp  &
                & -w*(-1.0d0)/dz/1.0d0*swad 
            amx1(row+1,row +7) = + dsp1(iz+1)*(1.0d0)*swsp
            amx1(row+1,row +3) = + qz2(iz)*(1.0d0)/L1/densef/p1*swdif
            amx1(row + 1, row) = +ms*k1*((1.0d0-rw1x(iz))-alfa1*(-1.0d0)*rw1x(iz))/(p1)*dense1/densef

            !-- box 2 ---

            amx1(row+2,row+2) =  &
                & -mw*k2*((1.0d0-rw2x(iz))-alfa2*rw2x(iz)*(-1.0d0))/(1.0d0-p2)*dense2/denses  &
                & - w/dz/1.0d0
            amx1(row+2,row-4) = -w*(-1.0d0)/dz/1.0d0
            ymx1(row+2) =  &
                & -w*(r2x(iz)-r2x(iz-1))/dz/1.0d0  &
                & -mw*k2*((1.0d0-rw2x(iz))*r2x(iz)-alfa2*rw2x(iz)*(1.0d0-r2x(iz)))/(1.0d0-p2)*dense2/denses
            amx1(row+2, row+3) = -mw*k2*((-1.0d0)*r2x(iz)-alfa2*(1.0d0-r2x(iz)))/(1.0d0-p2)*dense2/denses

            ymx1(row + 3) =  &
                & +qz2(iz)*(rw1x(iz)-rw2x(iz))/L2/densef/p2*swdif &
                & +qz3(iz)*(rw3x(iz)-rw2x(iz))/L2/densef/p2*swdif &
                & + dsp2(iz)*(rw2x(iz-1)-rw2x(iz))*swsp  &
                & + dsp2(iz+1)*(rw2x(iz+1)-rw2x(iz))*swsp &
                & -w*(rw2x(iz)-rw2x(iz-1))/dz/1.0d0*swad  &
                & +ms*k2*((1.0d0-rw2x(iz))*r2x(iz)-alfa2*(1.0d0-r2x(iz))*rw2x(iz))/(p2)*dense2/densef
            amx1(row+3,row + 3) =  &
                & +qz2(iz)*(-1.0d0)/L2/densef/p2*swdif  &
                & +qz3(iz)*(-1.0d0)/L2/densef/p2*swdif  &
                & + dsp2(iz)*(-1.0d0)*swsp  &
                & + dsp2(iz+1)*(-1.0d0)*swsp  &
                & -w*(1.0d0)/dz/1.0d0*swad   &
                & +ms*k2*((-1.0d0)*r2x(iz)-alfa2*(1.0d0-r2x(iz)))/(p2)*dense2/densef
            amx1(row+3,row - 3) = + dsp2(iz)*(1.0d0)*swsp  &
                & -w*(-1.0d0)/dz/1.0d0*swad 
            amx1(row+3,row +9) = + dsp2(iz+1)*(1.0d0)*swsp
            amx1(row+3,row +5) = +qz3(iz)*(1.0d0)/L2/densef/p2*swdif
            amx1(row+3,row +1) = +qz2(iz)*(1.0d0)/L2/densef/p2*swdif
            amx1(row + 3, row + 2) = +ms*k2*((1.0d0-rw2x(iz))-alfa2*(-1.0d0)*rw2x(iz))/(p2)*dense2/densef

            !-- box 3 ---

            amx1(row+4,row+4) =  &
                & -w*(1.0d0)/dz/1.0d0  &
                & -mw*k3*((1.0d0-rw3x(iz))-alfa3*rw3x(iz)*(-1.0d0))/(1.0d0-p3)*dense3/denses
            amx1(row+4,row-2) = -w*(-1.0d0)/dz/1.0d0
            ymx1(row+4) = &
                & -w*(r3x(iz)-r3x(iz-1))/dz/1.0d0 &
                & -mw*k3*((1.0d0-rw3x(iz))*r3x(iz)-alfa3*rw3x(iz)*(1.0d0-r3x(iz)))/(1.0d0-p3)*dense3/denses
            amx1(row+4, row+5) = -mw*k3*((-1.0d0)*r3x(iz)-alfa3*(1.0d0-r3x(iz)))/(1.0d0-p3)*dense3/denses

            ymx1(row + 5) = &
                & +qz3(iz)*(rw2x(iz)-rw3x(iz))/L3/densef/p3*swdif  &
                & + dsp3(iz)*(rw3x(iz-1)-rw3x(iz))*swsp  &
                & + dsp3(iz+1)*(rw3x(iz+1)-rw3x(iz))*swsp &
                & -w*(rw3x(iz)-rw3x(iz-1))/dz/1.0d0*swad  &
                & +ms*k3*((1.0d0-rw3x(iz))*r3x(iz)-alfa3*(1.0d0-r3x(iz))*rw3x(iz))/(p3)*dense3/densef
            amx1(row+5,row+5) =  &
                & +qz3(iz)*(-1.0d0)/L3/densef/p3*swdif  &
                & + dsp3(iz)*(-1.0d0)*swsp  &
                & + dsp3(iz+1)*(-1.0d0)*swsp &
                & -w*(1.0d0)/dz/1.0d0*swad   &
                & +ms*k3*((-1.0d0)*r3x(iz)-alfa3*(1.0d0-r3x(iz)))/(p3)*dense3/densef
            amx1(row+5,row - 1) =  &
            & + dsp3(iz)*(1.0d0)*swsp  &
            & -w*(-1.0d0)/dz/1.0d0*swad 
            amx1(row+5,row +11) = + dsp3(iz+1)*(1.0d0)*swsp
            amx1(row+5,row +3) = +qz3(iz)*(1.0d0)/L3/densef/p3*swdif
            amx1(row + 5, row + 4) = +ms*k3*((1.0d0-rw3x(iz))-alfa3*(-1.0d0)*rw3x(iz))/(p3)*dense3/densef

        end if

    end do

    ymx1 = -1.0d0*ymx1

    call DGESV(3*2*(nz),int(1),amx1,3*2*(nz),IPIV,ymx1,3*2*(nz),INFO)

    do iz = 1, nz
        row = 1 + 6*(iz -1)
        r1x(iz) = r1x(iz) + ymx1(row)
        rw1x(iz) = rw1x(iz) + ymx1(row+1)
        r2x(iz) = r2x(iz) + ymx1(row+2)
        rw2x(iz) =rw2x(iz) + ymx1(row+3)
        r3x(iz) = r3x(iz) + ymx1(row+4)
        rw3x(iz) = rw3x(iz) + ymx1(row+5)
        if (r1x(iz)/=0.0d0) ymx2(row) = abs(ymx1(row)/r1x(iz))
        if (rw1x(iz)/=0.0d0) ymx2(row+1) = abs(ymx1(row+1)/rw1x(iz))
        if (r2x(iz)/=0.0d0) ymx2(row+2) = abs(ymx1(row+2)/r2x(iz))
        if (rw2x(iz)/=0.0d0) ymx2(row+3) = abs(ymx1(row+3)/rw2x(iz))
        if (r3x(iz)/=0.0d0) ymx2(row+4) = abs(ymx1(row+4)/r3x(iz))
        if (rw3x(iz)/=0.0d0) ymx2(row+5) = abs(ymx1(row+5)/rw3x(iz))
    end do

    errit = maxval(ymx2)

    ! write (*,*) errit, info

end do

r1 = r1x
rw1 = rw1x
r2 = r2x
rw2 = rw2x
r3 = r3x
rw3 = rw3x

r1av = dz*(rri*0.5d0+r1(nz)*0.50d0+sum(r1(1:nz-1)))/L
rw1av = dz*(rwi*0.5d0+rw1(nz)*0.50d0+sum(rw1(1:nz-1)))/L
r2av = dz*(rri*0.5d0+r2(nz)*0.50d0+sum(r2(1:nz-1)))/L
rw2av = dz*(rwi*0.5d0+rw2(nz)*0.50d0+sum(rw2(1:nz-1)))/L
r3av = dz*(rri*0.5d0+r3(nz)*0.50d0+sum(r3(1:nz-1)))/L
rw3av = dz*(rwi*0.5d0+rw3(nz)*0.50d0+sum(rw3(1:nz-1)))/L

do iz = 1,  nz
    j1ts = j1ts &
        & + ms*L1*(1.0d0-p1)*w*denses*(merge(rri,r1(iz-1),iz==1)- r1(iz))  &
        & + ms*L1*dense1*dz*mw*k1*(r1(iz)*(1.0d0-rw1(iz)) &
        & -alfa1*rw1(iz)*(1.0d0-r1(iz)))

    j1t = j1t  &
        & +mw*L1*densef*w*p1*(rw1(iz-1)-rw1(iz)) *swad  &
        & - mw*L1*dense1*dz*ms*k1*(r1(iz)*(1.0d0-rw1(iz))  &
        & - alfa1*rw1(iz)*(1.0d0-r1(iz)))  &
        & +swdif*f1/rl*dz*mw/L*(rwi-rw1(iz))  &
        & +swdif*f2/rl*dz*mw/L*(rw2(iz)-rw1(iz))  &
        & + swsp * f1/rl*mw*((merge(rw1(iz+1),rwi,iz/=nz))-rw1(iz))

    j2ts = j2ts   &
        & + ms*L2*(1.0d0-p2)*w*denses*(r2(iz-1) - r2(iz))  &
        & + ms*L2*dense2*dz*mw*k2*(r2(iz)*(1.0d0-rw2(iz))   &
        & -alfa2*rw2(iz)*(1.0d0-r2(iz)))

    j2t = j2t  &
        & +mw*L2*densef*w*p2*(rw2(iz-1)-rw2(iz)) *swad  &
        & - mw*L2*dense2*dz*ms*k2*(r2(iz)*(1.0d0-rw2(iz))  &
        & - alfa2*rw2(iz)*(1.0d0-r2(iz)))  &
        & +swdif*f2/rl*dz*mw/L*(rw1(iz)-rw2(iz))  &
        & +swdif*f3/rl*dz*mw/L*(rw3(iz)-rw2(iz))  &
        & + swsp * f2/rl*mw*((merge(rw2(iz+1),  &
        & rw1(iz)*(1.0d0-swbndsw)+rwi*swbndsw,iz/=nz))-rw2(iz))

    j3ts = j3ts   &
        & + ms*L3*(1.0d0-p3)*w*denses*(r3(iz-1) - r3(iz))  &
        & + ms*L3*dense3*dz*mw*k3*(r3(iz)*(1.0d0-rw3(iz))  &
        & -alfa3*rw3(iz)*(1.0d0-r3(iz)))

    j3t = j3t  &
        & +mw*L3*densef*w*p3*(rw3(iz-1)-rw3(iz)) *swad  &
        & - mw*L3*dense3*dz*ms*k3*(r3(iz)*(1.0d0-rw3(iz))  & 
        & - alfa3*rw3(iz)*(1.0d0-r3(iz)))  &
        & +swdif*f3/rl*dz*mw/L*(rw2(iz)-rw3(iz))  &
        & + swsp * f3/rl*mw*((merge(rw3(iz+1),  &
        & rw2(iz)*(1.0d0-swbndsw)+rwi*swbndsw,iz/=nz))-rw3(iz))

    jocn = jocn  &
        & +swdif*qz1(iz)*dz*mw*(rwi-rw1(iz))  &
        & + swsp * dz*dsp1(iz)*mw*densef*L1*p1*  &
        & ((merge(rw1(iz-1),rwi,iz/=1))-rw1(iz))  &
        & + swsp * dz*merge(0.0d0,  &
        & dsp1(iz+1)*mw*densef*L1*p1*(rw1(iz+1)-rw1(iz)), iz .eq. nz)  &
        & + swsp * dz*dsp2(iz)*mw*densef*L2*p2*  &
        & ((merge(rw2(iz-1),rwi,iz/=1))-rw2(iz))  &
        & + swsp * dz*merge(0.0d0,  &
        & dsp2(iz+1)*mw*densef*L2*p2*(rw2(iz+1)-rw2(iz)), iz .eq. nz)  &
        & + swsp * dz*dsp3(iz)*mw*densef*L3*p3*  &
        & ((merge(rw3(iz-1),rwi,iz/=1))-rw3(iz)) &
        & + swsp * dz*merge(0.0d0,  &
        & dsp3(iz+1)*mw*densef*L3*p3*(rw3(iz+1)-rw3(iz)), iz .eq. nz)

    jocn3 = jocn3   &
        & + mw*L3*densef*w*p3*(merge(rwi,rw3(iz-1),iz==1)-rw3(iz)) *swad  &
        & + mw*L2*densef*w*p2*(merge(rwi,rw2(iz-1),iz==1)-rw2(iz)) *swad  &
        & + mw*L1*densef*w*p1*(merge(rwi,rw1(iz-1),iz==1)-rw1(iz)) *swad

    jocn2 = jocn2  &
        & + ms*L1*(1.0d0-p1)*w*denses*(merge(rri,r1(iz-1),iz==1) - r1(iz))  &
        & + ms*L2*(1.0d0-p2)*w*denses*(merge(rri,r2(iz-1),iz==1) - r2(iz))  &
        & + ms*L3*(1.0d0-p3)*w*denses*(merge(rri,r3(iz-1),iz==1) - r3(iz))
end do

j1ts = rl*j1ts
j1t = -rl*j1t
j2ts = rl*j2ts
j2t = -rl*j2t
j3ts = rl*j3ts
j3t = -rl*j3t
jocn = -rl*jocn
jocn2 = rl*jocn2
jocn3 = -rl*jocn3

jocnss = jocn + jocn3 ! porewater included

if (abs((jocn2-(jocn+jocn3))/jocn2) > 1d-3) then
    print *, 'error',abs((jocn2-(jocn+jocn3))/jocn2)     
    pause
endif 

endsubroutine o_iso_ha

endmodule o_iso_ha_mod