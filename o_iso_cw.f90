! program o_iso_cw_main
! implicit none

! call o_iso_cw_sense_dsw_prof  
! call o_iso_cw_sense_dsw_Phan
! call o_iso_cw_dr_Phan
! endprogram o_iso_cw_main

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine o_iso_cw_sense_dsw_prof  

implicit none

real(kind=8) :: w = 5d-5
real(kind=8) :: q = 3d-1
real(kind=8) :: poro = 0.20d0
real(kind=8) :: tau = 1d6
real(kind=8) :: tempcw = 15.0d0
real(kind=8) :: ztot = 20.0d0
real(kind=8) :: sat = 0.50d0
real(kind=8) :: e = 50.0d0

real(kind=8) :: j1,j2,j3,cdf,shf,dsf

integer,parameter :: n1 = 11
integer,parameter :: nz = 200
integer,parameter :: niso = 2

real(kind=8) beta_eq, beta_kin, lambda, gamma

integer iz, iiso, iiso2

real(kind=8),dimension(nz) :: z
real(kind=8),dimension(niso,nz) :: fr,fp
real(kind=8),dimension(niso) :: fsw,foc,frtop,dsw,doc,rstd,alfacw,k1,slp,itcpt

real(kind=8) f2r,r2d,d2r,r2f,r2dp,mwl_Luz,dp2d
real(kind=8) beta_eq_pack,beta_eq_sharp
real(kind=8) alfacw_savin
real(kind=8) t2alfacw,t2beta_eq

character (100) workdir
character (50) beta_eq_model, alfacw_model
!!!!!! ------------------------------------------------------------------------
workdir = '../oiso_output'
call system ('mkdir -p '//trim(adjustl(workdir)))
open (20, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof.txt',status = 'replace')

beta_eq = (1d0/16d0 - 1d0/17d0)/(1d0/16d0 - 1d0/18d0)
beta_eq = beta_eq_pack(tempcw+273.15d0) ! Pack and Herwartz (2014)
beta_eq = beta_eq_sharp(tempcw+273.15d0) ! Sharp et al. (2016)
beta_eq_model = 'sharp16'
beta_eq = t2beta_eq(tempcw+273.15d0,beta_eq_model) ! Sharp et al. (2016)
beta_kin = beta_eq
beta_kin = 1d0 ! this assume kinetics of 17O/16O exchange is not different from that of 18O/16O

lambda = 0.5305d0
gamma = 0d0

rstd(1) = 2.0052d-3 
rstd(2) = 3.799d-4 
dsw(1) = 0d0
dsw(2) = 0d0
dsw(2) = dp2d(mwl_Luz(dsw(1),rstd(1)),rstd(2))
doc(1) = 5.7d0 
doc(2) = 2.86d0 
! doc(2) = 3.0d0 
alfacw(1) = exp(25.0d0/1d3)
alfacw(1) = alfacw_savin(tempcw + 273.0d0)
alfacw_model = 'savin'
alfacw(1) = t2alfacw(tempcw + 273.0d0,alfacw_model)
alfacw(2) = alfacw(1)**beta_eq
k1(1) = 10.0d0**(-8.50d0)
k1(2) = k1(1)**beta_kin

do iiso = 1, niso    
    if (iiso == 1) iiso2 = 2
    if (iiso == 2) iiso2 = 1
    fsw(iiso) = r2f(d2r(dsw(iiso),rstd(iiso)),d2r(dsw(iiso2),rstd(iiso2)))
    foc(iiso) = r2f(d2r(doc(iiso),rstd(iiso)),d2r(doc(iiso2),rstd(iiso2)))

    call o_iso_cw(  &
        & w, k1(iiso),q,tau,poro,fsw(iiso),nz,foc(iiso),rstd(iiso),alfacw(iiso),sat,e,tempcw &! input
        & ,j1,j2,j3,cdf,shf,dsf,z,fr(iiso,:),fp(iiso,:),frtop(iiso)   &! output
        & ,ztot  &! inout
        & )
enddo 

write(20,*) 0d0,0d0,r2d(f2r(frtop(1),frtop(2)),rstd(1)), r2d(f2r(fsw(1),fsw(2)),rstd(1)) &
        & , r2d(f2r(frtop(2),frtop(1)),rstd(2)), r2d(f2r(fsw(2),fsw(1)),rstd(2))   &
        & , r2dp(f2r(frtop(2),frtop(1)),rstd(2))   &
            & - ( lambda* r2dp(f2r(frtop(1),frtop(2)),rstd(1)) + gamma)  &
        & , r2dp(f2r(fsw(2),fsw(1)),rstd(2))   &
            & - ( lambda* r2dp(f2r(fsw(1),fsw(2)),rstd(1)) + gamma)  
do iz = 1, nz
    write(20,*) z(iz),z(iz)/ztot, r2d(f2r(fr(1,iz),fr(2,iz)),rstd(1)), r2d(f2r(fp(1,iz),fp(2,iz)),rstd(1))  &
        & , r2d(f2r(fr(2,iz),fr(1,iz)),rstd(2)), r2d(f2r(fp(2,iz),fp(1,iz)),rstd(2))   &
        & , r2dp(f2r(fr(2,iz),fr(1,iz)),rstd(2))   &
            & - ( lambda* r2dp(f2r(fr(1,iz),fr(2,iz)),rstd(1)) + gamma)  &
        & , r2dp(f2r(fp(2,iz),fp(1,iz)),rstd(2))   &
            & - ( lambda* r2dp(f2r(fp(1,iz),fp(2,iz)),rstd(1)) + gamma)  
enddo 
close(20)

call o_iso_cw_sense_dsw( &
    & w,q,poro,k1,tau,n1,nz,doc,rstd,alfacw,sat,e,tempcw,niso,lambda,gamma  &!input 
    & ,slp,itcpt  &! output 
    )

print *, slp,itcpt

! call o_iso_cw_sense_Phan

endsubroutine o_iso_cw_sense_dsw_prof

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine o_iso_cw_sense_dsw( &
    & w,q,poro,k1,tau,n1,nz,doc,rstd,alfacw,sat,e,tempcw,niso,lambda,gamma  &!input 
    & ,slp,itcpt  &! output 
    )
implicit none

real(kind=8),intent(in) :: w,q,poro,tau,sat,e,tempcw,lambda,gamma

integer, intent(in) :: n1,nz,niso 
real(kind=8),dimension(niso),intent(in) :: k1,doc,rstd,alfacw
real(kind=8),dimension(niso),intent(out) :: slp,itcpt

real(kind=8),dimension(niso,n1) :: j1,j2,j3,dswl,fswl
real(kind=8),dimension(niso,n1) :: cdf,shf,dsf,frtop
real(kind=8),dimension(niso,n1,nz) :: fr,fp
real(kind=8),dimension(nz) :: z
real(kind=8),dimension(niso) :: avx,avy,sumx2,sumy2,sumxy,crr
real(kind=8),dimension(niso) :: dsw,fsw,foc
real(kind=8),dimension(niso) :: alfacw_loc,k1_loc
real(kind=8) ztot
integer i1, imin, iz, iiso, iiso2
character (100) workdir

real(kind=8) dp2d,mwl_Luz,r2f,d2r,r2dp,f2r,r2d,tc2d18orw,beta_eq_pack,beta_eq_sharp
real(kind=8) t2alfacw,t2beta_eq
real(kind=8) tempcw_loc,beta_eq, beta_kin
character (50) beta_eq_model, alfacw_model, d18orw_model
! logical :: temp_sense = .true.
logical :: temp_sense = .false.
!---------------------------------------------------      
workdir = '../oiso_output'

call system ('mkdir -p '//trim(adjustl(workdir)))

open (11, file =trim(adjustl(workdir))//'/'//'o_iso_cw_flux_18.txt',status = 'replace')
open (12, file =trim(adjustl(workdir))//'/'//'o_iso_cw_flux_17.txt',status = 'replace')
open (15, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof_sld_d18.txt',status = 'replace')
open (16, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof_pw_d18.txt',status = 'replace')
open (17, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof_sld_d17.txt',status = 'replace')
open (18, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof_pw_d17.txt',status = 'replace')
open (19, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof_sld_capD17.txt',status = 'replace')
open (20, file =trim(adjustl(workdir))//'/'//'o_iso_cw_prof_pw_capD17.txt',status = 'replace')

imin = 0
   
beta_eq_model = 'sharp16' 
alfacw_model = 'savin'
d18orw_model = 'bowen08'
   
tempcw_loc = tempcw 
k1_loc = k1
alfacw_loc = alfacw

do i1 = 1, n1

    dsw(1) = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
    dswl(1,i1) = dsw(1)
    if (dsw(1) <= -8d0) imin = i1 ! imin wants to satisfy dswl(imin) = -8 
    
    dsw(2) = dp2d(mwl_Luz(dsw(1),rstd(1)),rstd(2))
    dswl(2,i1) = dsw(2)
    
    if (temp_sense) then 
        tempcw_loc = -5.0d0 + 35.0d0*(i1-1.0d0)/(n1-1.0d0)
        dsw(1) = tc2d18orw(tempcw_loc,d18orw_model)
        dswl(1,i1) = dsw(1)
        if (dsw(1) <= -8d0) imin = i1 ! imin wants to satisfy dswl(imin) = -8     
        dsw(2) = dp2d(mwl_Luz(dsw(1),rstd(1)),rstd(2))
        dswl(2,i1) = dsw(2)
        beta_eq = (1d0/16d0 - 1d0/17d0)/(1d0/16d0 - 1d0/18d0)
        beta_eq = beta_eq_pack(tempcw_loc+273.15d0) ! Pack and Herwartz (2014)
        beta_eq = t2beta_eq(tempcw_loc+273.15d0,beta_eq_model) ! Sharp et al. (2016)
        beta_kin = beta_eq
        beta_kin = 1d0 ! this assume kinetics of 17O/16O exchange is not different from that of 18O/16O
        alfacw_loc(1) = t2alfacw(tempcw_loc + 273.15d0,alfacw_model)
        alfacw_loc(2) = alfacw_loc(1)**beta_eq
        k1_loc(1) = 10.0d0**(-8.50d0)
        k1_loc(2) = k1_loc(1)**beta_kin
    endif 
        
    
    do iiso = 1, niso
        if (iiso == 1) iiso2 = 2
        if (iiso == 2) iiso2 = 1
        fsw(iiso) = r2f(d2r(dsw(iiso),rstd(iiso)),d2r(dsw(iiso2),rstd(iiso2)))
        foc(iiso) = r2f(d2r(doc(iiso),rstd(iiso)),d2r(doc(iiso2),rstd(iiso2)))
        
        fswl(iiso,i1) = fsw(iiso)

        call o_iso_cw(  &
            & w, k1_loc(iiso),q,tau,poro,fsw(iiso),nz,foc(iiso),rstd(iiso),alfacw_loc(iiso),sat,e,tempcw_loc &! input
            & ,j1(iiso,i1),j2(iiso,i1),j3(iiso,i1),cdf(iiso,i1),shf(iiso,i1),dsf(iiso,i1)  &! output
            & ,z,fr(iiso,i1,:),fp(iiso,i1,:),frtop(iiso,i1)   &! output
            & ,ztot  &! inout
            & )
    enddo

end do

write(15,*) 0.0d0, 0.0d0  &
    & , (r2d(f2r(frtop(1,i1),frtop(2,i1)),rstd(1)), i1=1,n1)  
write(16,*) 0.0d0, 0.0d0  &
    & , (r2d(f2r(fswl(1,i1),fswl(2,i1)),rstd(1)), i1 = 1,n1)  

write(17,*) 0.0d0, 0.0d0  &
    & , (r2d(f2r(frtop(2,i1),frtop(1,i1)),rstd(2)), i1=1,n1)  
write(18,*) 0.0d0, 0.0d0  &
    & , (r2d(f2r(fswl(2,i1),fswl(1,i1)),rstd(2)), i1 = 1,n1)  

write(19,*) 0.0d0, 0.0d0  &
    & , (  &
        & r2dp(f2r(frtop(2,i1),frtop(1,i1)),rstd(2))   &
        & - ( lambda* r2dp(f2r(frtop(1,i1),frtop(2,i1)),rstd(1)) + gamma)  &
        & , i1=1,n1)  
write(20,*) 0.0d0, 0.0d0  &
    & , (  &
        & r2dp(f2r(fswl(2,i1),fswl(1,i1)),rstd(2))   &
        & - ( lambda* r2dp(f2r(fswl(1,i1),fswl(2,i1)),rstd(1)) + gamma)  &
        & , i1=1,n1)  

do iz = 1,nz
    write(15,*) z(iz), z(iz)/ztot &
    &  , (r2d(f2r(fr(1,i1,iz),fr(2,i1,iz)),rstd(1)),i1=1,n1) 
    write(16,*) z(iz), z(iz)/ztot &
    &  , (r2d(f2r(fp(1,i1,iz),fp(2,i1,iz)),rstd(1)),i1=1,n1) 
    
    write(17,*) z(iz), z(iz)/ztot &
    &  , (r2d(f2r(fr(2,i1,iz),fr(1,i1,iz)),rstd(2)),i1=1,n1) 
    write(18,*) z(iz), z(iz)/ztot &
    &  , (r2d(f2r(fp(2,i1,iz),fp(1,i1,iz)),rstd(2)),i1=1,n1) 
    
    write(19,*) z(iz), z(iz)/ztot &
    &  , (  &
        &  r2dp(f2r(fr(2,i1,iz),fr(1,i1,iz)),rstd(2))   &
        & - ( lambda* r2dp(f2r(fr(1,i1,iz),fr(2,i1,iz)),rstd(1)) + gamma)  &
        & ,i1=1,n1)
    write(20,*) z(iz), z(iz)/ztot &
    &  , (  &
        &  r2dp(f2r(fp(2,i1,iz),fp(1,i1,iz)),rstd(2))   &
        & - ( lambda* r2dp(f2r(fp(1,i1,iz),fp(2,i1,iz)),rstd(1)) + gamma)  &
        & ,i1=1,n1) 
end do

do iiso = 1,niso
    avx(iiso) = sum(dswl(iiso,:))/n1
    avy(iiso) = sum(j1(iiso,:))/n1
    sumx2(iiso) = sum((dswl(iiso,:)-avx(iiso))**2.0d0)
    sumy2(iiso) = sum((j1(iiso,:)-avy(iiso))**2.0d0)
    sumxy(iiso) = sum((dswl(iiso,:)-avx(iiso))*(j1(iiso,:)-avy(iiso)))
    slp(iiso) = sumxy(iiso)/sumx2(iiso)
    itcpt(iiso) = avy(iiso) - slp(iiso)*avx(iiso)
    crr(iiso) = sumxy(iiso)**2.0d0/sumx2(iiso)/sumy2(iiso)

    if (crr(iiso) <0.999d0) then
        print *, 'ERRRRRORRRRRE'
        print *, iiso,crr(iiso)
        pause
    end if 
    
    if (iiso == 1) then 
        write(11,*)  k1(iiso),tau,slp(iiso),itcpt(iiso),-itcpt(iiso)/slp(iiso) & 
        &  , merge(1.0d0,0.0d0, &
        &  (shf(iiso,n1)<=7.0d0).and. &
        &  (shf(iiso,imin)>=3.0d0))
    else
        write(12,*)  k1(iiso),tau,slp(iiso),itcpt(iiso),-itcpt(iiso)/slp(iiso) & 
        &  , merge(1.0d0,0.0d0, &
        &  (shf(iiso,n1)<=7.0d0).and. &
        &  (shf(iiso,imin)>=3.0d0))
    endif 
enddo 

close(11)
close(12)
close(15)
close(16)
close(17)
close(18)
close(19)
close(20)

endsubroutine o_iso_cw_sense_dsw

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine o_iso_cw_sense_dsw_Phan
implicit none

real(kind=8) :: w = 5d-5,q = 3d-1,poro = 0.2d0,tau = 1d6,sat = 0.5d0,e = 50d0, k1 = 10d0**(-8.5d0) 

integer, parameter :: n1 = 11,nz = 200, niso =2
real(kind=8),dimension(niso) :: doc,rstd,alfacw
real(kind=8),dimension(niso) :: slp,itcpt

real(kind=8),dimension(niso,n1) :: j1,j2,j3,dswl,fswl
real(kind=8),dimension(niso,n1) :: cdf,shf,dsf,frtop
real(kind=8),dimension(niso,n1,nz) :: fr,fp
real(kind=8),dimension(nz) :: z
real(kind=8),dimension(niso) :: avx,avy,sumx2,sumy2,sumxy,crr
real(kind=8),dimension(niso) :: dsw,fsw,foc
real(kind=8),dimension(niso) :: alfacw_loc,k1_loc
real(kind=8) ztot
integer i1, imin, iz, iiso, iiso2
character (100) workdir

real(kind=8) dp2d,mwl_Luz,r2f,d2r,r2dp,f2r,r2d,tc2d18orw,beta_eq_pack,beta_eq_sharp
real(kind=8) t2alfacw,t2beta_eq
real(kind=8) tempcw_loc,beta_eq, beta_kin
character (50) beta_eq_model, alfacw_model, d18orw_model(2)
! logical :: temp_sense = .true.
logical :: temp_sense = .false.
      
integer, parameter :: nry = 58, nrx = 25
real(kind=8) mry(nrx,nry), fl(nry), fa(nry), fsr(nry), fb, Drw(niso,2,nry)
real(kind=8) rco2(nry), temp(nry), fd(nry), tp(nry),frw(nry) 
integer iry
!---------------------------------------------------

open(30, file = './Royer2014.txt',status = 'old')
read(30,*) mry
close(30)
      
frw(:) = mry(5,:)
fl(:) = mry(6,:)
fa(:) = mry(7,:)
fsr(:) = mry(16,:)
rco2(:) = mry(18,:)
temp(:) = mry(19,:)
fd(:) = mry(11,:)
tp(:) = mry(1,:)
   
beta_eq_model = 'sharp16' 
alfacw_model = 'savin'
d18orw_model(1) = 'bowen08'
d18orw_model(2) = 'dansgaad64'

rstd(1) = 2.0052d-3 
rstd(2) = 3.799d-4 

doc(1) = 5.7d0
doc(2) = 2.86d0
! doc(2) = 3.0d0 

open (57,file='../oiso_output/GS-cw.txt',status='replace')

do iry = 1,nry
   
    tempcw_loc = temp(iry) 

    beta_eq = t2beta_eq(tempcw_loc+273.0d0,beta_eq_model)
    beta_kin = 1d0

    k1_loc(1) = k1
    k1_loc(2) = k1_loc(1)**beta_kin

    alfacw_loc(1) = t2alfacw(tempcw_loc+273.0d0,alfacw_model)
    alfacw_loc(2) = alfacw_loc(1)**beta_eq
    
    Drw(1,1,iry) = tc2d18orw(tempcw_loc,d18orw_model(1))
    Drw(1,2,iry) = tc2d18orw(tempcw_loc,d18orw_model(2))
    Drw(2,1,iry) = dp2d(mwl_Luz(Drw(1,1,iry),rstd(1)),rstd(2))
    Drw(2,2,iry) = dp2d(mwl_Luz(Drw(1,2,iry),rstd(1)),rstd(2))


    if (tp(iry) >= 380.0d0 ) then          !!! geocarb type assumption
        fb = rco2(iry)**0.5d0*fd(iry)/frw(iry)
    else if ((tp(iry) < 380.0d0).and.(tp(iry)>=350.0d0)) then
        fb = rco2(iry)**0.5d0*fd(iry)/frw(iry)*(tp(iry)-350.0d0)/30.0d0  &
            & +fd(iry)/frw(iry)*(2.0d0*rco2(iry)/(1.0d0+rco2(iry)))**0.4d0*(380.0d0-tp(iry))/30.0d0
    else if (tp(iry)<350.0d0) then
        fb = fd(iry)/frw(iry)*(2.0d0*rco2(iry)/(1.0d0+rco2(iry)))**0.4d0  
    end if

    w = 5d-5*frw(iry)
    tau = 10.0d0**(6.0d0)
    tau = tau*fb
    ztot = w*tau

    do i1 = 1, n1

        dsw(1) = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
        dswl(1,i1) = dsw(1)
        if (dsw(1) <= -8d0) imin = i1 ! imin wants to satisfy dswl(imin) = -8 
        
        dsw(2) = dp2d(mwl_Luz(dsw(1),rstd(1)),rstd(2))
        dswl(2,i1) = dsw(2)
            
        
        do iiso = 1, niso
            if (iiso == 1) iiso2 = 2
            if (iiso == 2) iiso2 = 1
            fsw(iiso) = r2f(d2r(dsw(iiso),rstd(iiso)),d2r(dsw(iiso2),rstd(iiso2)))
            foc(iiso) = r2f(d2r(doc(iiso),rstd(iiso)),d2r(doc(iiso2),rstd(iiso2)))
            
            fswl(iiso,i1) = fsw(iiso)

            call o_iso_cw(  &
                & w, k1_loc(iiso),q,tau,poro,fsw(iiso),nz,foc(iiso),rstd(iiso),alfacw_loc(iiso),sat,e,tempcw_loc &! input
                & ,j1(iiso,i1),j2(iiso,i1),j3(iiso,i1),cdf(iiso,i1),shf(iiso,i1),dsf(iiso,i1)  &! output
                & ,z,fr(iiso,i1,:),fp(iiso,i1,:),frtop(iiso,i1)   &! output
                & ,ztot  &! inout
                & )
        enddo

    end do

    j1 = j1*1.35d14*fa(iry)*max(0.0d0,(1.0d0-0.65d0*fl(iry)))

    do iiso = 1,niso
        avx(iiso) = sum(dswl(iiso,:))/n1
        avy(iiso) = sum(j1(iiso,:))/n1
        sumx2(iiso) = sum((dswl(iiso,:)-avx(iiso))**2.0d0)
        sumy2(iiso) = sum((j1(iiso,:)-avy(iiso))**2.0d0)
        sumxy(iiso) = sum((dswl(iiso,:)-avx(iiso))*(j1(iiso,:)-avy(iiso)))
        slp(iiso) = sumxy(iiso)/sumx2(iiso)
        itcpt(iiso) = avy(iiso) - slp(iiso)*avx(iiso)
        crr(iiso) = sumxy(iiso)**2.0d0/sumx2(iiso)/sumy2(iiso)

        if (crr(iiso) <0.999d0) then
            print *, 'ERRRRRORRRRRE'
            print *, iiso,crr(iiso)
            pause
        end if 
    enddo 

    write(57,*) slp(1), itcpt(1),Drw(1,1,iry), Drw(1,2,iry)  &
        & ,slp(2), itcpt(2),Drw(2,1,iry), Drw(2,2,iry)

enddo 

close(57)

endsubroutine o_iso_cw_sense_dsw_Phan

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine o_iso_cw_dr_Phan
implicit none

real(kind=8) :: w = 5d-5,q = 3d-1,poro = 0.2d0,tau = 1d6,sat = 0.5d0,e = 50d0, k1 = 10d0**(-8.5d0) 
real(kind=8) :: lambda = 0.5305d0, gamma = 0d0

integer, parameter :: n1 = 2,nz = 200, niso =2
real(kind=8),dimension(niso) :: doc,rstd,alfacw

real(kind=8),dimension(niso,n1) :: j1,j2,j3
real(kind=8),dimension(niso,n1) :: cdf,shf,dsf,frtop
real(kind=8),dimension(niso,n1,nz) :: fr,fp
real(kind=8),dimension(nz) :: z
real(kind=8),dimension(niso) :: dsw,fsw,foc
real(kind=8),dimension(niso) :: alfacw_loc,k1_loc
real(kind=8) ztot
integer i1, imin, iz, iiso, iiso2
character (100) workdir

real(kind=8) dp2d,mwl_Luz,r2f,d2r,r2dp,f2r,r2d,tc2d18orw,beta_eq_pack,beta_eq_sharp
real(kind=8) t2alfacw,t2beta_eq
real(kind=8) tempcw_loc,beta_eq, beta_kin
character (50) beta_eq_model, alfacw_model, d18orw_model(2)
! logical :: temp_sense = .true.
logical :: temp_sense = .false.
      
integer, parameter :: nry = 58, nrx = 25
real(kind=8) mry(nrx,nry), fl(nry), fa(nry), fsr(nry), fb, Drw(niso,2,nry)
real(kind=8) rco2(nry), temp(nry), fd(nry), tp(nry),frw(nry)
real(kind=8) d18o_max(4,nry),d18o_min(4,nry) 
integer iry
!---------------------------------------------------

open(30, file = './Royer2014.txt',status = 'old')
read(30,*) mry
close(30)

open(30, file = '../oiso_output/d18sw_dyn_max.txt',status = 'old')
read(30,*) d18o_max
close(30)

open(30, file = '../oiso_output/d18sw_dyn_min.txt',status = 'old')
read(30,*) d18o_min
close(30)
      
frw(:) = mry(5,:)
fl(:) = mry(6,:)
fa(:) = mry(7,:)
fsr(:) = mry(16,:)
rco2(:) = mry(18,:)
temp(:) = mry(19,:)
fd(:) = mry(11,:)
tp(:) = mry(1,:)
   
beta_eq_model = 'sharp16' 
alfacw_model = 'savin'
d18orw_model(1) = 'bowen08'
d18orw_model(2) = 'dansgaad64'

rstd(1) = 2.0052d-3 
rstd(2) = 3.799d-4 

doc(1) = 5.7d0
doc(2) = 2.86d0
! doc(2) = 3.0d0 

open (57,file='../oiso_output/d18r_cw_dyn.txt',status='replace')
open (58,file='../oiso_output/d17r_cw_dyn.txt',status='replace')
open (59,file='../oiso_output/capd17r_cw_dyn.txt',status='replace')

do iry = 1,nry
   
    tempcw_loc = temp(iry) 

    beta_eq = t2beta_eq(tempcw_loc+273.0d0,beta_eq_model)
    beta_kin = 1d0

    k1_loc(1) = k1
    k1_loc(2) = k1_loc(1)**beta_kin

    alfacw_loc(1) = t2alfacw(tempcw_loc+273.0d0,alfacw_model)
    alfacw_loc(2) = alfacw_loc(1)**beta_eq
    
    Drw(1,1,iry) = tc2d18orw(tempcw_loc,d18orw_model(1))
    Drw(1,2,iry) = tc2d18orw(tempcw_loc,d18orw_model(2))
    Drw(2,1,iry) = dp2d(mwl_Luz(Drw(1,1,iry),rstd(1)),rstd(2))
    Drw(2,2,iry) = dp2d(mwl_Luz(Drw(1,2,iry),rstd(1)),rstd(2))


    if (tp(iry) >= 380.0d0 ) then          !!! geocarb type assumption
        fb = rco2(iry)**0.5d0*fd(iry)/frw(iry)
    else if ((tp(iry) < 380.0d0).and.(tp(iry)>=350.0d0)) then
        fb = rco2(iry)**0.5d0*fd(iry)/frw(iry)*(tp(iry)-350.0d0)/30.0d0  &
            & +fd(iry)/frw(iry)*(2.0d0*rco2(iry)/(1.0d0+rco2(iry)))**0.4d0*(380.0d0-tp(iry))/30.0d0
    else if (tp(iry)<350.0d0) then
        fb = fd(iry)/frw(iry)*(2.0d0*rco2(iry)/(1.0d0+rco2(iry)))**0.4d0  
    end if

    w = 5d-5*frw(iry)
    tau = 10.0d0**(6.0d0)
    tau = tau*fb
    ztot = w*tau

    do i1 = 1, n1
        
        if (i1 == 1) dsw(1) = maxval(d18o_max(2:,iry))
        if (i1 == n1) dsw(1) = minval(d18o_min(2:,iry)) + minval(Drw(1,:,iry))
        
        dsw(2) = dp2d(mwl_Luz(dsw(1),rstd(1)),rstd(2))
        ! print *, dsw(1),dsw(2)
        
        do iiso = 1, niso
            if (iiso == 1) iiso2 = 2
            if (iiso == 2) iiso2 = 1
            fsw(iiso) = r2f(d2r(dsw(iiso),rstd(iiso)),d2r(dsw(iiso2),rstd(iiso2)))
            foc(iiso) = r2f(d2r(doc(iiso),rstd(iiso)),d2r(doc(iiso2),rstd(iiso2)))

            call o_iso_cw(  &
                & w, k1_loc(iiso),q,tau,poro,fsw(iiso),nz,foc(iiso),rstd(iiso),alfacw_loc(iiso),sat,e,tempcw_loc &! input
                & ,j1(iiso,i1),j2(iiso,i1),j3(iiso,i1),cdf(iiso,i1),shf(iiso,i1),dsf(iiso,i1)  &! output
                & ,z,fr(iiso,i1,:),fp(iiso,i1,:),frtop(iiso,i1)   &! output
                & ,ztot  &! inout
                & )
        enddo

    end do

    write(57,*) tp(iry),r2d(f2r(fr(1,1,1),fr(2,1,1)),rstd(1)),r2d(f2r(fr(1,2,1),fr(2,2,1)),rstd(1))
    write(58,*) tp(iry),r2d(f2r(fr(2,1,1),fr(1,1,1)),rstd(2)),r2d(f2r(fr(2,2,1),fr(1,2,1)),rstd(2))
    write(59,*) tp(iry),r2dp(f2r(fr(2,1,1),fr(1,1,1)),rstd(2))-(lambda*r2dp(f2r(fr(1,1,1),fr(2,1,1)),rstd(1)) +gamma) &
        & ,r2dp(f2r(fr(2,2,1),fr(1,2,1)),rstd(2))-(lambda*r2dp(f2r(fr(1,2,1),fr(2,2,1)),rstd(1)) +gamma)

enddo 

close(57)
close(58)
close(59)

endsubroutine o_iso_cw_dr_Phan

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine o_iso_cw_sense_all
implicit none

real(kind=8) :: w = 5d-5
real(kind=8) :: q = 3d-1
real(kind=8) :: poro = 0.20d0
real(kind=8) :: k1 = 10.0d0**(-8.50d0)
real(kind=8) :: dsw = 0.0d0
real(kind=8) :: tau = 1d6
real(kind=8) :: tempcw = 15.0d0
real(kind=8) :: alfacw =exp(25.0d0/1d3)
real(kind=8) :: doc = 5.7d0
real(kind=8) :: ztot = 20.0d0
real(kind=8) :: sat = 0.50d0
real(kind=8) :: e = 50.0d0
real(kind=8) :: rstd = 2.0052d-3 
real(kind=8) fsw,foc,frtop

integer, parameter :: nz = 200
integer, parameter :: n1 = 11
integer, parameter :: n2 = 101
integer, parameter :: n3 = 101

real(kind=8) z(nz),fr(nz),fp(nz)
real(kind=8) j1(n1,n2,n3),j2(n1,n2,n3),j3(n1,n2,n3),dswl(n1)
real(kind=8) cdf(n1,n2,n3),shf(n1,n2,n3),dsf(n1,n2,n3)

real(kind=8) avx,avy,sumx2,sumy2,sumxy,slp,itcpt,crr

integer iz, row, col, i1, i2, i3
integer rem
character (100) workdir
character (2) chr
!---------------------------------------------------      
workdir = '../oiso_output'

call system ('mkdir -p '//trim(adjustl(workdir)))

open (12, file =trim(adjustl(workdir))//'/'//'o_iso_cw_flux.csv',status = 'replace')
open (22, file =trim(adjustl(workdir))//'/'//'o_iso_cw_flux.dat',status = 'replace')
open (32, file =trim(adjustl(workdir))//'/'//'o_iso_cw_shift.dat',status = 'replace')

rem = 1

do i3 = 1, n3

    ! w = 5d-6*10.0d0**(2.0d0*(i3-1.0d0)/(n3-1.0d0))
    k1 = 1d-11*10.0d0**(5.0d0*(i3-1.0d0)/(n3-1.0d0))
    ! q = 3d-3*10.0d0**(4.0d0*(i3-1.0d0)/(n3-1.0d0))

    do i2 = 1, n2

        tau = 1d4*10.0d0**(2.0d0*(i2-1.0d0)/(n2-1.0d0))
        ! k1 = 1d-11*10.0d0**(5.0d0*(i2-1.0d0)/(n2-1.0d0))
        ! poro = 0.050d0+ (0.40d0*(i2-1.0d0)/(n2-1.0d0))

        do i1 = 1, n1

            dsw = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
            dswl(i1) = dsw

            fsw = (dsw/1d3+1.0d0)*rstd/((dsw/1d3+1.0d0)*rstd+1.0d0)
            foc = (doc/1d3+1.0d0)*rstd/((doc/1d3+1.0d0)*rstd+1.0d0)
            call o_iso_cw(  &
                & w, k1,q,tau,poro,fsw,nz,ztot,foc,rstd,alfacw,sat,e,tempcw &! input
                & ,j1(i1,i2,i3),j2(i1,i2,i3),j3(i1,i2,i3),cdf(i1,i2,i3),shf(i1,i2,i3),dsf(i1,i2,i3),z,fr,fp,frtop   &! output
                & ,ztot  &! inout
                & )
            
            
            write(chr,'(i2.2)') i1

            open (11, file=trim(adjustl(workdir))//'/'//'o_iso_cw'//chr//'.csv',status = 'replace')

            write(11,*) 0.0d0,',', 0.0d0 &
            &  ,',',  &
            &  1d3*(fsw/(1.0d0-fsw)/rstd-1.0d0) &
            &  ,',', &
            &  1d3*(frtop/(1.0d0-frtop)/rstd-1.0d0)

            do iz = 1,nz
                write(11,*) z(iz),',', z(iz)/ztot &
                &  ,',', &
                &  1d3*(fp(iz)/(1.0d0-fp(iz))/rstd-1.0d0) &
                &  ,',',  &
                &  1d3*(fr(iz)/(1.0d0-fr(iz))/rstd-1.0d0)
            end do

            close (11)

        end do

        avx = sum(dswl(:))/n1
        avy = sum(j1(:,i2,i3))/n1
        sumx2 = sum((dswl(:)-avx)**2.0d0)
        sumy2 = sum((j1(:,i2,i3)-avy)**2.0d0)
        sumxy = sum((dswl(:)-avx)*(j1(:,i2,i3)-avy))
        slp = sumxy/sumx2
        itcpt = avy - slp*avx
        crr = sumxy**2.0d0/sumx2/sumy2

        if (crr<0.999d0) then
            write(*,*) 'ERRRRRORRRRRE'
            pause
        end if 

        write(12,*)  &
        &  k1 &
        &  ,',', &  
        &  tau &
        &  ,',', &    
        &   slp &
        &  ,',',  &   
        &   itcpt &
        &  ,',',    & 
        &   -itcpt/slp

        if (rem == i3) then
            write(22,*)  k1, tau, slp,itcpt, -itcpt/slp
            write(32,*)  k1, tau  &
            &   , shf(7,i2,i3),dsf(7,i2,i3) &
            &   , shf(n1,i2,i3),dsf(n1,i2,i3) &
            &   , merge(1.0d0,0.0d0, &
            &    (shf(n1,i2,i3)<=7.0d0).and. &
            &    (shf(7,i2,i3)>=3.0d0))
        else if (rem /= i3) then
            write(22,*)  ''
            write(32,*)  ''
            write(22,*)  k1, tau, slp,itcpt, -itcpt/slp
            write(32,*)  k1, tau &
            &   , shf(7,i2,i3),dsf(7,i2,i3) &
            &   , shf(n1,i2,i3),dsf(n1,i2,i3) &
            &   , merge(1.0d0,0.0d0, &
            &    (shf(n1,i2,i3)<=7.0d0).and. &
            &    (shf(7,i2,i3)>=3.0d0))
            rem =i3
        end if 

    end do

end do 

close(12)
close(22)
close(32)

endsubroutine o_iso_cw_sense_all

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine o_iso_cw(  &
    & w, k1,q,tau,poro,fsw,nz,foc,rstd,alfacw,sat,e,tempcw &! input
    & ,j1,j2,j3,cdf,shf,dsf,z,fr,fp,frtop   &! output
    & ,ztot  &! inout
    & )
implicit none
real(kind=8),intent(in) :: w 
real(kind=8),intent(in) :: q 
real(kind=8),intent(in) :: poro 
real(kind=8),intent(in) :: fsw 
real(kind=8),intent(in) :: tau 
real(kind=8),intent(in) :: k1 
real(kind=8),intent(in) :: foc
real(kind=8),intent(in) :: rstd
real(kind=8),intent(in) :: alfacw
real(kind=8),intent(in) :: tempcw
real(kind=8),intent(in) :: sat
real(kind=8),intent(in) :: e
integer,intent(in) :: nz

real(kind=8),intent(out) :: j1,j2,j3
real(kind=8),intent(out) :: cdf,shf,dsf,frtop
real(kind=8),dimension(nz),intent(out) :: z,fr,fp

real(kind=8),intent(inout) :: ztot

real(kind=8) :: ms = 31.3d0
real(kind=8) :: mf = 55.6d0
real(kind=8) :: rg = 8.3d-3
real(kind=8) :: temp1 = 5.0d0

real(kind=8) kcw
real(kind=8) dz
real(kind=8) v
real(kind=8),dimension(nz) :: frx,fpx

real(kind=8) amx(2*nz,2*nz), ymx(2*nz), emx(2*nz)
integer info
real(kind=8) imbr
integer ipiv(2*nz)
external DGESV

real(kind=8) :: tol = 1d-6
real(kind=8) error 
integer iz, row, col

real(kind=8) densb
real(kind=8) :: densf = 1.0d3 ! kg/m3
real(kind=8) :: denss = 2.8d3

real(kind=8) error1, error2, error3

!---------------------------------------------------    

ztot = w*tau

kcw = k1*exp(-e*(1.0d0/(273.0d0+tempcw)-1.0d0/(273.0d0+temp1))/rg)

densb = (1.0d0-poro)*denss

do iz = 1,nz
    z(iz) = ztot*(iz) /(nz) 
end do

dz = z(2) - z(1)

fp = fsw
fr = foc

fpx = fp
frx = fr

error = 1d5

do while (error>tol)

    amx = 0.0d0
    ymx = 0.0d0

    do iz = 1, nz
        row = 2*(iz-1)+1
        if (iz/=nz) then
            amx(row, row) =  w*ms*(1.0d0-poro)*(-1.0d0)/dz*denss &
                -kcw*ms*mf*((1.0d0-fpx(iz))*1.0d0-1.0d0*alfacw*fpx(iz)*(-1.0d0))*densb

            amx(row,row + 2) = w*ms*(1.0d0-poro)*(1.0d0)/dz*denss

            ymx(row) = w*ms*(1.0d0-poro)*(frx(iz+1)-frx(iz))/dz*denss &
                -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)-1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb

        else if (iz==nz) then 
            amx(row, row) =  w*ms*(1.0d0-poro)*(-1.0d0)/dz*denss &
                -kcw*ms*mf*((1.0d0-fpx(iz))*1.0d0-1.0d0*alfacw*fpx(iz)*(-1.0d0))*densb


            ymx(row) = w*ms*(1.0d0-poro)*(foc-frx(iz))/dz*denss &
                -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)-1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb

        end if

        col = row  + 1
        amx(row, col) =  -kcw*ms*mf*((-1.0d0)*frx(iz)-1.0d0*alfacw*1.0d0*(1.0d0-frx(iz)))*densb

    end do

    do iz = 1, nz
        row = 2*(iz-1)+2
        if (iz/=1) then
            amx(row, row) =  -q*mf*(1.0d0)/dz*densf &
                -kcw*ms*mf*((-1.0d0)*frx(iz)-1.0d0*alfacw*1.0d0*(1.0d0-frx(iz)))*densb

            amx(row,row - 2) = -q*mf*(-1.0d0)/dz*densf

            ymx(row) = -q*mf*(fpx(iz)-fpx(iz-1))/dz*densf &
                -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)-1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb

        else if (iz==1) then 
            amx(row, row) =  -q*mf*(1.0d0)/dz*densf-kcw*ms*mf*((-1.0d0)*frx(iz) &
                -1.0d0*alfacw*1.0d0*(1.0d0-frx(iz)))*densb

            ymx(row) = -q*mf*(fpx(iz)-fsw)/dz*densf-kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz) &
                -1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb

        end if

        col = row  - 1
        amx(row, col) =  -kcw*ms*mf*((1.0d0-fpx(iz))*1.0d0-1.0d0*alfacw*fpx(iz)*(-1.0d0))*densb


    end do 

    ymx = -ymx

    call DGESV(2*nz,int(1),amx,2*nz,IPIV,ymx,2*nz,INFO) 

    emx = 0.0d0

    do iz = 1, nz
        row = 2*(iz-1)+1
        frx(iz) = frx(iz)+ymx(row)
        if (frx(iz)/=0.0d0) emx(row) = abs(ymx(row)/frx(iz))
    end do 

    do iz = 1, nz
        row = 2*(iz-1)+2
        fpx(iz) = fpx(iz)+ymx(row)
        if (fpx(iz)/=0.0d0) emx(row) = abs(ymx(row)/fpx(iz))
    end do 

    error = maxval(emx)

end do


frtop = foc-  q*densf*mf*(fsw - fpx(nz))/w/denss/ms/(1.0d0-poro)   ! calculation with mass balance

cdf=(1.0d0-fpx(1))*frx(1)/alfacw/fpx(1)/(1.0d0-frx(1))

j1=0.0d0
j2=0.0d0
j3=0.0d0
shf=0.0d0
dsf = 0.0d0

fr = frx
fp = fpx

j1 = w*ms*denss*(fr(nz)-fr(1))*(1.0d0-poro)
! j1 = w*ms*denss*(foc-frtop)*(1.0d0-poro)
j2 = j2+ 0.50d0*densb*ms*mf*kcw*dz*((1.0d0-fp(1))*fr(1) &
    -alfacw*(1.0d0-fr(1))*fp(1) +(1.0d0-fp(nz))*fr(nz)-alfacw*(1.0d0-fr(nz))*fp(nz))
j3 = -q*mf*densf*(fp(nz)-fp(1))
! j3 = -q*mf*densf*(fp(nz)-fsw)
shf= 1d3*(fr(1)/(1.0d0-fr(1))-fr(nz)/(1.0d0-fr(nz)))/rstd
dsf= 1d3*(fr(1)/(1.0d0-fr(1))/rstd-1.0d0)
do iz = 2, nz-1
    j2 = j2 + densb*ms*mf*kcw*dz*((1.0d0-fp(iz))*fr(iz)-alfacw*(1.0d0-fr(iz))*fp(iz))  
end do     

error1 = abs(j1-j2)/abs(j1)
error2 = abs(j3-j2)/abs(j2)
error3 = abs(j1-j3)/abs(j1)
if (any((/error1,error2,error3/)> tol)) then 
    write(*,*) 'error', error1, error2, error3
    ! write(*,*) j1,j2,j3
    ! pause
endif 

endsubroutine o_iso_cw
