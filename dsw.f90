program dsw_main

use o_iso_cw_mod
use o_iso_ha_mod

implicit none

call o_iso_cw_sense_dsw_prof 
call o_iso_ha_sense_dsw_prof

call o_iso_ha_sense_dsw_Phan
call o_iso_cw_sense_dsw_Phan

call dsw_ss
call dsw_dyn

call o_iso_cw_dr_Phan
call o_iso_ha_drdp_Phan

endprogram dsw_main

!###################################################################################################
!###################################################################################################
!###################################################################################################

subroutine dsw_ss
implicit none

integer, parameter :: nry = 58, nrx = 25, nrx_ha = 4, nrx_cw = 8
integer iry
real(kind=8) :: lambda = 0.5305, gamma = 0d0
real(kind=8) mry(nrx,nry),tp(nry), dsw(2,3,nry)
real(kind=8) cw_gs(nrx_cw,nry), ha_gs(nrx_ha,nry)
real(kind=8) rstd(2)

real(kind=8) d2dp
! ----------------------------------------------------------
rstd(1) = 2.0052d-3 
rstd(2) = 3.799d-4 

open(30, file = './Royer2014.txt',status = 'old')
read(30,*) mry
close(30)

tp(:) = mry(1,:)  ! age

open(30, file = '../oiso_output/GS-oca.txt',status = 'old')
read(30,*) ha_gs
close(30)

open(30, file = '../oiso_output/GS-cw.txt',status = 'old')
read(30,*) cw_gs
close(30)

dsw(1,1,:) = -(ha_gs(2,:)+cw_gs(2,:))/(ha_gs(1,:)+cw_gs(1,:))
dsw(1,2,:) = -(ha_gs(2,:)+cw_gs(2,:)+cw_gs(1,:)*cw_gs(3,:))/(ha_gs(1,:)+cw_gs(1,:))
dsw(1,3,:) = -(ha_gs(2,:)+cw_gs(2,:)+cw_gs(1,:)*cw_gs(4,:))/(ha_gs(1,:)+cw_gs(1,:))

dsw(2,1,:) = -(ha_gs(4,:)+cw_gs(6,:))/(ha_gs(3,:)+cw_gs(1,:))
dsw(2,2,:) = -(ha_gs(4,:)+cw_gs(6,:)+cw_gs(5,:)*cw_gs(7,:))/(ha_gs(3,:)+cw_gs(5,:))
dsw(2,3,:) = -(ha_gs(4,:)+cw_gs(6,:)+cw_gs(5,:)*cw_gs(8,:))/(ha_gs(3,:)+cw_gs(5,:))

open(30, file = '../oiso_output/d18sw_ss.txt',status = 'replace')
open(31, file = '../oiso_output/d17sw_ss.txt',status = 'replace')
open(32, file = '../oiso_output/capd17sw_ss.txt',status = 'replace')
do iry = 1,nry
    write(30,*) tp(iry), dsw(1,:,iry) 
    write(31,*) tp(iry), dsw(2,:,iry) 
    write(32,*) tp(iry), d2dp(dsw(2,1,iry),rstd(2))-(lambda*d2dp(dsw(1,1,iry),rstd(1))+gamma)   &
        & , d2dp(dsw(2,2,iry),rstd(2))-(lambda*d2dp(dsw(1,2,iry),rstd(1))+gamma) &
        & , d2dp(dsw(2,3,iry),rstd(2))-(lambda*d2dp(dsw(1,3,iry),rstd(1))+gamma) 
enddo
close(30)
close(31)
close(32)

endsubroutine dsw_ss

!###################################################################################################
!###################################################################################################
!###################################################################################################

subroutine dsw_dyn
implicit none 

real(kind=8) :: MO = 1.5d24/18.0d0

integer, parameter :: nry = 58, nrx_ha = 4, nrx_cw = 8, nrx_ss = 4
integer iry, idrw
real(kind=8) :: lambda = 0.5305, gamma = 0d0
real(kind=8) tp(nry), dsw(2,3,nry)
real(kind=8) cw_gs(nrx_cw,nry), ha_gs(nrx_ha,nry), d18o_ss(nrx_ss,nry), d17o_ss(nrx_ss,nry)
real(kind=8) rstd(2),dswi_min(2),dswi_max(2), MxO(2,3,nry)

real(kind=8) d2dp,dp2d,mwl_Luz,d2r,r2f,f2r,r2d
! ----------------------------------------------------------
rstd(1) = 2.0052d-3 
rstd(2) = 3.799d-4 

open(30, file = '../oiso_output/d17sw_ss.txt',status = 'old')
read(30,*) d17o_ss
close(30)

open(30, file = '../oiso_output/d18sw_ss.txt',status = 'old')
read(30,*) d18o_ss
close(30)

tp(:) = d18o_ss(1,:)  ! age

open(30, file = '../oiso_output/GS-oca.txt',status = 'old')
read(30,*) ha_gs
close(30)

open(30, file = '../oiso_output/GS-cw.txt',status = 'old')
read(30,*) cw_gs
close(30)

dswi_min(1) = min(-10d0,minval(d18o_ss(2:,1)))
dswi_max(1) = max(-10d0,maxval(d18o_ss(2:,1)))


dswi_min(2) = dp2d(mwl_Luz(dswi_min(1),rstd(1)),rstd(2))
dswi_max(2) = dp2d(mwl_Luz(dswi_max(1),rstd(1)),rstd(2))
! dswi_min(2) = min(dp2d(mwl_Luz(dswi_min(1),rstd(1)),rstd(2)),minval(d17o_ss(2:,1)))
! dswi_max(2) = max(dp2d(mwl_Luz(dswi_max(1),rstd(1)),rstd(2)),maxval(d17o_ss(2:,1)))

! start with min d18O

open(30, file = '../oiso_output/d18sw_dyn_min.txt',status = 'replace')
open(31, file = '../oiso_output/d17sw_dyn_min.txt',status = 'replace')
open(32, file = '../oiso_output/capd17sw_dyn_min.txt',status = 'replace')

MxO(1,:,1) = MO*r2f(d2r(dswi_min(1),rstd(1)),d2r(dswi_min(2),rstd(2)))
MxO(2,:,1) = MO*r2f(d2r(dswi_min(2),rstd(2)),d2r(dswi_min(1),rstd(1)))

dsw(1,:,1) = dswi_min(1)
dsw(2,:,1) = dswi_min(2)

do iry = 2, nry
    MxO(1,1,iry) = MxO(1,1,iry-1) + 10d6*( &
        & ha_gs(2,iry-1) + dsw(1,1,iry-1)*ha_gs(1,iry-1) &
        & + cw_gs(2,iry-1) + dsw(1,1,iry-1)*cw_gs(1,iry-1) &
        )
    MxO(1,2,iry) = MxO(1,2,iry-1) + 10d6*( &
        & ha_gs(2,iry-1) + dsw(1,2,iry-1)*ha_gs(1,iry-1) &
        & + cw_gs(2,iry-1) + (dsw(1,2,iry-1)+cw_gs(3,iry-1))*cw_gs(1,iry-1) &
        )
    MxO(1,3,iry) = MxO(1,3,iry-1) + 10d6*( &
        & ha_gs(2,iry-1) + dsw(1,3,iry-1)*ha_gs(1,iry-1) &
        & + cw_gs(2,iry-1) + (dsw(1,3,iry-1)+cw_gs(4,iry-1))*cw_gs(1,iry-1) &
        )
        
    MxO(2,1,iry) = MxO(2,1,iry-1) + 10d6*( &
        & ha_gs(4,iry-1) + dsw(2,1,iry-1)*ha_gs(3,iry-1) &
        & + cw_gs(6,iry-1) + dsw(2,1,iry-1)*cw_gs(5,iry-1) &
        )
    MxO(2,2,iry) = MxO(2,2,iry-1) + 10d6*( &
        & ha_gs(4,iry-1) + dsw(2,2,iry-1)*ha_gs(3,iry-1) &
        & + cw_gs(6,iry-1) + (dsw(2,2,iry-1)+cw_gs(7,iry-1))*cw_gs(5,iry-1) &
        )
    MxO(2,3,iry) = MxO(2,3,iry-1) + 10d6*( &
        & ha_gs(4,iry-1) + dsw(2,3,iry-1)*ha_gs(3,iry-1) &
        & + cw_gs(6,iry-1) + (dsw(2,3,iry-1)+cw_gs(8,iry-1))*cw_gs(5,iry-1) &
        )
    do idrw = 1, 3
        dsw(1,idrw,iry) = r2d(f2r(MxO(1,idrw,iry)/MO,MxO(2,idrw,iry)/MO),rstd(1))
        dsw(2,idrw,iry) = r2d(f2r(MxO(2,idrw,iry)/MO,MxO(1,idrw,iry)/MO),rstd(2))
    enddo
enddo

do iry = 1,nry
    write(30,*) tp(iry), dsw(1,:,iry) 
    write(31,*) tp(iry), dsw(2,:,iry) 
    write(32,*) tp(iry), d2dp(dsw(2,1,iry),rstd(2))-(lambda*d2dp(dsw(1,1,iry),rstd(1))+gamma)   &
        & , d2dp(dsw(2,2,iry),rstd(2))-(lambda*d2dp(dsw(1,2,iry),rstd(1))+gamma) &
        & , d2dp(dsw(2,3,iry),rstd(2))-(lambda*d2dp(dsw(1,3,iry),rstd(1))+gamma) 
enddo

close(30)
close(31)
close(32)


! then with max d18O

open(30, file = '../oiso_output/d18sw_dyn_max.txt',status = 'replace')
open(31, file = '../oiso_output/d17sw_dyn_max.txt',status = 'replace')
open(32, file = '../oiso_output/capd17sw_dyn_max.txt',status = 'replace')

MxO(1,:,1) = MO*r2f(d2r(dswi_max(1),rstd(1)),d2r(dswi_max(2),rstd(2)))
MxO(2,:,1) = MO*r2f(d2r(dswi_max(2),rstd(2)),d2r(dswi_max(1),rstd(1)))

dsw(1,:,1) = dswi_max(1)
dsw(2,:,1) = dswi_max(2)

do iry = 2, nry
    MxO(1,1,iry) = MxO(1,1,iry-1) + 10d6*( &
        & ha_gs(2,iry-1) + dsw(1,1,iry-1)*ha_gs(1,iry-1) &
        & + cw_gs(2,iry-1) + dsw(1,1,iry-1)*cw_gs(1,iry-1) &
        )
    MxO(1,2,iry) = MxO(1,2,iry-1) + 10d6*( &
        & ha_gs(2,iry-1) + dsw(1,2,iry-1)*ha_gs(1,iry-1) &
        & + cw_gs(2,iry-1) + (dsw(1,2,iry-1)+cw_gs(3,iry-1))*cw_gs(1,iry-1) &
        )
    MxO(1,3,iry) = MxO(1,3,iry-1) + 10d6*( &
        & ha_gs(2,iry-1) + dsw(1,3,iry-1)*ha_gs(1,iry-1) &
        & + cw_gs(2,iry-1) + (dsw(1,3,iry-1)+cw_gs(4,iry-1))*cw_gs(1,iry-1) &
        )
        
    MxO(2,1,iry) = MxO(2,1,iry-1) + 10d6*( &
        & ha_gs(4,iry-1) + dsw(2,1,iry-1)*ha_gs(3,iry-1) &
        & + cw_gs(6,iry-1) + dsw(2,1,iry-1)*cw_gs(5,iry-1) &
        )
    MxO(2,2,iry) = MxO(2,2,iry-1) + 10d6*( &
        & ha_gs(4,iry-1) + dsw(2,2,iry-1)*ha_gs(3,iry-1) &
        & + cw_gs(6,iry-1) + (dsw(2,2,iry-1)+cw_gs(7,iry-1))*cw_gs(5,iry-1) &
        )
    MxO(2,3,iry) = MxO(2,3,iry-1) + 10d6*( &
        & ha_gs(4,iry-1) + dsw(2,3,iry-1)*ha_gs(3,iry-1) &
        & + cw_gs(6,iry-1) + (dsw(2,3,iry-1)+cw_gs(8,iry-1))*cw_gs(5,iry-1) &
        )
    do idrw = 1, 3
        dsw(1,idrw,iry) = r2d(f2r(MxO(1,idrw,iry)/MO,MxO(2,idrw,iry)/MO),rstd(1))
        dsw(2,idrw,iry) = r2d(f2r(MxO(2,idrw,iry)/MO,MxO(1,idrw,iry)/MO),rstd(2))
    enddo
enddo

do iry = 1,nry
    write(30,*) tp(iry), dsw(1,:,iry) 
    write(31,*) tp(iry), dsw(2,:,iry) 
    write(32,*) tp(iry), d2dp(dsw(2,1,iry),rstd(2))-(lambda*d2dp(dsw(1,1,iry),rstd(1))+gamma)   &
        & , d2dp(dsw(2,2,iry),rstd(2))-(lambda*d2dp(dsw(1,2,iry),rstd(1))+gamma) &
        & , d2dp(dsw(2,3,iry),rstd(2))-(lambda*d2dp(dsw(1,3,iry),rstd(1))+gamma) 
enddo

close(30)
close(31)
close(32)

endsubroutine dsw_dyn