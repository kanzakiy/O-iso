
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++          FUNCTIONS        ++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function r2d(ratio,rstd)
implicit none
real(kind=8) r2d,ratio,rstd
r2d=(ratio/rstd-1d0)*1d3
endfunction r2d

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function d2r(delta,rstd)
implicit none
real(kind=8) d2r,delta,rstd
d2r=(delta*1d-3+1d0)*rstd
endfunction d2r

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function r2dp(ratio,rstd)
implicit none
real(kind=8) r2dp,ratio,rstd
r2dp=1d3*log(ratio/rstd)
endfunction r2dp

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dp2r(dp,rstd)
implicit none
real(kind=8) dp2r,dp,rstd
dp2r=rstd*exp(dp*1d-3)
endfunction dp2r

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dp2d(dp,rstd)
implicit none
real(kind=8) dp2d,dp,rstd,dp2r, r2d
dp2d = r2d(dp2r(dp,rstd),rstd)  
endfunction dp2d

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function d2dp(d,rstd)
implicit none
real(kind=8) d2dp,d,rstd,d2r, r2dp
d2dp = r2dp(d2r(d,rstd),rstd)  
endfunction d2dp

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function r2f(r1,r2)
implicit none
real(kind=8) r2f,r1,r2
r2f = r1/(1d0+r1+r2)
endfunction r2f

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function f2r(f1,f2)
implicit none
real(kind=8) f2r,f1,f2
f2r = f1/(1d0-f1-f2)
endfunction f2r

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function mwl_Luz(d18o,rstd)  ! returning d'17 based on d18O and meteoric water line by Luz and Barkan (2010) 
implicit none
real(kind=8) mwl_Luz,d18o,d2r,r2dp,rstd
mwl_Luz = 0.528d0*r2dp(d2r(d18o,rstd),rstd) + 0.033d0
endfunction mwl_Luz

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function beta_eq_pack(t)  ! equilibirum exponent for 17/16O vs. 18/16O fractionation suggested by Pack and Herwartz (2014) 
implicit none
real(kind=8) beta_eq_pack,t  ! t [K]
beta_eq_pack = -740d0/t/t + 0.5305d0
endfunction beta_eq_pack

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function beta_eq_sharp(t)  ! equilibirum exponent for 17/16O vs. 18/16O fractionation suggested by Sharp et al. (2016) 
implicit none
real(kind=8) beta_eq_sharp,t ! t [K]
beta_eq_sharp = -1.85d0/t + 0.5305d0
endfunction beta_eq_sharp

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function t2beta_eq(t,model)  ! equilibirum exponent for 17/16O vs. 18/16O fractionation suggested by Sharp et al. (2016) 
implicit none
real(kind=8) t2beta_eq,t ! t [K]
character(50) model

if (trim(adjustl(model)) == 'sharp16') then 
    t2beta_eq = -1.85d0/t + 0.5305d0
elseif (trim(adjustl(model)) == 'pack14') then 
    t2beta_eq = -740d0/t/t + 0.5305d0
endif 

endfunction t2beta_eq

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function alfacw_savin(t)  ! equilibirum exponent for 17/16O vs. 18/16O fractionation suggested by Sharp et al. (2016) 
implicit none
real(kind=8) alfacw_savin,t ! t [K]
alfacw_savin = exp((2.5d6*(t)**(-2.0d0)-2.87d0)/1d3)   ! (Savin&Lee)
endfunction alfacw_savin

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function tc2d18orw(tc,model)  ! converting temp [C] to rainwater d18O based on modern empirical relationship
implicit none
real(kind=8) tc2d18orw,tc ! tc [C]
character(50) model

if (trim(adjustl(model)) == 'bindeman19_lin') then 
    tc2d18orw = 0.297d0*tc-11.64d0
elseif (trim(adjustl(model)) == 'bindeman19_quod') then 
    tc2d18orw = -0.0123d0*tc**2d0+0.64d0*tc - 13.05d0
elseif (trim(adjustl(model)) == 'dansgaad64') then 
    tc2d18orw = 0.69d0*tc - 13.6d0
elseif (trim(adjustl(model)) == 'bowen08') then 
    tc2d18orw = 0.299d0*tc - 12.4d0
endif 

endfunction tc2d18orw

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function t2alfacw(t,model)  ! converting temp [K] to alpha 
implicit none
real(kind=8) t2alfacw,t ! t [K]
character(50) model

if (trim(adjustl(model)) == 'savin') then 
    t2alfacw = exp((2.5d6*(t)**(-2.0d0)-2.87d0)/1d3) 
endif 

endfunction t2alfacw

!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************

function alpha_zhao03(tc,beta)
implicit none
real(kind=8) alpha_zhao03, tc, beta
alpha_zhao03 = (6.673d6/(273.0d0+tc)**2.0d0+10.398d3/(273.0d0+tc)-4.78d0) &
    & *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+tc)))*beta &
    & -(2.194d6/(273.0d0+tc)**2.0d0+15.163d3/(273.0d0+tc)-4.72d0)+1.767d0*(2.0d0*beta-1.0d0)
alpha_zhao03 = exp(alpha_zhao03/1d3)
endfunction alpha_zhao03

!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************

function k_arrhenius(kref,E,tc)
implicit none
real(kind=8) k_arrhenius, E, tc,kref
real(kind=8) :: Rg = 8.3d-3
k_arrhenius = kref*exp(-E*(1.0d0/(273.0d0+tc)-1.0d0/(278.0d0))/Rg)
endfunction k_arrhenius

!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************
!**************************************************************************************************************************************************************