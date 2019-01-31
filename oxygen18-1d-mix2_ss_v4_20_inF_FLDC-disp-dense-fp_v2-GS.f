      program oxygenisotope
      !!
      ! file is made from v3
      !!  this file tries to use explicit for spreading of crust
      !! use version 3 for implicit method
      !! only box 1 is open to ocean
      !! calculation of steady state values alone
      !! file is derived from oxygen18-1d-mix2.f
      !! R is replaced by F (18O/18O+16O)
      implicit none
      
      double precision :: o18wi = -10.0d0
      double precision rwi
      
      double precision :: o18mi = 7.0d0 
      double precision rmi
      
      double precision rp1i,rp2i,rp3i
      
      double precision :: L1 = 0.5d0 ! km
      double precision :: L2 = 1.5d0     
      double precision :: L3 = 4.0d0  
      
      integer, parameter :: nz = 200
      
      double precision z(nz)
      double precision dz
      
      integer iz
      
      double precision :: fsp = 3.0d0 ! km^2 yr^-10
      
      double precision :: rl = 1d5 ! km  ridge length
      double precision  w, L, q
      
      double precision :: T1 = 5.0d0 ! deg C
      double precision :: T2 = 200.0d0
      double precision :: T3 = 350.0d0  
      
      double precision :: kref = 10.0d0**(-8.50d0) ! yr^-10
      double precision k1, k2, k3
      
      double precision :: E1 = 50.0d0 ! kj mol^-1
      double precision :: E2 = 50.0d0
      double precision :: E3 = 50.0d0
      
      double precision :: Rg = 8.3d-3 ! kj mol^-1 K^-1
      
      double precision :: dense = 2.7d12 ! kg km^-3
      
      double precision :: tres = 1.0d6! yr
      
      double precision :: rsmow = 0.0020052d0
      
      double precision :: o18ri = 5.7d0
      double precision rri 
      
      double precision :: ms = 0.5d0/16.0d0*1.0d3 ! mol kg^-1
      double precision :: mw = 0.89d0/16.0d0*1.0d3
      
      double precision :: alfa1 = exp(27.0d0/1d3)
      double precision :: alfa2 = exp(7.0d0/1d3)
      double precision :: alfa3 = exp(2.5d0/1d3)
      
      double precision :: f1 = 1.0d16  ! kg yr^-1
      double precision :: f2 = 1.0d13
      double precision :: f3 = 1.0d13
      
      double precision :: p1 = 0.25d0
      double precision :: p2 = 0.05d0
      double precision :: p3 = 0.02d0
            
      integer it, row
      integer, parameter :: nt = 100  ! 1e5
      double precision dt
      double precision time(nt)
      
      double precision r1(nz),r2(nz),r3(nz),rw1(nz),rw2(nz),rw3(nz)
      double precision r1x(nz),r2x(nz),r3x(nz)
      double precision rw1x(nz),rw2x(nz),rw3x(nz)
      
      double precision rt(nz),rwt(nz),rtav,rwtav
      
      double precision r1av,r2av,r3av
      double precision rw1av,rw2av,rw3av
      
      double precision j1t,j2t,j3t
      
      double precision j1ts,j2ts,j3ts
      
      double precision jocn, jocn2, jocn3
      
      integer, parameter :: nmx = 3*2*(nz)
      
      integer ipiv(nmx)
      integer info
      
      external DGESV
      
      double precision amx1(nmx,nmx), ymx1(nmx)
      
      double precision ymx2(nmx)
      
      double precision errit, itl
      double precision :: tol = 1d-8
      
      integer, parameter :: nw = 1000
      integer :: spc = int(nt/nw)
      
      integer, parameter :: ns = 101
      integer is1, is2, is3, is4, is5
      integer, parameter :: nps1 = 11  ! o18wi
      integer, parameter :: nps2 = 1
      integer, parameter :: nps3 = 1
      integer, parameter :: nps4 = 1
      integer, parameter :: nps5 = 11
      
      double precision r1ss(nps1,nps2,nps3,nps4)
      double precision r2ss(nps1,nps2,nps3,nps4)
      double precision r3ss(nps1,nps2,nps3,nps4)
      double precision rw1ss(nps1,nps2,nps3,nps4)
      double precision rw2ss(nps1,nps2,nps3,nps4)
      double precision rw3ss(nps1,nps2,nps3,nps4)
      double precision j1o18(nps1,nps2,nps3,nps4)
      double precision j2o18(nps1,nps2,nps3,nps4)
      double precision j3o18(nps1,nps2,nps3,nps4)
      double precision j1o18s(nps1,nps2,nps3,nps4)
      double precision j2o18s(nps1,nps2,nps3,nps4)
      double precision j3o18s(nps1,nps2,nps3,nps4)
      
      double precision jocnss(nps1,nps2,nps3,nps4)
      
      double precision error
      double precision :: factor1 = 2.0d0
      double precision :: factor2 = 2.0d0  

      double precision :: volc = 1.0d0     
      
      integer, parameter :: nrec = 10
      integer reclis(nrec)
      character(2) chr
      character(256) chr2
      character(1) chr3
      integer irec, iter  
      
      double precision :: swad = 1.0d0    ! 1.0 if porewater is carried with rock by spreading
      double precision :: swsp = 1.0d0    ! 1.0 if water flows single path (lateral mixing)
      double precision :: swmtl = 0.0d0  !  1.0 if initial porewater has mantle isotopic comp. 
                                         !    or that of seawater in eq. with rock
                                          !  0.0 if seawater comp.
      double precision :: swbndsw = 1.0d0  ! 1.0  if initial iso. comp. is that of seawater
                                            !  in single flow path     
      
      double precision :: swdif = 1.0d0  !  switch for diffuse flow (vertical mixing)      
      
      character linebuf*256
      
      double precision avx,avy,sumx2,sumy2,sumxy,slp,itcpt,crr
      double precision do18x(nps1)
      
      double precision qz1(nz),qz2(nz),qz3(nz)
      double precision :: a = -2.0d0
      double precision b1,b2,b3
      
      double precision f1c,f2c,f3c
      
      double precision dsp1(nz),dsp2(nz),dsp3(nz)
      
      double precision :: densef = 1.0d12 ! kg km^-3
      double precision :: denses = 3.0d12 ! kg km^-3
      double precision dense1, dense2, dense3
      
      double precision :: beta = 0.876002291d0 ! anddesite (Zhao and Zheng, 2003)
      
      double precision :: csv = 1.0d0
      
      integer, parameter :: nry = 58, nrx = 25
      double precision mry(nrx,nry), fl(nry), fa(nry), fsr(nry)
      double precision rco2(nry), temp(nry), tp(nry)
      integer iry

      !---------------------------
            
      alfa1 = (6.673d6/(273.0d0+T1)**2.0d0+10.398d3/(273.0d0+T1)-4.78d0)
     $  *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+T1)))*beta
     $  -(2.194d6/(273.0d0+T1)**2.0d0+15.163d3/(273.0d0+T1)-4.72d0)
     $   +1.767d0*(2.0d0*beta-1.0d0)
            
      alfa2 = (6.673d6/(273.0d0+T2)**2.0d0+10.398d3/(273.0d0+T2)-4.78d0)
     $  *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+T2)))*beta
     $  -(2.194d6/(273.0d0+T2)**2.0d0+15.163d3/(273.0d0+T2)-4.72d0)
     $   +1.767d0*(2.0d0*beta-1.0d0)
            
      alfa3 = (6.673d6/(273.0d0+T3)**2.0d0+10.398d3/(273.0d0+T3)-4.78d0)
     $  *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+T3)))*beta
     $  -(2.194d6/(273.0d0+T3)**2.0d0+15.163d3/(273.0d0+T3)-4.72d0)
     $   +1.767d0*(2.0d0*beta-1.0d0)
      
      write(*,*) alfa1, alfa2, alfa3
      
      alfa1 = exp(alfa1/1d3)
      alfa2 = exp(alfa2/1d3)
      alfa3 = exp(alfa3/1d3)
            
      write(chr2,*) '-(w_pwad(oc)-(dif-sp)-'//
     $ 'dense-fp_20FLDC_dsp_v2-200_chk-b2)'
      
      open(30, file ='C:/Users/YK/Desktop/Isotope'//
     $   '/Royer2014.txt',status = 'old')
      read(30,*) mry
      close(30)
      
      fl(:) = mry(6,:)
      fa(:) = mry(7,:)
      fsr(:) = mry(16,:)
      rco2(:) = mry(18,:)
!       rco2(:) = mry(20,:)
      temp(:) = mry(19,:)
      tp(:) = mry(1,:)
      
!       write(*,*) fl
!       write(*,*) "---------------"
!       
!       write(*,*) fa
!       write(*,*) "---------------"
!       
!       write(*,*) fsr
!       write(*,*) "---------------"
!       
!       write(*,*) rco2
!       write(*,*) "---------------"
!       
      ! write(*,*) temp
      ! write(*,*) "---------------"
!       
      ! write(*,*) tp
      ! write(*,*) "---------------"
      
      open (57, file='C:/Users/YK/Desktop/Isotope/'//
     $  'GS-oca'//trim(adjustl(chr2))
     $    //'.csv', 
     $                             status='replace')
      
      do iry = 1, nry
!       iry = nry
      
      kref = 10.0d0**(-8.50d0)!*fsr(iry)!*rco2(iry)**0.50d0
      fsp = 3.0d0*fsr(iry)
      tres = 1d6/fsr(iry)
      
      ! T1 = 5.0d0 + 3.0d0*(temp(iry) -temp(nry))
      T1 = 5.0d0 + 1.9d0*(temp(iry) -temp(nry))
        
        f1=1d16*(fsp/3.0d0)!**0.50d0
        f2=1d13*(fsp/3.0d0)!**0.50d0
        f3=1d13*(fsp/3.0d0)!**0.50d0
      
      write(*,*) log10(kref),fsp,log10(tres), T1
      
      ! do is2 = nps2,1,-1
      is2= 1
      
      ! write(*,*) is2
      
      ! write(chr3,'(i1.1)') nps2-is2+1
      ! write(linebuf,*) '<Graph'//chr3//'>'
      
      ! write(71,*) trim(adjustl(linebuf))
      ! write(74,*) trim(adjustl(linebuf))
      
      
!       a = -2.0d0*10.0d0*10.0d0**(-2.0d0*(is2-1.0d0)/(nps2-1.0d0))
      ! f1 = 1d16*10.d0**((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-2.0d0))
!       swsp = ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0))
!       tres = 100.0d6*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-4.0d0))
        ! fsp = 1.0d0*9.0d0/
      ! $      (3.0d0**(2.0d0*(1.0d0*is2 - 1.0d0)/(1.0d0*nps2 - 1.0d0)))
!         fsp = 3d2*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-2.0d0))
        ! fsp = 3.0d0*32.0d0/
      ! $      (2.0d0**(5.0d0*(1.0d0*is2 - 1.0d0)/(1.0d0*nps2 - 1.0d0)))
        ! fsp = 3d0*10.d0**
      ! $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-0.5d0))
!       k1 = 1d-6*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-4.0d0))
!       rl = 1d4*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(2.0d0))
        
!         f1=1d16*1d1*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-1.0d0))
!         f2=1d13*1d1*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-1.0d0))
!         f3=1d13*1d1*10.d0**
!      $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-1.0d0))
        
        ! f1=1d16*(fsp/3.0d0)**0.50d0
        ! f2=1d13*(fsp/3.0d0)**0.50d0
        ! f3=1d13*(fsp/3.0d0)**0.50d0
       
        ! if (is2 == 1) then
          ! L1 = 0.50d0
          ! L2 = 1.50d0
          ! L3 = 4.0d0
        ! else 
          ! L1 = 0.50d0*2.2967d0*(10.d0**
      ! $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-0.5d0)))**(-1.463d0)
          ! L2 = 1.50d0*2.2967d0*(10.d0**
      ! $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-0.5d0)))**(-1.463d0)
          ! L3 = 4.0d0*2.2967d0*(10.d0**
      ! $   ((1.0d0*is2-1.0d0)/(1.0d0*nps2-1.0d0)*(-0.5d0)))**(-1.463d0)
        ! end if
      
      ! do is3 = nps3,1,-1
      is3 = 1
      ! f2 = 1d14*10.d0**((1.0d0*is3-1.0d0)/(1.0d0*nps3-1.0d0)*(-2.0d0))
      ! kref = 1d-6*10.d0**
      ! $   ((1.0d0*is3-1.0d0)/(1.0d0*nps3-1.0d0)*(-5.0d0))
!       T1 = 105.0d0 + ((1.0d0*is3-1.0d0)/(1.0d0*nps3-1.0d0)*(-100.0d0))
!       k1 = 10.0d0**(-8.50d0)*
!      $    exp(-E1*(1.0d0/(273.0d0+T1)-1.0d0/(273.0d0+5.0d0))/Rg)
!       alfa1 = exp((2.58d6*(273.0d0+T1)**(-2.0d0)-6.38337d0)/1d3)
      ! f1 = 1d16*10.d0**((1.0d0*is3-1.0d0)/(1.0d0*nps3-1.0d0)*(-2.0d0))
      ! p1 = 0.50d0/(2.0d0**(4.0d0*(1.0d0*is3 - 1.0d0)
      ! $  /(1.0d0*nps3 - 1.0d0)))
      ! p2 = 0.05d0/(2.0d0**(5.0d0*(1.0d0*is3 - 1.0d0)
      ! $  /(1.0d0*nps3 - 1.0d0)))
       ! T1 = 5.0d0 + ((1.0d0*is3-1.0d0)/(1.0d0*nps3-1.0d0)*(55.0d0))
        ! T2 = 200.0d0 + ((1.0d0*is3-1.0d0)/(1.0d0*nps3-1.0d0)*(150.0d0))
      
      ! if ((mod(is3-1,4) == 0).and.(is2==3)) then
      ! write(chr3,'(i1.1)') 2-(is3-1)/4+1
      ! write(linebuf,*) '<Graph'//chr3//'>'
      
      ! write(72,*) ''!trim(adjustl(linebuf))
      ! write(75,*) ''!trim(adjustl(linebuf))
      ! end if
      
      ! do is4 = nps4,1,-1
      is4 = 1
      ! f3 = 1d14*10.d0**((1.0d0*is4-1.0d0)/(1.0d0*nps4-1.0d0)*(-2.0d0))
      ! tres = 100.0d6*10.d0**
      ! $   ((1.0d0*is4-1.0d0)/(1.0d0*nps4-1.0d0)*(-4.0d0))
      ! tres = 10.0d0**(7.50d0)*10.d0**
      ! $   ((1.0d0*is4-1.0d0)/(1.0d0*nps4-1.0d0)*(-2.0d0))
       ! a = -2.0d0*10.0d0*10.0d0**(-2.0d0*(is4-1.0d0)/(nps4-1.0d0))
       
      ! p2 = 0.05d0/(2.0d0**(5.0d0*(1.0d0*is4 - 1.0d0)
      ! $  /(1.0d0*nps4 - 1.0d0)))
       
      ! p3 = 0.02d0/(2.0d0**(7.0d0*(1.0d0*is4 - 1.0d0)
      ! $  /(1.0d0*nps4 - 1.0d0)))
       
       ! T2 = 200.0d0 + ((1.0d0*is4-1.0d0)/(1.0d0*nps4-1.0d0)*(150.0d0))
       ! T3 = 350.0d0 + ((1.0d0*is4-1.0d0)/(1.0d0*nps4-1.0d0)*(250.0d0))
      
      ! if ((mod(is4-1,4) == 0).and.(is2==3).and.(is3==5)) then
      ! write(chr3,'(i1.1)') 2-(is4-1)/4+1
      ! write(linebuf,*) '<Graph'//chr3//'>'
      
      ! write(73,*) ''!trim(adjustl(linebuf))
      ! write(76,*) ''!trim(adjustl(linebuf))
      ! end if
      
      do is1 = nps1,1,-1
      
      o18wi = -1.0d0*(1.0*is1-1.0d0)/(nps1 - 1.0d0)*20.0d0
      do18x(is1) = o18wi
      
      rri = rsmow*(1d-3*o18ri + 1.0d0)/
     $   (1.0d0+rsmow*(1d-3*o18ri + 1.0d0))
      rwi = rsmow*(1d-3*o18wi + 1.0d0)/
     $   (1.0d0+rsmow*(1d-3*o18wi + 1.0d0))
      rmi = rsmow*(1d-3*o18mi + 1.0d0)/
     $    (1.0d0+rsmow*(1d-3*o18mi + 1.0d0))
      
!       write(*,*) rri,rwi
      
      rp1i = rwi
      rp2i = rwi*alfa2
      rp3i = rwi*alfa3
      
      k1 = kref*exp(-E1*(1.0d0/(273.0d0+T1)-1.0d0/(278.0d0))/Rg)
      k2 = kref*exp(-E2*(1.0d0/(273.0d0+T2)-1.0d0/(278.0d0))/Rg)
      k3 = kref*exp(-E3*(1.0d0/(273.0d0+T3)-1.0d0/(278.0d0))/Rg)
      
            
      alfa1 = (6.673d6/(273.0d0+T1)**2.0d0+10.398d3/(273.0d0+T1)-4.78d0)
     $  *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+T1)))*beta
     $  -(2.194d6/(273.0d0+T1)**2.0d0+15.163d3/(273.0d0+T1)-4.72d0)
     $   +1.767d0*(2.0d0*beta-1.0d0)
            
      alfa2 = (6.673d6/(273.0d0+T2)**2.0d0+10.398d3/(273.0d0+T2)-4.78d0)
     $  *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+T2)))*beta
     $  -(2.194d6/(273.0d0+T2)**2.0d0+15.163d3/(273.0d0+T2)-4.72d0)
     $   +1.767d0*(2.0d0*beta-1.0d0)
            
      alfa3 = (6.673d6/(273.0d0+T3)**2.0d0+10.398d3/(273.0d0+T3)-4.78d0)
     $  *EXP((1.0d0-beta)/(8.31441d-3*(273.0d0+T3)))*beta
     $  -(2.194d6/(273.0d0+T3)**2.0d0+15.163d3/(273.0d0+T3)-4.72d0)
     $   +1.767d0*(2.0d0*beta-1.0d0)
      
      ! write(*,*) alfa1, alfa2, alfa3
      
      alfa1 = exp(alfa1/1d3)
      alfa2 = exp(alfa2/1d3)
      alfa3 = exp(alfa3/1d3)
      
      ! f1 = 5.0d0*1.0d16/T1
      ! f2 = 200.0d0*1.0d13/T2
      ! f3 = 350.0d0*1.0d13/T3
      
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
      
!       write (*,*) L,z
      
      dz = z(2) - z(1)
      
      f1c = dz*(0.50d0*(qz1(1)+qz1(nz))+sum(qz1(2:nz-1)))*rl
      f2c = dz*(0.50d0*(qz2(1)+qz2(nz))+sum(qz2(2:nz-1)))*rl
      f3c = dz*(0.50d0*(qz3(1)+qz3(nz))+sum(qz3(2:nz-1)))*rl
      
      if(abs((f1-f1c)/f1) > 1.0d-2) then
        write(*,*) log10(f1),log10(f1c),abs((f1-f1c)/f1) *1d2
      end if
      
      if(abs((f2-f2c)/f2) > 1.0d-2) then
        write(*,*) log10(f2),log10(f2c),abs((f2-f2c)/f2) *1d2
      end if
      
      if(abs((f3-f3c)/f3) > 1.0d-2) then
        write(*,*) log10(f3),log10(f3c),abs((f3-f3c)/f3) *1d2 
      end if
      
      do iz = 1, nz
        dsp1(iz) = 1d-3*0.15d0*(1d3*
     $   (-0.00088d0*(z(iz)/w*1d-6)**2.0d0
     $     +0.322d0*(z(iz)/w*1d-6)+2.0d0))**0.61d0    ! Logngitudial macrodispersitivity (km)
     $     *qz1(iz)/p1/L1/densef/dz/dz*(
     $   -0.00088d0*(z(iz)/w*1d-6)**2.0d0
     $     +0.322d0*(z(iz)/w*1d-6)+2.0d0)
        dsp2(iz) = 1d-3*0.15d0*(1d3*
     $   (-0.00088d0*(z(iz)/w*1d-6)**2.0d0
     $     +0.322d0*(z(iz)/w*1d-6)+2.0d0))**0.61d0    ! Logngitudial macrodispersitivity (km)
     $     *qz2(iz)/p2/L2/densef/dz/dz*(
     $   -0.00088d0*(z(iz)/w*1d-6)**2.0d0
     $     +0.322d0*(z(iz)/w*1d-6)+2.0d0)
        dsp3(iz) = 1d-3*0.15d0*(1d3*
     $   (-0.00088d0*(z(iz)/w*1d-6)**2.0d0
     $     +0.322d0*(z(iz)/w*1d-6)+2.0d0))**0.61d0    ! Logngitudial macrodispersitivity (km)
     $     *qz3(iz)/p3/L3/densef/dz/dz*(
     $  -0.00088d0*(z(iz)/w*1d-6)**2.0d0
     $     +0.322d0*(z(iz)/w*1d-6)+2.0d0)
      end do
      
!       dt = dz/w
      
      do irec=1,nrec
        reclis(irec)=irec*nt/nrec
      end do
      
!       write(*,*) nmx
      
!       do is2 = 1,3
!         if (is2 == 2) then
!            k1 = k2
!            p1 = p2
!            f1 = f2
!            L1 = L2
!            alfa1 = alfa2
!         else if (is2 == 3) then
!            k1 = k3
!            p1 = p3
!            f1 = f3
!            L1 = L3
!            alfa1 = alfa3
!         end if
        
      
!       r1 = rri
!       rw1 = rwi *(1.0d0-swmtl) + rp1i*swmtl
!       r2 = rri
!       rw2 = rwi *(1.0d0-swmtl) + rp2i*swmtl
!       r3 = rri
!       rw3 = rwi *(1.0d0-swmtl) + rp3i*swmtl
!       rw1 = (rri-(ms+p1*mw/(1.0d0-p1))/(mw*k1))/alfa1
      
      it = 1
      
      irec = 0
      
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
      
!       write(*,*) r1x, rw1x
      
      errit = 1d4
      
      do while (errit > tol)
      
      amx1 = 0.0d0
      ymx1 = 0.0d0
      ymx2 = 0.0d0
      
!       write(*,*) ymx1
      
      do iz = 1, nz
      
        row = 1 + 6*(iz-1)
        
!         write(*,*) row
        
        if (iz == 1) then
          !-- box 1 ---
          amx1(row,row) =  -mw*k1*((1.0d0-rw1x(iz))
     $     -alfa1*(-1.0d0)*rw1x(iz))/(1.0d0-p1)*dense1/denses
     $        - w/dz/1.0d0
          ymx1(row) = 
     $          -w*(r1x(iz)-rri)/dz/1.0d0
     $           -mw*k1*((1.0d0-rw1x(iz))*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz))*rw1x(iz))
     $          /(1.0d0-p1)*dense1/denses
!        $        + w*(r1x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row, row+1) = -mw*k1*((-1.0d0)*r1x(iz)
     $     -alfa1*(1.0d0-r1x(iz)))
     $          /(1.0d0-p1)*dense1/denses
          
          ymx1(row + 1) = 
     $                +qz1(iz)*(rwi-rw1x(iz))/L1/densef/p1*swdif
     $             +qz2(iz)*(rw2x(iz)-rw1x(iz))/L1/densef/p1*swdif
     $         + dsp1(iz)*(rwi-rw1x(iz))*swsp
     $         + dsp1(iz+1)*(rw1x(iz+1)-rw1x(iz))*swsp
     $     -w*(rw1x(iz)-rwi)/dz/1.0d0*swad *(1.0d0-swmtl)
!      $     +w*(-rmi)/dz/1.0d0*swad *swmtl
     $     -w*(rw1x(iz)-rp1i)/dz/1.0d0*swad *swmtl
     $           +ms*k1*((1.0d0-rw1x(iz))*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz))*rw1x(iz))
     $          /(p1)*dense1/densef

          amx1(row+1,row + 1) = 
     $                +qz1(iz)*(-1.0d0)/L1/densef/p1*swdif
     $              +qz2(iz)*(-1.0d0)/L1/densef/p1*swdif
     $         + dsp1(iz)*(-1.0d0)*swsp
     $         + dsp1(iz+1)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad *(1.0d0-swmtl)
!      $     +w*(-rmi)/dz/1.0d0*swad *swmtl
     $     -w*(1.0d0)/dz/1.0d0*swad *swmtl
     $           +ms*k1*((-1.0d0)*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz)))
     $          /(p1)*dense1/densef
          amx1(row+1,row +7) = 
     $         + dsp1(iz+1)*(1.0d0)*swsp
          amx1(row+1,row +3) = 
     $              +qz2(iz)*(1.0d0)/L1/densef/p1*swdif
!      $    + w*(1.0d0)/dz/2.0d0*swad 
          amx1(row + 1, row) = 
     $           +ms*k1*((1.0d0-rw1x(iz))
     $            -alfa1*(-1.0d0)*rw1x(iz))
     $          /(p1)*dense1/densef
          
          !-- box 2 ---
          
          amx1(row+2,row+2) =  
     $        -mw*k2*((1.0d0-rw2x(iz))
     $          -alfa2*rw2x(iz)*(-1.0d0))
     $           /(1.0d0-p2)*dense2/denses
     $        - w/dz/1.0d0
          ymx1(row+2) =  
     $          -w*(r2x(iz)-rri)/dz/1.0d0
     $        -mw*k2*((1.0d0-rw2x(iz))*r2x(iz)
     $          -alfa2*rw2x(iz)*(1.0d0-r2x(iz)))
     $           /(1.0d0-p2)*dense2/denses
!        $        + w*(r2x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row+2, row+3) =
     $        -mw*k2*((-1.0d0)*r2x(iz)
     $          -alfa2*(1.0d0-r2x(iz)))
     $           /(1.0d0-p2)*dense2/denses
          
          ymx1(row + 3) = 
     $         +qz2(iz)*(rw1x(iz)-rw2x(iz))/L2/densef/p2*swdif
     $         +qz3(iz)*(rw3x(iz)-rw2x(iz))/L2/densef/p2*swdif
     $         + dsp2(iz)*(rwi-rw2x(iz))*swsp
     $         + dsp2(iz+1)*(rw2x(iz+1)-rw2x(iz))*swsp
     $     -w*(rw2x(iz)-rwi)/dz/1.0d0*swad *(1.0d0-swmtl)
     $     -w*(rw2x(iz)-rp2i)/dz/1.0d0*swad *swmtl
     $           +ms*k2*((1.0d0-rw2x(iz))*r2x(iz)
     $            -alfa2*(1.0d0-r2x(iz))*rw2x(iz))
     $          /(p2)*dense2/densef
          amx1(row+3,row + 3) =
     $         +qz2(iz)*(-1.0d0)/L2/densef/p2*swdif
     $         +qz3(iz)*(-1.0d0)/L2/densef/p2*swdif
     $         + dsp2(iz)*(-1.0d0)*swsp
     $         + dsp2(iz+1)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad *(1.0d0-swmtl)
     $     -w*(1.0d0)/dz/1.0d0*swad *swmtl
     $           +ms*k2*((-1.0d0)*r2x(iz)
     $            -alfa2*(1.0d0-r2x(iz)))
     $          /(p2)*dense2/densef
          amx1(row+3,row +9) = 
     $         + dsp2(iz+1)*(1.0d0)*swsp
          amx1(row+3,row +5) = 
     $         +qz3(iz)*(1.0d0)/L2/densef/p2*swdif
          amx1(row+3,row +1) = 
     $         +qz2(iz)*(1.0d0)/L2/densef/p2*swdif
          amx1(row + 3, row + 2) = 
     $           +ms*k2*((1.0d0-rw2x(iz))
     $            -alfa2*(-1.0d0)*rw2x(iz))
     $          /(p2)*dense2/densef
          
          !-- box 3 ---
          
          amx1(row+4,row+4) = 
     $          -w*(1.0d0)/dz/1.0d0
     $        -mw*k3*((1.0d0-rw3x(iz))
     $          -alfa3*rw3x(iz)*(-1.0d0))
     $           /(1.0d0-p3)*dense3/denses
          ymx1(row+4) = 
     $          -w*(r3x(iz)-rri)/dz/1.0d0
     $        -mw*k3*((1.0d0-rw3x(iz))*r3x(iz)
     $          -alfa3*rw3x(iz)*(1.0d0-r3x(iz)))
     $           /(1.0d0-p3)*dense3/denses
!        $        + w*(r3x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row+4, row+5) = 
     $        -mw*k3*((-1.0d0)*r3x(iz)
     $          -alfa3*(1.0d0-r3x(iz)))
     $           /(1.0d0-p3)*dense3/denses
          
          ymx1(row + 5) = 
     $         +qz3(iz)*(rw2x(iz)-rw3x(iz))/L3/densef/p3*swdif
     $         + dsp3(iz)*(rwi-rw3x(iz))*swsp
     $         + dsp3(iz+1)*(rw3x(iz+1)-rw3x(iz))*swsp
     $     -w*(rw3x(iz)-rwi)/dz/1.0d0*swad *(1.0d0-swmtl)
     $     -w*(rw3x(iz)-rp3i)/dz/1.0d0*swad *swmtl
     $           +ms*k3*((1.0d0-rw3x(iz))*r3x(iz)
     $            -alfa3*(1.0d0-r3x(iz))*rw3x(iz))
     $          /(p3)*dense3/densef
          amx1(row+5,row + 5) = 
     $         +qz3(iz)*(-1.0d0)/L3/densef/p3*swdif
     $         + dsp3(iz)*(-1.0d0)*swsp
     $         + dsp3(iz+1)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad *(1.0d0-swmtl)
     $     -w*(1.0d0)/dz/1.0d0*swad *swmtl
     $           +ms*k3*((-1.0d0)*r3x(iz)
     $            -alfa3*(1.0d0-r3x(iz)))
     $          /(p3)*dense3/densef
          amx1(row+5,row +11) = 
     $         + dsp3(iz+1)*(1.0d0)*swsp
          amx1(row+5,row +3) = 
     $         +qz3(iz)*(1.0d0)/L3/densef/p3*swdif
          amx1(row + 5, row + 4) = 
     $           +ms*k3*((1.0d0-rw3x(iz))
     $            -alfa3*(-1.0d0)*rw3x(iz))
     $          /(p3)*dense3/densef
       
        else if ((iz == nz)) then
            !-- box 1 ---
          amx1(row,row) =  -mw*k1*((1.0d0-rw1x(iz))
     $     -alfa1*(-1.0d0)*rw1x(iz))/(1.0d0-p1)*dense1/denses
     $        - w/dz/1.0d0
          amx1(row,row-6) = 
     $          -w*(-1.0d0)/dz/1.0d0
          ymx1(row) = 
     $          -w*(r1x(iz)-r1x(iz-1))/dz/1.0d0
     $           -mw*k1*((1.0d0-rw1x(iz))*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz))*rw1x(iz))
     $          /(1.0d0-p1)*dense1/denses
!        $        + w*(r1x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row, row+1) = -mw*k1*((-1.0d0)*r1x(iz)
     $     -alfa1*(1.0d0-r1x(iz)))
     $          /(1.0d0-p1)*dense1/denses
          
          ymx1(row + 1) = 
     $                +qz1(iz)*(rwi-rw1x(iz))/L1/densef/p1*swdif
     $            +qz2(iz)*(rw2x(iz)-rw1x(iz))/L1/densef/p1*swdif
     $         + dsp1(iz)*(rw1x(iz-1)-rw1x(iz))*swsp
     $     -w*(rw1x(iz)-rw1x(iz-1))/dz/1.0d0*swad 
     $           +ms*k1*((1.0d0-rw1x(iz))*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz))*rw1x(iz))
     $          /(p1)*dense1/densef

          amx1(row+1,row + 1) = 
     $                +qz1(iz)*(-1.0d0)/L1/densef/p1*swdif
     $              +qz2(iz)*(-1.0d0)/L1/densef/p1*swdif
     $         + dsp1(iz)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad 
     $           +ms*k1*((-1.0d0)*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz)))
     $          /(p1)*dense1/densef
          amx1(row+1,row - 5) = 
     $         + dsp1(iz)*(1.0d0)*swsp
     $     -w*(-1.0d0)/dz/1.0d0*swad 
!           amx1(row+1,row +7) = 
!      $         + f1/rl*(1.0d0)/L1/dense/p1/dz/1.0d0*swsp
          amx1(row+1,row +3) = 
     $              +qz2(iz)*(1.0d0)/L1/densef/p1*swdif
!      $    + w*(1.0d0)/dz/2.0d0*swad 
          amx1(row + 1, row) = 
     $           +ms*k1*((1.0d0-rw1x(iz))
     $            -alfa1*(-1.0d0)*rw1x(iz))
     $          /(p1)*dense1/densef
          
          !-- box 2 ---
          
          amx1(row+2,row+2) =  
     $        -mw*k2*((1.0d0-rw2x(iz))
     $          -alfa2*rw2x(iz)*(-1.0d0))
     $           /(1.0d0-p2)*dense2/denses
     $        - w/dz/1.0d0
          amx1(row+2,row-4) = 
     $          -w*(-1.0d0)/dz/1.0d0
          ymx1(row+2) =  
     $          -w*(r2x(iz)-r2x(iz-1))/dz/1.0d0
     $        -mw*k2*((1.0d0-rw2x(iz))*r2x(iz)
     $          -alfa2*rw2x(iz)*(1.0d0-r2x(iz)))
     $           /(1.0d0-p2)*dense2/denses
!        $        + w*(r2x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row+2, row+3) =
     $        -mw*k2*((-1.0d0)*r2x(iz)
     $          -alfa2*(1.0d0-r2x(iz)))
     $           /(1.0d0-p2)*dense2/denses
          
          ymx1(row + 3) = 
     $         +qz2(iz)*(rw1x(iz)-rw2x(iz))/L2/densef/p2*swdif
     $         +qz3(iz)*(rw3x(iz)-rw2x(iz))/L2/densef/p2*swdif
     $    + dsp2(iz)*(rw2x(iz-1)-rw2x(iz))*swsp
     $     -w*(rw2x(iz)-rw2x(iz-1))/dz/1.0d0*swad 
     $           +ms*k2*((1.0d0-rw2x(iz))*r2x(iz)
     $            -alfa2*(1.0d0-r2x(iz))*rw2x(iz))
     $          /(p2)*dense2/densef
          amx1(row+3,row + 3) =
     $         +qz2(iz)*(-1.0d0)/L2/densef/p2*swdif
     $         +qz3(iz)*(-1.0d0)/L2/densef/p2*swdif
     $         + dsp2(iz)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad 
     $           +ms*k2*((-1.0d0)*r2x(iz)
     $            -alfa2*(1.0d0-r2x(iz)))
     $          /(p2)*dense2/densef
          amx1(row+3,row - 3) =
     $         + dsp2(iz)*(1.0d0)*swsp
     $     -w*(-1.0d0)/dz/1.0d0*swad 
!           amx1(row+3,row +9) = 
!      $         + f2/rl*(1.0d0)/L2/dense/p2/dz/1.0d0*swsp
          amx1(row+3,row +5) = 
     $         +qz3(iz)*(1.0d0)/L2/densef/p2*swdif
          amx1(row+3,row +1) = 
     $         +qz2(iz)*(1.0d0)/L2/densef/p2*swdif
          amx1(row + 3, row + 2) = 
     $           +ms*k2*((1.0d0-rw2x(iz))
     $            -alfa2*(-1.0d0)*rw2x(iz))
     $          /(p2)*dense2/densef
          
          !-- box 3 ---
          
          amx1(row+4,row+4) = 
     $          -w*(1.0d0)/dz/1.0d0
     $        -mw*k3*((1.0d0-rw3x(iz))
     $          -alfa3*rw3x(iz)*(-1.0d0))
     $           /(1.0d0-p3)*dense3/denses
          amx1(row+4,row-2) = 
     $          -w*(-1.0d0)/dz/1.0d0
          ymx1(row+4) = 
     $          -w*(r3x(iz)-r3x(iz-1))/dz/1.0d0
     $        -mw*k3*((1.0d0-rw3x(iz))*r3x(iz)
     $          -alfa3*rw3x(iz)*(1.0d0-r3x(iz)))
     $           /(1.0d0-p3)*dense3/denses
!        $        + w*(r3x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row+4, row+5) = 
     $        -mw*k3*((-1.0d0)*r3x(iz)
     $          -alfa3*(1.0d0-r3x(iz)))
     $           /(1.0d0-p3)*dense3/denses
          
          ymx1(row + 5) = 
     $         +qz3(iz)*(rw2x(iz)-rw3x(iz))/L3/densef/p3*swdif
     $   + dsp3(iz)*(rw3x(iz-1)-rw3x(iz))*swsp
     $     -w*(rw3x(iz)-rw3x(iz-1))/dz/1.0d0*swad 
     $           +ms*k3*((1.0d0-rw3x(iz))*r3x(iz)
     $            -alfa3*(1.0d0-r3x(iz))*rw3x(iz))
     $          /(p3)*dense3/densef
          amx1(row+5,row + 5) = 
     $         +qz3(iz)*(-1.0d0)/L3/densef/p3*swdif
     $         + dsp3(iz)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad 
     $           +ms*k3*((-1.0d0)*r3x(iz)
     $            -alfa3*(1.0d0-r3x(iz)))
     $          /(p3)*dense3/densef
          amx1(row+5,row - 1) = 
     $         + dsp3(iz)*(1.0d0)*swsp
     $     -w*(-1.0d0)/dz/1.0d0*swad 
!           amx1(row+5,row +11) = 
!      $         + f3/rl*(1.0d0)/L3/dense/p3/dz/1.0d0*swsp
          amx1(row+5,row +3) = 
     $         +qz3(iz)*(1.0d0)/L3/densef/p3*swdif
          amx1(row + 5, row + 4) = 
     $           +ms*k3*((1.0d0-rw3x(iz))
     $            -alfa3*(-1.0d0)*rw3x(iz))
     $          /(p3)*dense3/densef
     
     
        
        else 
            !-- box 1 ---
          amx1(row,row) =  -mw*k1*((1.0d0-rw1x(iz))
     $     -alfa1*(-1.0d0)*rw1x(iz))/(1.0d0-p1)*dense1/denses
     $        - w/dz/1.0d0
          amx1(row,row-6) = 
     $          -w*(-1.0d0)/dz/1.0d0
          ymx1(row) = 
     $          -w*(r1x(iz)-r1x(iz-1))/dz/1.0d0
     $           -mw*k1*((1.0d0-rw1x(iz))*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz))*rw1x(iz))
     $          /(1.0d0-p1)*dense1/denses
!        $        + w*(r1x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row, row+1) = -mw*k1*((-1.0d0)*r1x(iz)
     $     -alfa1*(1.0d0-r1x(iz)))
     $          /(1.0d0-p1)*dense1/denses
          
          ymx1(row + 1) = 
     $                +qz1(iz)*(rwi-rw1x(iz))/L1/densef/p1*swdif
     $            +qz2(iz)*(rw2x(iz)-rw1x(iz))/L1/densef/p1*swdif
     $         + dsp1(iz)*(rw1x(iz-1)-rw1x(iz))*swsp
     $         + dsp1(iz+1)*(rw1x(iz+1)-rw1x(iz))*swsp
     $     -w*(rw1x(iz)-rw1x(iz-1))/dz/1.0d0*swad 
     $           +ms*k1*((1.0d0-rw1x(iz))*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz))*rw1x(iz))
     $          /(p1)*dense1/densef

          amx1(row+1,row + 1) = 
     $                +qz1(iz)*(-1.0d0)/L1/densef/p1*swdif
     $              +qz2(iz)*(-1.0d0)/L1/densef/p1*swdif
     $         + dsp1(iz)*(-1.0d0)*swsp
     $         + dsp1(iz+1)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad 
     $           +ms*k1*((-1.0d0)*r1x(iz)
     $            -alfa1*(1.0d0-r1x(iz)))
     $          /(p1)*dense1/densef
          amx1(row+1,row - 5) = 
     $         + dsp1(iz)*(1.0d0)*swsp
     $     -w*(-1.0d0)/dz/1.0d0*swad 
          amx1(row+1,row +7) = 
     $         + dsp1(iz+1)*(1.0d0)*swsp
          amx1(row+1,row +3) = 
     $         + qz2(iz)*(1.0d0)/L1/densef/p1*swdif
!      $    + w*(1.0d0)/dz/2.0d0*swad 
          amx1(row + 1, row) = 
     $           +ms*k1*((1.0d0-rw1x(iz))
     $            -alfa1*(-1.0d0)*rw1x(iz))
     $          /(p1)*dense1/densef
          
          !-- box 2 ---
          
          amx1(row+2,row+2) =  
     $        -mw*k2*((1.0d0-rw2x(iz))
     $          -alfa2*rw2x(iz)*(-1.0d0))
     $           /(1.0d0-p2)*dense2/denses
     $        - w/dz/1.0d0
          amx1(row+2,row-4) = 
     $          -w*(-1.0d0)/dz/1.0d0
          ymx1(row+2) =  
     $          -w*(r2x(iz)-r2x(iz-1))/dz/1.0d0
     $        -mw*k2*((1.0d0-rw2x(iz))*r2x(iz)
     $          -alfa2*rw2x(iz)*(1.0d0-r2x(iz)))
     $           /(1.0d0-p2)*dense2/denses
!        $        + w*(r2x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row+2, row+3) =
     $        -mw*k2*((-1.0d0)*r2x(iz)
     $          -alfa2*(1.0d0-r2x(iz)))
     $           /(1.0d0-p2)*dense2/denses
          
          ymx1(row + 3) = 
     $         +qz2(iz)*(rw1x(iz)-rw2x(iz))/L2/densef/p2*swdif
     $         +qz3(iz)*(rw3x(iz)-rw2x(iz))/L2/densef/p2*swdif
     $    + dsp2(iz)*(rw2x(iz-1)-rw2x(iz))*swsp
     $    + dsp2(iz+1)*(rw2x(iz+1)-rw2x(iz))*swsp
     $     -w*(rw2x(iz)-rw2x(iz-1))/dz/1.0d0*swad 
     $           +ms*k2*((1.0d0-rw2x(iz))*r2x(iz)
     $            -alfa2*(1.0d0-r2x(iz))*rw2x(iz))
     $          /(p2)*dense2/densef
          amx1(row+3,row + 3) =
     $         +qz2(iz)*(-1.0d0)/L2/densef/p2*swdif
     $         +qz3(iz)*(-1.0d0)/L2/densef/p2*swdif
     $         + dsp2(iz)*(-1.0d0)*swsp
     $         + dsp2(iz+1)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad 
     $           +ms*k2*((-1.0d0)*r2x(iz)
     $            -alfa2*(1.0d0-r2x(iz)))
     $          /(p2)*dense2/densef
          amx1(row+3,row - 3) =
     $         + dsp2(iz)*(1.0d0)*swsp
     $     -w*(-1.0d0)/dz/1.0d0*swad 
          amx1(row+3,row +9) = 
     $         + dsp2(iz+1)*(1.0d0)*swsp
          amx1(row+3,row +5) = 
     $         +qz3(iz)*(1.0d0)/L2/densef/p2*swdif
          amx1(row+3,row +1) = 
     $         +qz2(iz)*(1.0d0)/L2/densef/p2*swdif
          amx1(row + 3, row + 2) = 
     $           +ms*k2*((1.0d0-rw2x(iz))
     $            -alfa2*(-1.0d0)*rw2x(iz))
     $          /(p2)*dense2/densef
          
          !-- box 3 ---
          
          amx1(row+4,row+4) = 
     $          -w*(1.0d0)/dz/1.0d0
     $        -mw*k3*((1.0d0-rw3x(iz))
     $          -alfa3*rw3x(iz)*(-1.0d0))
     $           /(1.0d0-p3)*dense3/denses
          amx1(row+4,row-2) = 
     $          -w*(-1.0d0)/dz/1.0d0
          ymx1(row+4) = 
     $          -w*(r3x(iz)-r3x(iz-1))/dz/1.0d0
     $        -mw*k3*((1.0d0-rw3x(iz))*r3x(iz)
     $          -alfa3*rw3x(iz)*(1.0d0-r3x(iz)))
     $           /(1.0d0-p3)*dense3/denses
!        $        + w*(r3x(iz)-rri)/dz/1.0d0
!           amx1(row,row+2) = w/dz/1.0d0
          amx1(row+4, row+5) = 
     $        -mw*k3*((-1.0d0)*r3x(iz)
     $          -alfa3*(1.0d0-r3x(iz)))
     $           /(1.0d0-p3)*dense3/denses
          
          ymx1(row + 5) = 
     $         +qz3(iz)*(rw2x(iz)-rw3x(iz))/L3/densef/p3*swdif
     $   + dsp3(iz)*(rw3x(iz-1)-rw3x(iz))*swsp
     $   + dsp3(iz+1)*(rw3x(iz+1)-rw3x(iz))*swsp
     $     -w*(rw3x(iz)-rw3x(iz-1))/dz/1.0d0*swad 
     $           +ms*k3*((1.0d0-rw3x(iz))*r3x(iz)
     $            -alfa3*(1.0d0-r3x(iz))*rw3x(iz))
     $          /(p3)*dense3/densef
          amx1(row+5,row+5) = 
     $         +qz3(iz)*(-1.0d0)/L3/densef/p3*swdif
     $         + dsp3(iz)*(-1.0d0)*swsp
     $         + dsp3(iz+1)*(-1.0d0)*swsp
     $     -w*(1.0d0)/dz/1.0d0*swad 
     $           +ms*k3*((-1.0d0)*r3x(iz)
     $            -alfa3*(1.0d0-r3x(iz)))
     $          /(p3)*dense3/densef
          amx1(row+5,row - 1) = 
     $         + dsp3(iz)*(1.0d0)*swsp
     $     -w*(-1.0d0)/dz/1.0d0*swad 
          amx1(row+5,row +11) = 
     $         + dsp3(iz+1)*(1.0d0)*swsp
          amx1(row+5,row +3) = 
     $         +qz3(iz)*(1.0d0)/L3/densef/p3*swdif
          amx1(row + 5, row + 4) = 
     $           +ms*k3*((1.0d0-rw3x(iz))
     $            -alfa3*(-1.0d0)*rw3x(iz))
     $          /(p3)*dense3/densef
     
     
     
        
        end if
        
      end do
      
!       write(*,*) ymx1
      
      ymx1 = -1.0d0*ymx1
      
!       open (22, file='mx1.txt', 
!      $                             status='replace')
!       open (23, file='mx2.txt', 
!      $                             status='replace')
!         
!       do iz = 1, 2*(nz-1)
!         write (22,*) (amx1(iz,is), is = 1,2*(nz-1))
!         write (23,*) ymx1(iz)
!       end do
!         
!       close(22)
!       close(23)
      
      call DGESV(nmx,int(1),amx1,nmx,IPIV,ymx1,nmx,INFO)

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
      
      
!       write (*,*) 
!      $  (1d3*(rw1(iz)/rsmow/(1.0d0-rw1(iz))-1.0d0),iz=1,nz,40)
!       write (*,*) 
!      $  (1d3*(r1(iz)/rsmow/(1.0d0-r1(iz))-1.0d0),iz=1,nz,40)
!       write (*,*) 
!      $  (1d3*(rw2(iz)/rsmow/(1.0d0-rw2(iz))-1.0d0),iz=1,nz,40)
!       write (*,*) 
!      $  (1d3*(r2(iz)/rsmow/(1.0d0-r2(iz))-1.0d0),iz=1,nz,40)
!       write (*,*) 
!      $  (1d3*(rw3(iz)/rsmow/(1.0d0-rw3(iz))-1.0d0),iz=1,nz,40)
!       write (*,*) 
!      $  (1d3*(r3(iz)/rsmow/(1.0d0-r3(iz))-1.0d0),iz=1,nz,40)
      
      errit = maxval(ymx2)
      
      write (*,*) errit, info
      
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
      
!       write(*,*) r1av(it),r2av(it),r3av(it)
      
      do iz = 1,  nz
        j1ts = j1ts 
     $   + ms*L1*(1.0d0-p1)*w*denses*(merge(rri,r1(iz-1),iz==1)- r1(iz))
     $   + ms*L1*dense1*dz*mw*k1*(r1(iz)*(1.0d0-rw1(iz))
     $   -alfa1*rw1(iz)*(1.0d0-r1(iz)))
        
        j1t = j1t
     $   +mw*L1*densef*w*p1*(rw1(iz-1)-rw1(iz)) *swad
     $  - mw*L1*dense1*dz*ms*k1*(r1(iz)*(1.0d0-rw1(iz))
     $     - alfa1*rw1(iz)*(1.0d0-r1(iz)))
     $  +swdif*f1/rl*dz*mw/L*(rwi-rw1(iz))
     $  +swdif*f2/rl*dz*mw/L*(rw2(iz)-rw1(iz))
     $  + swsp * f1/rl*mw*((merge(rw1(iz+1),rwi,iz/=nz))-rw1(iz))
     
        j2ts = j2ts 
     $   + ms*L2*(1.0d0-p2)*w*denses*(r2(iz-1) - r2(iz))
     $   + ms*L2*dense2*dz*mw*k2*(r2(iz)*(1.0d0-rw2(iz))
     $   -alfa2*rw2(iz)*(1.0d0-r2(iz)))
        
        j2t = j2t
     $   +mw*L2*densef*w*p2*(rw2(iz-1)-rw2(iz)) *swad
     $  - mw*L2*dense2*dz*ms*k2*(r2(iz)*(1.0d0-rw2(iz))
     $     - alfa2*rw2(iz)*(1.0d0-r2(iz)))
     $  +swdif*f2/rl*dz*mw/L*(rw1(iz)-rw2(iz))
     $  +swdif*f3/rl*dz*mw/L*(rw3(iz)-rw2(iz))
     $  + swsp * f2/rl*mw*((merge(rw2(iz+1),
     $   rw1(iz)*(1.0d0-swbndsw)+rwi*swbndsw,iz/=nz))-rw2(iz))
     
        j3ts = j3ts 
     $   + ms*L3*(1.0d0-p3)*w*denses*(r3(iz-1) - r3(iz))
     $   + ms*L3*dense3*dz*mw*k3*(r3(iz)*(1.0d0-rw3(iz))
     $   -alfa3*rw3(iz)*(1.0d0-r3(iz)))
        
        j3t = j3t
     $   +mw*L3*densef*w*p3*(rw3(iz-1)-rw3(iz)) *swad
     $  - mw*L3*dense3*dz*ms*k3*(r3(iz)*(1.0d0-rw3(iz)) 
     $  - alfa3*rw3(iz)*(1.0d0-r3(iz)))
     $  +swdif*f3/rl*dz*mw/L*(rw2(iz)-rw3(iz))
     $  + swsp * f3/rl*mw*((merge(rw3(iz+1),
     $       rw2(iz)*(1.0d0-swbndsw)+rwi*swbndsw,iz/=nz))-rw3(iz))
        
        jocn = jocn 
     $  +swdif*qz1(iz)*dz*mw*(rwi-rw1(iz))
     $  + swsp * dz*dsp1(iz)*mw*densef*L1*p1*
     $           ((merge(rw1(iz-1),rwi,iz/=1))-rw1(iz))
     $  + swsp * dz*merge(0.0d0,
     $       dsp1(iz+1)*mw*densef*L1*p1*(rw1(iz+1)-rw1(iz))
     $  , iz .eq. nz)
     $  + swsp * dz*dsp2(iz)*mw*densef*L2*p2*
     $        ((merge(rw2(iz-1),rwi,iz/=1))-rw2(iz))
     $  + swsp * dz*merge(0.0d0,
     $       dsp2(iz+1)*mw*densef*L2*p2*(rw2(iz+1)-rw2(iz))
     $  , iz .eq. nz)
     $  + swsp * dz*dsp3(iz)*mw*densef*L3*p3*
     $        ((merge(rw3(iz-1),rwi,iz/=1))-rw3(iz))
     $  + swsp * dz*merge(0.0d0,
     $       dsp3(iz+1)*mw*densef*L3*p3*(rw3(iz+1)-rw3(iz))
     $  , iz .eq. nz)
        
        jocn3 = jocn3 
     $  + mw*L3*densef*w*p3*(merge(rwi,rw3(iz-1),iz==1)-rw3(iz)) *swad
     $  + mw*L2*densef*w*p2*(merge(rwi,rw2(iz-1),iz==1)-rw2(iz)) *swad
     $  + mw*L1*densef*w*p1*(merge(rwi,rw1(iz-1),iz==1)-rw1(iz)) *swad
        
        jocn2 = jocn2 
     $  + ms*L1*(1.0d0-p1)*w*denses*(merge(rri,r1(iz-1),iz==1) - r1(iz))
     $  + ms*L2*(1.0d0-p2)*w*denses*(merge(rri,r2(iz-1),iz==1) - r2(iz))
     $  + ms*L3*(1.0d0-p3)*w*denses*(merge(rri,r3(iz-1),iz==1) - r3(iz))
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
      
      r1ss(is1,is2,is3,is4) = r1av
      r2ss(is1,is2,is3,is4) = r2av
      r3ss(is1,is2,is3,is4) = r3av
      rw1ss(is1,is2,is3,is4) = rw1av
      rw2ss(is1,is2,is3,is4) = rw2av
      rw3ss(is1,is2,is3,is4) = rw3av
      jocnss(is1,is2,is3,is4) = jocn + jocn3 ! porewater included
      
      if (abs((jocn2-(jocn+jocn3))/jocn2) > 1d-3) then
        write (*,*) 'error',abs((jocn2-(jocn+jocn3))/jocn2)
     $  ,'@','(',is1,',',is2,',',is3,',',is4,')'     
      end if 
      
      end do
      
      avx = sum(do18x(:))/nps1
      avy = sum(jocnss(:,is2,is3,is4))/nps1
      sumx2 = sum((do18x(:)-avx)**2.0d0)
      sumy2 = sum((jocnss(:,is2,is3,is4)-avy)**2.0d0)
      sumxy = sum((do18x(:)-avx)*(jocnss(:,is2,is3,is4)-avy))
      slp = sumxy/sumx2
      itcpt = avy - slp*avx
      crr = sumxy**2.0d0/sumx2/sumy2
      
      if (crr < 0.999d0) then 
        write (*,*) 'error', crr,is2,is3,is4
      end if
      
      write(57,*) 
     $     slp
     $    , ',', 
     $     itcpt
     $    , ',', 
     $     -itcpt/slp
      
      
      end do
     
      close(57)
      
      end
      
      
      subroutine del_spaces(s)
      character (*), intent (inout) :: s
      character (len=len(s)) tmp
      integer i, j
      j = 1
      do i = 1, len(s)
        if (s(i:i)==' ') cycle
          tmp(j:j) = s(i:i)
          j = j + 1
      end do
      s = tmp(1:j-1)
      end subroutine del_spaces