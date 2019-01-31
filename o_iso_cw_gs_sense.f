      program o_iso_cw
      
      double precision :: w = 5d-5
      double precision :: q = 3d-1
      double precision v
      double precision :: ms = 31.3d0
      double precision :: mf = 55.6d0
      double precision :: poro = 0.20d0
      double precision :: sat = 0.50d0
      
      double precision :: k1 = 10.0d0**(-8.50d0)
      double precision :: e = 50.0d0
      double precision :: rg = 8.3d-3
      
      double precision kcw
      double precision :: temp1 = 5.0d0
      double precision :: tempcw = 15.0d0
      
      double precision :: alfacw =exp(25.0d0/1d3)
      
      double precision :: doc = 5.7d0
      
      double precision :: dsw = 0.0d0
      
      double precision :: ztot = 20.0d0
      double precision :: tau = 1d6
      integer, parameter :: nz = 200
      double precision z(nz),fr(nz),fp(nz),frx(nz),fpx(nz)
      double precision dz
      
      integer, parameter :: nmx = 2*(nz)
      
      double precision amx(nmx,nmx), ymx(nmx), emx(nmx)
      integer info
      double precision imbr
      integer ipiv(nmx)
      external DGESV
      
      double precision :: tol = 1d-6
      double precision error 
      
      integer iz, row, col, i1, i2, i3
      
      double precision :: rstd = 2.0052d-3 
      double precision fsw,foc
      
      double precision densb
      double precision :: densf = 1.0d3 ! kg/m3
      double precision :: denss = 2.8d3
      
      integer, parameter :: n1 = 11
      integer, parameter :: n2 = 1!9
      integer, parameter :: n3 = 1!9
      
      double precision j1(n1,n2,n3),j2(n1,n2,n3),j3(n1,n2,n3),dswl(n1)
      double precision cdf(n1,n2,n3)
      
      character (2) chr
      
      double precision avx,avy,sumx2,sumy2,sumxy,slp,itcpt,crr
      
      integer, parameter :: nry = 58, nrx = 25
      double precision mry(nrx,nry), fl(nry), fa(nry), fsr(nry)
      double precision rco2(nry), temp(nry), fd(nry), tp(nry),frw(nry) 
      double precision fsv(nry),fv(nry)
      integer iry
      double precision rwt(3,nry), d18Oc(3,nry)
      
      integer is
C        integer, parameter :: ns = 8
C       integer, parameter :: ns = 3  !  Fb
C       integer, parameter :: ns = 2  !  fl
C       integer, parameter :: ns = 5 !  kref
      integer, parameter :: ns = 6 !  kref
      
      double precision :: kref = 10.0d0**(-8.50d0)
      
      integer :: flg_ff = 0
      integer :: flg_kref = 0
      integer :: flg_std = 1
      integer :: flg_uplift = 0
      integer :: flg_fb = 0
      integer :: flg_fl = 0
      
      double precision :: wstd = 5d-5
      double precision :: fb
!---------------------------------------------------   
      
      open(30, file = 'C:/Users/YK/Desktop/Isotope/'//
     $ 'Royer2014.txt',status = 'old')
      read(30,*) mry
      close(30)
      
      open (80, file ='C:/Users/YK/Desktop/Isotope/'//
     $ 'o_iso_cw_flux_sense.txt',status = 'replace')
      
      do is = 1, ns
      
      
      frw(:) = mry(5,:)
      fl(:) = mry(6,:)
      fa(:) = mry(7,:)
      fsr(:) = mry(16,:)
      rco2(:) = mry(18,:)
      temp(:) = mry(19,:)
      fd(:) = mry(11,:)
      tp(:) = mry(1,:)
      
      fsv(:) = mry(24,:)
      fv(:) = mry(25,:)
      
      ! write(*,*) fv
      ! write(*,*) fsv
      
      if (flg_ff == 1) then 
      
      if (is ==1) then
      
      write(*,*) is
      
      ! frw(:) = 1.0d0
      ! fl(:) = 1.0d0
      ! fa(:) = 1.0d0
      ! fd(:) = 1.0d0
      ! fsr(:) = 1.0d0
      ! rco2(:) = 1.0d0
      ! temp(:) = mry(19,nry)
      
      else if (is ==2) then
      
      write(*,*) is
      
      ! frw(:) = 1.0d0
      fl(:) = 1.0d0
      fa(:) = 1.0d0
      fd(:) = 1.0d0
      fsr(:) = 1.0d0
      rco2(:) = 1.0d0
      temp(:) = mry(19,nry)
      
      else if (is ==3) then
      
      write(*,*) is
      
      frw(:) = 1.0d0
      ! fl(:) = 1.0d0
      fa(:) = 1.0d0
      fd(:) = 1.0d0
      fsr(:) = 1.0d0
      rco2(:) = 1.0d0
      temp(:) = mry(19,nry)
      
      else if (is ==4) then
      
      write(*,*) is
      
      frw(:) = 1.0d0
      fl(:) = 1.0d0
      ! fa(:) = 1.0d0
      fd(:) = 1.0d0
      fsr(:) = 1.0d0
      rco2(:) = 1.0d0
      temp(:) = mry(19,nry)
      
      else if (is ==5) then
      
      write(*,*) is
      
      frw(:) = 1.0d0
      fl(:) = 1.0d0
      fa(:) = 1.0d0
      ! fd(:) = 1.0d0
      fsr(:) = 1.0d0
      rco2(:) = 1.0d0
      temp(:) = mry(19,nry)
      
      else if (is ==6) then
      
      write(*,*) is
      
      frw(:) = 1.0d0
      fl(:) = 1.0d0
      fa(:) = 1.0d0
      fd(:) = 1.0d0
      ! fsr(:) = 1.0d0
      rco2(:) = 1.0d0
      temp(:) = mry(19,nry)
      
      else if (is ==7) then
      
      write(*,*) is
      
      frw(:) = 1.0d0
      fl(:) = 1.0d0
      fa(:) = 1.0d0
      fd(:) = 1.0d0
      fsr(:) = 1.0d0
      ! rco2(:) = 1.0d0
      temp(:) = mry(19,nry)
      
      else if (is ==8) then
      
      write(*,*) is
      
      frw(:) = 1.0d0
      fl(:) = 1.0d0
      fa(:) = 1.0d0
      fd(:) = 1.0d0
      fsr(:) = 1.0d0
      rco2(:) = 1.0d0
      ! temp(:) = mry(19,nry)
      
      end if
      
      end if
      
      if (flg_kref == 1) then 
      
      write(*,*) is
      kref = 10.0d0**(-9.0d0)*10.0d0**(1.0d0*(is-1.0d0)/(ns-1.0d0))
      
      end if
      
      if (flg_uplift == 1) then 
      
      write(*,*) is
      wstd = 5d-6*10.0d0**(2.0d0*(is-1.0d0)/(ns-1.0d0))
      
      end if
      
      if (flg_std == 1) then 
      
      write(*,*) is
      
      end if
      
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
!       write(*,*) temp
!       write(*,*) "---------------"
!       
      ! write(*,*) tp
      ! write(*,*) "---------------"
      
      
!       do i3 = 1, n3
      i3 = 1
      
      do iry = 1, nry
      
      
!       w = 5d-6*10.0d0**(2.0d0*(i3-1.0d0)/(n3-1.0d0))
!       k1 = 1d-11*10.0d0**(5.0d0*(i3-1.0d0)/(n3-1.0d0))
      
      k1 = kref
      
      if (tp(iry) >= 380.0d0 ) then          !!! geocarb type assumption
        fb = rco2(iry)**0.5d0*fd(iry)/frw(iry)
      else if ((tp(iry) < 380.0d0).and.(tp(iry)>=350.0d0)) then
        fb = rco2(iry)**0.5d0*fd(iry)/frw(iry)
     $    *(tp(iry)-350.0d0)/30.0d0
     $    +fd(iry)/frw(iry)
     $   *(2.0d0*rco2(iry)/(1.0d0+rco2(iry)))**0.4d0  
     $    *(380.0d0-tp(iry))/30.0d0
      else if (tp(iry)<350.0d0) then
        fb = fd(iry)/frw(iry)
     $   *(2.0d0*rco2(iry)/(1.0d0+rco2(iry)))**0.4d0 
      end if
      
      if (flg_fb == 1) then
      if (is == 2) then 
      fb = rco2(iry)**0.5d0*fd(iry) /frw(iry)  
      else if (is == 3) then 
      fb = rco2(iry)**0.25d0*fd(iry) /frw(iry)
      end if 
      end if
      
!       do i2 = 1, n2
      i2 = 1
      
!       tau = 1d4*10.0d0**(2.0d0*(i2-1.0d0)/(n2-1.0d0))
!       k1 = 1d-11*10.0d0**(5.0d0*(i2-1.0d0)/(n2-1.0d0))

      
      ! w = 5d-5!*fsr(iry)
      w = wstd*frw(iry)
      
      ! ztot = 30.0d0
      ! tau = ztot/w
      
      ! tau = 10.0d0**(5.75d0)
      ! tau = 10.0d0**(5.5d0)
      tau = 10.0d0**(6.0d0)
      ztot = w*tau
      ztot = ztot*fb
      
      alfacw =exp((2.5d6*(temp(iry)+273.0d0)**(-2.0d0)-2.87d0)/1d3)   ! (Savin&Lee)
      
      kcw = k1*exp(-e*(1.0d0/(273.0d0+temp(iry))
     $    -1.0d0/(273.0d0+temp1))/rg)
      
!       write(*,*) iry, w, tau, alfacw, kcw
      
!       ztot = w*tau
      
      do i1 = 1, n1
      
      dsw = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
      
      densb = (1.0d0-poro)*denss
      
      do iz = 1,nz
        z(iz) = ztot*(iz) /(nz) 
      end do
      
      dz = z(2) - z(1)
      
      fsw = (dsw/1d3+1.0d0)*rstd/((dsw/1d3+1.0d0)*rstd+1.0d0)
      foc = (doc/1d3+1.0d0)*rstd/((doc/1d3+1.0d0)*rstd+1.0d0)
      
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
        amx(row, row) =  w*ms*(1.0d0-poro)*(-1.0d0)/dz*denss
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*1.0d0
     $  -1.0d0*alfacw*fpx(iz)*(-1.0d0))*densb
        
        amx(row,row + 2) = w*ms*(1.0d0-poro)*(1.0d0)/dz*denss
        
        ymx(row) = w*ms*(1.0d0-poro)*(frx(iz+1)-frx(iz))/dz*denss
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)
     $  -1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb
        
        else if (iz==nz) then 
        amx(row, row) =  w*ms*(1.0d0-poro)*(-1.0d0)/dz*denss
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*1.0d0
     $  -1.0d0*alfacw*fpx(iz)*(-1.0d0))*densb
        
        
        ymx(row) = w*ms*(1.0d0-poro)*(foc-frx(iz))/dz*denss
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)
     $  -1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb
       
        end if
        
!         if (iz/=1) then
          col = row  + 1
          amx(row, col) =  
     $  -kcw*ms*mf*((-1.0d0)*frx(iz)
     $  -1.0d0*alfacw*1.0d0*(1.0d0-frx(iz)))*densb
        
!         end if
          
        
      end do
      
      do iz = 1, nz
        row = 2*(iz-1)+2
        if (iz/=1) then
        amx(row, row) =  -q*mf*(1.0d0)/dz*densf
     $  -kcw*ms*mf*((-1.0d0)*frx(iz)
     $  -1.0d0*alfacw*1.0d0*(1.0d0-frx(iz)))*densb
        
        amx(row,row - 2) = -q*mf*(-1.0d0)/dz*densf
        
        ymx(row) = -q*mf*(fpx(iz)-fpx(iz-1))/dz*densf
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)
     $  -1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb
        
        else if (iz==1) then 
        amx(row, row) =  -q*mf*(1.0d0)/dz*densf
     $  -kcw*ms*mf*((-1.0d0)*frx(iz)
     $  -1.0d0*alfacw*1.0d0*(1.0d0-frx(iz)))*densb
        
        ymx(row) = -q*mf*(fpx(iz)-fsw)/dz*densf
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*frx(iz)
     $  -1.0d0*alfacw*fpx(iz)*(1.0d0-frx(iz)))*densb
       
        end if
        
        
!         if (iz/=nz) then
          col = row  - 1
          amx(row, col) =  
     $  -kcw*ms*mf*((1.0d0-fpx(iz))*1.0d0
     $  -1.0d0*alfacw*fpx(iz)*(-1.0d0))*densb
        
!         end if
      
      end do 
      
      ymx = -ymx
      
      call DGESV(nmx,int(1),amx,nmx,IPIV,ymx,nmx,INFO) 
      
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
      
!       write(*,*) 'error',error
      
      end do
      
!       write(chr,'(i2.2)') i1
!       if ((i2==9).and.(i3==5)) then
!       open (11, file='o_iso_cw'//chr//'.csv',status = 'replace')
!       
!       do iz = 1,nz
!         write(11,*) z(iz)/ztot
!      $  ,',',
!      $  1d3*(fpx(iz)/(1.0d0-fpx(iz))/rstd-1.0d0)
!      $  ,',',
!      $  1d3*(frx(iz)/(1.0d0-frx(iz))/rstd-1.0d0)
!       end do
!       
!       close (11)
!       
!       end if
      
      cdf(i1,i2,i3)=(1.0d0-fpx(1))*frx(1)/
     $  alfacw/fpx(1)/(1.0d0-frx(1))
      
      j1(i1,i2,i3)=0.0d0
      j2(i1,i2,i3)=0.0d0
      j3(i1,i2,i3)=0.0d0
      
      fr = frx
      fp = fpx
      
      j1(i1,i2,i3) = w*ms*denss*(fr(nz)-fr(1))*(1.0d0-poro)
      j2(i1,i2,i3) = j2(i1,i2,i3)+ 0.50d0*densb*ms*mf*kcw*dz*
     $ ((1.0d0-fp(1))*fr(1)-alfacw*(1.0d0-fr(1))*fp(1) +
     $  (1.0d0-fp(nz))*fr(nz)-alfacw*(1.0d0-fr(nz))*fp(nz))
      j3(i1,i2,i3) = -q*mf*densf*(fp(nz)-fp(1))
      do iz = 2, nz-1
        j2(i1,i2,i3) = j2(i1,i2,i3) + densb*ms*mf*kcw*dz*
     $  ((1.0d0-fp(iz))*fr(iz)-alfacw*(1.0d0-fr(iz))*fp(iz))  
      end do     
      
!       write(*,*) j1(i1,i2,i3),j2(i1,i2,i3),j3(i1,i2,i3)
      
      ! error1 = abs(j1(i1,i2,i3)-j2(i1,i2,i3))/abs(j1(i1,i2,i3))
      ! error2 = abs(j3(i1,i2,i3)-j2(i1,i2,i3))/abs(j2(i1,i2,i3))
      ! error3 = abs(j1(i1,i2,i3)-j3(i1,i2,i3))/abs(j1(i1,i2,i3))
      ! write(*,*) 'error', error1, error2, error3
      
      end do
      
      do i1 = 1, n1
        dsw = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
        dswl(i1) = dsw
      
      end do
      
      ! j1 = j1*1.5d14*fa(iry)*(1.0d0-0.75d0*fl(iry))
      ! j1 = j1*1.35d14*fa(iry)*max(0.0d0,(1.0d0-0.75d0*fl(iry)))
      ! j1 = j1*1.35d14*fa(iry)*max(0.0d0,(1.0d0-0.65d0*fl(iry)))
      ! j1 = j1*1.5d14*0.250d0*fa(iry)/fl(iry)
      ! j1 = j1*1.35d14*fa(iry)*fsv(iry)
      ! j1 = j1*1.35d14*fa(iry)*fv(iry)
      
      
      if (flg_fl == 1) then 
      
      if (is == 1) then 
      j1 = j1*1.35d14*fa(iry)*max(0.0d0,(1.0d0-0.65d0*fl(iry)))
      else if (is == 2) then 
      j1 = j1*1.35d14*fa(iry)*fsv(iry)
      end if 
      
      else 
      j1 = j1*1.35d14*fa(iry)*max(0.0d0,(1.0d0-0.65d0*fl(iry)))
      
      end if 
      
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
      end if 
      
!       write(*,*) fa(iry), fl(iry)
      
!       slp = slp*1.5d14*0.250d0*fa(iry)*(1.0d0-fl(iry))
!      $   /(1.0d0-fl(nry))
      
!       itcpt = itcpt*1.5d14*0.250d0*fa(iry)*(1.0d0-fl(iry))
!      $   /(1.0d0-fl(nry))
      
      write(80,*) slp, itcpt, -itcpt/slp
      
!       end do
      
      end do 
      
      write(80,*) ""
      
      end do
      
      close(80)
      
      end 
      
      