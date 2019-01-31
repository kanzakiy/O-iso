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
      double precision :: denss = 3.0d3
      
      integer, parameter :: n1 = 11
      integer, parameter :: n2 = 101
      integer, parameter :: n3 = 101
      
      double precision j1(n1,n2,n3),j2(n1,n2,n3),j3(n1,n2,n3),dswl(n1)
      double precision cdf(n1,n2,n3),shf(n1,n2,n3),dsf(n1,n2,n3)
      
      character (2) chr
      
      double precision avx,avy,sumx2,sumy2,sumxy,slp,itcpt,crr
      
      double precision frtop, error1, error2, error3
      integer rem
!---------------------------------------------------      
      
      open (12, file ='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw_flux.csv',status = 'replace')
      open (13, file ='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw_cdf.csv',status = 'replace')
      open (21, file='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw.dat',status = 'replace')
      open (22, file ='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw_flux.dat',status = 'replace')
      open (23, file ='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw_cdf.dat',status = 'replace')
      open (32, file ='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw_shift.dat',status = 'replace')
      
      rem = 1
      
      do i3 = 1, n3
      
      ! w = 5d-6*10.0d0**(2.0d0*(i3-1.0d0)/(n3-1.0d0))
      k1 = 1d-11*10.0d0**(5.0d0*(i3-1.0d0)/(n3-1.0d0))
      ! q = 3d-3*10.0d0**(4.0d0*(i3-1.0d0)/(n3-1.0d0))
      
      do i2 = 1, n2
      
      tau = 1d4*10.0d0**(2.0d0*(i2-1.0d0)/(n2-1.0d0))
      ! k1 = 1d-11*10.0d0**(5.0d0*(i2-1.0d0)/(n2-1.0d0))
      ! poro = 0.050d0+ (0.40d0*(i2-1.0d0)/(n2-1.0d0))
      
      ztot = w*tau
      
      do i1 = 1, n1
      
      dsw = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
      
      kcw = k1*exp(-e*(1.0d0/(273.0d0+tempcw)-1.0d0/(273.0d0+temp1))/rg)
      
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
      
      write(chr,'(i2.2)') i1
      if ((i2==8).and.(i3==5)) then
      
      ! write(*,*) tau
      
      open (11, file='C:/users/user/desktop/isotope-res/'//
     $    'o_iso_cw'//chr//'.csv',status = 'replace')
      
      frtop = foc- 
     $   q*densf*mf*(fsw - fpx(nz))/w/denss/ms/(1.0d0-poro)   ! calculation with mass balance
      
      
      ! frtop = (w*ms*(1.0d0-poro)*(frx(1))/dz*denss
      ! $   +kcw*ms*mf*alfacw*fsw*densb)/
      ! $   (w*ms*(1.0d0-poro)/dz*denss
      ! $   +kcw*ms*mf*((1.0d0-fsw)*frtop+alfacw*fsw)*densb)
      
      ! w*ms*(1.0d0-poro)*(frx(0)-frtop)/dz*denss
      ! $  -kcw*ms*mf*((1.0d0-fsw)*frtop
      ! $  -1.0d0*alfacw*fsw*(1.0d0-frtop))*densb = 0
      
      write(11,*) 0.0d0
     $  ,',',
     $  1d3*(fsw/(1.0d0-fsw)/rstd-1.0d0)
     $  ,',',
     $  1d3*(frtop/(1.0d0-frtop)/rstd-1.0d0)
      
      do iz = 1,nz
        write(11,*) z(iz)/ztot
     $  ,',',
     $  1d3*(fpx(iz)/(1.0d0-fpx(iz))/rstd-1.0d0)
     $  ,',',
     $  1d3*(frx(iz)/(1.0d0-frx(iz))/rstd-1.0d0)
      end do
      
      close (11)
      
      write(21,*) 0.0d0
     $  ,
     $  1d3*(fsw/(1.0d0-fsw)/rstd-1.0d0)
     $  ,
     $  1d3*(frtop/(1.0d0-frtop)/rstd-1.0d0)
     
      do iz = 1,nz
        write(21,*) z(iz)/ztot
     $  ,
     $  1d3*(fpx(iz)/(1.0d0-fpx(iz))/rstd-1.0d0)
     $  ,
     $  1d3*(frx(iz)/(1.0d0-frx(iz))/rstd-1.0d0)
      end do
      
      write(21,*) ''
      
      end if
      
      cdf(i1,i2,i3)=(1.0d0-fpx(1))*frx(1)/
     $  alfacw/fpx(1)/(1.0d0-frx(1))
      
      j1(i1,i2,i3)=0.0d0
      j2(i1,i2,i3)=0.0d0
      j3(i1,i2,i3)=0.0d0
      shf(i1,i2,i3)=0.0d0
      dsf(i1,i2,i3) = 0.0d0
      
      fr = frx
      fp = fpx
      
      j1(i1,i2,i3) = w*ms*denss*(fr(nz)-fr(1))*(1.0d0-poro)
      ! j1(i1,i2,i3) = w*ms*denss*(foc-frtop)*(1.0d0-poro)
      j2(i1,i2,i3) = j2(i1,i2,i3)+ 0.50d0*densb*ms*mf*kcw*dz*
     $ ((1.0d0-fp(1))*fr(1)-alfacw*(1.0d0-fr(1))*fp(1) +
     $  (1.0d0-fp(nz))*fr(nz)-alfacw*(1.0d0-fr(nz))*fp(nz))
      j3(i1,i2,i3) = -q*mf*densf*(fp(nz)-fp(1))
      ! j3(i1,i2,i3) = -q*mf*densf*(fp(nz)-fsw)
      shf(i1,i2,i3)= 1d3*(fr(1)/(1.0d0-fr(1))-fr(nz)/(1.0d0-fr(nz)))
     $  /rstd
      dsf(i1,i2,i3)= 1d3*(fr(1)/(1.0d0-fr(1))
     $  /rstd-1.0d0)
      do iz = 2, nz-1
        j2(i1,i2,i3) = j2(i1,i2,i3) + densb*ms*mf*kcw*dz*
     $  ((1.0d0-fp(iz))*fr(iz)-alfacw*(1.0d0-fr(iz))*fp(iz))  
      end do     
      
      error1 = abs(j1(i1,i2,i3)-j2(i1,i2,i3))/abs(j1(i1,i2,i3))
      error2 = abs(j3(i1,i2,i3)-j2(i1,i2,i3))/abs(j2(i1,i2,i3))
      error3 = abs(j1(i1,i2,i3)-j3(i1,i2,i3))/abs(j1(i1,i2,i3))
      write(*,*) 'error', error1, error2, error3
      ! write(*,*) j1(i1,i2,i3),j2(i1,i2,i3),j3(i1,i2,i3)
      
      end do
      
      do i1 = 1, n1
        dsw = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
        dswl(i1) = dsw
      
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
      end if 
      
      write(12,*) 
     $  k1
     $  ,',',  
     $  tau
     $  ,',',     
     $   slp
     $  ,',',    
     $   itcpt
     $  ,',',    
     $   -itcpt/slp
      
      if (rem == i3) then
        write(22,*)  k1, tau, slp,itcpt, -itcpt/slp
        write(32,*)  k1, tau
     $   , shf(7,i2,i3),dsf(7,i2,i3)
     $   , shf(n1,i2,i3),dsf(n1,i2,i3)
     $   , merge(1.0d0,0.0d0,
     $    (shf(n1,i2,i3)<=7.0d0).and.
     $    (shf(7,i2,i3)>=3.0d0))
      else if (rem /= i3) then
        write(22,*)  ''
        write(32,*)  ''
        write(22,*)  k1, tau, slp,itcpt, -itcpt/slp
        write(32,*)  k1, tau
     $   , shf(7,i2,i3),dsf(7,i2,i3)
     $   , shf(n1,i2,i3),dsf(n1,i2,i3)
     $   , merge(1.0d0,0.0d0,
     $    (shf(n1,i2,i3)<=7.0d0).and.
     $    (shf(7,i2,i3)>=3.0d0))
        rem =i3
      end if 
      
      end do
      
      end do 
      
      close(12)
      
      rem = 1 
      
      do i3 = 1, n3
      
!       w = 5d-6*10.0d0**(2.0d0*(i3-1.0d0)/(n3-1.0d0))
      k1 = 1d-11*10.0d0**(5.0d0*(i3-1.0d0)/(n3-1.0d0))
      ! q = 3d-3*10.0d0**(4.0d0*(i3-1.0d0)/(n3-1.0d0))
      
      do i2 = 1, n2
      
      tau = 1d4*10.0d0**(2.0d0*(i2-1.0d0)/(n2-1.0d0))
!       k1 = 1d-11*10.0d0**(5.0d0*(i2-1.0d0)/(n2-1.0d0))
      ! poro = 0.050d0+ (0.40d0*(i2-1.0d0)/(n2-1.0d0))
      
      do i1 = 1, n1
      
      dsw = -20.0d0 + 20.0d0*(i1-1.0d0)/(n1-1.0d0)
      
      write(13,*)
     $  k1
     $  ,',',
     $  tau
     $  ,',',
     $  dsw
     $  ,',',
     $  cdf(i1,i2,i3)
      
      if (rem == i3) then
        write(23,*)  k1, tau, dsw, cdf(i1,i2,i3)
      else if (rem /= i3) then
        write(23,*)  ''
        write(23,*)  k1, tau, dsw, cdf(i1,i2,i3)
        rem =i3
      end if 
      
      end do
      end do 
      end do
      
      close(13)
      close(21)
      close(22)
      close(23)
      close(32)
      
      end 
      
      