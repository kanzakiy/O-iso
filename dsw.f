      program seawater
      
      double precision :: MO = 1.5d24/18.0d0
      
      character(256) chr2
      
      integer, parameter :: nry = 58, nrx = 22
      double precision mry(nrx,nry), fl(nry), fa(nry), fsr(nry)
      double precision rco2(nry), temp(nry), tp(nry), dans(nry)
      double precision bow(nry)
      
      double precision moc(3,nry)
      double precision mha(nry), bha(nry), ssha(nry)
      
      double precision mxcw(3,nry)
      double precision mcw(nry), bcw(nry), sscw(nry)
      
      double precision dsw(6,nry)
      
      integer i1, i2
      integer :: n1 = 13
      
      double precision dswi, M18O(3), Fsw(3)
      double precision :: rstd = 2.0052d-3
      integer rem
      
      double precision maxdy(nry), mindy(nry)
      ! ----------------------------------
      
      write(chr2,*) '-(w_pwad(oc)-(dif-sp)-'//
     $ 'dense-fp_20FLDC_dsp_v2-200_chk-b2)'
      
      open(30, file ='C:/Users/YK/Desktop/Isotope'//
     $   '/royer-forcing.txt',status = 'old')
      read(30,*) mry
      close(30)
      
      fl(:) = mry(6,:)
      fa(:) = mry(7,:)
      fsr(:) = mry(16,:)
      rco2(:) = mry(18,:)
!       rco2(:) = mry(20,:)
      temp(:) = mry(19,:)
      tp(:) = mry(1,:)
      dans(:) = mry(21,:)
      bow(:) = mry(22,:)
      
      ! write(*,*) bow
      
      open (57, file='C:/Users/YK/Desktop/Isotope/'//
     $  'GS-oca'//trim(adjustl(chr2))
     $    //'.csv', 
     $                             status='old')
      
      read (57,*) moc
      close(57)
      
      mha(:) = moc(1,:)
      bha(:) = moc(2,:)
      ssha(:) = moc(3,:)
      
      ! write(*,*) mha
      
      open (12, file ='C:/Users/YK/Desktop/Isotope/'//
     $ 'o_iso_cw_flux.csv',status = 'old')
      read(12, *) mxcw
      close(12)
      
      mcw(:) = mxcw(1,:)
      bcw(:) = mxcw(2,:)
      sscw(:) = mxcw(3,:)
      
      ! write(*,*) mcw
      
      open(11,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'o-iso-cw-dynamic.txt', status = 'replace')
      
      open(10,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'o-iso-cw-steadystate.txt', status = 'replace')
     
      open(9,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'rainwater-dy.txt', status = 'replace')
     
      open(8,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'rainwater-ss.txt', status = 'replace')
     
      open(7,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'likely-dsw-dy.txt', status = 'replace')
     
      open(6,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'likely-dsw-ss.txt', status = 'replace')
     
      open(5,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'rainwater.txt', status = 'replace')
      
      dsw = 0.0d0
      
      do i2 = 1, nry
      
      dsw(4,i2) = -(bha(i2)+bcw(i2))/(mha(i2)+mcw(i2))
      dsw(5,i2) = -(bha(i2)+bcw(i2)+mcw(i2)*dans(i2))
     $   /(mha(i2)+mcw(i2))
      dsw(6,i2) = -(bha(i2)+bcw(i2)+mcw(i2)*bow(i2))
     $   /(mha(i2)+mcw(i2))
      
      write(10,*) tp(i2), dsw(4,i2), dsw(5,i2), dsw(6,i2)
      write(8,*) tp(i2), dsw(4,i2), dsw(5,i2)+dans(i2),
     $  dsw(6,i2)+bow(i2)
      write(6,*) tp(i2), maxval(dsw(4:6,i2))
      
      end do
      
      do i2 = nry, 1, -1
      write(6,*) tp(i2), minval(dsw(4:6,i2))
      
      end do
      
      write(6,*) tp(1), maxval(dsw(4:6,1))
      
      rem = 1
      
      do i1 = 1, n1
      
      dswi = -10.0d0 + 12.0d0*(i1 - 1.0d0)/(n1 - 1.0d0)
      
      dsw(:3,1) = dswi
      
      Fsw(:) = rstd*(dswi/1d3+1.0d0)/(rstd*(dswi/1d3+1.0d0)+1.0d0)
      
      M18O(:) = MO*Fsw(:)
      
      do i2 = 2, nry
      
      M18O(1) = M18O(1) + 10.0d6*(mcw(i2-1)*dsw(1,i2-1) + bcw(i2-1)
     $                + mha(i2-1)*dsw(1,i2-1) + bha(i2-1))
      
      M18O(2) = M18O(2) + 10.0d6*(mcw(i2-1)
     $    *(dsw(2,i2-1) +dans(i2-1)) + bcw(i2-1)
     $      +mha(i2-1)*dsw(2,i2-1) + bha(i2-1))
      
      M18O(3) = M18O(3) + 10.0d6*(mcw(i2-1)
     $    *(dsw(3,i2-1) +bow(i2-1)) + bcw(i2-1)
     $      +mha(i2-1)*dsw(3,i2-1) + bha(i2-1))
     
      Fsw(:) = M18O(:)/MO
      
      dsw(:3,i2) = 1d3*(Fsw(:)/(1.0d0-Fsw(:))/rstd - 1.0d0)
      
      end do
      
      if (rem == i1) then
      do i2 = 1, nry
      write(11,*) tp(i2), dsw(1,i2), dsw(2,i2), dsw(3,i2)
      write(9,*) dsw(1,i2), dsw(2,i2)+dans(i2),
     $   dsw(3,i2)+bow(i2)
      end do 
      
      else if (rem/=i1) then
      write(11,*) ''
      write(9,*) ''
      do i2 = 1, nry
      write(11,*) tp(i2), dsw(1,i2), dsw(2,i2), dsw(3,i2)
      write(9,*) dsw(1,i2), dsw(2,i2)+dans(i2),
     $   dsw(3,i2)+bow(i2)
      end do 
      
      rem = i1
      
      end if
      
      if (i1 == 1) then
      mindy(:) = min(dsw(3,:),min(dsw(1,:),dsw(2,:)))
      do i2 = 1, nry
      write(5,*) dsw(1,i2), dsw(2,i2)+dans(i2),
     $   dsw(3,i2)+bow(i2)
      end do
      else if (i1 == 6) then
      maxdy(:) = max(max(dsw(1,:),dsw(2,:)),dsw(3,:))
      ! else if (i1 == 2) then
      ! do i2 = 1, nry
      ! write(5,*) dsw(1,i2), dsw(2,i2)+dans(i2),
      ! $   dsw(3,i2)+bow(i2)
      ! end do
      end if
      
      end do
      
      do i2 = 1, nry
      write(7,*) tp(i2), maxdy(i2)
      end do
      
      do i2 = nry, 1, -1
      write(7,*) tp(i2), mindy(i2)
      end do
      
      write(7,*) tp(1), maxdy(1)
      
      close (11)
      close (10)
      close (9)
      close (8)
      close (7)
      close (6)
      close (5)
      
      end