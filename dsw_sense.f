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
      integer :: n1 = 2
      
      double precision dswi, M18O(3), Fsw(3)
      double precision :: rstd = 2.0052d-3
      integer rem
      
      double precision maxdy(nry), mindy(nry)
      
      integer is, ismx
C        integer, parameter :: ns = 8   ! ff
C       integer, parameter :: ns = 3   ! fb
C       integer, parameter :: ns = 2   ! fl
C       integer, parameter :: ns = 11
C       integer, parameter :: ns = 5
      integer, parameter :: ns = 6
      
      double precision dswimx, dswimn
      
      double precision frw(nry), fd(nry)
      double precision rain(6,nry)
      ! ----------------------------------
     
      open(7,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'likely-dsw-dy_sense.txt', status = 'replace')
     
      open(6,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'likely-dsw-ss_sense.txt', status = 'replace')
     
      open(5,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'rainwater-list.txt', status = 'replace')
     
      open(9,file = 'C:/Users/YK/Desktop/isotope/'//
     $  'seawater-list.txt', status = 'replace')
      
      write(chr2,*) '-(w_pwad(oc)-(dif-sp)-'//
     $ 'dense-fp_20FLDC_dsp_v2-200_chk-b2)'
      
      open(30, file ='C:/Users/YK/Desktop/Isotope'//
     $   '/royer-forcing.txt',status = 'old')
      read(30,*) mry
      close(30)
      
      frw(:) = mry(5,:)
      fl(:) = mry(6,:)
      fa(:) = mry(7,:)
      fd(:) = mry(11,:)
      fsr(:) = mry(16,:)
      rco2(:) = mry(18,:)
!       rco2(:) = mry(20,:)
      temp(:) = mry(19,:)
      tp(:) = mry(1,:)
      dans(:) = mry(21,:)
      bow(:) = mry(22,:)
      
      ! write(*,*) bow
      
      open (57, file='C:/Users/YK/Desktop/Isotope/'//
     $  'GS-oca_sense'//trim(adjustl(chr2))
     $    //'.txt', 
     $                             status='old')
      
      open (12, file ='C:/Users/YK/Desktop/Isotope/'//
     $ 'o_iso_cw_flux_sense.txt',status = 'old')
      
      do is = 1, ns
      
      read (57,*) moc
      
      mha(:) = moc(1,:)
      bha(:) = moc(2,:)
      ssha(:) = moc(3,:)
      
      ! write(*,*) mha
      read(12, *) mxcw
      
      mcw(:) = mxcw(1,:)
      bcw(:) = mxcw(2,:)
      sscw(:) = mxcw(3,:)
      
      
      dsw = 0.0d0
      
      do i2 = 1, nry
      
      dsw(4,i2) = -(bha(i2)+bcw(i2))/(mha(i2)+mcw(i2))
      dsw(5,i2) = -(bha(i2)+bcw(i2)+mcw(i2)*dans(i2))
     $   /(mha(i2)+mcw(i2))
      dsw(6,i2) = -(bha(i2)+bcw(i2)+mcw(i2)*bow(i2))
     $   /(mha(i2)+mcw(i2))
      
      ! write(10,*) tp(i2), dsw(4,i2), dsw(5,i2), dsw(6,i2)
      ! write(8,*) tp(i2), dsw(4,i2), dsw(5,i2)+dans(i2),
      ! $  dsw(6,i2)+bow(i2)
      write(6,*) tp(i2), maxval(dsw(4:6,i2))
      
      end do
      
      do i2 = nry, 1, -1
      write(6,*) tp(i2), minval(dsw(4:6,i2))
      
      end do
      
      write(6,*) tp(1), maxval(dsw(4:6,1))
      write(6,*) ''
      
      dswimx = maxval(dsw(4:6,1))
      dswimn = minval(dsw(4:6,1))
      
      
      rem = 1
      
      do i1 = 1, n1
      
      ! dswi = -10.0d0 + 12.0d0*(i1 - 1.0d0)/(n1 - 1.0d0)
      
      if (i1==1) dswi = min(-10.0d0,dswimn)
      if (i1==2) dswi = max(-7.0d0,dswimx)
      
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
      
      ! if (rem == i1) then
      ! do i2 = 1, nry
      ! write(11,*) tp(i2), dsw(1,i2), dsw(2,i2), dsw(3,i2)
      ! write(9,*) dsw(1,i2), dsw(2,i2)+dans(i2),
      ! $   dsw(3,i2)+bow(i2)
      ! end do 
      
      ! else if (rem/=i1) then
      ! write(11,*) ''
      ! write(9,*) ''
      ! do i2 = 1, nry
      ! write(11,*) tp(i2), dsw(1,i2), dsw(2,i2), dsw(3,i2)
      ! write(9,*) dsw(1,i2), dsw(2,i2)+dans(i2),
      ! $   dsw(3,i2)+bow(i2)
      ! end do 
      
      ! rem = i1
      
      ! end if
      
      if (i1 == 1) then
      mindy(:) = min(dsw(3,:),min(dsw(1,:),dsw(2,:)))
      rain(1,:) = dsw(1,:)
      rain(2,:) = dsw(2,:) + dans(:)
      rain(3,:) = dsw(3,:) + bow(:)
      ! do i2 = 1, nry
      ! write(5,*) dsw(1,i2), dsw(2,i2)+dans(i2),
      ! $   dsw(3,i2)+bow(i2)
      ! end do
      else if (i1 == 2) then
      maxdy(:) = max(max(dsw(1,:),dsw(2,:)),dsw(3,:))
      rain(4,:) = dsw(1,:)
      rain(5,:) = dsw(2,:) + dans(:)
      rain(6,:) = dsw(3,:) + bow(:)
      ! else if (i1 == 2) then
      ! do i2 = 1, nry
      ! write(5,*) dsw(1,i2), dsw(2,i2)+dans(i2),
      ! $   dsw(3,i2)+bow(i2)
      ! end do
      end if
      
      end do
      
      do i2 = 1, nry
      write(7,*) tp(i2), maxdy(i2)
      write(5,*) rain(1,i2), rain(2,i2), rain(3,i2),
     $   rain(4,i2), rain(5,i2), rain(6,i2)
      write(9,*) mindy(i2),maxdy(i2)
      end do
      
      do i2 = nry, 1, -1
      write(7,*) tp(i2), mindy(i2)
      end do
      
      write(7,*) tp(1), maxdy(1)
      write(7,*) ''
      write(5,*) ''
      write(9,*) ''
      
      
      
      
      end do
      
      close (5)
      close (9)
      close (7)
      close (6)
      close(12)
      close(57)
      ! close (5)
      
      end