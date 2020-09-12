module o_iso_ref_mod
implicit none

real(kind=8),parameter :: rsmow_18 = 2.0052d-3 
real(kind=8),parameter :: rsmow_17 = 3.799d-4 
real(kind=8),parameter :: d18_mc = 5.7d0 
! real(kind=8),parameter :: d17_mc = 2.86d0   
real(kind=8),parameter :: d17_mc = 3.074982726259945d0   
! real(kind=8),parameter :: d18_mc = 5.3d0 
! real(kind=8),parameter :: d17_mc = 2.7d0  

! logical,parameter :: switch_KIE = .true.
logical,parameter :: switch_KIE = .false.

endmodule o_iso_ref_mod