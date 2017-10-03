module param
   implicit none

   character(len=*), parameter :: runname = "Eta1_12L"
   real(8), parameter :: dt = 0.00025d0*2.0d0/3.0d0
   real(8), parameter :: nu = 0.0000102d0
   real(8), parameter :: wp = 2.0d0
   real(8), parameter :: u0 = 1.0d0

   integer, parameter :: n = 2**11*3
   integer, parameter :: np1 = 96
   integer, parameter :: np2 = 64

   integer, parameter :: nh = n/2
   real(8), parameter :: dn3 = dble(n)**3
   integer, parameter :: np = np1*np2
   integer, parameter :: np0 = np*2
   integer, parameter :: nb1 = n/np1
   integer, parameter :: nb2 = n/np2
   integer, parameter :: nbh1 = nb1/2
   integer, parameter :: nbh2 = nb2/2
   integer, parameter :: nk = nh*7/4
   real(8), parameter :: kmax = n*sqrt(2.0d0)/3.0d0
   real(8), parameter :: fkmin = 0.5d0
   real(8), parameter :: fkmax = 2.5d0
   real(8), parameter :: kmax2 = kmax**2
   real(8), parameter :: fkmin2 = fkmin**2
   real(8), parameter :: fkmax2 = fkmax**2

   character(len=*), parameter :: idir = "./in/"
   character(len=*), parameter :: odir = "./out/"
   character(len=*), parameter :: ddir = "./data/"

end module param
