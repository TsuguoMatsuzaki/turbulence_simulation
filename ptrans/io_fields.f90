module io_fields
   use param
   use mpi_common
   implicit none

   integer, parameter :: npi1 = 96
   integer, parameter :: npi2 = 64
   integer, parameter :: npo1 = 96
   integer, parameter :: npo2 = 64

contains

   subroutine read_fields(u, v, w, it)
      real(8), dimension(2, n, nb1, nbh2) :: u, v, w
      integer :: it
      integer :: m1, m2
      integer :: js, je, ks, ke, ich5
      character(len=17) :: filename

      do m2 = 0, npi2/np2-1
         do m1 = 0, npi1/np1-1
            js = n/npi1*m1+1
            je = js+n/npi1-1
            ks = nh/npi2*m2+1
            ke = ks+nh/npi2-1
            ich5 = (npi1*npi2/np2*me2+npi1/np1*me1)+npi1*m2+m1
            write(filename, '("_", i5.5, "_p", i5.5, ".bin")') it, ich5
            open(11, file = idir//runname//filename, form = "unformatted")
            read(11) u(1:2, 1:n, js:je, ks:ke), &
            &        v(1:2, 1:n, js:je, ks:ke), &
            &        w(1:2, 1:n, js:je, ks:ke)
            close(11)
         end do
      end do

   end subroutine read_fields

   subroutine write_fields(u, v, w, it)
      real(8), dimension(2, n, nb1, nbh2) :: u, v, w
      integer :: it
      integer :: m1, m2
      integer :: js, je, ks, ke, ich5
      character(len=17) :: filename

      do m2 = 0, npo2/np2-1
         do m1 = 0, npo1/np1-1
            js = n/npo1*m1+1
            je = js+n/npo1-1
            ks = nh/npo2*m2+1
            ke = ks+nh/npo2-1
            ich5 = (npo1*npo2/np2*me2+npo1/np1*me1)+npo1*m2+m1
            write(filename, '("_", i5.5, "_p", i5.5, ".bin")') it, ich5
            open(11, file = odir//runname//filename, form = "unformatted")
            write(11) u(1:2, 1:n, js:je, ks:ke), &
            &         v(1:2, 1:n, js:je, ks:ke), &
            &         w(1:2, 1:n, js:je, ks:ke)
            close(11)
         end do
      end do

   end subroutine write_fields

end module io_fields
