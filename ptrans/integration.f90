module integration
   use param
   use mpi_common
   use common
   use forcing
   use nonlinear
   use timer
   implicit none

   real(8), parameter :: dt2 = dt*0.5d0
   real(8), parameter :: cz1  = 1.0d0 - sqrt(0.5d0)
   real(8), parameter :: ct11 = 1.0d0 - 3.0d0*cz1
   real(8), parameter :: ct12 = 2.0d0*cz1
   real(8), parameter :: cz2  = 1.0d0 + sqrt(0.5d0)
   real(8), parameter :: ct21 = 1.0d0 - 3.0d0*cz2
   real(8), parameter :: ct22 = 2.0d0*cz2
   real(8), parameter :: cr3  = 1.0d0/3.0d0

contains

   subroutine integrate(u, v, w, hu, hv, hw)
      real(8), dimension(2, n, nb1, nbh2) :: u, v, w
      real(8), dimension(2, n, nb1, nbh2) :: hu, hv, hw
      real(8), dimension(2, n, nb1, nbh2) :: qu, qv, qw
      real(8), dimension(n, nb1, nbh2) :: vis
      real(8) :: dur, dui, dvr, dvi, dwr, dwi
      integer :: i, j, k

      call apply_forcing(u, v, w, vis)

!$omp parallel do private(dur,dui,dvr,dvi,dwr,dwi)
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               dur = (hu(1, i, j, k) + vis(i, j, k)*u(1, i, j, k))*dt
               dui = (hu(2, i, j, k) + vis(i, j, k)*u(2, i, j, k))*dt
               u(1, i, j, k) = u(1, i, j, k) + dur*0.5d0
               u(2, i, j, k) = u(2, i, j, k) + dui*0.5d0
               qu(1, i, j, k) = dur
               qu(2, i, j, k) = dui

               dvr = (hv(1, i, j, k) + vis(i, j, k)*v(1, i, j, k))*dt
               dvi = (hv(2, i, j, k) + vis(i, j, k)*v(2, i, j, k))*dt
               v(1, i, j, k) = v(1, i, j, k) + dvr*0.5d0
               v(2, i, j, k) = v(2, i, j, k) + dvi*0.5d0
               qv(1, i, j, k) = dvr
               qv(2, i, j, k) = dvi

               dwr = (hw(1, i, j, k) + vis(i, j, k)*w(1, i, j, k))*dt
               dwi = (hw(2, i, j, k) + vis(i, j, k)*w(2, i, j, k))*dt
               w(1, i, j, k) = w(1, i, j, k) + dwr*0.5d0
               w(2, i, j, k) = w(2, i, j, k) + dwi*0.5d0
               qw(1, i, j, k) = dwr
               qw(2, i, j, k) = dwi
            end do
         end do
      end do

      call calc_nonlinear(u, v, w, hu, hv, hw)

!$omp parallel do private(dur,dui,dvr,dvi,dwr,dwi)
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               dur = (hu(1, i, j, k) + vis(i, j, k)*u(1, i, j, k))*dt
               dui = (hu(2, i, j, k) + vis(i, j, k)*u(2, i, j, k))*dt
               u(1, i, j, k) = u(1, i, j, k) + cz1*(dur - qu(1, i, j, k))
               u(2, i, j, k) = u(2, i, j, k) + cz1*(dui - qu(2, i, j, k))
               qu(1, i, j, k) = ct11*qu(1, i, j, k) + ct12*dur
               qu(2, i, j, k) = ct11*qu(2, i, j, k) + ct12*dui

               dvr = (hv(1, i, j, k) + vis(i, j, k)*v(1, i, j, k))*dt
               dvi = (hv(2, i, j, k) + vis(i, j, k)*v(2, i, j, k))*dt
               v(1, i, j, k) = v(1, i, j, k) + cz1*(dvr - qv(1, i, j, k))
               v(2, i, j, k) = v(2, i, j, k) + cz1*(dvi - qv(2, i, j, k))
               qv(1, i, j, k) = ct11*qv(1, i, j, k) + ct12*dvr
               qv(2, i, j, k) = ct11*qv(2, i, j, k) + ct12*dvi

               dwr = (hw(1, i, j, k) + vis(i, j, k)*w(1, i, j, k))*dt
               dwi = (hw(2, i, j, k) + vis(i, j, k)*w(2, i, j, k))*dt
               w(1, i, j, k) = w(1, i, j, k) + cz1*(dwr - qw(1, i, j, k))
               w(2, i, j, k) = w(2, i, j, k) + cz1*(dwi - qw(2, i, j, k))
               qw(1, i, j, k) = ct11*qw(1, i, j, k) + ct12*dwr
               qw(2, i, j, k) = ct11*qw(2, i, j, k) + ct12*dwi
            end do
         end do
      end do

      call calc_nonlinear(u, v, w, hu, hv, hw)

!$omp parallel do private(dur,dui,dvr,dvi,dwr,dwi)
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               dur = (hu(1, i, j, k) + vis(i, j, k)*u(1, i, j, k))*dt
               dui = (hu(2, i, j, k) + vis(i, j, k)*u(2, i, j, k))*dt
               u(1, i, j, k) = u(1, i, j, k) + cz2*(dur - qu(1, i, j, k))
               u(2, i, j, k) = u(2, i, j, k) + cz2*(dui - qu(2, i, j, k))
               qu(1, i, j, k) = ct21*qu(1, i, j, k) + ct22*dur
               qu(2, i, j, k) = ct21*qu(2, i, j, k) + ct22*dui

               dvr = (hv(1, i, j, k) + vis(i, j, k)*v(1, i, j, k))*dt
               dvi = (hv(2, i, j, k) + vis(i, j, k)*v(2, i, j, k))*dt
               v(1, i, j, k) = v(1, i, j, k) + cz2*(dvr - qv(1, i, j, k))
               v(2, i, j, k) = v(2, i, j, k) + cz2*(dvi - qv(2, i, j, k))
               qv(1, i, j, k) = ct21*qv(1, i, j, k) + ct22*dvr
               qv(2, i, j, k) = ct21*qv(2, i, j, k) + ct22*dvi

               dwr = (hw(1, i, j, k) + vis(i, j, k)*w(1, i, j, k))*dt
               dwi = (hw(2, i, j, k) + vis(i, j, k)*w(2, i, j, k))*dt
               w(1, i, j, k) = w(1, i, j, k) + cz2*(dwr - qw(1, i, j, k))
               w(2, i, j, k) = w(2, i, j, k) + cz2*(dwi - qw(2, i, j, k))
               qw(1, i, j, k) = ct21*qw(1, i, j, k) + ct22*dwr
               qw(2, i, j, k) = ct21*qw(2, i, j, k) + ct22*dwi
            end do
         end do
      end do

      call calc_nonlinear(u, v, w, hu, hv, hw)

!$omp parallel do private(dur,dui,dvr,dvi,dwr,dwi)
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               dur = (hu(1, i, j, k) + vis(i, j, k)*u(1, i, j, k))*dt2
               dui = (hu(2, i, j, k) + vis(i, j, k)*u(2, i, j, k))*dt2
               u(1, i, j, k) = u(1, i, j, k) + (dur - qu(1, i, j, k))*cr3
               u(2, i, j, k) = u(2, i, j, k) + (dui - qu(2, i, j, k))*cr3

               dvr = (hv(1, i, j, k) + vis(i, j, k)*v(1, i, j, k))*dt2
               dvi = (hv(2, i, j, k) + vis(i, j, k)*v(2, i, j, k))*dt2
               v(1, i, j, k) = v(1, i, j, k) + (dvr - qv(1, i, j, k))*cr3
               v(2, i, j, k) = v(2, i, j, k) + (dvi - qv(2, i, j, k))*cr3

               dwr = (hw(1, i, j, k) + vis(i, j, k)*w(1, i, j, k))*dt2
               dwi = (hw(2, i, j, k) + vis(i, j, k)*w(2, i, j, k))*dt2
               w(1, i, j, k) = w(1, i, j, k) + (dwr - qw(1, i, j, k))*cr3
               w(2, i, j, k) = w(2, i, j, k) + (dwi - qw(2, i, j, k))*cr3
            end do
         end do
      end do

      call calc_nonlinear(u, v, w, hu, hv, hw)

   end subroutine integrate

end module integration
