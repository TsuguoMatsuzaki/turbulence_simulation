module spectrum
   use param
   use mpi_common
   use common
   implicit none

contains

   subroutine calc_spectrum(u, v, w, hu, hv, hw, sp)
      real(8), dimension(2, n, nb1, nbh2) :: u, v, w
      real(8), dimension(2, n, nb1, nbh2) :: hu, hv, hw
      real(8), dimension(4, 0:nk) :: sp
      real(8) :: uu, uh
      integer :: i, j, k, kk

      sp = 0.0d0
!$omp parallel do reduction(+:sp) private(kk,uu,uh)
      do k = 1, nbh2
         if (me2 == 0 .and. k == 1) then
            do j = 1, nb1
               do i = 1, n
                  kk = nint(sqrt(wk(i, j, k)))

                  uu = (u(1, i, j, k)**2+u(2, i, j, k)**2 &
                  &    +v(1, i, j, k)**2+v(2, i, j, k)**2 &
                  &    +w(1, i, j, k)**2+w(2, i, j, k)**2)*0.5d0
                  uh = (u(1, i, j, k)*hu(1, i, j, k)+u(2, i, j, k)*hu(2, i, j, k) &
                  &    +v(1, i, j, k)*hv(1, i, j, k)+v(2, i, j, k)*hv(2, i, j, k) &
                  &    +w(1, i, j, k)*hw(1, i, j, k)+w(2, i, j, k)*hw(2, i, j, k))
                  sp(1, kk) = sp(1, kk) + uu
                  sp(2, kk) = sp(2, kk) + uu*wk(i, j, k)
                  sp(3, kk) = sp(3, kk) + uh
                  sp(4, kk) = sp(4, kk) + uh*wk(i, j, k)
               end do
            end do
         else
            do j = 1, nb1
               do i = 1, n
                  kk = nint(sqrt(wk(i, j, k)))

                  uu = (u(1, i, j, k)**2+u(2, i, j, k)**2 &
                  &    +v(1, i, j, k)**2+v(2, i, j, k)**2 &
                  &    +w(1, i, j, k)**2+w(2, i, j, k)**2)
                  uh = (u(1, i, j, k)*hu(1, i, j, k)+u(2, i, j, k)*hu(2, i, j, k) &
                  &    +v(1, i, j, k)*hv(1, i, j, k)+v(2, i, j, k)*hv(2, i, j, k) &
                  &    +w(1, i, j, k)*hw(1, i, j, k)+w(2, i, j, k)*hw(2, i, j, k))*2.0d0
                  sp(1, kk) = sp(1, kk) + uu
                  sp(2, kk) = sp(2, kk) + uu*wk(i, j, k)
                  sp(3, kk) = sp(3, kk) + uh
                  sp(4, kk) = sp(4, kk) + uh*wk(i, j, k)
               end do
            end do
         end if
      end do

      call mpi_allreduce(MPI_IN_PLACE, sp, 4*(nk+1), MPI_REAL8, &
      &                  MPI_SUM, comm, ierr)

   end subroutine calc_spectrum

end module spectrum
