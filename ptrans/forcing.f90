module forcing
   use param
   use mpi_common
   use common
   implicit none

contains

   subroutine apply_forcing(u, v, w, vis)
      real(8), dimension(2, n, nb1, nbh2) :: u, v, w
      real(8), dimension(n, nb1, nbh2) :: vis
      real(8) :: ene, ene0, ens, ens0, enef, enef0, nrm
      real(8), dimension(3) :: en
      real(8) :: alpha
      integer :: i, j, k

      ene = 0.0d0
      ene0 = 0.0d0
      ens = 0.0d0
      ens0 = 0.0d0
      enef = 0.0d0
      enef0 = 0.0d0
!$omp parallel do reduction(+:ene,ene0,ens,ens0,enef,enef0) private(nrm)
      do k = 1, nbh2
         if (me2 == 0 .and. k == 1) then
            do j = 1, nb1
               do i = 1, n
                  nrm = (u(1, i, j, k)**2+u(2, i, j, k)**2 &
                  &     +v(1, i, j, k)**2+v(2, i, j, k)**2 &
                  &     +w(1, i, j, k)**2+w(2, i, j, k)**2)*msk(i, j, k)
                  ene0 = ene0 + nrm
                  ens0 = ens0 + nrm*wk(i, j, k)
                  enef0 = enef0 + nrm*mskf(i, j, k)
               end do
            end do
         else
            do j = 1, nb1
               do i = 1, n
                  nrm = (u(1, i, j, k)**2+u(2, i, j, k)**2 &
                  &     +v(1, i, j, k)**2+v(2, i, j, k)**2 &
                  &     +w(1, i, j, k)**2+w(2, i, j, k)**2)*msk(i, j, k)
                  ene = ene + nrm
                  ens = ens + nrm*wk(i, j, k)
                  enef = enef + nrm*mskf(i, j, k)
               end do
            end do
         end if
      end do

      en(1) = ene0*0.5d0 + ene
      en(2) = ens0*0.5d0 + ens
      en(3) = enef0*0.5d0 + enef

      call mpi_allreduce(MPI_IN_PLACE, en, 3, MPI_REAL8, &
      &                  MPI_SUM, comm, ierr)

      alpha = nu*en(2)/en(3)

      if (myrank == 0 ) write(*, *) en, alpha

!$omp parallel do
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               vis(i, j, k) = alpha*mskf(i, j, k) - nu*wk(i, j, k)
            end do
         end do
      end do

   end subroutine apply_forcing

end module forcing
