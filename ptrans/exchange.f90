module exchange
   use param
   use mpi_common
   implicit none

contains

   subroutine ex1(s, r)
      real(8), dimension((nb1+1)*nb2*nb1, np1) :: s, r

      call mpi_alltoall(s, (nb1+1)*nb2*nb1, MPI_REAL8, &
      &                 r, (nb1+1)*nb2*nb1, MPI_REAL8, comm1, ierr)

   end subroutine ex1

   subroutine ex2(s, r)
      real(8), dimension((nb2+2)*nb1*nb2, np2) :: s, r

      call mpi_alltoall(s, (nb2+2)*nb1*nb2, MPI_REAL8, &
      &                 r, (nb2+2)*nb1*nb2, MPI_REAL8, comm2, ierr)

   end subroutine ex2

end module exchange
