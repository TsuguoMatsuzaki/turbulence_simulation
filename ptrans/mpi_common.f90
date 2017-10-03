module mpi_common
   use mpi
   use param
   implicit none

   integer :: ierr, nprocs, myrank
   integer :: comm, comm1, comm2
   integer :: me, me1, me2

contains

   subroutine init_mpi
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

      if (np0 /= nprocs) then
         print *, "Error!!! : np0 /= nprocs"
         stop
      end if

      call mpi_comm_split(MPI_COMM_WORLD, myrank/np, myrank, comm, ierr)
      call mpi_comm_rank(comm, me, ierr)

      call mpi_comm_split(comm, me/np1, me, comm1, ierr)
      call mpi_comm_split(comm, mod(me,np1), me, comm2, ierr)
      call mpi_comm_rank(comm1, me1, ierr)
      call mpi_comm_rank(comm2, me2, ierr)
   end subroutine init_mpi

   subroutine finalize_mpi
      call mpi_comm_free(comm1, ierr)
      call mpi_comm_free(comm2, ierr)
      call mpi_comm_free(comm, ierr)
      call mpi_finalize(ierr)
   end subroutine finalize_mpi

end module mpi_common
