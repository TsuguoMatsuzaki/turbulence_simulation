module timer
   use mpi_common
   implicit none

   private
   integer, parameter :: max_ntimer = 20
   integer :: ntimer = 0
   real(8), dimension(max_ntimer) :: times = 0.0d0
   character(len=20), dimension(max_ntimer) :: names = "*"
   logical, dimension(max_ntimer) :: running = .false.

   public :: start_timer
   public :: stop_timer
   public :: print_times

contains

   subroutine start_timer(name)
      character(len=*) :: name
      integer :: itimer

      do itimer = 1, ntimer+1
         if (itimer == ntimer+1) then
            if (itimer > max_ntimer) then
               call print_error(name, 1, "start_timer")
               return
            end if
            names(itimer) = name
            ntimer = ntimer + 1
            exit
         end if
         if (trim(names(itimer)) == name) exit
      end do

      if (running(itimer)) then
         call print_error(name, 3, "start_timer")
         return
      end if

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      times(itimer) = times(itimer) - mpi_wtime()

      running(itimer) = .true.

   end subroutine start_timer

   subroutine stop_timer(name)
      character(len=*) :: name
      integer :: itimer

      do itimer = 1, ntimer+1
         if (itimer == ntimer+1) then
            call print_error(name, 2, "stop_timer")
            return
         end if
         if (trim(names(itimer)) == name) exit
      end do

      if (.not.running(itimer)) then
         call print_error(name, 4, "stop_timer")
         return
      end if

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      times(itimer) = times(itimer) + mpi_wtime()

      running(itimer) = .false.

   end subroutine stop_timer

   subroutine print_times
      integer :: itimer

      if (myrank == 0) then
         do itimer = 1, ntimer
            if (running(itimer)) then
               call print_error(trim(names(itimer)), 3, "print_times")
               cycle
            end if
            write(*, '("Time of timer: ", a, " is ", f10.4)') trim(names(itimer)), times(itimer)
         end do
      end if

   end subroutine print_times

   subroutine print_error(timer_name, err_num, sub_name)
      character(len=*) :: timer_name, sub_name
      integer :: err_num

      select case (err_num)
         case (1)
            write(0, '("[Error in ", a, " at rank ", i5, "] Exceeded the maximum number of timers.")') &
            &    sub_name, myrank
         case (2)
            write(0, '("[Error in ", a, " at rank ", i5, "] Timer: ", a, " dose not exist.")') &
            &    sub_name, myrank, timer_name
         case (3)
            write(0, '("[Error in ", a, " at rank ", i5, "] Timer: ", a, " was not stoped.")') &
            &    sub_name, myrank, timer_name
         case (4)
            write(0, '("[Error in ", a, " at rank ", i5, "] Timer: ", a, " was not started.")') &
            &    sub_name, myrank, timer_name
      end select

   end subroutine

end module timer
