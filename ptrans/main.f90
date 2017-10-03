program main
   use param
   use mpi_common
   use common
   use io_fields
   use output
   use fft
   use spectrum
   use nonlinear
   use integration
   use timer
   implicit none

   real(8), dimension(2, n, nb1, nbh2) :: u, v, w
   real(8), dimension(2, n, nb1, nbh2) :: hu, hv, hw
   real(8), dimension(4, 0:nk) :: sp
   integer :: it, it0, itn, itsta, itspe, itfld

   it0 = 6450
   itn = 10
   itsta = 1
   itspe = 1
   itfld = 11

   call init_mpi
   call init_common
   call init_fft
   call open_file(it0, itn)

   call read_fields(u, v, w, it0)

   call start_timer("main")

   call calc_nonlinear(u, v, w, hu, hv, hw)

   call calc_spectrum(u, v, w, hu, hv, hw, sp)
   call output_spectrum(sp, it0)
   call output_statistics(sp, it0)

   do it = it0+1, it0+itn

      call integrate(u, v, w, hu, hv, hw)

      if (mod(it-it0,itsta) == 0) then
         call calc_spectrum(u, v, w, hu, hv, hw, sp)

         if (mod(it-it0,itspe) == 0) then
            call output_spectrum(sp, it)
         end if

         call output_statistics(sp, it)
      end if

      if (mod(it-it0,itfld) == 0) then
         call write_fields(u, v, w, it)
      end if

   end do

   call stop_timer("main")

   call close_file

   call print_times

   call finalize_mpi

end program main
