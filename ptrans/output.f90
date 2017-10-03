module output
   use param
   use mpi_common
   use common
   implicit none

contains

   subroutine open_file(it0, itn)
      integer :: it0, itn
      character(len=11) :: ch11

      if (myrank == 0) then
         write(ch11, '(i5.5,"-",i5.5)') it0, it0+itn
         open(100, file = ddir//runname//"_"//ch11//"_stat.txt")
         write(100, '("#", a9, a10, 8a18)') "t","t/T","E","Om","ep","Skew",&
         &                                  "Re_lm","Re_L","lm","L"
      end if
   end subroutine open_file

   subroutine output_statistics(sp, it)
      real(8), dimension(4, 0:nk) :: sp
      integer :: it, kk
      real(8) :: ene, ens, skew0, skew20, eps, lam, L, skew, uu0, Re_lam, Re_L, t, tot

      ene = 0.0d0
      ens = 0.0d0
      skew0 = 0.0d0
      skew20 = 0.0d0
      L = 0.0d0
      do kk = nh, 1, -1
         ene = ene + sp(1, kk)
         ens = ens + sp(2, kk)
         skew0 = skew0 + sp(3, kk)*kk**2
         skew20 = skew20 + sp(4, kk)
         L = L + sp(1, kk)/dble(kk)
      end do

      eps = 2.0d0*nu*ens
      lam = sqrt(5.0d0*ene/ens)
      L = L*(3.0d0*pi/4.0d0/ene)
!     skew = skew0*(3.0d0/7.0d0)*sqrt(15.0d0/2.0d0/ens**3)
      skew = skew20*(3.0d0/7.0d0)*sqrt(15.0d0/2.0d0/ens**3)
      uu0 = sqrt(2.0d0*ene/3.0d0)
      Re_lam = uu0*lam/nu
      Re_L = uu0*L/nu
      t = it*dt
      tot = wp*u0*t

      if (myrank == 0) then
         write(100, '(2f10.6,8e18.10)') t, tot, ene, ens, eps, skew, Re_lam, Re_L, lam, L
      end if

   end subroutine output_statistics

   subroutine output_spectrum(sp, it)
      real(8), dimension(4, 0:nk) :: sp
      integer :: it, kk
      character(len=5) :: ch5

      if (myrank == 0) then
         write(ch5, '(i5.5)') it
         open(200, file = ddir//runname//"_"//ch5//"_sp.txt")
         do kk = 1, nk
            write(200, '(f10.3,4e18.8)') real(kk), sp(1:4, kk)
         end do
         close(200)
      end if

   end subroutine output_spectrum

   subroutine close_file
      if (myrank == 0) close(100)
   end subroutine close_file

end module output
