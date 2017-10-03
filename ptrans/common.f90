module common
   use param
   use mpi_common
   implicit none

   real(8), dimension(n) :: akx
   real(8), dimension(nb1) :: aky
   real(8), dimension(nbh2) :: akz
   real(8), dimension(n, nb1, nbh2) :: wk, owk, msk, mskf
   real(8) :: pi

contains

   subroutine init_common
      integer :: i, j, k

      pi = 4.0d0*atan(1.0d0)

!$omp parallel do
      do i = 1, nh
         akx(i) = dble(i-1)
         akx(i+nh) = dble(i-1-nh)
      end do

!$omp parallel do
      do j = 1, nb1
         aky(j) = akx(nb1*me1+j)
      end do

!$omp parallel do
      do k = 1, nbh2
         akz(k) = akx(nbh2*me2+k)
      end do

!$omp parallel do
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               wk(i, j, k) = akx(i)**2 + aky(j)**2 + akz(k)**2
               if (wk(i, j, k) /= 0.0d0) then
                  owk(i, j, k) = 1.0d0/wk(i, j, k)
               else
                  owk(i, j, k) = 0.0d0
               end if
               if (wk(i, j, k) < kmax2) then
                  msk(i, j, k) = 1.0d0
               else
                  msk(i, j, k) = 0.0d0
               end if
               if (fkmin2 <= wk(i, j, k) .and . wk(i, j, k) < fkmax2) then
                  mskf(i, j, k) = 1.0d0
               else
                  mskf(i, j, k) = 0.0d0
               end if
            end do
         end do
      end do

   end subroutine init_common

end module common
