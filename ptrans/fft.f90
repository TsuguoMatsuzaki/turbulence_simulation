module fft
   use param
   use exchange
   use timer
   implicit none

   private
   integer, parameter :: nblk = 16
   integer(8) :: pcf, pcb, prf, prb

   public :: init_fft, rft3df, rft3db

contains

   subroutine init_fft
      complex(8), dimension(n) :: dumc
      real(8), dimension(n) :: dumr

      dumc = dcmplx(0.0d0, 0.0d0)
      dumr = 0.0d0

      call dfftw_plan_dft_1d(pcf, n, dumc, dumc, -1, 1)
      call dfftw_plan_dft_1d(pcb, n, dumc, dumc,  1, 1)
      call dfftw_plan_dft_r2c_1d(prf, n, dumr, dumc, 1)
      call dfftw_plan_dft_c2r_1d(prb, n, dumc, dumr, 1)

   end subroutine init_fft

   subroutine rft3df(a)
      real(8), dimension(n, nb1, nb2) :: a
      real(8), dimension(max((nb1+1)*nb2,(nb2+2)*nb1), n) :: b, c

      call start_timer("forward_fft")

      call mrft1df(a, b)

      call start_timer("forward_ex2")
      call ex2(b, c)
      call stop_timer("forward_ex2")

      call mcft1d2f(c, a, b)

      call start_timer("forward_ex1")
      call ex1(b, c)
      call stop_timer("forward_ex1")

      call mcft1d1f(c, a)

      call stop_timer("forward_fft")

   end subroutine rft3df

   subroutine rft3db(a)
      real(8), dimension(n, nb1, nb2) :: a
      real(8), dimension(max((nb1+1)*nb2,(nb2+2)*nb1), n) :: b, c

      call start_timer("backward_fft")

      call mcft1d1b(a, b)

      call start_timer("backward_ex1")
      call ex1(b, c)
      call stop_timer("backward_ex1")

      call mcft1d2b(c, a, b)

      call start_timer("backward_ex2")
      call ex2(b, c)
      call stop_timer("backward_ex2")

      call mrft1db(c, a)

      call stop_timer("backward_fft")

   end subroutine rft3db

   subroutine mrft1df(a, b)
      real(8), dimension(n, nb1*nb2) :: a
      real(8), dimension(2, nbh2+1, nb1*nb2, np2) :: b
      real(8), dimension(2, nh+1) :: c
      integer :: ij, k, k1, ip

!$omp parallel do private(c,k1)
      do ij = 1, nb1*nb2
         call dfftw_execute_dft_r2c(prf, a(1,ij), c)
         do ip = 1, np2
            do k = 1, nbh2
               k1 = k+nbh2*(ip-1)
               b(1:2, k, ij, ip) = c(1:2, k1)/dn3
            end do
         end do
      end do

   end subroutine mrft1df

   subroutine mrft1db(a, b)
      real(8), dimension(2, nbh2+1, nb1*nb2, np2) :: a
      real(8), dimension(n, nb1*nb2) :: b
      real(8), dimension(2, nh+1) :: c
      integer :: ij, k, k1, ip

!$omp parallel do private(c,k1)
      do ij = 1, nb1*nb2
         do ip = 1, np2
            do k = 1, nbh2
               k1 = k+nbh2*(ip-1)
               c(1:2, k1) = a(1:2, k, ij, ip)
            end do
         end do
         c(1, nh+1) = 0.0d0
         c(2, nh+1) = 0.0d0
         call dfftw_execute_dft_c2r(prb, c, b(1,ij))
      end do

   end subroutine mrft1db

   subroutine mcft1d2f(a, b, c)
      real(8), dimension(2, nbh2+1, nb1, n) :: a
      real(8), dimension(2, n, nbh2, nb1) :: b
      real(8), dimension(2, nb1+1, nbh2, nb1, np1) :: c
      integer :: i, j, jj, k, kk, ip

!$omp parallel do
      do i = 1, nb1
         do kk = 1, nbh2, nblk
            do jj = 1, n, nblk
               do k = kk, min(kk+nblk-1,nbh2)
                  do j = jj, min(jj+nblk-1,n)
                     b(1:2, j, k, i) = a(1:2, k, i, j)
                  end do
               end do
            end do
            do k = kk, min(kk+nblk-1,nbh2)
               call dfftw_execute_dft(pcf, b(1,1,k,i), b(1,1,k,i))
               do ip = 1, np1
                  do j = 1, nb1
                     c(1:2, j, k, i, ip) = b(1:2, j+(ip-1)*nb1, k, i)
                  end do
               end do
            end do
         end do
      end do

   end subroutine mcft1d2f

   subroutine mcft1d1f(a, b)
      real(8), dimension(2, nb1+1, nbh2, n) :: a
      real(8), dimension(2, n, nb1, nbh2) :: b
      integer :: i, ii, j, jj, k

!$omp parallel do
      do k = 1, nbh2
         do jj = 1, nb1, nblk
            do ii = 1, n, nblk
               do j = jj, min(jj+nblk-1,nb1)
                  do i = ii, min(ii+nblk-1,n)
                     b(1:2, i, j, k) = a(1:2, j, k, i)
                  end do
               end do
            end do
            do j = jj, min(jj+nblk-1,nb1)
               call dfftw_execute_dft(pcf, b(1,1,j,k), b(1,1,j,k))
            end do
         end do
      end do

   end subroutine mcft1d1f

   subroutine mcft1d1b(a, b)
      real(8), dimension(2, n, nb1, nbh2) :: a
      real(8), dimension(2, nb1+1, nbh2, n) :: b
      integer :: i, ii, j, jj, k

!$omp parallel do
      do k = 1, nbh2
         do jj = 1, nb1, nblk
            do j = jj, min(jj+nblk-1,nb1)
               call dfftw_execute_dft(pcb, a(1,1,j,k), a(1,1,j,k))
            end do
            do ii = 1, n, nblk
               do j = jj, min(jj+nblk-1,nb1)
                  do i = ii, min(ii+nblk-1,n)
                     b(1:2, j, k, i) = a(1:2, i, j, k)
                  end do
               end do
            end do
         end do
      end do

   end subroutine mcft1d1b

   subroutine mcft1d2b(a, b, c)
      real(8), dimension(2, nb1+1, nbh2, nb1, np1) :: a
      real(8), dimension(2, n, nbh2, nb1) :: b
      real(8), dimension(2, nbh2+1, nb1, n) :: c
      integer :: i, j, jj, k, kk, ip

!$omp parallel do
      do i = 1, nb1
         do kk = 1, nbh2, nblk
            do k = kk, min(kk+nblk-1,nbh2)
               do ip = 1, np1
                  do j = 1, nb1
                     b(1:2, j+(ip-1)*nb1, k, i) = a(1:2, j, k, i, ip)
                  end do
               end do
               call dfftw_execute_dft(pcb, b(1,1,k,i), b(1,1,k,i))
            end do
            do jj = 1, n, nblk
               do k = kk, min(kk+nblk-1,nbh2)
                  do j = jj, min(jj+nblk-1,n)
                     c(1:2, k, i, j) = b(1:2, j, k, i)
                  end do
               end do
            end do
         end do
      end do

   end subroutine mcft1d2b

end module fft
