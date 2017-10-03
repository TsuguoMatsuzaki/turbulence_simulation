module nonlinear
   use param
   use mpi_common
   use common
   use fft
   use timer
   implicit none

contains

   subroutine calc_nonlinear(u, v, w, hu, hv, hw)
      real(8), dimension(2, n, nb1, nbh2) :: u, v, w
      real(8), dimension(2, n, nb1, nbh2) :: hu, hv, hw
      real(8), dimension(2, n, nb1, nbh2) :: uu, vv, ww
      real(8), dimension(2, n, nb1, nbh2) :: uv, uw, vw
      real(8) :: sk, ssk, csk, akhr, akhi
      integer :: i, j, k
      integer, dimension(MPI_STATUS_SIZE) :: ist

      if (myrank < np) then
         call start_timer("nonlinear1")
!$omp parallel do
         do k = 1, nbh2
            do j = 1, nb1
               do i = 1, n
                  uu(1:2, i, j, k) = u(1:2, i, j, k)/sqrt(2.0d0)
                  vv(1:2, i, j, k) = v(1:2, i, j, k)/sqrt(2.0d0)
                  ww(1:2, i, j, k) = w(1:2, i, j, k)/sqrt(2.0d0)
               end do
            end do
         end do

         call calc_convolution(uu, vv, ww, uv, uw, vw)

!$omp parallel do
         do k = 1, nbh2
            do j = 1, nb1
               do i = 1, n
                  hu(1, i, j, k) = akx(i)*uu(2, i, j, k) + aky(j)*uv(2, i, j, k) &
                  &                                      + akz(k)*uw(2, i, j, k)
                  hu(2, i, j, k) = akx(i)*uu(1, i, j, k) + aky(j)*uv(1, i, j, k) &
                  &                                      + akz(k)*uw(1, i, j, k)
                  hv(1, i, j, k) = akx(i)*uv(2, i, j, k) + aky(j)*vv(2, i, j, k) &
                  &                                      + akz(k)*vw(2, i, j, k)
                  hv(2, i, j, k) = akx(i)*uv(1, i, j, k) + aky(j)*vv(1, i, j, k) &
                  &                                      + akz(k)*vw(1, i, j, k)
                  hw(1, i, j, k) = akx(i)*uw(2, i, j, k) + aky(j)*vw(2, i, j, k) &
                  &                                      + akz(k)*ww(2, i, j, k)
                  hw(2, i, j, k) = akx(i)*uw(1, i, j, k) + aky(j)*vw(1, i, j, k) &
                  &                                      + akz(k)*ww(1, i, j, k)
               end do
            end do
         end do

         call start_timer("sendrecv1")
         call mpi_sendrecv(hu, n*nb1*nb2, MPI_REAL8, myrank+np, 0, &
         &                 uu, n*nb1*nb2, MPI_REAL8, myrank+np, 0, MPI_COMM_WORLD, ist, ierr)
         call mpi_sendrecv(hv, n*nb1*nb2, MPI_REAL8, myrank+np, 0, &
         &                 vv, n*nb1*nb2, MPI_REAL8, myrank+np, 0, MPI_COMM_WORLD, ist, ierr)
         call mpi_sendrecv(hw, n*nb1*nb2, MPI_REAL8, myrank+np, 0, &
         &                 ww, n*nb1*nb2, MPI_REAL8, myrank+np, 0, MPI_COMM_WORLD, ist, ierr)
         call stop_timer("sendrecv1")

         call stop_timer("nonlinear1")
      else
         call start_timer("nonlinear2")
!$omp parallel do private(sk,ssk,csk)
         do k = 1, nbh2
            do j = 1, nb1
               do i = 1, n
                  sk = (akx(i)+aky(j)+akz(k))*pi/n
                  ssk = sin(sk)/sqrt(2.0d0)
                  csk = cos(sk)/sqrt(2.0d0)

                  uu(1, i, j, k) = csk*u(1, i, j, k) - ssk*u(2, i, j, k)
                  uu(2, i, j, k) = ssk*u(1, i, j, k) + csk*u(2, i, j, k)
                  vv(1, i, j, k) = csk*v(1, i, j, k) - ssk*v(2, i, j, k)
                  vv(2, i, j, k) = ssk*v(1, i, j, k) + csk*v(2, i, j, k)
                  ww(1, i, j, k) = csk*w(1, i, j, k) - ssk*w(2, i, j, k)
                  ww(2, i, j, k) = ssk*w(1, i, j, k) + csk*w(2, i, j, k)
               end do
            end do
         end do

         call calc_convolution(uu, vv, ww, uv, uw, vw)

!$omp parallel do private(sk,ssk,csk,akhr,akhi)
         do k = 1, nbh2
            do j = 1, nb1
               do i = 1, n
                  sk = (akx(i)+aky(j)+akz(k))*pi/n
                  ssk = sin(sk)
                  csk = cos(sk)

                  akhr = akx(i)*uu(1,i,j,k)+aky(j)*uv(1,i,j,k)+akz(k)*uw(1,i,j,k)
                  akhi = akx(i)*uu(2,i,j,k)+aky(j)*uv(2,i,j,k)+akz(k)*uw(2,i,j,k)
                  uu(1,i,j,k) = -ssk*akhr + csk*akhi
                  uu(2,i,j,k) =  csk*akhr + ssk*akhi
                  akhr = akx(i)*uv(1,i,j,k)+aky(j)*vv(1,i,j,k)+akz(k)*vw(1,i,j,k)
                  akhi = akx(i)*uv(2,i,j,k)+aky(j)*vv(2,i,j,k)+akz(k)*vw(2,i,j,k)
                  vv(1,i,j,k) = -ssk*akhr + csk*akhi
                  vv(2,i,j,k) =  csk*akhr + ssk*akhi
                  akhr = akx(i)*uw(1,i,j,k)+aky(j)*vw(1,i,j,k)+akz(k)*ww(1,i,j,k)
                  akhi = akx(i)*uw(2,i,j,k)+aky(j)*vw(2,i,j,k)+akz(k)*ww(2,i,j,k)
                  ww(1,i,j,k) = -ssk*akhr + csk*akhi
                  ww(2,i,j,k) =  csk*akhr + ssk*akhi
               end do
            end do
         end do

         call start_timer("sendrecv2")
         call mpi_sendrecv(uu, n*nb1*nb2, MPI_REAL8, myrank-np, 0, &
         &                 hu, n*nb1*nb2, MPI_REAL8, myrank-np, 0, MPI_COMM_WORLD, ist, ierr)
         call mpi_sendrecv(vv, n*nb1*nb2, MPI_REAL8, myrank-np, 0, &
         &                 hv, n*nb1*nb2, MPI_REAL8, myrank-np, 0, MPI_COMM_WORLD, ist, ierr)
         call mpi_sendrecv(ww, n*nb1*nb2, MPI_REAL8, myrank-np, 0, &
         &                 hw, n*nb1*nb2, MPI_REAL8, myrank-np, 0, MPI_COMM_WORLD, ist, ierr)
         call stop_timer("sendrecv2")
         call stop_timer("nonlinear2")

      end if

      call start_timer("nonlinear1")
      call start_timer("nonlinear2")
!$omp parallel do private(sk,ssk,csk,akhr,akhi)
      do k = 1, nbh2
         do j = 1, nb1
            do i = 1, n
               hu(1:2, i, j, k) = hu(1:2, i, j, k) + uu(1:2, i, j, k)
               hv(1:2, i, j, k) = hv(1:2, i, j, k) + vv(1:2, i, j, k)
               hw(1:2, i, j, k) = hw(1:2, i, j, k) + ww(1:2, i, j, k)

               akhr = (akx(i)*hu(1,i,j,k)+aky(j)*hv(1,i,j,k) &
               &                         +akz(k)*hw(1,i,j,k))*owk(i,j,k)
               akhi = (akx(i)*hu(2,i,j,k)+aky(j)*hv(2,i,j,k) &
               &                         +akz(k)*hw(2,i,j,k))*owk(i,j,k)

               hu(1, i, j, k) = ( hu(1, i, j, k) - akhr*akx(i))*msk(i, j, k)
               hu(2, i, j, k) = (-hu(2, i, j, k) + akhi*akx(i))*msk(i, j, k)
               hv(1, i, j, k) = ( hv(1, i, j, k) - akhr*aky(j))*msk(i, j, k)
               hv(2, i, j, k) = (-hv(2, i, j, k) + akhi*aky(j))*msk(i, j, k)
               hw(1, i, j, k) = ( hw(1, i, j, k) - akhr*akz(k))*msk(i, j, k)
               hw(2, i, j, k) = (-hw(2, i, j, k) + akhi*akz(k))*msk(i, j, k)
            end do
         end do
      end do
      call stop_timer("nonlinear1")
      call stop_timer("nonlinear2")

   end subroutine calc_nonlinear

   subroutine calc_convolution(uu, vv, ww, uv, uw, vw)
      real(8), dimension(n, nb1, nb2) :: uu, vv, ww
      real(8), dimension(n, nb1, nb2) :: uv, uw, vw
      integer :: i, j, k

      call rft3db(uu)
      call rft3db(vv)
      call rft3db(ww)

!$omp parallel do
      do k = 1, nb2
         do j = 1, nb1
            do i = 1, n
               uv(i, j, k) = uu(i, j, k)*vv(i, j, k)
               uw(i, j, k) = uu(i, j, k)*ww(i, j, k)
               vw(i, j, k) = vv(i, j, k)*ww(i, j, k)
               uu(i, j, k) = uu(i, j, k)*uu(i, j, k)
               vv(i, j, k) = vv(i, j, k)*vv(i, j, k)
               ww(i, j, k) = ww(i, j, k)*ww(i, j, k)
            end do
         end do
      end do

      call rft3df(uu)
      call rft3df(uv)
      call rft3df(uw)
      call rft3df(vv)
      call rft3df(vw)
      call rft3df(ww)

   end subroutine calc_convolution

end module nonlinear
