module Variables
  implicit none

! Matrices
  integer :: m, m1, m2, n, n1, n2, mn

  double precision, allocatable :: u(:,:), v(:,:)      ! velocities
  double precision, allocatable :: p(:,:)              ! pressure

  double precision, allocatable :: fu(:,:), fv(:,:)    ! work varialbes
  double precision, allocatable :: a(:,:), x(:), b(:)  ! linear system

! Computational Region  
  double precision :: xlen, ylen
  double precision :: dx, odx, odx2, dtodx, odx2re, odtdx
  double precision :: dy, ody, ody2, dtody, ody2re, odtdy
  double precision :: u_init,  dt,   eps
  double precision :: re, cs, cw, cc, ce, cn

! Control variables
  integer :: i, j, it, ip

! Work variables						     

end module Variables
program Cavity2D
  use Variables


  call Set_Params( 64, 64, 100.0d0 )  ! # of grids, Reynolds number
  call Alloc_Variables
  call Make_CoefMatrix

  do it = 1, 1000

    call Set_Boundary
    call Calc_RHS
    call Solve_Poisson
    call UV_Updates

    if( mod(it,100) == 0 ) call Output( 0 )

  enddo

  call Output( 20 )
  
  stop

end program Cavity2D
