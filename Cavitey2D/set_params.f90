subroutine  Set_Params( arg_M, arg_N, arg_RE )
  use Variables

  integer :: arg_M, arg_N
  double precision :: arg_RE
  
  m  = arg_M
  n  = arg_N
  re = arg_RE

  m1 = m+1
  m2 = m+2

  n1 = n+1
  n2 = n+2

  mn = m*n

  xlen = 1.0d0
  ylen = 1.0d0

  u_init = 1.0d0


! Initialization of variables  
  dx   = xlen/m
  dy   = ylen/n
  
!  dt   = 0.01
  dt   = 0.25*re*dx*dx

  print '(a,e15.6)', "dt = ", 0.25*re*dx*dx

  odx    = 1.0d0/dx
  ody    = 1.0d0/dy
  odx2   = 1.0d0/(dx*dx)
  ody2   = 1.0d0/(dy*dy)
  odx2re = odx2/re
  ody2re = ody2/re

  dtodx  = dt*odx
  dtody  = dt*ody

  odtdx  = 1.0d0/dt/dx
  odtdy  = 1.0d0/dt/dy


  eps    = 1.0e-5      ! Tolerance of a linear system

  write(6,'(a,3i5)') "# of meshes: ", m, n


  return

end subroutine Set_Params
