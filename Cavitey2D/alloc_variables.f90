subroutine  Alloc_Variables
  use Variables

  allocate( u (0:m, 0:n1) )  ! velovity in x-direction
  allocate( fu(0:m, 0:n1) )  ! RHS in x-direction

  allocate( v (0:m1,0:n ) )  ! velocity in y-direction
  allocate( fv(0:m1,0:n ) )  ! RHS in y-direction

  allocate( p(0:m1,0:n1) )   ! pressure

  allocate( a(5,mn) )   ! Coefficient matrix
  allocate( x(-m+1:mn+m) ) ! solution vector, part of x is for work area
  allocate( b(mn)   )   ! RHS of linear system


! zero clear
  u = 0.0d0
  v = 0.0d0
  p = 0.0d0

  a = 0.0d0
  x = 0.0d0
  b = 0.0d0


  return

end subroutine Alloc_Variables
