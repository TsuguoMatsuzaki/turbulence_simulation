subroutine  Solve_Poisson 
  use Variables
  integer ::  iter
  double precision :: resb, tmp, omega, omg1

  
! Right hand side of Poisson Eq.
  ip = 0
  do j = 1, n
    do i = 1, m
      ip    = ip + 1
      b(ip) = - (fu(i,j)-fu(i-1,j))*odtdx - (fv(i,j)-fv(i,j-1))*odtdy
    end do
  end do
  b(1) = 0.0d0    ! Base value of pressure p=0.0


! Initial values of an iterative method of Poisson Eq.
  ip = 0
  do j = 1, n
    do i = 1, m
      ip    = ip + 1
      x(ip) = p(i,j)
    end do
  end do


  call Lin_SOR( m, mn, a, x, b, eps )


  ip = 0
  do j = 1, n
    do i = 1, m
      ip     = ip + 1
      p(i,j) = x(ip)
    end do
  end do


  return

end subroutine Solve_Poisson
