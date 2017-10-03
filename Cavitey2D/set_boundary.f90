subroutine Set_Boundary
  use Variables

! velocities

  do j = 1, n                  ! u=0 at both sides
    u(0,j) = 0.0d0
    u(m,j) = 0.0d0
  end do

  do i = 0, m                  ! u=0 at bottom, u=1 at top
    u(i,0)  = -u(i,1)
    u(i,n1) = 2.0d0*u_init - u(i,n)
  end do


  do i = 1, m                  ! v=0 at top and bottom
    v(i,0) = 0.0d0
    v(i,n) = 0.0d0
  enddo

  do j = 0, n                  ! v=0 at both sides
    v(0, j) = -v(1,j)
    v(m1,j) = -v(m,j)
  end do


  return

end subroutine Set_Boundary
