subroutine  UV_Updates
  use Variables

! update velocity in x-direction
  do j = 1, n
    do i = 1, m-1

      u(i,j) = fu(i,j) - dtodx*( p(i+1,j  ) - p(i,j) )

    end do
  end do

! update velocity in y-direction
  do j = 1, n-1
    do i = 1, m

      v(i,j) = fv(i,j) - dtody*( p(i  ,j+1) - p(i,j) )

    end do
  end do
    

  return

end subroutine UV_Updates
