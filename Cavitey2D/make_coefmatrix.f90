subroutine  Make_CoefMatrix 
  use Variables

  ip = 0

  do j = 1, n

    if( j == 1 )  then    ! at the bottom
       cs = 0.0d0
    else 
       cs = -ody2
    endif

    if( j == n )  then    ! at the top
       cn = 0.0d0
    else
       cn = -ody2
    endif

    do i = 1, m

      if( i == 1 ) then   ! at the left boundary
         cw = 0.0d0
      else
         cw = -odx2
      endif

      if( i == m ) then   ! at the right boundary
         ce = 0.0d0
      else 
         ce = -odx2
      endif

    
      ip       = ip + 1
      a(1,ip)  = cs
      a(2,ip)  = cw
      a(3,ip)  = -(cs+cw+ce+cn)
      a(4,ip)  = ce
      a(5,ip)  = cn

      if( ip == 1 ) then 
         a(1,ip)  = 0.0
         a(2,ip)  = 0.0
         a(3,ip)  = 1.0d0
         a(4,ip)  = 0.0
         a(5,ip)  = 0.0
      endif

    enddo

  enddo


  return

end subroutine Make_CoefMatrix
