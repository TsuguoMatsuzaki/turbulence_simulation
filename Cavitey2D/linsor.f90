subroutine  Lin_SOR( m, mn, a, x, b, eps ) 
  
  double precision ::  a(5,mn), x(-m+1:mn+m), b(mn)

  integer ::  iter
  double precision :: eps
  double precision :: res, resb, tmp, omega, omg1

  
  omega = 1.8d0
  omg1  = 1.0d0 - omega


  resb = 0.0d0
  do ip = 1, mn
    resb = max( resb, abs(b(ip)) )
  enddo

  res = 0.0d0
  do ip = 1, mn
    res = max( res, abs( b(ip) -  a(1,ip)*x(ip-m)  &
                               -  a(2,ip)*x(ip-1)  &
                               -  a(3,ip)*x(ip  )  &
                               -  a(4,ip)*x(ip+1)  &
                               -  a(5,ip)*x(ip+m) )  )
  enddo
  res = res/resb 


  iter = 0

  do while( res > eps ) 

    iter = iter + 1

    do ip = 1, mn
      tmp  = ( b(ip) -  a(1,ip)*x(ip-m)  &
                       -  a(2,ip)*x(ip-1)  &
                       -  a(4,ip)*x(ip+1)  &
                       -  a(5,ip)*x(ip+m) ) / a(3,ip)
      x(ip) = omg1*x(ip) + omega*tmp
    enddo


    res = 0.0d0
    do ip = 1, mn
      res = max( res, abs( b(ip) -  a(1,ip)*x(ip-m)  &
                                 -  a(2,ip)*x(ip-1)  &
                                 -  a(3,ip)*x(ip  )  &
                                 -  a(4,ip)*x(ip+1)  &
                                 -  a(5,ip)*x(ip+m) )  )
    enddo
    res = res/resb

  enddo

  write(6,*)  iter, res

  return

end subroutine Lin_SOR
