subroutine  Calc_RHS 
  use Variables


! Convection & diffusion terms of u-eq.
  do j = 1, n
    do i = 1, m-1

      wuu_w = 0.25d0*(u(i  ,j)+u(i-1,j))*(u(i  ,j)+u(i-1,j))
      wuu_e = 0.25d0*(u(i+1,j)+u(i  ,j))*(u(i+1,j)+u(i  ,j))

      wuv_s = 0.25d0*(u(i,j-1)+u(i,j  ))*(v(i,j-1)+v(i+1,j-1))
      wuv_n = 0.25d0*(u(i,j  )+u(i,j+1))*(v(i,j  )+v(i+1,j  ))
  
      conv  = odx*( wuu_e - wuu_w )  + ody*( wuv_n - wuv_s )

      diff  = odx2re*( u(i+1,j  ) - 2.0d0*u(i,j) + u(i-1,j  ) ) &
            + ody2re*( u(i  ,j+1) - 2.0d0*u(i,j) + u(i  ,j-1) )

      fu(i,j) = u(i,j) + dt*( diff - conv )

    end do
  enddo


! Convection & diffusion terms of u-eq.
  do j = 1, n-1
    do i = 1, m

      wvu_w = 0.25d0*(v(i-1,j)+v(i  ,j))*(u(i-1,j)+u(i-1,j+1))
      wvu_e = 0.25d0*(v(i  ,j)+v(i+1,j))*(u(i  ,j)+u(i  ,j+1))

      wvv_s = 0.25d0*(v(i,j-1)+v(i,j  ))*(v(i,j-1)+v(i,j  ))
      wvv_n = 0.25d0*(v(i,j  )+v(i,j+1))*(v(i,j  )+v(i,j+1))
  
      conv  = odx*( wvu_e - wvu_w )  + ody*( wvv_n - wvv_s )

      diff  = odx2re*( v(i+1,j  ) - 2.0d0*v(i,j) + v(i-1,j  ) ) &
            + ody2re*( v(i  ,j+1) - 2.0d0*v(i,j) + v(i  ,j-1) )

      fv(i,j) = v(i,j) + dt*( diff - conv )

    end do
  end do


  return

end subroutine Calc_RHS
