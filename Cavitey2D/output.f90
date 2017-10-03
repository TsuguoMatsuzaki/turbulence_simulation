subroutine  Output( FNO )
  use Variables
 
  integer :: FNO, UNIT

  write(6,*) "F_Unit = ", FNO

  if( FNO == 0 ) then
     UNIT = 6
  else
     UNIT = FNO
     open( UNIT, file="./foo" )
  endif

  if( mod(m,2) == 0 )  then
     mh = m/2
     write(UNIT,'(1p,2e15.4)') ( dy*real(j), (u(mh,j)+u(mh,j+1))/2,j=n,0,-1)
  else
     mh = m/2
     write(UNIT,'(1p,2e15.4)') &
         ( dy*real(j),        &
           (u(mh,j)+u(mh,j+1)+u(mh+1,j)+u(mh+1,j+1))/4, j=n,0,-1)
  endif


  return

end subroutine Output
