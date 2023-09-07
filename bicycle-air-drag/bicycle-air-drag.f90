program bicicleta
  implicit none
  integer i
  real(8) vi, P, m, t, D, ro, A, vv, va, N
  vi=4.0d0;P=400.0d0;m=70.0d0;t=300.0d0;D=0.1d0;ro=1.3d0;A=1.0d0/3.0d0;N=t/d
  vv=0.0d0;va=0.0d0;

  open(10,file='vacuo.dat')
  open(20,file='arrasto.dat')

  write(10,*)0.0d0, vi
  write(20,*)0.0d0, vi

  do i=1, int(N)

    vv = vv + D * xlr8(P,m,vv+vi)
    va = va + D * (xlr8(P,m,va+vi) - arrasto(ro,A,m,va+vi))

    write(10,*)i*D, vv+vi
    write(20,*)i*D, va+vi

  end do

  close(10)
  close(20)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8) function xlr8(P,m,v)
    real(8) P, m, v

    xlr8 = P/(m*v)
  end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  real(8) function arrasto(ro,A,m,v)
    real(8) ro, A, m, v

    arrasto = (ro*(v**2.0d0)*A)/(m*2.0d0)
  end function

end program
