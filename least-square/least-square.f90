program minq
  implicit none

  real(8) xi, yi, Sx, Sy, Sxx, Sxy, a, b, D
  integer ios, N
  character(40) arg

  n = 0; Sx = 0.0d0; Sy = 0.0d0; Sxx = 0.0d0; Sxy = 0.0d0;
  call getarg(1, arg)
  open(10, file= arg)

  do
    read(10,*,iostat=ios) xi, yi
    if (ios/=0) EXIT
    N = N + 1
    Sx = Sx + xi
    Sy = Sy + yi
    Sxx = Sxx + xi**2.0d0
    Sxy = Sxy + xi*yi
  end do

   D = (real(N,8) * Sxx) - (Sx)**2.0d0
   a = (real(N,8) * Sxy - Sx*Sy) * D**(-1)
   b = (Sxx*Sy - Sx*Sxy) * D**(-1)

   close(10)

   write(*,'(/,A,F19.16,A,F19.16,/)')"       ", a,"   ", b

end program
