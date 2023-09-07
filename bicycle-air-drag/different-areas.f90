program Avariado
  implicit none
  integer i, j
  real(8) vi, P, m, t, D, ro, A, va, N
  P=400.0d0;m=70.0d0;t=300.0d0;D=0.1d0;ro=1.3d0;N=t/d;j=10

  open(10,file="a339.dat")
  open(20,file="a374.dat")
  open(30,file="a460.dat")
  open(40,file="Froome.dat")
  open(50,file="Pantani.dat")

  do
    vi=4.0d0;va=0.0d0;
    write(j,*)0.0d0, vi
    if(j.eq.10)then
      A=0.339d0
    else if(j.eq.20)then
      A=0.374d0
    else if(j.eq.30)then
      A=0.460d0
    else if(j.eq.40)then
      A=0.258d0
    else if(j.eq.50)then
      A=0.228d0
    end if

    do i=1, int(N)

      va = va + D * (xlr8(P,m,va+vi) - arrasto(ro,A,m,va+vi))

      write(j,*)i*D, va+vi

    end do
    if(j.eq.50)exit
    j = j+10
  end do

  close(10)
  close(20)
  close(30)
  close(40)
  close(50)

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
