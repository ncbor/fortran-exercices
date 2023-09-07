program graficar_diagrama
  implicit none
  real(8)r,h,x,xi,q,xs
  integer i
  h=1.0d-4;r=1.0d0+h

  open(10,file='db.dat')
  write(10,*)0.0d0,0.0d0
  write(10,*)1.0d0,0.0d0
  do while(r.le.4.0d0)
    if(r.gt.1.0d0.and.r.le.3.0d0)then
      x=1.0d0-1.0d0/r
      write(10,*)r,x
    elseif(r.gt.3.0d0.and.r.le.1.0d0+sqrt(6.0d0))then
      q = sqrt(-3.0d0-2.0d0*r+r**2.0d0)
      x =1.0d0/(2.0d0*r)*(1.0d0+r+q)
      write(10,*)r,x
      x=1.0d0/(2.0d0*r)*(1.0d0+r-q)
      write(10,*)r,x
    elseif(r.gt.1.0d0+sqrt(6.0d0).and.r.le.4.0d0)then
      x=0.5
      do i=1,2000
        x = mapa(x,r)
      end do
      xs = x
      do i =1,100
        write(10,*)r,x
        x=mapa(x,r)
        if(abs(x-xs).lt.1.0d-3)exit
      end do
    end if
    r=r+h
  end do
  close(10)

contains
  real(8)function mapa(x,r)
    real(8)x,r
    mapa = r*x*(1-x)
  end function
end program
