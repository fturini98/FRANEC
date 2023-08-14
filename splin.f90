! questa routine costruisce una spline cubica "naturale",
! ossia con derivata seconda nulla sul bordo.
! in input accetta:
! le coordinate dei punti del grigliato: (x, y)
! il numero di punti: n
! in output ritorna i valori delle derivate seconde: y2 

! La funzione e' strutturata in modo da vettorializzare 
! la maggior parte delle istruzioni. A questo scopo i tre loop
! iniziali sono stati separati.
subroutine makespline(x,y,n,y2)
  implicit none

  integer :: n
  real,dimension(n) :: x,y,y2
  integer,parameter :: NMAX=500
  integer :: i,k
  real,parameter :: qn = 0., un = 0.
  real :: sig,p
  real,dimension(NMAX) :: u, dx1, dx2, dy1, r

  y2(1) = 0.
  u(1) = 0.

  do i=2,n-1
     dx2(i) = x(i+1)-x(i-1)
  end do

  do i=2,n
     dx1(i) = x(i)-x(i-1)
     dy1(i) = y(i)-y(i-1)
     r(i) = dy1(i)/dx1(i)
  end do

  do i=2,n-1
     sig = dx1(i)/dx2(i)
     p = 1./(sig*y2(i-1)+2.)
     y2(i) = (sig-1.)*p
     u(i) = (6.*(r(i+1)-r(i))/dx2(i)-sig*u(i-1))*p
  end do

  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
     y2(k) = y2(k)*y2(k+1)+u(k)
  end do
  return
end subroutine makespline

! questa routine utilizza i valori delle derivate calcolate da makespline
! per calcolare il valore della spline nel punto di coordinata x
! in input accetta: 
! le coordinate dei punti del grigliato: (xa, ya)
! il numero di punti: n
! i valori delle derivate seconde: y2a
! in output ritorna il valore calcolato: y 
subroutine usesplint(xa,ya,y2a,n,x,y)
  implicit none
  
  integer :: n
  real,dimension(n) :: xa,y2a,ya
  integer :: k,khi,klo
  real ::  a,b,h, x,y

  klo = 1
  khi = n
  do while(khi-klo > 1)
     k = (khi+klo)/2
     if(xa(k) > x) then
        khi = k
     else
        klo = k
     endif
  end do

  h = xa(khi)-xa(klo)
!  if (h == 0.) stop 'bad xa input in splint'
  a = (xa(khi)-x)/h
  b = (x-xa(klo))/h
  y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  return
end subroutine usesplint
