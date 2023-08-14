subroutine CUB(A,B,ALFA,BETA,GAMMA,DELTA) 

  real,dimension(4) :: A,B
  real :: ALFA,BETA,GAMMA,DELTA
  real :: D1,C1,C2,D2,C3,C4,C5,C6,D3,B12

  B12 = B(1)*B(1)   
  D1 = B(2)*(B(2)+B(1)) 
  C1 = B(3)*(B(3)+B(1))-D1 
  C2 = B(3)-B(2) 
  D2 = (A(2)-A(1))/(B(2)-B(1)) 
  C3 = (A(3)-A(1))/(B(3)-B(1))-D2 
  C4 = B(4)*(B(4)+B(1))-D1 
  C5 = B(4)-B(2) 
  C6 = (A(4)-A(1))/(B(4)-B(1))-D2 
  D3 = C1*C5-C2*C4 
  ALFA = (C3*C5-C2*C6)/D3 
  BETA = (C1*C6-C4*C3)/D3 
  GAMMA = D2-ALFA*(D1+B12)-BETA*(B(2)+B(1)) 
  DELTA = A(1)-ALFA*B12*B(1)-BETA*B12-GAMMA*B(1)

  return 
end subroutine CUB


subroutine CUB2dr(A, B, idx, xABGD) 

  implicit none

  real,dimension(4) :: B
  real,dimension(4,4) :: A
  real,dimension(320) :: xABGD
  integer :: idx
  real :: D1,C1,C2,D2,C3,C4,C5,C6,D3,B12, tmp1, tmp2, tmp3, tmp4

  integer :: i, ofs

  B12 = B(1)*B(1)   
  tmp4 = B(2)+B(1)
  D1 = B(2)*tmp4 
  C1 = B(3)*(B(3)+B(1))-D1 
  C2 = B(3)-B(2) 
  C4 = B(4)*(B(4)+B(1))-D1 
  tmp1 = B(2)-B(1)
  tmp2 = B(3)-B(1)
  C5 = B(4)-B(2) 
  tmp3 = B(4)-B(1)
  D3 = C1*C5-C2*C4 

  do i =1,4
     ofs = idx+(i-1)*4
     D2 = (A(2,i)-A(1,i))/tmp1 
     C3 = (A(3,i)-A(1,i))/tmp2-D2 
     C6 = (A(4,i)-A(1,i))/tmp3-D2 
     
     xABGD(ofs+1) = (C3*C5-C2*C6)/D3 
     xABGD(ofs+2) = (C1*C6-C4*C3)/D3 
     xABGD(ofs+3) = D2-xABGD(ofs+1)*(D1+B12)-xABGD(ofs+2)*tmp4
     xABGD(ofs+4) = A(1,i)-xABGD(ofs+1)*B12*B(1)- &
          xABGD(ofs+2)*B12- xABGD(ofs+3)*B(1) 
  end do

  return 
end subroutine CUB2dr


subroutine CUB2dt(A, B, xtABGD, is, ie) 

  implicit none

  real,dimension(4) :: B
  real,dimension(80) :: xtABGD, A
  integer :: is, ie

  integer :: idx, inizio, fine
  real :: D1,C1,C2,D2,C3,C4,C5,C6,D3,B12, tmp1, tmp2, tmp3, tmp4

  integer :: i, ofs

  B12 = B(1)*B(1)   
  tmp4 = B(2)+B(1)
  D1 = B(2)*tmp4 
  C1 = B(3)*(B(3)+B(1))-D1 
  C2 = B(3)-B(2) 
  C4 = B(4)*(B(4)+B(1))-D1 
  tmp1 = B(2)-B(1)
  tmp2 = B(3)-B(1)
  C5 = B(4)-B(2) 
  tmp3 = B(4)-B(1)
  D3 = C1*C5-C2*C4 

  inizio = 2*is -1
  fine = 2*ie

  do i =inizio,fine
     ofs = 4*(i-1) + 1
     D2 = (A(ofs+1)-A(ofs))/tmp1 
     C3 = (A(ofs+2)-A(ofs))/tmp2-D2 
     C6 = (A(ofs+3)-A(ofs))/tmp3-D2 
     
     xtABGD(ofs) = (C3*C5-C2*C6)/D3 
     xtABGD(ofs+1) = (C1*C6-C4*C3)/D3 
     xtABGD(ofs+2) = D2-xtABGD(ofs)*(D1+B12)-xtABGD(ofs+1)*tmp4
     xtABGD(ofs+3) = A(ofs)-xtABGD(ofs)*B12*B(1)- &
          xtABGD(ofs+1)*B12- xtABGD(ofs+2)*B(1) 
  end do

  return 
end subroutine CUB2dt
