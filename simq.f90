subroutine SIMQ(A,B,N,KS) 
  use interfaccia
  
  use nummod

  implicit none

  real,dimension(361) :: A
  real,dimension(19) :: B
  integer :: N, KS

  integer :: JJ, J, JY, IT, I, IJ, IMAX, I1, I2, IQS, IX, IXJ, IA, IB, IC
  integer :: IXJX, K, JX, JJX, NY
  real :: BIGA, savev
  real,parameter :: TOL = 0.

  KS = 0 
  JJ = -N 

  do J=1,N 
     JY = J+1 
     JJ = JJ+N+1 
     BIGA = 0. 
     IT = JJ-J 
     do I=J,N 
        IJ = IT+I 
        if(abs(BIGA)-abs(A(IJ)) < 0) then 
           BIGA = A(IJ) 
           IMAX = I
        endif
     end do
     if(abs(BIGA)-TOL <= 0) then 
        KS=1 
        return
     endif
     I1 = J+N*(J-2) 
     IT = IMAX-J 
     do K=J,N 
        I1 = I1+N 
        I2 = I1+IT 
        savev = A(I1) 
        A(I1) = A(I2) 
        A(I2) = savev 
        A(I1) = A(I1)/BIGA 
     enddo
     savev = B(IMAX) 
     B(IMAX) = B(J) 
     B(J) = savev/BIGA 
     if(J-N == 0) then 
        exit 
     endif
     IQS = N*(J-1) 
     do IX=JY,N 
        IXJ = IQS+IX 
        IT = J-IX 
        do JX=JY,N 
           IXJX = N*(JX-1)+IX 
           JJX = IXJX+IT 
           A(IXJX) = A(IXJX)-(A(IXJ)*A(JJX)) 
        end do
        B(IX) = B(IX)-(B(J)*A(IXJ)) 
     end do
  end do
  NY = N-1 
  IT = N*N 
  do J=1,NY 
     IA = IT-J 
     IB = N-J 
     IC = N 
     do K=1,J 
        B(IB) = B(IB)-A(IA)*B(IC) 
        IA = IA-N 
        IC = IC-1 
     end do
  end do
  return 
end subroutine SIMQ
