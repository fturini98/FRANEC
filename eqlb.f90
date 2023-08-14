subroutine EQLB(MAXME) 
  use fisica
  use equil
  use chimic
  use overshoot
  use nummod

  implicit none

  integer :: maxme

  real :: SCN = 0., SCNO = 0.
  integer :: inin, ifin, j, k
  real :: a, b, c, x2equi, x4equi, x5equi, x6equi 

  ININ = 1 
  IFIN = 1 
  do J=2,MAXME-1 
     if(G(6,J) >= 0. .and. G(6,J-1) < 0.) ININ = J 
     if(G(6,J) >= 0. .and. G(6,J+1) < 0.) then 
        IFIN = J 
        ! Inserimento overshooting
        if(KOVER == 1 .and. G(6,1) >= 0.0 .and. NMD > 1) IFIN = L3
        exit
     endif
  end do
  if(IFIN == 1) return 
  A = COEFF(1) 
  ! controllo una possibile divisione per 0
  if(abs(A) < 1.0d-60) then    ! non dovrei essere qui...
     write(67,*) "In EQLB ho i coeff nulli. Torno al chiamante."
     return
  endif
  B = XXX(3,1)/4.*COEFF(2) 
  C = -XXX(1,1)**2*COEFF(3)/2. 
  X2EQUI = 3.*(-B+sqrt(B**2-4.*A*C))/(2.*A) 
  if( .not. (XXV(2,1) <= XXX(2,1) .and. XXX(2,1) < X2EQUI) ) then
     if( .not. (XXV(2,1) >= XXX(2,1) .and. XXX(2,1) > X2EQUI) ) then
        do K=ININ,IFIN 
           XXX(2,K) = X2EQUI 
        end do
     endif
  endif

  A = COEFF(4) 
  B = COEFF(5) 
  SCNO = XXX(4,1)/12.+XXX(5,1)/14.+XXX(6,1)/16. 
  SCN = SCNO-XXX(6,1)/16.
  X5EQUI = 14.*SCN/(1.+B/A) 
  X4EQUI = 12./14.*X5EQUI*B/A 
  X6EQUI = 16.*(SCNO - X5EQUI/14. - X4EQUI/12.)
  if( XXX(4,1) > X4EQUI .and. XXV(4,1) < X4EQUI ) then 
     do K=ININ,IFIN 
        XXX(5,K) = X5EQUI 
        XXX(4,K) = X4EQUI 
        xxx(6,k) = X6EQUI
     end do
  endif
  if( XXX(4,1) < X4EQUI .and. XXV(4,1) > X4EQUI ) then 
     do K=ININ,IFIN 
        XXX(5,K) = X5EQUI 
        XXX(4,K) = X4EQUI
        xxx(6,k) = X6EQUI
     end do
  endif
  if( XXX(5,1) > X5EQUI .and. XXV(5,1) < X5EQUI ) then 
     do K=ININ,IFIN 
        XXX(5,K) = X5EQUI 
        XXX(4,K) = X4EQUI 
        xxx(6,k) = X6EQUI
     end do
  endif
  if( XXX(5,1) < X5EQUI .and. XXV(5,1) > X5EQUI ) then 
     do K=ININ,IFIN 
        XXX(5,K) = X5EQUI 
        XXX(4,K) = X4EQUI 
        xxx(6,k) = X6EQUI
     end do
  endif
  
  return 

end subroutine EQLB
