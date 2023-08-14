subroutine EVOLUT(INI, IFI) 
  use interfaccia
  use fisica
  use strut
  use serchi
  use chimic
  use chim
  use nummod
  use zone_conv
  use overshoot
  use mesh
  use sceltachim
  use equil
  use costanti
  implicit none

  integer :: INI, IFI

  integer :: K, IFIN, J, LCAR, JF, M
  real :: UNO,DUE,TRE,quattro
  real :: P, RO, T

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="evolut        "
  

  IFIN = IFI-1 
 
  ! ciclo di aggiornamento abbondanze elementi
  do K=INI,IFIN 
     LCAR = 0 
     do J=1,MELE 
        XX(J) = XXX(J,K) 
     end do

     if(XX(1) <= 0. .and. XX(3) <= 0.)then 
        XX(1) = XSERV(1,K) 
        XX(3) = XSERV(2,K) 
        LCAR = 1 
     endif
!!$     T = 1.0d6*G(4,K) 
!!$     P = 1.0d17*G(3,K) - arad_3 * (T**4) 
        
     T =  1.0d6* (G(4,K) +  G(4,K+1)) *0.5
     P = 1.0d17* (G(3,K) +  G(3,K+1))*0.5 - arad_3 * (T**4) 
 
     call STATE(caller,P,T,RO,dum,dum,dum) 
     call EPSI(RO,T,K,1,HT1,UNO,DUE,TRE,quattro) 

     do J=1,MELE 
        XXX(J,K) = XX(J) 
     end do
    
     if(LCAR == 1) then 
        XSERV(1,K) = XXX(1,K) 
        XSERV(2,K) = XXX(3,K) 
        XXX(1,K) = 0. 
        XXX(3,K) = 0. 
     endif
  end do
  return 


end subroutine EVOLUT
