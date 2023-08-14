subroutine HBREZO (AMCO, YCO, CCO, NCO, OCO)
  use fisica
  use chimic
  use indi
  use chim
  use numer
  use costanti
  implicit none

  real :: AMCO, YCO, CCO, NCO, OCO

  integer,parameter,dimension(6) :: NMESH = (/20,10,5,5,3,2/)
  integer :: k, kcor, j, ind, iter, i, nme, kk, niter, n, k1, k2, l, num
  real :: core, rat, dm
  

  CORE = AMCO * 1.d33 * Msun 
  do K=1,180 
     if(G(5,K) >= CORE) exit 
  end do
  if(K > (LAXME-3)) then
     write(2,111) 
     write(66,*)'60 - hbrezo'
     write(66,111) 
     stop
  endif

  RAT = CORE/G(5,K) 
  KCOR = K 
  do J=1,K 
     G(5,J) = G(5,J)*RAT 
  end do
  IND = 0 
  do ITER=1,2 
     do I=1,3 
        IND = IND+1 
        DM = G(5,K+1)-G(5,K) 
        NME = NMESH(IND)
        DM = DM/real(NME+1) 
        KK = K+1 
        NITER = LAXME-K 
        do N=KK,LAXME 
           K1 = N+NME 
           GG(5,K1) = G(5,N) 
        end do
        do N=1,NME 
           K1 = K+N 
           G(5,K1) = G(5,K)+DM*real(N) 
        end do
        LAXME = LAXME+NME 
        do N=1,NITER 
           K1 = K+NME+N 
           G(5,K1) = GG(5,K1) 
        end do
        K = K+NME+1 
     end do
     K = KCOR 
  end do
  do L=1,K 
     XXX(1,L) = 0. 
     XXX(2,L) = 0. 
     XXX(3,L) = YCO 
     XXX(4,L) = CCO 
     XXX(5,L) = NCO 
     XXX(6,L) = OCO 
  end do
  K2 = K+1 
  do L=K2,LAXME 
     do I=1,MELE 
        XXX(I,L) = XX(I)
     end do
  end do

  if(IPRALL /= 1) return 

  write(2,222)(G(5,NUM),NUM=1,LAXME) 
  write(2,222)(XXX(1,NUM),NUM=1,LAXME) 
222 format(1P,13E10.3) 
111  format(12X,'MASSA DI CORE TROPPO GRANDE, MUOIO DISPERATO',//) 

  return 

end subroutine HBREZO
