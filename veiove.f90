subroutine VEIOVE(EM,R,EL,P,T,RO,NUMERO,IPP) 
  use interfaccia
  use fisica
  use strut
  use fitt
  use chimic
  use indi
  use chim
  use numer
  use costanti
  implicit none

  real :: EM, R, EL, P, T, RO
  integer :: NUMERO, IPP

  integer :: num, k, k1, k2, it, salta 
  real :: cpp, dm, dr, dad, cappa, tt, gravi, ppp, eps, epsa
  real :: prova1, qt, dl, dp, dt, drad, grad, ecar

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="veiove        "
  
  if(NUMERO-1 <= 0) then 
     EM = 0. 
     EL = 0. 
     R = 0. 
     P = PCEN 
     T = TCEN 
     call STATE(caller,PCEN,TCEN,ROCEN,dum,CPP,dum) 
     NUM = LAST-1
  else 
     NUM = LAXME-LAST 
  endif

  do K=1,NUM 
     K2 = K 
     K1 = K2+1 
     if(NUMERO == 2) then
        K2 = LAXME-K+1 
        K1 = K2-1 
     endif
     DM = G(5,K1)-G(5,K2) 
     EM = EM+DM 
     do IT=1,MELE 
        XX(IT) = XXX(IT,K2) 
     end do
     call STATE(caller,P,T,RO,DAD,CPP,dum) 
     if(K > 1 .or. NUMERO > 1) then
        !!DR = DM/(12.566*RO*(R**2))
        DR = DM / (4.d0 * pigre * RO * (R**2)) 
     else
        !!DR = (3.*DM/(12.566*RO))**(1./3.)
        DR = (3. * DM / (4.d0 * pigre * RO))**(1./3.) 
     endif
     R = R+DR 
     if(R <= 0.) return
     call KAPPA(RO,T,CAPPA,XX(3)) 
     TT = T*1.d-6 
     GRAVI = 0. 
     PPP = P*1.d-17 
     if(IPP == 1) call EPSIG(PPP,T,RO,CPP,HT1,K1,GRAVI,LAXME,DAD,IPP) 

     call EPSI(RO,T,K1,0,0.0,EPS,EPSA,ecar,PROVA1) 
     call NEUTR(RO,T,QT) 
     EPS = EPS+QT+GRAVI 
     DL = DM*EPS 
     EL = EL+DL 
     if(EL <= 0. .and. NUMERO == 2) return
     !!DP = -6.668d-8*EM*RO*DR/(R**2)
     DP = -Ggrav * EM * RO * DR / (R**2) 
     P = P+DP 
     if(P <= 0.) return
     !!DRAD = 3.955d9*CAPPA*(P/(T**4))*(EL/EM)
     DRAD = 1.d8 * cte_grad * CAPPA * (P / (T**4)) * (EL / EM) 
     salta = 0
     if(IREAD == 2 .and. IPP == 0) then
        GRAD = DAD
        salta = 1
     endif
     if(salta == 0) then
        if(DRAD-DAD < 0) then 
           GRAD = DRAD
        else 
           GRAD = DAD
        endif
     endif
     DT = GRAD*T*DP/P 
     T = T+DT 
     if(T <= 0.) return

     if(NABLA ==  0) then 
        G(1,K1) = R 
        G(2,K1) = EL 
        G(3,K1) = P 
        G(4,K1) = T 
        G(6,K1) = DRAD-DAD 
        GG(6,K1) = DRAD-DAD
     else 
        cycle
     endif
  end do

  return 
end subroutine VEIOVE
