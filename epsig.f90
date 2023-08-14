subroutine EPSIG(P,T,RO,CP,HT,JF,GRAVI,MAXMV,DAD,IPP) 
  use interfaccia
  use fisica
  use chimic
  use chim
  use numer
  use nummod

  implicit none

  real :: P, T, RO, CP, HT, GRAVI, DAD
  integer :: JF, MAXMV, IPP

  real :: prad, pr, emas, rat, dp, dt, cpma, pmuin, pmuiv, pmuen, pmuev
  real :: pmun, pmuv, ppunto, tpunto, vpunto, pmu_new, p_old
  real :: t_old, pmu_old, vpunto_new, prad_old
  real,dimension(MELE) :: xxsave
  integer :: l, k1, i, isup, iinf

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="epsig         "

  integer :: IBLKG = 0  !! mettere a 1 per bloccare i pulsi con Egrav

  if(IBLKG == 1 .and. XXX(1,1) < 1e-6 .and. XXX(3,1) < 0.15) return 0 

  PRAD = 2.52d-15*(T**4) 
  PR = P*1.d17 
  EMAS = (G(5,JF)+G(5,JF+1))/2. 
  if(IPP == 1) EMAS = EMAS*1.d-33 

  if(GG(5,MAXMV) < EMAS) then 
     L = MAXMV
  else if (GG(5,2) > EMAS) then 
     L = 2 
  else 
     ISUP = MAXMV 
     IINF = 2 
     do while (ISUP - IINF > 1)
        L = (ISUP + IINF)/2 
        if(GG(5,L) > EMAS) then 
           ISUP = L 
        else 
           IINF = L 
        endif
     end do
     L = ISUP
  end if

  K1 = L-1 
  RAT = (EMAS-GG(5,K1))/(GG(5,L)-GG(5,K1)) 
  DP = PR-(GG(3,K1)+RAT*(GG(3,L)-GG(3,K1)))*1.d17 
  DT = T-(GG(4,K1)+RAT*(GG(4,L)-GG(4,K1)))*1.d6 

  PPUNTO = DP*T*CP*DAD/PR 
  TPUNTO = -CP*DT

  GRAVI = (PPUNTO+TPUNTO)/(3.1558d7*HT)
  return

end subroutine EPSIG
