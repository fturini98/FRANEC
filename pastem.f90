subroutine PASTEM(MAXME,MAXMV,PROV,PROVV,ECNO,ECNOV,FG6,FG6V, TMAX, &
     div_pastem, max_pastem)
  use interfaccia
  use fisica
  use strut
  use nuconv
  use chimic
  use tempe
  use chim
  use numer
  use nummod

  implicit none
  
  integer :: MAXME, MAXMV
  real :: PROV, PROVV, ECNO, ECNOV, FG6, FG6V, TMAX, div_pastem, max_pastem
  
  integer :: maxmo, k, l, k1, ip, it
  real :: tim, timr, timl, timma, timx, timp, timt, timhel, timcno, dmdt
  real :: dldtm, drdtm, dxdtm, dipv, ditv, rat, pp, tt, dip, dit, value
  real :: valp, valx, dpdtm, dtdtm, vall, ppro, dehel, timxn
  real :: cohel, pfg6, decar, ecn, decno, ht2, ht3, cocar, timcar

  ! questo blocco conta il numero di passi fatti
  ! serve per far decadere il divisore del passo temporale
  integer, save :: contatore, firsttime = 1
  logical :: precision ! accresciuta precisione H burn

  if(firsttime == 1) then
     contatore = nmd
     firsttime = 0
  endif

  HT1V = HT1 

  if(NMD > 10 .and. ETAFPR > 0 .and. iread == 4) then
     HT1 = 1.d3
     return
  endif

  TIM = 1.d10/(EMTOT*EMTOT*G(4,1)) 
  if(XXX(3,1) <= 0.) TIM = 1.d8/(EMTOT*EMTOT*EMTOT) 
  TIMR = TIM 
  TIML = TIM 
  TIMMA = TIM 
  TIMX = TIM 
  TIMP = TIM 
  TIMT = TIM 
  TIMHEL = TIM 
  TIMCNO = TIM 
  DMDT = (G(5,MAXME)-GG(5,MAXMV))/HT1V 
  DLDTM = log10(G(2,MAXME)/GG(2,MAXMV))/HT1V 
  DRDTM = (TEFF-TEFFV)/HT1V 
  !******* CALCOLO DI DXDTM A SECONDA DELLA COMBUSTIONE CENTRALE *********
  DXDTM = (XCE(1,2)-XCE(1,1))/HT1V 
  if(XXX(1,1) <= 0.) DXDTM = (XCE(3,2)-XCE(3,1))/HT1V 
  if(XXX(1,1) <= 0. .and. XXX(3,1) <= 0.) DXDTM = (XCE(4,2)-XCE(4,1))/HT1V
  !******* CONTROLLA VARIAZIONI MASSIME DI R - L - P - T *****************
  DIPV = 0. 
  DITV = 0. 
  MAXMO = MAXME-3 
  do K=3,MAXMO 
     do L=2,MAXMV 
        if(G(5,K)-GG(5,L) < 0) then 
           exit
        else if(G(5,K)-GG(5,L) == 0) then 
           K1 = L-1 
           if(XXX(1,1) <= 0. .and. XXX(3,1) <= 0. .and. G(4,K) < 20.) exit
           if(XXX(1,1) <= 0. .and. XXX(3,1) > 0.9 .and. G(4,K) < 2.) exit
           if( (G(6,K-2)*G(6,K)) <= 0. ) exit
           if( (G(6,K-1)*G(6,K)) <= 0. ) exit
           if( (G(6,K+1)*G(6,K)) <= 0. ) exit 
           if( (G(6,K+2)*G(6,K)) <= 0. ) exit 
           if( (G(6,K)*GG(6,L)) <= 0. ) exit 
           if( G(6,K) > 0.9 .and. G(6,K) < 1.1 ) exit
           RAT = (G(5,K)-GG(5,K1))/(GG(5,L)-GG(5,K1)) 
           PP = GG(3,K1)+RAT*(GG(3,L)-GG(3,K1)) 
           TT = GG(4,K1)+RAT*(GG(4,L)-GG(4,K1)) 
           DIP = abs( ( PP - G(3,K) ) / PP ) * 100. 
           DIT = abs( ( TT - G(4,K) ) / TT ) * 100. 
           if(DIP-DIPV > 0) then 
              DIPV = DIP 
              IP = K 
           endif
           if(DIT-DITV <= 0) then 
              exit
           else 
              DITV = DIT 
              IT = K 
           endif
        else 
           cycle
        endif
     end do
  end do
  !***********************************************************************
  VALUE = 3. 
  if( TMAX >= 550. ) VALUE = 2. 
  !***********************************************************************
  if( XXX(1,1) > 0. .or. XXX(3,1) > 0. ) then 
     VALP = 10. 
  else 
     VALP = 7. 
  endif
  if(XXX(1,1) <= 0. .and. EMTOT > 24. .and. XXX(3,1) > 8.d-1) then 
     VALP = 5. 
     VALUE = 2. 
  endif
  !***********************************************************************
  if(XXX(1,1) > 0.) then 
     if(G(6,1) >= 0.) then 
        VALX = 0.4-(1. - XXX(1,1))/5. 
     else 
        VALX = 0.15-(1. - XXX(1,1))/10. 
     endif
  endif
  !***********************************************************************
  if(XXX(1,1) <= 0. .and. XXX(3,1) > 0.) then 
     VALX = 0.4-(1. - XXX(3,1))/5. 
  endif
  !***********************************************************************
  if(XXX(1,1) <= 0. .and. XXX(3,1) <= 0.) then 
     if(G(6,1) >= 0.) then 
        VALX = 0.4-(1.-XXX(4,1))/5. 
     else 
        VALX = 0.15-(1.-XXX(4,1))/10. 
     endif
  endif
  DPDTM = DIPV/HT1V 
  DTDTM = DITV/HT1V 
  !*******L=.03-R=.002-X=.01-P=.1-T=.03********************************** 
  VALL = .002 
  if(G(4,1) < 8.) VALL = .02 
  if(XXX(1,1) <= 0.) VALL = .01 
  !********************************************************************** 
  if(abs(DMDT) > 0.) TIMMA = abs(.0005*G(5,MAXME)/DMDT) 
  if(abs(DLDTM) > 0.) TIML = abs(VALL/DLDTM) 
  if(abs(DRDTM) > 0.) TIMR = abs(.002*TEFF/DRDTM) 
  !********************************************************************** 
  if(abs(DXDTM) > 0. .and. XXX(1,1) > 0.) then 
     if(XXX(1,1) > 1.d-3) then
        TIMX = abs(VALX*XCE(1,2)/DXDTM) 
     else
        TIMX = abs(VALX*1.d-3/DXDTM) 
     endif

     !! ema: limito il passo temporale in modo da avere una variazione su Xc
     !! inferiore al 1%
     precision = .false.
     if(precision) then
        TIMXn = TIMX
        if(XXX(1,1) >= 1.0d-10) then
           TIMXn = 0.05 * XXX(1,1) / abs(dXdtm)
        else if(XXX(1,1) >= 1.0d-20) then
           TIMXn = 0.10 * XXX(1,1) / abs(dXdtm)
        else if(XXX(1,1) >= 1.0d-30) then
           TIMXn = 0.20 * XXX(1,1) / abs(dXdtm)
        endif
        !
        if(TIMX > TIMXn .and. nmd >= 10) then
           TIMX = TIMXn
           !   write(*,'("  Tempo-H: ", 2(1x,1pe13.6),0p)') TIMX, TIMXn
        endif
     endif

  endif
  if(abs(DXDTM) > 0. .and. XXX(1,1) <= 0. .and. XXX(3,1) > 0.) then 
     if(XXX(3,1) > 1.0d-3) then
        TIMX = abs(VALX*XCE(3,2)/DXDTM) 
     else
        TIMX = abs(VALX*1.d-3/DXDTM) 
     endif
  endif
  if(abs(DXDTM) > 0. .and. XXX(1,1) <= 0. .and. XXX(3,1) <= 0.) then 
     if(XXX(4,1) > 1.0d-3) then
        TIMX = abs(VALX*XCE(4,2)/DXDTM) 
     else
        TIMX = abs(VALX*1.d-3/DXDTM) 
     endif
  endif
  !***********************************************************************
  if(abs(DPDTM) > 0.) then
     TIMP = abs(VALP/DPDTM) 
     TIMT = abs(VALUE/DTDTM) 
  endif
  !***********************************************************************
  PPRO = (PROV+PROVV)/2. 
  if(PROV > 1.d-2 .and. PPRO > 0.) then 
     DEHEL = abs((PROV-PROVV)/(HT1V*PPRO)) 
     COHEL = .10 
     if(XXX(3,1) <= 0.) COHEL = .03 
     if(DEHEL > 0.) TIMHEL = COHEL/DEHEL 
  endif
  if(FG6 > 1.d-4) then 
     PFG6 = (FG6V+FG6)/2. 
     if(PFG6 > 0.) then 
        DECAR = abs((FG6V-FG6)/(HT1V*PFG6)) 
        COCAR = .02 
        if(XXX(4,1) <= 0.) COCAR = .01 
        if(DECAR > 0.) TIMCAR = COCAR/DECAR
     endif
  endif

  if(ECNO > 2.) then 
     ECN = (ECNO+ECNOV)/2. 
     DECNO = abs((ECNO-ECNOV)/(HT1V*ECN)) 
     if(DECNO > 0.) TIMCNO = 0.02/DECNO 
  endif
  if( FG6 >= 1.d-4 ) TIMCNO = TIMCAR 
  HT1 = min(TIM,TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCNO,TIML,TIMMA) 

  ! faccio decadare il divisore a 1 dpo 400 passi,
  ! ma non nel caso che abbia passato un valore negativo
  if(nmd - contatore > 400 .and. div_pastem > 0) div_pastem = 1.0

  HT1 = HT1/abs(div_pastem)

  if(HT1 > max_pastem) HT1 = max_pastem

  if( FG6 < 1.d-4 ) then 
     write(2,110) 
  else 
     write(2,120) 
  endif

  write(2,101)TIM,TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCNO,TIML,TIMMA 
101 format(1X,1P,9D8.1) 
110 format(2X,'MASSA',3X,'RAGGIO',2X,'X-CEN',3X,'PRESS',3X,'TEMP',4X, &
       'ELIO',4X,'CNO',5X,'LUM',5X,'MPUNTO')                            
120 format(2X,'MASSA',3X,'RAGGIO',2X,'X-CEN',3X,'PRESS',3X,'TEMP',4X, &
       'ELIO',4X,'CAR',5X,'LUM',5X,'MPUNTO')                             

  HT2 = HT1V * 1.5 
  if(XXX(3,1) <= 0. .and. PPRO > 2.) HT2 = HT1V*1.05 
  if(HT1 > HT2) HT1 = HT2 
  HT3 = HT1V / 5.0 
  if(XXX(3,1) <= 0. .and. PPRO > 2.) HT3 = HT1V/1.05 
  if(HT1 < HT3) HT1 = HT3 
  return 

end subroutine PASTEM
