subroutine HENYEY(MAXME, MAXMV, ERROR, INUM, IERR) 
  use interfaccia
  use fisica
  use nummod
  use second
  use sistema
  use chimic
  use chim
  use numer
  use costanti
  implicit none
  
  real :: ERROR
  integer :: MAXME, MAXMV, INUM, IERR

  real,dimension(3) :: EPSN, EPNEU, GRAVI, RHO, PMAU, CSPA, GRSAD
  real,dimension(3) :: DELTA, ESP, CAP, DRAD, EPSE, DAD
  real,dimension(2) :: TMPVAR

  integer :: k, kk, j, max, jf, jff, i, jj, n, nn, ip
  real :: rr, zl, pp, tt, em, delp, delt, derl, derr, xh, xhe, tm, pr
  real :: p1, el1, em1, gi, acca, drada, ddelt, ddelp, ddell, dderr
  real :: dept, depp, drot, drop, dgrrt, dgrrp, dgrrl, vagra, grad, ecar
  real :: dgrp, dgrl, dgrr, dgrt, ddrop, ddrot, rr1, dg1, dg3, dg10p, dg10t 

  character(len=15) :: caller="henyey        "

  ERROR  = 0. 
  IERR = 0 
  INUM = 0 
  do K=1,4 
     ALF(K,1) = 0. 
     BET(K,1) = 0.
     GAM(K,1) = 0.
  end do
  BET(3,1) = 1.
  GAM(4,1) = 1.
  do K = 1,MAXME 
     KK = K + 1 
     do J = 1,5 
        DG(J,K) = G(J,KK) - G(J,K) 
     end do
  end do
  ! ************************* LOOP PRINCIPALE ****************************
  MAX = MAXME - 1 
  do JF = 1,MAX 
     JFF = JF + 1 
     ! ****** QUANTITA' FISICHE DEFINITI NEL PUNTO JF + 1/2 *****************
     RR = ( G(1,JF) + G(1,JFF) )/2. 
     ZL = ( G(2,JF) + G(2,JFF) )/2. 
     PP = ( G(3,JF) + G(3,JFF) )/2. 
     TT = ( G(4,JF) + G(4,JFF) )/2. 
     EM = ( G(5,JF) + G(5,JFF) )/2. 
     ! * CHIMICA DEL MESH JF E' QUELLA DELLA REGIONE COMPRESA TRA JF E JF+1 *
     do K = 1,MELE
!        XX(K) = (XXX(K,JF) + XXX(K, JFF)) / 2.
        XX(K) = XXX(K,JF)
     end do
     ! ******* DEFINIZIONE DELLE DERIVATE DI PRESSIONE E TEMPERATURA ********
     DELP = 1.d-5 * PP 
     DELT = 1.d-5 * TT 
     ! ******** TEST SUPERADIBATICO *****************************************
     DERL = 1.d-5 * ZL 
     DERR = 1.d-5 * RR 
     ! *********** FINE TEST SUPERADIABATICO ********************************
     XH  = XX(1) 
     XHE = XX(3) + XX(2)
     ! ************ CALCOLO DERIVATE DELLE VARIE QUANTITA' FISICHE **********
     do I = 1,3 
        if(I == 2) then 
           TT = TT + DELT 
        else if(I == 3) then 
           TT = TT - DELT 
           PP = PP + DELP 
        endif
        EPSN(  I ) = 0. 
        EPNEU( I ) = 0. 
        GRAVI( I ) = 0. 
        DELTA( I ) = 0. 
        TM = TT * 1.d6 
        PR = PP * 1.d17 - arad_3 * TM**4 
      
        call STATE(caller, PR,TM,RHO(I),DAD(I),CSPA(I),PMAU(I))

        call EPSI(RHO(I),TM,JF,0,HT1,EPSN(I),EPSE(I),ecar,DELTA(I))
        
        call NEUTR( RHO(I) , TM , EPNEU(I) ) 
        if ( iread /= 4 .and. (NMD >= 10 .or. IREAD == 2) ) then 
           call EPSIG(PP,TM,RHO(I),CSPA(I),HT1,JF,GRAVI(I),MAXMV,DAD(I),0) 
        endif
        
        call KAPPA(RHO(I), TM, CAP(I), XHE)
        
        
        DRAD(I) = cte_grad * CAP(I) * PP * ZL / (EM * TT**4) 
        !!DRAD(I) = 39.46 * CAP(I) * PP * ZL / (EM * TT**4) 
        ESP(I) = EPSN(I) + EPNEU(I) + GRAVI(I) 

        !********** CALCOLA GRADIENTE SUPERADIABATICO ***********************
        GRSAD(I) = DAD(I) 
        if(DRAD(I) >= DAD(I)) then 
           P1 = PP*100. 
           EL1 = ZL*100.
           EM1 = EM*1000.
           GI = 1.d13 * Ggrav * EM / (RR * RR) 
           !!GI = 6.668d5*EM/(RR*RR) 
           call SUPERA(RHO(I),P1,TT,CAP(I),PMAU(I),CSPA(I),DAD(I),GI,      &
                GRSAD(I),ACCA,DRAD(I))
        endif
        !****************************************************************
     end do
     PP = PP - DELP 
     ! ************ TEST SUPERADIABATICO ********************************
     TMPVAR(1) = 0. 
     TMPVAR(2) = 0.
     if(DRAD(1) >= DAD(1)) then 
        do JJ=1,2 
           if(JJ == 1) then 
              DRADA = DRAD(1)*(1.0+DERL/ZL) 
           else 
              DRADA = DRAD(1) 
              RR = RR + DERR 
           endif
           P1 = PP * 100. 
           EL1 = ZL * 100.
           EM1 = EM * 1000. 
           !!GI = 6.668d5 * EM/(RR*RR)
           GI = 1.d13 * Ggrav * EM / (RR * RR) 
           call SUPERA(RHO(1),P1,TT,CAP(1),PMAU(1),CSPA(1),DAD(1),GI,        &
                TMPVAR(JJ),ACCA,DRADA) 
        end do
        RR = RR-DERR 
     end if
     ! ************* FINE TEST SUPERADIABATICO **************************
     DDELT = 2.0 * DELT 
     DDELP = 2.0 * DELP 
     ! ************* TEST SUPERADIABATICO ********************************
     DDELL = 2.0 * DERL 
     DDERR = 2.0 * DERR 
     ! ************* FINE TEST SUPERADIABATICO ***************************
     DEPT  = (ESP(2) - ESP(1) )/DDELT 
     DEPP  = (ESP(3) - ESP(1) )/DDELP 
     DROT  = (RHO(2) - RHO(1) )/DDELT 
     DROP  = (RHO(3) - RHO(1) )/DDELP 
     DGRRT = (DRAD(2) - DRAD(1) )/DDELT 
     DGRRP = (DRAD(3) - DRAD(1) )/DDELP
     DGRRL = (cte_grad / 2.) * CAP(1) * PP / (EM * TT**4) 
     !!DGRRL = 19.73 * CAP(1) * PP / (EM * TT**4) 
     ! ** SCEGLIE IL GRADIENTE IN BASE AL CRITERIO DI SCHWARTZCHILD *********
     G(6,JF) = DRAD(1) - DAD(1) 
     if(MAIS <= 100) then 
        GG(6,JF) = G(6,JF) 
        VAGRA = G(6,JF) 
     else 
        VAGRA = GG(6,JF) 
     endif
     if( VAGRA < 0. ) then 
        GRAD = DRAD(1) 
        DGRT = DGRRT 
        DGRP = DGRRP 
        DGRL = DGRRL 
        DGRR = 0. 
     else 
        GRAD = GRSAD(1) 
        DGRT = (GRSAD(2) - GRSAD(1))/DDELT 
        DGRP = (GRSAD(3) - GRSAD(1))/DDELP 
        DGRL = 0.
        DGRR = 0.
        ! *************** TEST SUPERADIABATICO ******************************
        DGRL = (TMPVAR(1) - GRSAD(1))/DDELL 
        DGRR = (TMPVAR(2) - GRSAD(1))/DDERR 
        ! *************** FINE TEST SUPERADIABATICO *************************
     endif
     ! ********************** CALCOLA LE DERIVATE ***************************
     DDROP =  DROP / RHO(1) 
     DDROT =  DROT / RHO(1) 
     RR1   =  1. / RR 
     DG1   =  1. / DG(1,JF) 
     DG3   =  1. / DG(3,JF) 
     DG10P = -10. * DG(5,JF) * DEPP 
     DG10T = -10. * DG(5,JF) * DEPT 
     B(1,1) =  G(1,JF) * (RR1 + DG1) 
     C(1,1) =  G(1,JFF) * (RR1 - DG1) 
     B(2,1) =  0. 
     C(2,1) =  0. 
     B(3,1) =  G(3,JF) * (-DG3 - DDROP) 
     C(3,1) =  G(3,JFF) * (DG3 - DDROP) 
     B(4,1) =  G(4,JF) * (-DDROT) 
     C(4,1) =  G(4,JFF) * (-DDROT) 
     B(1,2) =  G(1,JF) * (RR1 - DG1) 
     C(1,2) =  G(1,JFF) * (RR1 + DG1) 
     B(2,2) =  0. 
     C(2,2) =  0. 
     B(3,2) =  G(3,JF) * DDROP 
     C(3,2) =  G(3,JFF) * DDROP 
     B(4,2) =  G(4,JF) * DDROT 
     C(4,2) =  G(4,JFF) * DDROT 
     B(1,3) =  0.
     C(1,3) =  0. 
     B(2,3) = -G(2,JF) 
     C(2,3) =  G(2,JFF) 
     B(3,3) =  G(3,JF) * DG10P 
     C(3,3) =  G(3,JFF) * DG10P 
     B(4,3) =  G(4,JF) * DG10T 
     C(4,3) =  G(4,JFF) * DG10T 
     B(1,4) = -G(1,JF) * DGRR 
     C(1,4) = -G(1,JFF) * DGRR 
     B(2,4) = -G(2,JF) * DGRL 
     C(2,4) = -G(2,JFF) * DGRL 
     B(3,4) = G(3,JF)*(DG(4,JF)*(.5*DG(3,JF)+PP)/(TT*DG(3,JF)**2)-DGRP) 
     C(3,4) = G(3,JFF)*(DG(4,JF)*(.5*DG(3,JF)-PP)/(TT*DG(3,JF)**2)-DGRP) 
     B(4,4) = G(4,JF )*(PP*(-.5*DG(4,JF)-TT)/(DG(3,JF)*TT*TT)-DGRT) 
     C(4,4) = G(4,JFF)*(PP*(-.5*DG(4,JF)+TT)/(DG(3,JF)*TT*TT)-DGRT) 
     ! ******************* VERIFICA EQUAZIONI *******************************
     !!E(1) = log(abs(DG(3,JF)*RR*RR/(DG(1,JF)*6.668d-2*EM*RHO(1) ))) 
     !!E(2) = log(abs(12.56637061d-3*RR*RR*RHO(1)*DG(1,JF)/DG(5,JF)))
     E(1) = log(abs(DG(3,JF)*RR*RR/(1.d6*DG(1,JF)*Ggrav*EM*RHO(1) ))) 
     E(2) = log(abs(4.d-3*pigre*RR*RR*RHO(1)*DG(1,JF)/DG(5,JF)))
     E(3) = DG(2,JF) - ESP(1) * 10.0 * DG( 5 , JF ) 
     E(4) = DG(4,JF) * PP / ( DG(3,JF) * TT ) - GRAD 
     do N = 1,4 
        ALF1(N) = -E(N) - B(1,N) * ALF(1,JF) - B(2,N) * ALF(2,JF) 
        BET1(N) = -B(1,N) * BET(1,JF) - B(2,N) * BET(2,JF) - B(3,N) 
        GAM1(N) = -B(1,N) * GAM(1,JF) - B(2,N) * GAM(2,JF) - B(4,N) 
     end do
     E(3) = E(3)/sqrt(DG(2,JF)**2+100.0*ESP(1)*ESP(1)*DG(5,JF)**2) 
     E(4) = E(4)/sqrt(PP*PP*DG(4,JF)**2/(TT*TT*DG(3,JF)**2)+GRAD*GRAD) 
 
     ! *********** CALCOLA QUANTO SONO VERIFICATE LE EQUAZIONI **************
     do N = 1,4 
        if( abs( E(N) ) <= abs( ERROR ) ) cycle
        ERROR  = E(N) 
        IERR = N 
        INUM = JF 
     end do
     ! ************** IN CASO DI DEBUG STAMPA VARIE QUANTITA' ***************
     if( IPRALL /= 0 ) then 
        write(2,110) JF 
        write(2,90) RR,ZL,PP,TT,EM,DEPP,DEPT,DROP,DROT,DGRL,DGRP,DGRT 
        write(2,90) (EPSN(NN),NN=1,3),(EPNEU(NN),NN=1,3),(GRAVI(NN),NN=1,3) 
        write(2,90)(RHO(NN),NN=1,3),(CAP(NN),NN=1,3),(DRAD(NN),NN=1,3),   &
             (ESP(NN),NN=1,3) 
        write(2,90) G(1,JF),G(1,JFF),DG(1,JF),G(2,JF),G(2,JFF),DG(2,JF),   &
             G(3,JF),G(3,JFF),DG(3,JF),G(4,JF),G(4,JFF),DG(4,JF)
        write(2,90) (XX(IP),IP=1,6),G(5,JF),G(5,JFF),DG(5,JF) 
        write(2,100) C(1,1),C(3,1),C(4,1),B(1,1),B(3,1),B(4,1),E(1) 
        write(2,100) C(1,2),C(3,2),C(4,2),B(1,2),B(3,2),B(4,2),E(2) 
        write(2,100) C(2,3),C(3,3),C(4,3),B(2,3),B(3,3),B(4,3),E(3) 
        write(2,100) C(2,4),C(3,4),C(4,4),B(2,4),B(3,4),B(4,4),E(4) 
     endif
     call RESNUC(JFF) 
 end do
 G(6,MAXME) = G(6,MAXME - 1) 
  ! ********************** FINE LOOP PRINCIPALE **************************
 return 

90 format(1X,1P,12E10.3) 
100 format(1X,1P,E10.3,'DX',E10.3,'DP',E10.3,'DT=', E10.3,'DXV',      &
         E10.3,'DPV',E10.3,'DTV',E10.3)         
110 format(/,1X,'MESH N.',I4) 
end subroutine HENYEY
