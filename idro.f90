! La routine trova le nuove abbondanze degli elementi (1,2,...,6)
! dopo una evoluzione di durata temporale PAS.
! In input riceve:
! RO: densita'
! T: usato solo per stampe di debug
! PAS: il passo temporale
! JF: il numero di mesh su cui lavorare
! INI,IFI: numero di elemento iniziale e finale su cui lavorare
! N: numero totale di elementi da processare
! YV: abbondanze iniziali degli elementi
! YF: abbondanze nuove da determinare
! NEQUI(): vettore che indica se gli elementi sono o no all'equilibrio

! La routine calcola il valore delle nuove abbondanze linearizzando la 
! dipendenza da T:
! YF = YV + (dY/dT) * PAS
! dove (dY/dT) e' la derivata temporale delle abbondanze. Tale valore e' 
! calcolato dalla routine nel vettore DX(). 
! In particolare si sceglie di valutare il valore di DX a meta' 
! dell'intervallo temporale, o meglio nel punto Y posto a meta' 
! dell'intervallo di abbondanza:  Y = 1/2 (YF+YV).
! L'algoritmo cerca quindi lo zero del sistema DB() = -(YF-YV)/PAS + DX
! mediante metodo di Raphson-Newton.
! Nel caso in cui l'elemento 1 (He3) sia all'equilibrio (NEQUI(1)=0) 
! l'equazione relativa viene modificata ricercando lo zero di
! DB(1) = DX(1) valutata in YF(1), ossia un punto stazionario per YF(1).

subroutine IDRO (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S, oishel) 
  use interfaccia
  use fisica
  use chequ
  use equil
  use overshoot
  use nummod
  use chimic
  use varie

  implicit none

  real :: RO, T, PAS
  integer :: JF, INI, IFI, N, oishel
  real,dimension(MELE) :: YV, YF 
  integer,dimension(MELE) :: NEQUI
  real,dimension(nsezi) :: S

  real,parameter :: CR = 0.5
  integer :: ni, ih, kk, j, ib, ig, lk, l, ks, ir, jl, m
  real :: scno, the3, beta, gamma, ef, dd, yyo, yyn, yyc, rpox, rpnz
  real :: rapp1, rapp2, rpca

  real,dimension(MELE,MELE) :: DXX
  real,dimension(MELE) :: Y, DX, DB, B
  real,dimension(MELE*MELE) :: A 
  logical :: linearCONV = .true.

  SCNO = YV(4)+YV(5)+YV(6) 

  ! Efficienza relativa del ramo CN rispetto a CN+NO
  if( S(7)+S(32) > 0 ) then
     RAPP1 = S(7)/(S(7)+S(32))
     RAPP2 = 1.0 -  RAPP1
  else
     RAPP1 = 0.0
     RAPP2 = 0.0
  endif
  
  !********************* TEMPO SCALA EQUILIBRIO HE3 **********************
  if( (G(6,JF) < 0. .and. KOVER == 0) .or. (KOVER == 1 .and. JF > L3 .and. &
       L3 > 0) .or. oISHEL == 1 ) then 
     THE3 = 1./((YV(1)*S(30)+YV(3)*S(31))*RO)
     if(THE3 < PAS) NEQUI(1) = 0 
     ISHEL = 1 
  endif
  
  NI = 0

  ! azzero la matrice jacobiana del sistema da risolvere
  do KK=ini,ifi 
     do IH=ini,ifi
        DXX(IH,KK) = 0.
     end do
  end do

  ! inizio le iterazioni (max 20)
  do NI=0,20
     ! calcolo punto centrale dell'intervallo di abbondanze
     do IH=INI,IFI 
        Y(IH) = .5*(YF(IH)+YV(IH)) 
     end do
     ! ******* CALCOLO DERIVATE DY(K)/DT IN T+PAS/2 SE RADIATIVO   **********
     ! ******* - - - - - - - - - - - -   IN    T    SE CONVETTIVO  **********
     if(NEQUI(1) > 0) then 
        DX(1) = RO*(Y(2)**2*S(29)/2.-Y(1)**2*S(30)-Y(3)*Y(1)*S(31)) 
     else 
        DX(1) = RO*(YF(2)**2*S(29)/2.-YF(1)**2*S(30)-YF(3)*YF(1)*S(31)) 
     endif
     DX(2) = RO*(Y(1)**2*S(30)-3.*Y(2)**2*S(29)/2.-Y(3)*Y(1)*S(31)       &
          -2.*Y(2)*(Y(4)*S(1)+Y(5)*S(5)+Y(6)*S(9)))
     DX(3) = RO*(Y(1)**2*S(30)/2.+Y(3)*Y(1)*S(31)+RAPP1*                 &
          Y(5)*Y(2)*S(5)+Y(6)*Y(2)*S(9))
     if(KECUIL(JF) == 0) then 
        DX(4) = RO*(Y(2)*(-Y(4)*S(1)+RAPP1*Y(5)*S(5))) 
        DX(5) = RO*(Y(2)*(-Y(5)*S(5)+Y(6)*S(9)+Y(4)*S(1))) 
        DX(6) = RO*(Y(2)*(-Y(6)*S(9)+RAPP2*Y(5)*S(5))) 
     endif

     ! SE CONVETTIVO INTEGRA LINEARMENTE se linearCONV = .true.
     ! CONTROLLA EQUILIBRI INTEGRATI DI HE3,C,N,O NELLA EQLB    
     if( ISHEL == 0 .and. NI == 0) then 
        COEFF(1) = COEFF(1)+RO*S(30)*(G(5,JF+1)-G(5,JF)) 
        COEFF(2) = COEFF(2)+RO*S(31)*(G(5,JF+1)-G(5,JF)) 
        COEFF(3) = COEFF(3)+RO*S(29)*(G(5,JF+1)-G(5,JF)) 
        COEFF(4) = COEFF(4)+RO*S(1)*(G(5,JF+1)-G(5,JF)) 
        COEFF(5) = COEFF(5)+RO*S(5)*RAPP1*(G(5,JF+1)-G(5,JF)) 

        if(linearCONV .or. xxx(1,1) > 0.99*XH ) then
           do KK=INI,IFI 
              YF(KK) = Y(KK)+DX(KK)*PAS 
           end do
           
           return 
        endif
     endif

     !************************* RAHPSON - NEWTON *************************** 
     !************* CALCOLO TERMINI NOTI E CONTROLLO CONVERGENZA *********** 

     BETA = 0.
     do J=INI,IFI 
        ! DB() e' il sistema di cui trovare lo 0 
        ! YF() sono le variabili indipendenti
        DB(J) = -(YF(J)-YV(J))/PAS*NEQUI(J)+DX(J) 
        if(.not. (YV(J) < 1.d-9 .or. NI == 0 .or. NEQUI(J) == 0)) then
           EF = abs(1. -(YF(J)-YV(J))/(PAS*DX(J))) 
           if( EF > BETA ) then 
              IB = J 
              BETA = EF 
           endif
        endif
     end do
     
     if( NI /= 0 ) then
        if( NI > 5 ) then 
           write(2,444) IB,BETA,IG,GAMMA,NI,JF 
        endif
        ! verifico se ho convergenza, nel caso esco
        if( BETA < 1.d-4 .or. GAMMA < 1.d-4) then
           if( nequi(1) > 0) then
              if(Y(1) >= 0.) return
              Y(1) = 0.
           else
              if(YF(1) >= 0.) return
              YF(1) = 0.
           endif
        endif
     endif
     ! ************************  DERIVATE -DB(J)/DYF  ***********************
     if(NEQUI(1) > 0) then 
        DXX(1,1) = 1./PAS+RO*(2.*Y(1)*S(30)+Y(3)*S(31))*CR 
        DXX(1,2) = -RO*Y(2)*S(29)*CR 
        DXX(1,3) = RO*Y(1)*S(31)*CR 
     else 
        DXX(1,1) = RO*(2.*YF(1)*S(30)+YF(3)*S(31)) 
        DXX(1,2) = -RO*YF(2)*S(29) 
        DXX(1,3) = RO*YF(1)*S(31) 
     endif
     DXX(2,1) = RO*(-2.*Y(1)*S(30)+Y(3)*S(31))*CR 
     DXX(2,2) = 1./PAS+RO*(3.*Y(2)*S(29)+2.*(Y(4)*S(1)+Y(5)*S(5)+Y(6)*   &
          S(9)))*CR
     DXX(2,3) = RO*Y(1)*S(31)*CR 
     DXX(3,1) = RO*(-Y(1)*S(30)-Y(3)*S(31))*CR 
     DXX(3,2) = RO*(-RAPP1*Y(5)*S(5)-Y(6)*S(9))*CR 
     DXX(3,3) = 1./PAS-RO*Y(1)*S(31)*CR 
     if(KECUIL(JF) == 0) then 
        DXX(2,4) = RO*2.*Y(2)*S(1)*CR 
        DXX(2,5) = RO*2.*Y(2)*S(5)*CR 
        DXX(2,6) = RO*2.*Y(2)*S(9)*CR 
        DXX(3,5) = -RO*RAPP1*Y(2)*S(5)*CR 
        DXX(3,6) = -RO*Y(2)*S(9)*CR 
        DXX(4,2) = RO*(-RAPP1*Y(5)*S(5)+Y(4)*S(1))*CR 
        DXX(4,4) = 1./PAS+RO*Y(2)*S(1)*CR 
        DXX(4,5) = -RO*RAPP1*Y(2)*S(5)*CR 
        DXX(5,2) = RO*(-Y(4)*S(1)+Y(5)*S(5)-Y(6)*S(6))*CR 
        DXX(5,4) = -RO*Y(2)*S(1)*CR 
        DXX(5,5) = 1./PAS+RO*Y(2)*S(5)*CR 
        DXX(5,6) = -RO*Y(2)*S(9)*CR 
        DXX(6,2) = RO*(-RAPP2*Y(5)*S(5)+Y(6)*S(9))*CR 
        DXX(6,5) = -RO*RAPP2*Y(2)*S(5)*CR 
        DXX(6,6) = 1./PAS+RO*Y(2)*S(9)*CR 
     endif
     ! *************  CALCOLO DELLE DIFFERENZE FINITE  ********************* 
     LK = 0 
     do J=INI,IFI 
        KK = 0 
        LK = LK+1 
        B(LK) = DB(J) 
        do L=INI,IFI 
           KK = KK+1 
           JL = (J-INI)*N+KK 
           A(JL) = DXX(L,J) 
        end do
     end do
     call SIMQ(A,B,N,KS) 
     if( KS == 1 ) then
        write(*,130) 
        write(2,130) 
        write(66,*) '81 - idro'
        write(66,130) 
        stop 
     endif
     GAMMA = 0.
     do J=INI,IFI 
        KK = J-INI+1
        YF(J) = YF(J)+B(KK)
        !*********** CALCOLA VARIAZIONI PERCENTUALI ELEMENTI IMPORTANTI ***
        if( YV(J) > 1.d-9 ) then 
           DD = abs(B(KK)/YF(J)) 
           if( DD > GAMMA ) then 
              IG = J 
              GAMMA = DD 
           endif
        endif
     end do

     !******************* CONTROLLO EQUILIBRI CNO ************************
     YYO = SCNO*RAPP2/(S(9)/S(5)+S(9)/S(1)*RAPP1+RAPP2) 
     RPOX = (YF(6)-YV(6))/(YYO-YV(6)) 
     ! matt 
     ! controllo di non avere abbondanze negative.
     ! Altrimenti evito il ciclo di controllo equilibri
     if(RPOX > 1.) then 
        KECUIL(JF) = 1 
        IFI = 3 
        N = 3 
        YF(6) = YYO 
        YF(5) = (SCNO-YF(6))/(1.+S(5)/(S(1)*RAPP1)) 
        YF(4) = (SCNO-YF(6))/(1.+S(1)/(S(5)*RAPP1)) 
        Y(6) = YF(6) 
        Y(5) = YF(5) 
        Y(4) = YF(4) 
     else if(YF(6) < SCNO) then
        YYN = (SCNO-YF(6))/(1.+S(5)/(S(1)*RAPP1)) 
        YYC = (SCNO-YF(6))/(1.+S(1)/(S(5)*RAPP1)) 
        RPNZ = (YF(5)-YV(5))/(YYN-YV(5)) 
        RPCA = (YF(4)-YV(4))/(YYC-YV(4)) 
        if(RPNZ > 1 .or. RPCA > 1.) then 
           KECUIL(JF) = 2 
           IR = 5 
           IFI = 3 
           N = 3 
           YF(5) = YYN 
           YF(4) = YYC 
           Y(6) = YF(6) 
           Y(5) = YF(5) 
           Y(4) = YF(4) 
        endif
     endif
  end do

  ! se sono qui vuol dire che sono uscito dal ciclo per troppe iterazioni
  write(*,250) 
  write(2,250) 
  write(2,113)T,RO,JF 
  write(66,*)'80 - idro'
  write(66,250) 
  write(66,113)T,RO,JF 
  stop 

113 format(1X,'T =',1P,E10.3,2X,'RO =',1P,E10.3,2X,'MESH =',I4) 
130 format(1X,'ATTENZIONE QUALCOSA NON VA NELLA EPSI. IL DET = 0 NEL  &
       &RAHPSON/NEWTON - vedi idro')
250 format(1X,'MORTO: LA EPSI ITERA TROPPO NEL RAPHSON/NEWTON-        &
       &vedi idro')
444 format(1X,'ITERAZIONI EPSI: DF =',I3,E10.3,' DX =',I3,E10.3,      &
       ' N. ITER.=',I3,' N. MESH=',I4) 

end subroutine IDRO
