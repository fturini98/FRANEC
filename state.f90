! Questa routine interpola le tabelle dell'equazione di stato
! In input riceve:
! PR: la pressione per cui valutare le quantita'
! TE: la temperatura per cui valutare le quantita'
! In output ritorna:
! RHO, GRAD, CSPE, PMOL

! L'interpolazione in T e P sulle tabelle e' effettuata linearmente.
! La routine implementa un metodo per calcolare solo le quantita' necessarie:
! All'ingresso controlla se le quantita' RHO, GRAD, CSPE, PMOL sono uguali 
! a 1.d99. In tal caso disabilita il calcolo della quantita' relativa.

! Le tabelle implementate sono:
! EOS NUOVA OPAL 2006          
! Straniero

! la routine prevede una "banda di raccordo" tra le due tabelle in relazione
! al valore di C nel mesh. Al di sotto di Clow usa OPAL, sopra 0.03 usa 
! Straniero, all'interno della banda interpola tra i valori 

subroutine STATE(caller, PR,TE,RHO,GRAD,CSPE,PMOL) 
  use interfaccia
  use chim
  use nummod
  use fasievolut

  implicit none

  character(len=15) :: caller
  real :: PR, TE, RHO, GRAD, CSPE, PMOL

  character(len=1) :: TESTO(72) 
  integer,parameter :: mx=5, nt=197, np=325
  !***********dimension per livermore 2006*****************               
  real,save,dimension(mx,nt,np) :: RO, CP, AD, GAMMA1, CHIRO, CHIT, amu_M
  real,save,dimension(nt) :: T6LA
  real,dimension(2,2) :: DE, DA, CS, XGAM, CHR, CHT, AMUM 
  real,dimension(2) :: DEN, DAD, CSP, PMOLM, XXGAM, CHIR, CHITI
  real,parameter,dimension(5) :: XXL =(/0.8,0.6,0.4,0.2,0.0/)
  real,parameter,dimension(6) :: XXL1 =(/1.0,0.8,0.7,0.6,0.3,0.0/) 
  !*************dimension per OSCAR***********************                
  real,save,dimension(9,58,33) :: RO1, AD1, CP1, PM
  real,dimension(2,5) :: PU, DE1, DA1, CS1
  real,dimension(5) :: DAD1, CSP1, PMU, DEN1
  real,dimension(58) :: TT, PP 

  integer :: MT1=0, MT2=0, MP1=0, MP2=0 
  real :: RATT = 0., RATP = 0., RATCC = 0.
  ! IES=1 se estrapolo in temperatura
  integer,save :: IS = 0, ISS = 0
  real,parameter :: UM = 2.302585
  
  ! pressioni iniziali per OPAL, lette da file
  real,dimension(nt) :: pp1

  integer :: ntfree, k, j, n, l, ll, kinf, ksup, jj, nci, ncf, ig
  real :: tt6, xt, tso, pnew, for1, for2, uf1, uf2, uf3, uf4
  real :: diff, t, p, ter
  real :: XNGAM, CAIRO, CAIT, cspetmp
  
  integer,parameter :: AGGIUNTE = 0
  integer :: calcolarho, calcolacspe, calcolagrad, calcolapmol
  integer, save :: SvK = 1
  integer,dimension(2) :: ind
  character(len=76) :: int1, int2
  integer :: raccordo
  real :: rhos, pmols, grads, cspes
  ! valore di xx(4) sotto cui non chiamare Oscar 
  ! e sopra cui chiamarlo sempre
  ! Cup = 0.03 perche' durante il flash produco 0.03 di carbonio  
  real,parameter :: Clow = 0.007, Cup = 0.03

  ! controllo quali sono le quantita' da calcolare
  ! un valore in input pari a 1.d99 disabilita il calcolo della 
  ! quantita' relativa
  calcolarho = 1 
  calcolacspe = 1 
  calcolagrad = 1 
  calcolapmol = 1
  if(RHO == 1.d99) calcolarho = 0
  if(GRAD == 1.d99) calcolagrad = 0 
  if(CSPE == 1.d99) calcolacspe = 0 
  if(PMOL == 1.d99) calcolapmol = 0

  ! controllo valori in input
  if (PR <= 0. .or. TE <= 0) then 
     write(*,*)
     write(*,321) PR, TE 
     write(*,*) "caller = ", caller
     write(66,*)'150 - state'
     write(66,321) PR, TE 
     write(66,*) "caller = ", caller
     stop 
  endif

  ! controllo T per Livermore
  if(TE <= 1870.0 .or. TE >= 200.0d6 ) goto 90 

  !********************************************                           
  !  se ho prodotto carbonio vado da Oscar                                
  raccordo = 0
  if(XX(4) > Cup) then
     ! vado direttamente nelle tabelle di Oscar
     goto 90 
  else if(XX(4) > Clow .and. xx(1) < 1d-35 .and. idfase <= fase_innHe) then
     ! imposto il calcolo a banda di raccordo
     raccordo = 1
  endif

  ! DEFINISCO PARAMETRI INIZIALI
  TT6 = TE*1.d-6 
  XT = log10(TT6) 

  ! SE HO 1870 < T < 2*10e8 LEGGO LE LIVERMORE                        
  ! lettura tabella nuova OPAL 
  if(ISS == 0) then 
     write(*,*) 'Lettura Livermore'
     write(67,*) "================================"
     write(67,*) 'Lettura OPAL + FreeEOS. Modello: ', nmd
     write(67,*) "================================"
     write(67,*)
     ISS = 1 

     read(5,'(a76)') int1
     read(5,'(a76)') int2

     ! scrittura referenze su run.log
     write(67,*) "@ ================================"
     write(67,*) "@         EOS + FreeEOS           "
     write(67,*) "@ ================================"
     write(67,'(" @ ",a76)') int1
     write(67,'(" @ ",a76)') int2
     write(67,*) "@ ================================"
     write(67,*)

     ! lettura grigliato pressioni inziali
     read(5,*)
     read(5,*) (PP1(ig),ig=1,nt)

     do K=1,mx 
        do J=1,nt 
           read(5,119)(RO(K,J,N),N=1,np),(CP(K,J,N),N=1,np),             &
                (AD(K,J,N),N=1,np),(GAMMA1(K,J,N),N=1,np), &
                (CHIRO(K,J,N),N=1,np),(CHIT(K,J,N),N=1,np), &
                (amu_M(K,J,N),N=1,np),T6LA(J)
        end do
     end do
  endif
  ! ******************************************************************    

  if(XX(1) == 0.) XX(1) = 1.d-40 
  if(XX(3) == 0.) XX(3) = 1.d-40 
  if(XX(4) == 0.) XX(4) = 1.d-40 
  if(XX(5) == 0.) XX(5) = 1.d-40 
  if(XX(6) == 0.) XX(6) = 1.d-40 
  
  !************SCELGO X T ED P  ***********************                   
  ! Nota: le tabelle sono in Pgas: Pnew=Pgas                              
  !****************************************************                   
  Pnew = PR 
  Pnew = log10(Pnew) 
  !    X(K-1) e' maggiore di X(k)                                         
  do K=2,5 
     if( XXL(K) <= XX(1) ) then 
        L = K 
        LL = K-1 
        RATCC = (XX(1)-XXL(K))/(XXL(K-1)-XXL(K)) 
        exit
     endif
  end do

  !      T6LA(MT2=K-1) e' minore di T6LA(MT1=K)                           
  if(svk < 2) svk = 2
  if(.not. (T6LA(SvK-1) < XT .and. T6LA(SvK) > XT)) then
     Kinf = 2 
     Ksup = nt 
     ! attenzione a XT NaN
     if(T6LA(2) > XT) then 
        K = 2 
     else if(T6LA(nt) < XT) then 
        ! aggiunta gestione bordo superiore
        K = nt 
     else
        do while (Ksup-Kinf > 1)
           K = (Kinf+Ksup)/2
           if(T6LA(K) >= XT) then 
              Ksup = K 
           else if(T6LA(K) < XT) then 
              Kinf = K 
           end if
        end do
        K = Ksup
     endif
  else
     k = Svk
  endif
  SvK = K

  MT1 = K 
  MT2 = MT1-1 
  RATT = (XT-T6LA(K-1))/(T6LA(K)-T6LA(K-1)) 

  FOR1 = Pnew-PP1(MT1) 
  FOR2 = Pnew-PP1(MT2) 
  if(FOR1 < 0. .or. FOR2 < 0.) goto 50 
  if(FOR1 > 16. .or. FOR2 > 16.) goto 50 
  MP1 = int(FOR1*20.0)+1 
  MP2 = int(FOR2*20.0)+1 
  RATP = (Pnew-(PP1(MT1)+real(MP1-1)/20.0))*20.0 

  ! ATTENZIONE!!!!!!  Qui faccio il controllo se ho valori nella          
  ! tabella di Livermore altrimenti vado alle tabelle di Oscar            
  do JJ=1,2 
     J = L-JJ+1

     UF1 = abs(RO(J,MT2,MP2)+99.9) 
     if(UF1 < 0.01) goto 90 

     UF2 = abs(RO(J,MT2,MP2+1)+99.9) 
     if(UF2 < 0.01) goto 90 

     UF3 = abs(RO(J,MT1,MP1)+99.9) 
     if(UF3 < 0.01) goto 90 

     UF4 = abs(RO(J,MT1,MP1+1)+99.9) 
     if(UF4 < 0.01) goto 90 

     ! ATTENZIONE!!!!!!  Faccio anche il controllo che non ci siano          
     ! Cv negativi altrimenti vado alle tabelle di Oscar                     
     UF1 = abs(CP(J,MT2,MP2)+99.9) 
     if(UF1 < 0.01) goto 90 

     UF2 = abs(CP(J,MT2,MP2+1)+99.9) 
     if(UF2 < 0.01) goto 90 

     UF3 = abs(CP(J,MT1,MP1)+99.9) 
     if(UF3 < 0.01) goto 90 

     UF4 = abs(CP(J,MT1,MP1+1)+99.9) 
     if(UF4 < 0.01) goto 90 
  end do

  !************************** CALCOLO DENSITA' ************************** 
  if(calcolarho == 1) then
     do JJ=1,2 
        J = L-JJ+1
        DE(1,JJ) = RO(J,MT2,MP2)+(RO(J,MT2,MP2+1)-RO(J,MT2,MP2))*RATP 
        DE(2,JJ) = RO(J,MT1,MP1)+(RO(J,MT1,MP1+1)-RO(J,MT1,MP1))*RATP 
        !************estrapolazione in T******************************       
        DEN(JJ) = exp(UM*(DE(1,JJ)+(DE(2,JJ)-DE(1,JJ))*RATT)) 
     end do
     RHO = DEN(1)+RATCC*(DEN(2)-DEN(1))
  endif

  !********************* CALCOLO CALORE SPECIFICO ************************
  if(calcolacspe == 1) then
     do JJ=1,2 
        J = L-JJ+1
        CS(1,JJ) = CP(J,MT2,MP2)+(CP(J,MT2,MP2+1)-CP(J,MT2,MP2))*RATP 
        CS(2,JJ) = CP(J,MT1,MP1)+(CP(J,MT1,MP1+1)-CP(J,MT1,MP1))*RATP 
        !************estrapolazione in T******************************        
        CSP(JJ) = exp(UM*(CS(1,JJ)+(CS(2,JJ)-CS(1,JJ))*RATT)) 
     end do
     CSPE = CSP(1)+RATCC*(CSP(2)-CSP(1)) 
  endif

  !********************* CALCOLO GRADIENTE ADIAB *************************
  if(calcolagrad == 1) then
     do JJ=1,2 
        J = L-JJ+1
        DA(1,JJ) = AD(J,MT2,MP2)+(AD(J,MT2,MP2+1)-AD(J,MT2,MP2))*RATP 
        DA(2,JJ) = AD(J,MT1,MP1)+(AD(J,MT1,MP1+1)-AD(J,MT1,MP1))*RATP 
        !************estrapolazione in T******************************        
        DAD(JJ) = DA(1,JJ)+(DA(2,JJ)-DA(1,JJ))*RATT 
     end do
     GRAD = DAD(1)+RATCC*(DAD(2)-DAD(1)) 
  endif

  ! le variabili nel blocco AGGIUNTE non sono usate
  ! posso evitarne il calcolo
  if(AGGIUNTE == 1) then
     !********************** CALCOLO GAMMA1*******************************
     do JJ=1,2 
        J = L-JJ+1
        XGAM(1,JJ) = GAMMA1(J,MT2,MP2)+(GAMMA1(J,MT2,MP2+1)             &
             -GAMMA1(J,MT2,MP2))*RATP
        XGAM(2,JJ) = GAMMA1(J,MT1,MP1)+(GAMMA1(J,MT1,MP1+1)             &
             -GAMMA1(J,MT1,MP1))*RATP  
        !************estrapolazione in T*****************************       
        XXGAM(JJ) = XGAM(1,JJ)+(XGAM(2,JJ)-XGAM(1,JJ))*RATT 
     end do
     XNGAM = XXGAM(1)+RATCC*(XXGAM(2)-XXGAM(1)) 

     !********************* CALCOLO CHIRO *************************          
     do JJ=1,2 
        J = L-JJ+1
        CHR(1,JJ) = CHIRO(J,MT2,MP2)+(CHIRO(J,MT2,MP2+1)-              &
             CHIRO(J,MT2,MP2))*RATP   
        CHR(2,JJ) = CHIRO(J,MT1,MP1)+(CHIRO(J,MT1,MP1+1)-              &
             CHIRO(J,MT1,MP1))*RATP                                        
        !************estrapolazione in T*****************************        
        CHIR(JJ) = CHR(1,JJ)+(CHR(2,JJ)-CHR(1,JJ))*RATT 
     end do
     CAIRO = CHIR(1)+RATCC*(CHIR(2)-CHIR(1)) 

     !********************* CALCOLO CHIT *************************           
     do JJ=1,2 
        J = L-JJ+1
        CHT(1,JJ) = CHIT(J,MT2,MP2)+(CHIT(J,MT2,MP2+1)-                &
             CHIT(J,MT2,MP2))*RATP 
        CHT(2,JJ) = CHIT(J,MT1,MP1)+(CHIT(J,MT1,MP1+1)-                &
             CHIT(J,MT1,MP1))*RATP                          
        !************estrapolazione in T******************************       
        CHITI(JJ) = CHT(1,JJ)+(CHT(2,JJ)-CHT(1,JJ))*RATT 
     end do
     CAIT = CHITI(1)+RATCC*(CHITI(2)-CHITI(1)) 
  endif

  !************** calcolo Peso molecolare medio ***********************
  if(calcolapmol == 1) then
     do JJ=1,2 
        J = L-JJ+1
        AMUM(1,JJ) = amu_M(J,MT2,MP2)+(amu_M(J,MT2,MP2+1)-              &
             amu_M(J,MT2,MP2))*RATP  
        AMUM(2,JJ) = amu_M(J,MT1,MP1)+(amu_M(J,MT1,MP1+1)-              &
             amu_M(J,MT1,MP1))*RATP        
        !************estrapolazione in T******************************       
        PMOLM(JJ) = AMUM(1,JJ)+(AMUM(2,JJ)-AMUM(1,JJ))*RATT 
     end do
     PMOL = PMOLM(1)+RATCC*(PMOLM(2)-PMOLM(1)) 
  endif

  !***********************************************************************
  if(XX(1) < 1.d-39) XX(1) = 0. 
  if(XX(3) < 1.d-39) XX(3) = 0. 
  if(XX(4) < 1.d-39) XX(4) = 0. 
  if(XX(5) < 1.d-39) XX(5) = 0. 
  if(XX(6) < 1.d-39) XX(6) = 0. 
  
  if(raccordo == 1) then
     ! sono nella banda di raccordo, chiamo anche Oscar e poi interpolo
     rhos = RHO
     grads = grad
     cspes = cspe
     pmols = pmol
     goto 90
  endif
  
  return 

  !*********************************************************************  
  !                        COMINCIA OSCAR                                 
  !********************************************************************   
90 continue 
  ! LETTURA TABELLE DI EOS Straniero
  if(IS == 0) then 
     IS = 1 
     read(8,100) TESTO 
     write(2,101) TESTO 
     write(*,101) TESTO 
     write(67,*) "================================"
     write(67,*), "Lettura EOS Straniero. Modello :", nmd
     write(67,101) TESTO 
     write(67,*) "================================"
     do K=1,9 
        do J=1,58 
           if(K == 3 .and. J < 40) cycle
           if(K == 4 .and. J < 40) cycle
           if(K == 5 .and. J < 40) cycle
           if(K >= 6 .and. J > 40) cycle
           read(8,102)(RO1(K,J,N),N=1,33),(AD1(K,J,N),N=1,33),               &
                (CP1(K,J,N),N=1,33),(PM(K,J,N),N=1,33),TT(J),PP(J) 
        end do
     end do
  endif
  ! IN ATTESA DELLE TABELLE METALLARE                                  
  DIFF = 1.-(XX(1)+XX(2)+XX(3)+XX(4)+XX(5)+XX(6)) 
  XX(6) = XX(6)+DIFF 

  if(XX(1) == 0.) XX(1) = 1.d-40 
  if(XX(3) == 0.) XX(3) = 1.d-40 
  if(XX(4) == 0.) XX(4) = 1.d-40 
  if(XX(5) == 0.) XX(5) = 1.d-40 
  if(XX(6) == 0.) XX(6) = 1.d-40 

  P = log10(PR) 
  T = log10(TE) 

  if(T < 3.3 .or. T > 9.6) goto 50 
  if( T < 6. ) then 
     do K=2,6 
        if( XXL1(K) < XX(1) ) then 
           RATCC = (XX(1)-XXL1(K))/(XXL1(K-1)-XXL1(K)) 
           L = K+4 
           LL = L-1 
           if( L == 6 ) LL = 1 
           if( K == 6 ) then 
              L = 2 
              LL = 9 
           endif
           exit
        endif
     end do
     if(T < 5) then 
        TER = T-3.3 
        MT1 = int(TER*20)+1 
     else 
        TER = T-5.0 
        MT1 = int(TER*5)+35 
     endif
     MT2 = MT1+1 
     FOR1 = P-PP(MT1) 
     FOR2 = P-PP(MT2) 
     if(FOR1 < 0. .or. FOR2 < 0.) goto 50 
     if(FOR1 > 16. .or. FOR2 > 16.) goto 50 
     MP1 = int(FOR1*2.)+1 
     MP2 = int(FOR2*2.)+1 
     RATT = (T-TT(MT1))/(TT(MT2)-TT(MT1)) 
     RATP = (P-(PP(MT1)+real(MP1-1)/2.))*2. 

     ind(1) = l
     ind(2) = ll

     !************************** CALCOLO DENSITA' ************************** 
     if(calcolarho == 1) then
        do JJ=1,2 
           J = ind(JJ)
           DE1(1,JJ) = RO1(J,MT2,MP2)+(RO1(J,MT2,MP2+1)-RO1(J,MT2,MP2))*RATP 
           DE1(2,JJ) = RO1(J,MT1,MP1)+(RO1(J,MT1,MP1+1)-RO1(J,MT1,MP1))*RATP 
           DEN1(JJ) = exp(UM*(DE1(2,JJ)+(DE1(1,JJ)-DE1(2,JJ))*RATT)) 
        end do
        RHO = DEN1(1)+RATCC*(DEN1(2)-DEN1(1)) 
     endif

     !********************* CALCOLO CALORE SPECIFICO ************************
     if(calcolacspe == 1) then
        do JJ=1,2 
           J = ind(JJ)
           CS1(1,JJ) = CP1(J,MT2,MP2)+(CP1(J,MT2,MP2+1)-CP1(J,MT2,MP2))*RATP 
           CS1(2,JJ) = CP1(J,MT1,MP1)+(CP1(J,MT1,MP1+1)-CP1(J,MT1,MP1))*RATP 
           CSP1(JJ) = exp(UM*(CS1(2,JJ)+(CS1(1,JJ)-CS1(2,JJ))*RATT)) 
        end do
        CSPE = CSP1(1)+RATCC*(CSP1(2)-CSP1(1)) 
     endif

     !********************* CALCOLO GRADIENTE ADIAB *************************
     if(calcolagrad == 1) then
        do JJ=1,2 
           J = ind(JJ)
           DA1(1,JJ) = AD1(J,MT2,MP2)+(AD1(J,MT2,MP2+1)-AD1(J,MT2,MP2))*RATP 
           DA1(2,JJ) = AD1(J,MT1,MP1)+(AD1(J,MT1,MP1+1)-AD1(J,MT1,MP1))*RATP 
           DAD1(JJ) = DA1(2,JJ)+(DA1(1,JJ)-DA1(2,JJ))*RATT 
        end do
        GRAD = DAD1(1)+RATCC*(DAD1(2)-DAD1(1)) 
     endif

     !********************* CALCOLO PESO MOLECOLARE *************************
     if(calcolapmol == 1) then
        do JJ=1,2 
            J = ind(JJ)
           PU(1,JJ) = PM(J,MT2,MP2)+(PM(J,MT2,MP2+1)-PM(J,MT2,MP2))*RATP 
           PU(2,JJ) = PM(J,MT1,MP1)+(PM(J,MT1,MP1+1)-PM(J,MT1,MP1))*RATP 
           PMU(JJ) = PU(2,JJ)+(PU(1,JJ)-PU(2,JJ))*RATT 
        end do
        PMOL = PMU(1)+RATCC*(PMU(2)-PMU(1)) 
     endif
  else 
     !********************************* T>6 *********************************
     TER = T-6. 
     MT1 = int(TER*5.)+40 
     MT2 = MT1+1 
     NCI = 1 
     NCF = 5 
     FOR1 = P-PP(MT1) 
     FOR2 = P-PP(MT2) 
     if(FOR1 < 0. .or. FOR2 < 0.) goto 50 
     if(FOR1 > 16. .or. FOR2 > 16.) goto 50 
     MP1 = int(FOR1*2.)+1  
     MP2 = int(FOR2*2.)+1 

     RATT = (T-TT(MT1))/(TT(MT2)-TT(MT1)) 
     RATP = (P-(PP(MT1)+real(MP1-1)/2.))*2. 

     !************************** CALCOLO DENSITA' ************************** 
     if(calcolarho == 1) then
        do JJ=NCI,NCF 
           J = JJ 
           DE1(1,JJ) = RO1(J,MT2,MP2)+(RO1(J,MT2,MP2+1)-RO1(J,MT2,MP2))*RATP 
           DE1(2,JJ) = RO1(J,MT1,MP1)+(RO1(J,MT1,MP1+1)-RO1(J,MT1,MP1))*RATP 
           DEN1(JJ) = DE1(2,JJ)+(DE1(1,JJ)-DE1(2,JJ))*RATT 
           DEN1(JJ) = exp(UM*DEN1(JJ)) 
        end do
        RHO = 1./(XX(1)/DEN1(1)+XX(3)/DEN1(2)+XX(4)/DEN1(3)+XX(5)/DEN1(4)+ &
             XX(6)/DEN1(5)) 
     endif
     

     !********************* CALCOLO CALORE SPECIFICO ************************
     if(calcolacspe == 1 .or. calcolagrad == 1) then
        do JJ=NCI,NCF 
           J = JJ 
           CS1(1,JJ) = CP1(J,MT2,MP2)+(CP1(J,MT2,MP2+1)-CP1(J,MT2,MP2))*RATP 
           CS1(2,JJ) = CP1(J,MT1,MP1)+(CP1(J,MT1,MP1+1)-CP1(J,MT1,MP1))*RATP 
           CSP1(JJ) = CS1(2,JJ)+(CS1(1,JJ)-CS1(2,JJ))*RATT 
           CSP1(JJ) = exp(CSP1(JJ)*UM) 
        end do
        CSPEtmp = XX(1)*CSP1(1)+XX(3)*CSP1(2)+XX(4)*CSP1(3)+XX(5)*CSP1(4)+ &
             XX(6)*CSP1(5) 
     endif
     if(calcolacspe == 1) CSPE = CSPEtmp

     !********************* CALCOLO GRADIENTE ADIAB *************************
     if(calcolagrad == 1) then
        do JJ=NCI,NCF 
           J = JJ 
           DA1(1,JJ) = AD1(J,MT2,MP2)+(AD1(J,MT2,MP2+1)-AD1(J,MT2,MP2))*RATP 
           DA1(2,JJ) = AD1(J,MT1,MP1)+(AD1(J,MT1,MP1+1)-AD1(J,MT1,MP1))*RATP 
           DAD1(JJ) = DA1(2,JJ)+(DA1(1,JJ)-DA1(2,JJ))*RATT 
        end do
        GRAD = (CSP1(1)*XX(1)*DAD1(1)+CSP1(2)*XX(3)*DAD1(2)+CSP1(3)* &
             XX(4)*DAD1(3)+CSP1(4)*XX(5)*DAD1(4)+CSP1(5)*XX(6)*DAD1(5))/CSPEtmp 
     endif

     !********************* CALCOLO PESO MOLECOLARE *************************
     if(calcolapmol == 1) then
        PMOL = 1./(XX(1)*2.+XX(3)*3./4.+XX(4)*7./12.+XX(5)*4./7.+XX(6)*9./16.) 
     endif
  endif
  !***********************************************************************

  XX(6) = XX(6)-DIFF 
  if(XX(1) < 1.d-39) XX(1) = 0. 
  if(XX(3) < 1.d-39) XX(3) = 0. 
  if(XX(4) < 1.d-39) XX(4) = 0. 
  if(XX(5) < 1.d-39) XX(5) = 0. 
  if(XX(6) < 1.d-39) XX(6) = 0. 

  if(raccordo == 1) then
     ! calcolo un raccordo a banda tra xx(4) = Clow e xx(4) = Cup
     ! sopra Cup uso solo Oscar
     if(calcolarho == 1) then
        rho = rhos + (rho-rhos)/(Cup-Clow)*(XX(4) - Clow)
     endif
     if(calcolagrad == 1) then
        grad = grads + (grad-grads)/(Cup-Clow)*(XX(4) - Clow)
     endif
     if(calcolacspe == 1) then
        cspe = cspes + (cspe-cspes)/(Cup-Clow)*(XX(4) - Clow)
     endif
     if(calcolapmol == 1) then
        pmol = pmols + (pmol-pmols)/(Cup-Clow)*(XX(4) - Clow)
     endif

  endif

  return 

50 continue 
  write(*,103) 
  write(2,103) 
  write(2,104)TE,XT,PR,Pnew 
  write(*,104)TE,XT,PR,10.**Pnew,FOR1,FOR2 
  write(66,*)'152 - state'
  write(66,103)
  write(66,*) "caller = ", caller
  stop 

100 format(72A1) 
101 format(2X,72A1) 
102 format(11F7.3) 
103 format('  USCITO FUORI TABELLE FISICHE TE, XT,') 
104 format(1X,E10.3) 
321 format(1X,'STATE: (P, T): ',1p, 2(e12.4,1x)) 
119 format(13F10.5) 

end subroutine STATE
