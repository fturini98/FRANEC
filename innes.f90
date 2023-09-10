subroutine INNES(ITCC,IPP,MAXMEin,ILEG,LBA, YB,ZB, He, Zeta, Alpha) 
  use interfaccia
  use mistura
  use fisica
  use strut
  use fitt
  use atm
  use chimic
  use indi
  use tempe
  use chim
  use numer
  use varie
  use costanti
  implicit none

  real :: YB, ZB, He, Zeta, Alpha
  integer :: ITCC, IPP, MAXMEin, ILEG, LBA

  real,parameter,dimension(180) :: ZN = (/0.0,2.d-7,5.d-7,1.d-6,5.d-6, &
       1.d-5,5.d-5,.0001,.0005,.001, &
       .0015,.002,.0025,.003,.004,.005,.006,.007,.008,.009,.01,.011,.012, &
       .013,.014,.015,.016,.017,.018,.019,.02,.021,.022,.023,.024,.025, &
       .026,.027,.028,.029,.03,.032,.034,.036,.038,.04,.042,.044,.046, &
       .048,.05,.052,.054,.056,.058,.06,.062,.064,.066,.068,.072,.076,.08, &
       .084,.088,.092,.096,.1,.104,.108,.112,.116,.12,.124,.128,.132,.136, &
       .14,.144,.148,.152,.156,.16,.164,.168,.172,.176,.18,.184,.188,.192, &
       .196,.2,.205,.21,.215,.22,.225,.23,.234,.242,.25,.26,.27,.28,.29, &
       .30,.3075,.315,.33,.345,.36,.375,.39,.405,.42,.435,.45,.465,.48, &
       .49,.5,.515,.53,.54,.55,.56,.57,.58,.59,.6,.61,.62,.63,.64,.65,.66, &
       .67,.68,.69,.7,.71,.72,.73,.74,.75,.76,.77,.78,.79,.8,.81,.82,.83, &
       .84,.85,.855,.86,.865,.87,.875,.88,.885,.89,.895,.9,.905,.91,.915, &
       .92,.9225,.925,.9275,.93,.932,.934,.936,.938,.939,.94/)

  real,dimension(180) :: ZNV
  real :: ELLPC = 0., ELLTC = 0., EMMU = 0. 

  integer :: k, ini, indu, counteras, counterat, stat
  ! regolare i parametri sottostanti che servono per impedire loop infiniti
  integer,parameter :: CASstop = 50, CATstop = 50

  real :: y, farz, emmn, emtot1, rtot, em, r, el, p, t, ro, rr 
  real :: as, at, au, av, toll, toll1, tceno, NCO, OCO, XMEin
  integer,save :: letto = 0
  real,save :: AMCO, YCO, CCO

  ! per fato
  real,save :: EL0, EL1, EL2, P0, P1, P2, R0, R1, R2, T0, T1, T2, DPC, DTC, &
       DELL, DTEF, EL0P, EL1P, EL2P, P0P, P1P, P2P, R0P, R1P, R2P, &
       T0P, T1P, T2P, CORZ

  counteras = 0
  counterat = 0
  if(ISUB == 1 .and. IBAT == 0.) write(*,444) 

  ILEG = 0 
  if(ITCC /= 1 .and. letto /= 1) then
     EMAXH = 0. 
     EMAHE = 0. 
     letto = 1
     read(LBA,*,end=122) Y
     read(LBA,*) XME !Legge la metallicita dal modstart
     read(LBA,*) FARZ
     read(LBA,*) ALFA
     read(LBA,*) BMAG
     read(LBA,*) AMCO
     read(LBA,*) ETAFPR
     read(LBA,*) CCO 
     YB = Y 
     ZB = XME 
     FRAZ = 1.-FARZ/100. 
     read(LBA,*) EMMU

     read(LBA,*) ELLOG
     read(LBA,*) TEFF
     read(LBA,*) ELLPC
     read(LBA,*) ELLTC
     read(LBA,*) LAST
     read(LBA,*) CORZ 
     if(CORZ <= .09) CORZ = 1. 

!!$     ! stat = 0  significa che non sono alla fine del file
!!$     read(LBA,*,iostat=stat) He
!!$     if(stat == 0) then
!!$        read(LBA,*) Zeta
!!$        read(LBA,*) Alpha
!!$        do_relax = .true.
!!$        write(67,*) "@ ================================"
!!$        write(67,*) "@      Rilassamento modstart          "
!!$        write(67,*) "@ ================================"
!!$        write(67,'(" @ ","He =",f7.4,", z =",e11.4,", Alpha =",f6.3)') &
!!$             He, Zeta, Alpha
!!$        write(67,*) "@ ================================"
!!$        write(67,*)
!!$     endif

     ZNV(1:maxmein) = ZN(1:maxmein)*FRAZ/.94 

     LAXME = MAXMEin 
     EMMN = EMMU 
     if(IREAD /= 4) then
        if(IREAD /= 2) EMAXH = AMCO 
        if(ISUB == 1) EMAHE = ETAFPR 
     endif

     ! aggiusto le abbondanze se parto da HB (pepper)
     if(iread == 4 .and. etafpr > 0.) call inout_hb(NCO,OCO,XMEin)

     ! vettore della chimica XX:
     !  1=XH,     2=XHE3,   3=XHE4,    4=XC12,   5=XN14,
     !  6=XO16,   7=XO18,   8=XNE20,   9=XNE22, 10=XMG24,
     ! 11=XMG25, 12=XMG26, 13=XSI28,  14=XNEU,  15=XNA23,
     ! 16=XC13,  17=XN15,  18=XO17,   19=XF19,  20=XNE21,
     ! 21=XFE,   22=XLI6,  23=XLI7,   24=XBE,   25=XB
     ! 26=XD

     
     XX(2) = XME*DEFAUHe3
     XX(3)= Y 
     XX(4) = XME*DEFAUC 
     XX(5) = XME*DEFAUN 
     XX(6) = XME*DEFAUO
     XX(7) = 0.
     XX(8) = 0.
     XX(9) = 0.
     XX(10) = 0.
     XX(11) = 0.
     XX(12) = 0.
     XX(13) = 0.
     XX(14) = 0. 
     XX(15) = 0.
     XX(16) = 0.
     XX(17) = 0.
     XX(18) = 0.
     XX(19) = 0. 
     XX(20) = 0.
     XX(21) = XME*DEFAUFe
     XX(22) = XME*DEFAULi6 
     XX(23) = XME*DEFAULi7 
     XX(24) = XME*DEFAUBe 
     XX(25) = XME*DEFAUB 
     XX(26) = DEFAUD 

     XX(1) = 1. - Y - XME - XX(2) - sum(XX(22:26))
     XH = XX(1)
     
     if(IREAD == 4) YCO = 1.-CCO-XMEin 
     INI = 0 
  endif

  TEF = 10.**TEFF 
  ELL = 1.d33 * (10.**ELLOG) * Lsun
  EMTOT = EMMU * Msun 
  EMTOT1 = 1.d33 * EMMU * Msun 
  PCEN = 10.**ELLPC 
  TCEN = 10.**ELLTC 
  
  do K=1,MAXMEin 
     xxx(1:26,k) = xx(1:26)
     G(5,K) = EMTOT1*ZNV(K) 
  end do

  if(ISUB == 1) goto 26 
  if(IREAD == 4) call HBREZO(AMCO, YCO, CCO, NCO, OCO)
  MAXMEin = LAXME 
  if(ITCC /= 0) then 
     if(IPP == 0) TCEN = TCEN+CCO 
     ELLTC = log10(TCEN) 
  endif
  write(2,666) EMMU,ELLOG,ELLPC,ELLTC,TEFF 

  toll = .0005   ! valore originario: .0005
  DELL = toll *ELL 
  DTEF = toll *TEF 
  DPC = toll *PCEN 
  DTC = toll *TCEN 
  NABLA = 0 
  IOTA = 0 
15 RTOT = sqrt(ELL/(7.125d-4*(TEF**4))) 
  ! integra dal centro verso il punto di fit
  call VEIOVE(EM,R,EL,P,T,RO,1,IPP)
  G(6,1) = G(6,2) 
  GG(6,1) = GG(6,2) 
  if(P <= 0. .or. T <= 0.) then      ! impossibile convergere
     write(67,*) 'INNES (P,T)', P, T
     ILEG = 1 
     return 
  endif

  if(NABLA < 0) then 
     ! Tc +DTc, Pc
     EL1 = EL 
     R1 = R 
     P1 = P 
     T1 = T 

     NABLA = 1 
     TCEN = TCEN-DTC 
     PCEN = PCEN+DPC 
     goto 15 
  else if(NABLA == 0) then 
     EL0 = EL 
     R0 = R 
     P0 = P 
     T0 = T 
     ! valori interni presi da veiove
     write(2,111) R,EL,P,T 
  else 
     ! Tc, Pc + DPc
     EL2 = EL 
     R2 = R 
     P2 = P 
     T2 = T 

     NABLA = 0 
     IOTA = -1 
     PCEN = PCEN-DPC 

     ! esterno con Te = Te + DTe
     TEF = TEF+DTEF 
     TEFF = log10(TEF) 
  endif

26 RTOT = sqrt( ELL / (4.0d0 * pigre * ssb * (TEF**4)) ) 
  !!RTOT = sqrt(ELL/(7.125d-4*(TEF**4))) 
  ELLOG = log10(ELL/(1.d33 * Lsun)) 
  TEFF = log10(TEF) 
  ! chiamata alla quatm, che a sua volta chiama atmos 3 volte con
  ! (L,Te) (L,Te+Dte) (L+DL,Te) e calcola le derivate alla base della 
  ! subatmosfera
  call QUATM(-1,INI,MAXMEin,LBA,EMTOT) 
  if(ISUB == 1) then
     write(66,*) '92 - INNES'
     write(66,*) 'ISUB = 1'
     stop
  endif
  R = RB*1.d10 
  RR = R 
  P = PB*1.d17 
  T = TB*1.d6 
  EL = ELL 
  EM = EMTOT1*FRAZ 
  G(1,MAXMEin) = R 
  G(2,MAXMEin) = ELL 
  G(3,MAXMEin) = P 
  G(4,MAXMEin) = T 

  ! integra dall'esterno fino al punto di fit
  call VEIOVE(EM,R,EL,P,T,RO,2,IPP) 
  G(6,MAXMEin) = G(6,MAXMEin-1) 
  GG(6,MAXMEin) = GG(6,MAXMEin-1) 
  if(R <= 0. .or. EL <= 0.) then        ! impossibile convergere
     write(67,*) 'INNES (R,EL) (ITCC,IPP)', R, EL, ITCC, IPP
     ILEG = 1 
     return 
  endif

  if(IOTA < 0) then
     ! Te + DTe, L
     EL1P = EL 
     R1P = R 
     P1P = P 
     T1P = T 

     IOTA = 1 
     TEF = TEF-DTEF 
     TEFF = log10(TEF) 
     ELL = ELL+DELL 
     ELLOG = log10(ELL/(1.d33*Lsun)) 
     goto 26 
  else if(IOTA == 0) then 
     EL0P = EL 
     R0P = R 
     P0P = P 
     T0P = T 
     write(2,222) R,EL,P,T 
     write(2,333) 
  else 
     ! Te, L + DL
     EL2P = EL 
     R2P = R 
     P2P = P 
     T2P = T 

     ELL = ELL-DELL 
     ELLOG = log10(ELL/(1.d33*Lsun)) 
     IOTA = 0 
     call FATO(INDU,AS,AT,AU,AV,IPP, EL0,EL1,EL2,P0,P1,P2,R0,R1,R2,T0,T1,T2, &
          DPC,DTC,DELL,DTEF,EL0P,EL1P,EL2P,P0P,P1P,P2P,R0P,R1P,R2P, &
          T0P,T1P,T2P,CORZ) 
     if(INDU == 1) then
        write(67,*) 'INNES: INDU = 1'
        ILEG = 1 
        return 
     endif
     goto 15 
  endif

  ! indicatori della convergenza del fit
  AS = EL0-EL0P 
  AT = R0-R0P 
  AU = P0-P0P 
  AV = T0-T0P 
  toll1 = 0.001   ! valore iniziale 0.001

  if(.not. (IREAD == 2 .and. IPP == 0)) then
     if(abs(AS) - toll1 *EL0 > 0) then
        NABLA = -1 
        TCEN = TCEN+DTC 
        ! per verificare se va in loop
        counteras = counteras + 1
        if(counteras > CASstop) then
           write(66, *) '90 - INNES'
           write(66, *) 'massimo indice ecceduto per AS'
           stop
        endif
        goto 15 
     endif
  endif
  if(abs(AT)- toll1*R0 > 0 .or. abs(AU)- toll1*P0 > 0 .or. &
       abs(AV)- toll1*T0 > 0) then
     NABLA = -1 
     TCEN = TCEN+DTC 
     ! per verificare se va in loop
     counterat = counterat + 1
     if(counterat > CATstop) then
        write(66, *) '91 - INNES'
        write(66, *) 'massimo indice ecceduto per AT'
        stop
     endif
     goto 15 
  endif
  write(2,888) 
  ELLOG = log10(G(2,MAXMEin)/(1.d33*Lsun)) 
  TEFF = log10(TEF) 
  ELLPC = log10(PCEN) 
  ELLTC = log10(TCEN)
  G(1,1) = 0. 
  G(2,1) = 0. 
  G(3,1) = PCEN 
  G(4,1) = TCEN 
  return 

122 stop 

888 format(/,1X,'CONVERGENZA RAGGIUNTA',/) 
666 format(5F8.4,I4,F4.1) 
333 format(///) 
111 format(1X,'INTERNO',6X,'R',E12.5,3X,'L',E12.5,3X,'P',E12.5,3X,'T' &
       ,E12.5)   
222 format(1X,'ESTERNO',6X,'R',E12.5,3X,'L',E12.5,3X,'P',E12.5,3X,'T' &
       ,E12.5) 
444 format(5X,'M',7X,'L',6X,'TE',5X,'Y',7X,'Z',8X,'P',7X,'T',         &
       6X,'R',5X,'RO',/)            
end subroutine INNES


