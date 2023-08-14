subroutine EPSI (RO,T,JF,NABLA,HT,EPS,EPSA,ECAR,DELTA) 
  use interfaccia
  use mistura
  use fisica
  use chequ
  use nummod
  use equil
  use neutri
  use chim
  use varie
  use chimic
  use overshoot
  use mesh
  use tempi_scala
  use accrescimento
  implicit none

  real :: RO, T, HT, EPS, EPSA, ECAR, DELTA
  integer :: JF, NABLA

  real :: XZORA 
  real,dimension(35) :: E
  integer,parameter :: IDEBUG = 0

  real,parameter,dimension(MELE) ::  Z = (/2.,1.,2., 6., 7., 8., 8.,10., &
       10.,12.,12.,12.,14.,0.,11., 6., 7., 8., 9.,10.,26.,3.,3.,4., 5.,1./)
  real,parameter,dimension(MELE) :: AT = (/3.,1.,4.,12.,14.,16.,18.,20., &
       22.,24.,25.,26.,28.,1.,23.,13.,15.,17.,19.,21.,56.,6.,7.,9.,11.,2./)

  ! controllare Q13 e Q32
  real,dimension(35) :: Q = (/ &
       1.943,   7.162,  7.551,  2.216,  7.297, &
       4.415,   4.966,  4.014,  0.600,  4.730, &
       1.192,   0.586,  7.993,  9.667,  8.114, &
       1.675,   2.431,  9.316,  6.739,  8.794, &
       -0.478, 10.615,  2.377, 11.693,  2.242, &
       4.621,   9.986,  7.274,  6.671, 12.859, &
       1.0,    11.004,  14.031, 3.616,  5.493 /)

  real,parameter :: AN = 6.022d23
  integer,parameter :: ISHORT = 1

  integer :: kk, i, ini, ifi, ir, n, jj, k, ni, j
  real :: xfe, pas, epsp, epsc, xhe3
  real :: dmass, rpp, zeta, shep, schep, xhi, dfe_fe, xtest, xymax
  real :: rd, rli6, rli7, rbe, rb 
  real :: t6, t612, t9, t92, t93, t912, t913, t923, t932
  real :: xx1old, roan, tmp, rapp1

  real :: schermo, sbep, sbee, cd, cli6,&
       cli7, cbe, cbe1, cbe2, cb, skli
  integer :: z1,z2, ncache, mc, nel

  real,dimension(nsezi) :: S
  real,dimension(nsezi+2) :: Srho
  real,dimension(6) :: xxorig

  real,dimension(MELE) :: YV, YF, Yv1 
  integer,dimension(MELE) :: NEQUI
  real,save,dimension(MELE) :: yfs

  type(pow_of_T) :: powT

  real,save :: SSum(nsezi+2)
  integer,save :: conv, isconv, stopconv, convettivo, jfstart, jfstop
  integer :: newconv

  real :: ADMM_M, admm_MP1
  logical :: printout = .false.
  real :: rr, pp, em, pm, rm, emm, hpmm, hpm
  integer :: oishel, m
  real :: GICOVER, GICOVERR, ERRO, DIFOV, GICO

  real :: tau_li6, tau_li7, tau_be9, tau_b11, tau_d

  XFE = XME*DEFAUFe 
  !***** ISHORT SETTA COMBUSTIONE ELIO CORTA (7 ELEMENTI) SE =1) *********
  PAS = HT*3.1558d7 
  EPSP = 0. 
  EPSC = 0. 
  EPSA = 0. 
  ECAR = 0. 
  EPS = 0. 
  DELTA = 1. 
  if( JF == 1 ) then 
     do KK=1,5 
        COEFF(KK) = 0.
     end do
     do KK=1,8 
        PPNEU(KK) = 0.
     end do
     ISHEL = 0 
     ! overshooting
     L3 = 0
  endif
!  if(T < 1.d5) then 
!     return 
!  end if
!  if(T < 1.d6 .and. NMD >= 40) then 
!     return 
!  end if
  if(XX(1) <= 0. .and. XX(3) > 0. .and. T <= 4.d7) then 
     return 
  end if
  if(XX(1) <= 0. .and. XX(3) <=  0. .and. T <= 4.d8) then 
     return 
  end if
  if(XX(4) <= 0. .and. T > 3.d9) then 
     write(*,150) JF
     write(2,150) JF
     write(66,*) '50 - epsi'
     write(66,150) JF
     stop
  end if
  ! ******** IL RESTO DEL PROGRAMMA INVERTE IDROGENO E ELIO3 ************ 
  XHE3 = XX(2) 
  XX(2) = XX(1) 
  XX(1) = XHE3 
  !************************ INIZIALIZATION RAP-NEW ***********************

  ! controllo di non avere abbondanze di elementi leggeri negative
  do I=22,MELE 
     if(xx(i) < 0.) then
        xx(i) = 1.d-40
     endif
  end do

  ! **** NEQUI(K) = 0  SETTA L' ELEMENTO K ALL'EQUILIBRIO (SE RADIAT) **** 
  do I=1,MELE 
     NEQUI(I) = 1 
  end do

  ! ** CONTROLLO FASE EVOLUTIVA:                                        **
  ! ** IR=1 --> COMBUST. IDROGENO     IR=2 --> COMBUST. IDROGENO + ELIO **
  ! ** IR=3 --> COMBUST. ELIO         IR=4 --> COMBUST. CARBONIO        **
  if(XX(2) > 0. .and. T <= 5.d8) then 
     INI = 1       ! indice elemento iniziale
     if(T < 1.d8) then 
        IFI = 6     ! indice elemento finale
        IR = 1 
     else 
        IFI = 9     ! indice elemento finale
        IR = 2 
     endif
  endif
  if(XX(2) <= 0. .and. T <= 5.d8) then 
     INI = 3 
     IFI = 14 
     if( ISHORT == 1 ) IFI = 9 
     IR = 3 
  endif
  if(T > 5.d8) then 
     INI = 2 
     IFI = 20 
     IR = 4 
     NEQUI(2) = 0 
     NEQUI(3) = 0 
  endif

  N = IFI-INI+1    ! conta il numero di elementi su cui lavorare
  do JJ=1,31 
     E(JJ) = 0. 
  end do

  do I=1,MELE 
     YV(I) = XX(I)/AT(I) 
     YF(I) = YV(I) 
  end do

  ! calcolo potenze di temperatura utili in cross
  T6 = T*1.d-6   
  T612 = sqrt(T6)
  T9 = T*1.d-9
  T92 = T9*T9
  T93 = T92*T9
  T912 = sqrt(T9)
  T913 = T9**(1./3.)
  T923 = T913*T913
  T932 = T9*T912

  powt%T9 = T9
  powt%T92 = T92
  powt%T93 = T93
  powt%T912 = T912
  powt%T913 = T913
  powt%T923 = T923
  powt%T932 = T932
  powt%T6 = T6
  powt%T612 = T612

  call CROSS(RO,T, YV, IR, ISHORT, S, powt, nabla)
  ncache = 1
  ROAN = RO*AN
  ! ************ CALCOLO COEFFICIENTE ENERGETICA NUCLEARE  ************** 
  if(IR <=  2) then 
     E(1) = Q(29)*S(29)*ROAN*YV(2)**2/2. 
     E(2) = Q(30)*S(30)*ROAN*YV(1)**2/2. 
     Q(31) = 1.586+(S(33)*11.407+17.397)/(1.+S(33))

     E(3) = Q(31)*S(31)*ROAN*YV(1)*YV(3)

     ! energia da deuterio S(34)
     E(35) = Q(35)*S(34)*ROAN*YV(2)*YV(26)

     ! ********* se c'e' deuterio devo de-commentare questa *****
     EPSP = 1.6022d-6*(E(1)+E(2)+E(3)+E(35))

     !******************** CALCOLO NEUTRINI PROTONE PROTONE ******
     DMASS = (G(5,JF+1)-G(5,JF))*1.d33
     PPNEU(1) = PPNEU(1)+S(29)*ROAN*YV(2)**2/2.*DMASS
     PPNEU(2) = PPNEU(2)+S(31)*ROAN*YV(1)*YV(3)/(S(33)+1.0)*DMASS
     PPNEU(3) = PPNEU(3)+S(31)*ROAN*YV(1)*YV(3)*S(33)/(S(33)+1.0)*DMASS
     E(4) = Q(32)*S(1)*ROAN*YV(2)*YV(4)
     E(5) = Q(33)*S(5)*ROAN*YV(2)*YV(5)
     E(6) = Q(34)*S(9)*ROAN*YV(2)*YV(6)
     EPSC = 1.6022d-6*(E(4)+E(5)+E(6))
     !****************** CALCOLO NEUTRINI pep ******************************
     RPP = S(29)*ROAN*YV(2)**2/2.
     ! giada
     PPNEU(7) = PPNEU(7)+1.102d-4*RO*(1.+XX(2))/2./T612*(1.+0.02 &
          *T6)*RPP*DMASS
     !******** CALCOLO SCHERMAGGIO ELETTRONICO PER Hep *******************
     ZETA = sqrt(sum(YV(1:MELE)*Z(1:MELE)*(Z(1:MELE) + 1.)))
     SCHep = exp(0.188*2*ZETA*sqrt(RO)/(T6*T612))
     !*************** CALCOLO NEUTRINI Hep ********************************
     SHep = SCHep*1.41d-13/T923*exp(-6.141/T913)
     PPNEU(8) = PPNEU(8)+SHep*ROAN*YV(2)*YV(1)*DMASS
     !************************* CALCOLO NEUTRINI CNO ***********************
     PPNEU(4) = PPNEU(4)+S(1)*ROAN*YV(2)*YV(4)*DMASS
     PPNEU(5) = PPNEU(5)+S(5)*ROAN*YV(2)*YV(5)*DMASS
     PPNEU(6) = PPNEU(6)+S(9)*ROAN*YV(2)*YV(6)*DMASS
  endif
  if(IR > 1) then
     E(7) = Q(28)*S(28)*RO**2*AN*YV(3)**3/6.
     E(8) = Q(2)*S(2)*ROAN*YV(3)*YV(4)
     E(9) = Q(6)*S(6)*ROAN*YV(3)*YV(5)
     E(10) = Q(10)*S(10)*ROAN*YV(3)*YV(6)
     E(11) = Q(14)*S(14)*ROAN*YV(3)*YV(7)
     E(12) = Q(18)*S(18)*ROAN*YV(3)*YV(8)
     E(13) = Q(21)*S(21)*ROAN*YV(3)*YV(9)
     E(14) = Q(22)*S(22)*ROAN*YV(3)*YV(9)
     E(15) = Q(27)*S(27)*ROAN*YV(3)*YV(10)
     do K=7,15
        EPSA = EPSA+1.6022d-6*E(K)
     end do
  endif
  if(IR == 4) then
     E(16) = Q(3)*S(3)*ROAN*YV(2)*YV(16)
     E(17) = Q(4)*S(4)*ROAN*YV(3)*YV(16)
     E(18) = Q(7)*S(7)*ROAN*YV(2)*YV(17)
     E(19) = Q(8)*S(8)*ROAN*YV(3)*YV(17)
     E(20) = Q(11)*S(11)*ROAN*YV(2)*YV(18)
     E(21) = Q(12)*S(12)*ROAN*YV(3)*YV(18)
     E(22) = Q(13)*S(13)*ROAN*YV(2)*YV(7)
     E(23) = Q(15)*S(15)*ROAN*YV(2)*YV(19)
     E(24) = Q(16)*S(16)*ROAN*YV(3)*YV(19)
     E(25) = Q(17)*S(17)*ROAN*YV(2)*YV(8)
     E(26) = Q(19)*S(19)*ROAN*YV(2)*YV(20)
     E(27) = Q(20)*S(20)*ROAN*YV(2)*YV(9)
     E(28) = Q(23)*S(23)*ROAN*YV(2)*YV(15)
     E(29) = Q(24)*S(24)*ROAN*YV(2)*YV(15)
     E(30) = Q(25)*S(25)*ROAN*YV(4)**2/2.
     E(31) = Q(26)*S(25)*ROAN*YV(4)**2/2.
     E(32) = Q(1)*S(1)*ROAN*YV(2)*YV(4)
     E(33) = Q(5)*S(5)*ROAN*YV(2)*YV(5)
     E(34) = Q(9)*S(9)*ROAN*YV(2)*YV(6)
     do K=16,34
        ECAR = ECAR+1.6022d-6*E(K)
     end do
  endif
  if(EPSP > 0. .and. EPSC > 0.) then
     DELTA = EPSP/(EPSP+EPSC)
  else
     DELTA = 1.
  endif
  EPS = ECAR+EPSA+EPSC+EPSP
  if( NABLA == 0 ) then
     XHI = XX(2)
     XX(2) = XX(1)
     XX(1) = XHI
     return
  endif

  ! metto da parte le abbondanze di H,He,CNO in ingresso ... 
  do j=1,6
     xxorig(j) = xx(j)
  end do
  xx1old = xx(2)

  ! calcolo S moltiplicate per rho (e per rapp1/2)
  RAPP1 = S(7)/(S(7)+S(32))
  Srho(1:nsezi) = RO * S
  Srho(28) = Srho(28) * RO
  Srho(nsezi+1) = RO * S(5) * rapp1
  Srho(nsezi+2) = RO * S(5) * (1.0 - rapp1)

  KECUIL(JF) = 0
  if(IR == 1 .and. KOVER == 1 .and. G(6,1) > 0.0 .and. ishel /= 1) then
     ! Inserimento overshooting
     OISHEL = ISHEL
     if(G(6,JF) < 0.0 .or. G(6,JF+1) < 0.0) then
        if(G(6,JF) >= 0.0 .and. G(6,JF+1) < 0.0) then
           PM = G(3,JF)*1.d17
           RM = G(1,JF)*1.d10
           EMM = G(5,JF)*1.d33
           HPMM = 1.884d8*PM*(RM**4)/EMM
           HPM = HPMM*1.d-33

           do M = JF, MAXME-1
              ADMM_m = G(5,M)-(G(5,JF)+par_OVER*HPM)
              ADMM_mp1 = G(5,M+1)-(G(5,JF)+par_OVER*HPM)
              if( ADMM_m <= 0.0 .and. ADMM_MP1 > 0.0 ) exit
           end do
           L3 = M
           if(printout) then
              GICOVER = G(5,L3+1)
              GICOVERR = G(5,JF)+par_OVER*HPM
              ERRO = G(5,L3+1)-(G(5,JF)+par_OVER*HPM)
              DIFOV = G(5,JF)+par_OVER*HPM-G(5,JF+1)
              GICO = G(5,JF+1)
              write(50,7362) NMD,GICOVER,GICOVERR,ERRO,GICO, &
                   DIFOV,HPM,JF, L3
           endif
        endif
        if(JF > L3) ISHEL = 1
     endif
  endif

  ! trattamento chimica standard
  if(IR == 1) then
     call IDRO(RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S, oishel)
  else if(IR == 2) then
     call IDRELI(RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S)
  else if(IR == 3) then
     call ELIO(RO,T,PAS,JF,INI,IFI,N,ISHORT, YV, YF, NEQUI, S)
  else if(IR == 4) then
     call CARBON(RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S)
  endif
  
  ! aggiorno le abbondanze in massa
  XX = YF*AT

  ! controllo la chimica
  call check_epsi(ir, nmd, jf, xx, xxorig, YV, YF, AT, nequi, &
       pas, RO, T, S, xfe)

  XHI = XX(2)
  XX(2) = XX(1)
  XX(1) = XHI

  !*****************************************************************
  !    calcoliamo le abbondanze di D, Li, Be, B
  !****************************************************************
  if(IR <= 2) then
   
     rd = YV(2)*YV(26)*Srho(35)
     rli6 = YV(2)*YV(22)*Srho(36)
     rli7 = YV(2)*YV(23)*Srho(37)
     rbe = YV(2)*YV(24)*Srho(38)
     rb = YV(2)*YV(25)*Srho(39)

!! ema: calcolo tempo scala  di distruzione degli 
!! elementi leggeri (in realta' e' l'inverso di tau) !!
   if(NMD >= 1)then
     Tau_nuc(22,JF) = (Srho(36)*YV(2)) !! Li6
     Tau_nuc(23,JF) = (Srho(37)*YV(2)) !! Li7
     Tau_nuc(24,JF) = (Srho(38)*YV(2)) !! Be9
     Tau_nuc(25,JF) = (Srho(39)*YV(2)) !! B11
     Tau_nuc(26,JF) = (Srho(35)*YV(2)) !! D

    ! *********VARIAZIONI DELLE ABBONDANZE NUOVO STILE ******
!     Tau_li6 = Tau_nuc(22,JF)*3.15d7
!     Tau_li7 = Tau_nuc(23,JF)*3.15d7
!     Tau_be9 = Tau_nuc(24,JF)*3.15d7
!     Tau_b11 = Tau_nuc(25,JF)*3.15d7
!     Tau_d = Tau_nuc(26,JF)*3.15d7

!     YF(22) = YV(22)*exp(-HT*Tau_li6) !! Li6
!     YF(23) = YV(23)*exp(-HT*Tau_li7) !! Li7
!     YF(24) = YV(24)*exp(-HT*Tau_be9) !! Be9
!     YF(25) = YV(25)*exp(-HT*Tau_b11) !! B11
!     YF(26) = YV(26)*exp(-HT*Tau_d) !! D
  end if

     ! *********VARIAZIONI DELLE ABBONDANZE STILE ******
     YF(22) = YV(22) - RLI6*PAS !! Li6
     YF(23) = YV(23) - RLI7*PAS !! Li7
     YF(24) = YV(24) - RBE*PAS !! Be9
     YF(25) = YV(25) - RB*PAS !! B11
     YF(26) = YV(26) - RD*PAS  !! D

     do JJ=22,26
        ! calcolo abbondanze in massa
        if(yf(jj) > 0.) then
           xx(jj) = yf(jj)*at(jj)
        else
           ! se negative, impostate a 0.
           xx(jj) = 0.
        endif
     end do

  endif

  !************************* DEBUGGING **********************************
  if(IDEBUG == 1 .and. JF == 1 ) then
     write(2,121) JF,NI,G(6,JF),T,RO,PAS
     do JJ=INI,IFI
        YV(JJ) = YV(JJ)*AT(JJ)
     end do
     write(2,124)
     write(2,122)YV(2),YV(1),(YV(J),J=3,MELE)
     write(2,125)
     write(2,122)(XX(J),J=1,MELE)
  endif
  return

121 format(1X,'MESH:',I4,'  N ITERA.:',I3,'  DG, T, RO, PASSO:' &
       ,1P,4E15.8)
141 format(1X,'MESH:',I4,'  N ITERA.:',I3,'  EUILIBR.:',I3)
122 format(1X,1P,10E12.5)
142 format(1X,1P,8E9.2)
124 format(1X,'COMPOSIZIONE CHIMICA INIZIALE')
125 format(1X,'COMPOSIZIONE CHIMICA FINALE')
150 format(1X,'SUPIN RICADDE E PIU NON PARVE FORA ALLO MESCIO',I4)
3459 format(2I4, 5D12.5)
7362 format(I5,6e13.5,2I5)

end subroutine EPSI




! ********************************************************************
! ****************** ROUTINE DI CONTROLLO CHIMICA ********************
! ********************************************************************

subroutine check_epsi(ir, nmd, jf, xx, xxorig, YV, YF, AT, nequi, &
     pas, RO, T, S, xfe)
  use parametri
  use chimic
  use varie
  
  implicit none
  
  integer :: ir, nmd, jf
  real,dimension(mele) :: xx, AT, Yv, YF, YV1
  real,dimension(6) :: xxorig
  real ::  pas, RO, T, xfe
  real,dimension(nsezi) :: S
  integer,dimension(MELE) :: NEQUI

  integer :: jj, INI, IFI, N, i
  integer,parameter :: idrocheck = 0
  real :: dfe_fe, xzora, XTEST, xymax

  ! controllo che l'idrogeno non sia fuori controllo
  ! se e' aumentato del 10% rispetto all'ingresso nella idro
  ! allora correggo la situazione usando le abbondanze del mesh precedente
  if(ir == 1 .and. xx(2) > 1.1*xxorig(2)) then
     ! imposto le abbondanze uguali a quelle del mesh precedente
     write(67,'("mod ",i5," XH mesh ",i5," old = ",e9.4," new = ",e9.4)') &
          nmd, jf, xxorig(2), xx(2)
     write(67,'(2(i5,1x),1x,6(e9.3,1x))') nmd, jf, (xx(JJ),JJ=1,6) 

     xx(2) = xxx(1,JF-1)
     xx(1) = xxx(2,JF-1)
     do JJ=3,6
        xx(JJ) = xxx(jj,jf-1)
     end do
     write(67,'(2(i5,1x),1x,6(e9.3,1x))') nmd, jf, (xx(JJ),JJ=1,6)
  endif

  ! controllo che la idro non abbia dato soluzioni errate per He3,CNO
  if(ir == 1 .and. (xx(1) < 0. .or. xx(4) < 0. .or. xx(5) < 0. &
       .or. xx(6) < 0.) ) then
     if(idrocheck == 1 ) then
        ! solo print
        write(64,'(2I5,6e14.5)') nmd, jf, (xxorig(jj),jj=1,6)
        write(64,'(2I5,6e14.5)') nmd, jf, (xx(jj),jj=1,6)
     endif
  endif
  
  ! QUESTO BLOCCHETTO CONTROLLA EVENTUALI SBORDI DELLA SOMMA H+HE 
  dfe_fe = (XX(21)-xfe)/xfe
  xzora = (1.0+dfe_fe)*XME
  if(IR < 3) then
     XTEST = XX(1)+XX(3)+XX(2)
     XYMAX = 1.0-XZORA
     if(XTEST > XYMAX) then
        !  evito oscillazioni della chimica 
        XX(3) = 1.0 - XX(2) - XX(1) - xzora
     endif
  endif

  ! controllo che He non sbordi, nel caso correggo l'abbondanza
  if(ir == 3 .and. XX(3) > 1.) then
     XX(3) = 1.0 - XX(2) - XX(1) - xzora
  endif

end subroutine check_epsi
