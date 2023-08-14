subroutine STAMPA(TEMPO,DMAS,NMOD,NABLA,KSD,MAXMEin,MAXMV,IFVA,IFMO, &
     lumin, TMAX, nsmorza, fase, logCNO, kover, nsole)
  use interfaccia
  use fisica
  use mesh
  use strut
  use chequ
  use nummod
  use serchi
  use nuconv
  use atmosfere
  use neutri
  use second
  use chimic
  use tempe
  use chim
  use numer
  use varie
  use print1
  use fasievolut
  use mistura
  use costanti
  use def_io

  implicit none

  real :: TEMPO, DMAS, TMAX
  integer :: NMOD, NABLA, KSD, MAXMEin, MAXMV, IFVA, IFMO, nsmorza, nsole
  integer :: fase, logCNO, kover
  type(luminosity) :: lumin

  real,dimension(8) :: A, FF, HH
  real,dimension(9) :: CV, RPQ
  real,dimension(4) :: AA, DFIS, DFIX
  real,dimension(mele) :: ES
  real,dimension(11,10) :: CC
  real,dimension(11,lim) :: STAR
  real,dimension(6,10) :: DD
  real,dimension(13) :: EE

  integer,dimension(4) :: IVME
  integer,dimension(10) :: IMINT, IMEST
  real,dimension(lim) :: DRAD
  real,dimension(lim) :: E1, E2, E3, E4, E5, E6, E7, E13
  real,dimension(8) :: CL_seg, GA_seg, flussi
  real :: AGE

  real,parameter,dimension(8) :: GAMATRIX = (/1.18d-9, 7.32d-9, 2.43d-6, &
       6.18d-9, 11.6d-9, 11.7d-9, 21.5d-9, 0.73d-5/)
  real,parameter,dimension(8) :: CLMATRIX = (/0., 0.24d-9, 1.09d-6, 0.17d-9, &
       0.68d-9, 0.69d-9, 1.6d-9, 0.39d-5/)
  character(len=8),parameter :: ATM_USATA(3) = (/'BH05)   ','CK03)   ', &
       'VECCHIA)'/)

  integer :: k, j, il, l1, n1, icov, ienv, kcore, icore, kkcore, ka2
  integer :: kshhe, kshhy, ntmax, nhmax, nhemax, ngrmax, nnumax, l, l3
  integer :: n, n3, kk, iwd, jf, jff, jfa, jrona, jrana 
  integer :: i, lma, lmb, lea, leb, km, ll, jj,  ii, idx_atm 
  real :: rag, rtot, emsol, genpp, gencar, gencno, genhe, gengra
  real :: gennu, corco, corhe, elsu, temp, wa, wb, wc, rr, zl, pp, tt 
  real :: em, xhe, tm, pr, epsn, gravi, rho, dad, csp, pmau
  real :: epse, delta, epneu, cap, grsad, p1, el1, em1, gi, acca, esp
  real :: epsh, roc, delce, rohemx, rohmx, xha, xhb, xea, xeb, rat, fisv
  real :: gm, xi, xhe3, xhe4, xc, xn, xo, xo18, xfe, xd, xli6, xli7
  real :: xbe, xb, xmg26, xsi28, xne20, xne22, xmg24, xmg25, xna23, xs32
  real :: vari, dnorm, ga_tot, cl_tot, xcen, zx_sup
  real :: summa, ecar, HorHe, z, dmt, dm, xfe1, dfe_fe, xzora, cht, logg
  integer,save :: neverset = 1
  integer,parameter :: printover = 1
  real :: envconv, feh

  character(len=15) :: caller="stampa        "
  
  ! produzione energia nucleare nelle zone CNO, He, H
  lumin%pn = 0.0

  ! per ripartenza dopo blocco
  if(lumin%M_bordo_conv >= 0.0) neverset = 0
  
  MAXME = MAXMEin
  MAXMIV = MAXMV

  imint = 0
  imest = 0
  cc = 0.
  dd = 0.

  bb = 0.
  cv = 0.
  rpq = 0.
  ff = 0.
  hh = 0.
  aa = 0.

  ! QUANTITA' VARIE
  AGE = log10(TEMPO)
  A(1) = log10(HT1)
  A(2) = ELLOG
  A(3) = TEFF
  RAG  = (1. / Rsun) * sqrt(1.d13 * Lsun / (4.d0 * pigre * ssb) ) * &
       & sqrt((10.**ELLOG)) / (10.**(2.*TEFF))
  RTOT = RAG * Rsun
  EMSOL = EMTOT / Msun
  A(4) = log10(RAG) 
  LOGG = 13.0 + log10(Ggrav * EMTOT / (RTOT * RTOT))
  A(5) = LOGG
  A(6) = EMSOL
  A(7) = ETAFPR
  A(8) = DMAS/(MSUN * HT1)

  ! selezione l'etichetta di atmosfera: 3 vecchia, 2 CK, 1 BH
  idx_atm = 3 - N_CK03 - 2*N_BH05   

  write(*,100) AGE , A , MAIS , MAXMEin , NMOD, ATM_USATA(idx_atm)
  write(2,100) AGE , A , MAIS , MAXMEin , NMOD, ATM_USATA(idx_atm)

  ! STAMPA LA CHIMICA SUPERFICIALE
  write(ioCHIMICASUP,1981) NMOD,AGE,XXX(1,MAXMEin-1),XXX(26,MAXMEin-1), &
       (XXX(il,MAXMEin-1),il=2,6),XXX(21,MAXMEin-1), & 
       (XXX(il,MAXMEin-1),il=22,25)

  ! ho il Sole: output
  if(nmd == nsole) then
     zx_sup = (1.d0 - XXX(1,MAXMEin-1) - XXX(3,MAXMEin-1) ) / XXX(1,MAXMEin-1)
     open(457,file='Modello_Sole.DAT')
       write(457,'(5(1x,1PE17.10),0P)') a(2), a(3), rag, a(5), zx_sup
     close(457)
  end if

  ! CERCA ZONE CONVETTIVE
  nsmorza = 0

  L1    = 0
  N1    = 0
  ICOV  = 0
  IENV  = 0
  KCORE = 0
  ICORE = 0
  if( G(6, 1) >= 0. ) then
     ICOV =  1
     KCORE = MAXMEin
     do K = 2,MAXMEin
        if( G(6, K) <= 0. ) then
           KKCORE = K
           KCORE = K
           exit
        endif
     end do
  endif

  KA2 = 10
  if( KCORE > 0 .and. KCORE < MAXMEin ) ICORE = 1
  if( ICORE == 1 ) KA2 = KCORE + 2
  do K = KA2 , MAXMEin
     if( G(6, K)  >=  0. .and. G(6, K-1) < 0. ) then
        L1 = L1 + 1
        IMINT(L1) = K
     else if( G(6, K) < 0. .and. G(6, K-1) >= 0. ) then
        if(K-IMINT(L1) > 4) then
           N1 = N1 + 1
           IMEST(N1) = K
        else
           IMINT(L1) = 0
           L1 = L1-1 
        endif
     endif
  end do
  if( L1 == N1 + 1 ) IENV = 1

  ! STAMPA ANDAMENTI FISICA ALL'INTERNO
  KSHHE  = 0
  KSHHY  = 0
  NTMAX  = 1
  NHMAX  = 1
  NHEMAX = 1
  NGRMAX = 1
  NNUMAX = 1
  L      = 1
  L3     = 1
  N      = 1
  N3     = 1
  TMAX   = G(4, 1)
  AHMAX  = 0.
  AHEMAX = 0.
  AGRMAX = 0.
  ANUMAX = 0.
  GENPP  = 0.
  GENCAR = 0.
  GENCNO = 0.
  GENHE  = 0.
  GENGRA = 0.
  GENNU  = 0.
  CORCO  = 0.
  CORHE  = 0.
  do J = 1,6
     G(J, MAXMEin + 1) = G(J, MAXMEin)
  end do

  do K = 1,MAXMEin
     KK = K + 1
     do J = 1 , 5
        DG(J, K) = G(J, KK) - G(J, K)
     end do
  end do

  ELSU = G(2, MAXMEin)
  RTOT = (1.d6 / sqrt(4.d0 * pigre * ssb) ) * sqrt(ELSU)/ (10**(2.*TEFF))

  ! LOOP PRINCIPALE
  IWD = 0

  if(NABLA == 3) then
     call ATMOS(ELLOG,TEFF,WA,WB,WC,EMTOT,1,IWD)
     ! intestazione FISICA.DAT
     write(ioFISICA,'("#   MOD   AGE       MESH")')
     write(ioFISICA,'("# ", i5,1x,f14.10,1x,i5)') nmd, age, maxmein
     write(ioFISICA,108)
  endif

  do JF = 1,MAXMEin
     JFF = JF + 1
     RR = ( G(1, JF) + G(1, JFF) )/2.
     ZL = ( G(2, JF) + G(2, JFF) )/2.
     PP = ( G(3, JF) + G(3, JFF) )/2.
     TT = ( G(4, JF) + G(4, JFF) )/2.
     EM = ( G(5, JF) + G(5, JFF) )/2.
     do K = 1,6
        XX(K) = XXX(K, JF)
     end do
     XHE = XX(2) + XX(3)
     TM = TT * 1.d6
     PR = PP * 1.d17 - arad_3 * TM**4
     EPSN  = 0.
     GRAVI = 0.

     call STATE(caller, PR, TM, RHO, DAD, CSP, PMAU)

     call EPSI(RHO, TM, JF, 0, HT1, EPSN, EPSE, ECAR, DELTA)

     call NEUTR(RHO, TM, EPNEU)
     if ( NMD > 10 .or. IREAD == 2 ) then
        call EPSIG(PP,TM,RHO,CSP,HT1,JF,GRAVI,MAXMV,DAD,0)
     endif

     call KAPPA(RHO, TM, CAP, XHE)
     DRAD(JF) = cte_grad * CAP * PP * ZL/(EM * TT**4)

     GRSAD = 0.
     if(DRAD(JF) >= DAD) then
        P1 = PP*100.
        EL1 = ZL*100.
        EM1 = EM*1000.
        GI = 1.d13 * Ggrav * EM / (RR * RR)
        call SUPERA(RHO,P1,TT,CAP,PMAU,CSP,DAD,GI,GRSAD,ACCA,DRAD(JF))
     endif

     ESP    = EPSN + EPNEU + GRAVI
     EPSH   = EPSN - EPSE - ECAR
     GENPP  = GENPP  + DELTA * EPSH * 10. * DG(5, JF)
     GENCNO = GENCNO + (1. - DELTA ) * EPSH * 10. * DG(5, JF)
     GENHE  = GENHE  + EPSE  * 10. * DG(5, JF)
     GENCAR = GENCAR + ECAR  * 10. * DG(5, JF)
     GENGRA = GENGRA + GRAVI * 10. * DG(5, JF)
     GENNU  = GENNU  + EPNEU * 10. * DG(5, JF)
     if( JF == 1 ) ROC   = RHO
     if( JF == 1 ) DELCE = DELTA
     ! CERCA MASSIMO DELLE 3A 
     if( EPSE > AHEMAX  ) then
        AHEMAX = EPSE
        NHEMAX = JF
        ROHEMX = RHO
     endif
     ! CERCA MASSIMO DI H
     if( EPSH > AHMAX  ) then
        AHMAX = EPSH
        NHMAX = JF
        ROHMX = RHO
     endif
     ! CERCA MASSIMO DEI NEUTRINI 
     if( EPNEU > ANUMAX ) then
        ANUMAX = EPNEU
        NNUMAX = JF
     endif
     ! CERCA MASSIMO GRAVITAZIONALE 
     if( abs(GRAVI) > abs(AGRMAX) ) then
        AGRMAX = GRAVI
        NGRMAX = JF
     endif
     ! CERCA TEMPERATURA MASSIMA  
     if( TT > TMAX ) then
        TMAX    = TT
        NTMAX   = JF
        AA(1) = EM / EMTOT
        AA(2) = log10(PP) + 17.
        AA(3) = log10(TT) +  6.
        AA(4) = log10(RHO)
     endif

     ! CALCOLA QUANTITA' AL BORDO DEL CORE CONVETTIVO
     if( .not. (JF /= KCORE .or. ICOV /= 1 .or. KCORE == MAXMEin )) then
        BB(1) = EM / EMTOT
        BB(2) = RR / RTOT
        BB(3) = ZL / G(2, MAXMEin)
        BB(4) = log10(PP) + 17.
        BB(5) = log10(TT) +  6.
        BB(6) = log10(RHO)
        BB(7) = XXX(1, JF + 1)
        BB(8) = XXX(2, JF + 1)
        BB(9) = XXX(3, JF + 1)
        BB(10) = XXX(4, JF + 1)
        BB(11) = XXX(5, JF + 1)
        BB(12) = XXX(6, JF + 1)
     endif

     ! CALCOLA GRANDEZZE RELATIVE A BORDO INTERNO SHELL CONVETTIVE
     ! OPPURE A CONVEZIONE SUPERFICIALE
     if( .not. (L1 == 0 .or. JF /= IMINT( L )) ) then
        CC(1, L) = EM / EMTOT
        CC(2, L) = RR / RTOT
        CC(3, L) = ZL / G(2, MAXMEin)
        CC(4, L) = log10(PP) + 17.
        CC(5, L) = log10(TT) +  6.
        CC(6, L) = log10(RHO)
        JFA = JF
        if( IENV == 1 .and. JF == IMINT(L1) ) JFA = MAXMEin-1
        CC(7, L) = XXX(1, JFA)
        CC(8, L) = XXX(3, JFA)
        CC(9, L) = XXX(4, JFA)
        CC(10, L) = XXX(5, JFA)
        CC(11, L) = XXX(6, JFA)
        L = L + 1
     endif
     ! CALCOLA GRANDEZZE RELATIVE A BORDO ESTERNO SHELL CONVETTIVE
     if( .not. (L1 == 0 .or. JF /= IMEST(N)) ) then
        DD(1, N) = EM / EMTOT
        DD(2, N) = RR / RTOT
        DD(3, N) = ZL / G(2, MAXMEin)
        DD(4, N) = log10(PP) + 17.
        DD(5, N) = log10(TT) +  6.
        DD(6, N) = log10(RHO)
        N = N + 1
     endif

     ! controllo per core He, in cui calcolare abbondanze medie CNO
     ! Il controllo e' ripreso da sotto.
     if(JF > 1) then
        if(XXX(1,JF-1) <= 0. .and. XXX(1,JF) > 0. .and. &
          logCNO == 1) call mix_cno(jf-1,MAXMEin)
     endif
 
     ! produzione nucleare nelle zone
     If(xxx(3,jf) == 0) then
        lumin%pn(1) = lumin%pn(1) + epsn
     else if( xxx(1,jf) == 0) then
        lumin%pn(2) = lumin%pn(2) + epsn
     else
        lumin%pn(3) = lumin%pn(3) + epsn
     endif
     

     ! CONTROLLO DELLA CIACIO PER NSMORZA
     if( L1 > 0 ) then
        if( IENV == 1 .and. JF == IMINT(L1) ) then
           nsmorza = JF
        endif
     endif

     ! controllo la produzione di e. nuc. alla discontinuita'
     ! lasciata dalla convezione superficiale
     if(jf > 1) then
        if(fase <= fase_exH .and. G(5, jf)/ EMTOT >= lumin%M_bordo_conv .and. &
             G(5, jf-1)/ EMTOT < lumin%M_bordo_conv) then
           lumin%nuclear = EPSN
        endif
     endif

     ! STAMPA VARIE QUANTITA' SE NABLA=3 OGNI KSD MESHPOINTS
     if(NABLA == 3) then
        ! scrive FISICA.DAT
        if( ICOV == 1 .and. JF == 1 ) write(ioFISICA,101)
        if( ICOV == 1 .and. JF == KCORE ) write(ioFISICA,102)
        if( .not. (L1 == 0 .or. L1 == 1 .and. IENV == 1) ) then
           if( .not. (IENV == 1 .and. JF == IMINT(L1)) ) then
              if( JF == IMINT(L3) ) write(ioFISICA,103)
              if( JF == IMINT(L3) ) L3 = L3 + 1
              if( JF == IMEST(N3) ) write(ioFISICA,104)
              if( JF == IMEST(N3) ) N3 = N3 + 1
           endif
        endif
        if( JF /= 1 ) then
           if( XXX(3,JF-1) <= 0. .and. XXX(3,JF) > 0. ) write(ioFISICA,105)
           if( XXX(1,JF-1) <= 0. .and. XXX(1,JF) > 0. ) write(ioFISICA,106)
        endif

        if(L1 > 0) then
           if( IENV == 1 .and. JF == IMINT( L1 )) then
              write(ioFISICA,107)
              if(fase <= fase_exH .and. (G(5, jf)/ EMTOT < &
                   lumin%M_bordo_conv .or. neverset == 1)) then
                 ! trovo il minimo (in M) del bordo convezione sup.
                 ! nel corso dell'evoluzione
                 lumin%M_bordo_conv = G(5, jf)/ EMTOT
                 neverset = 0
              endif
              
           endif
        endif

        JRONA = 0
        JRANA = JF - ( JF / KSD ) * KSD
        if( JF == 1 .or. JF == MAXMEin ) JRONA = 1
        if( JRANA /= 0 .and. JRONA == 0 ) cycle

        EE(1) = EM / EMTOT
        E1(JF) = EE(1)
        EE(2)  = RR / RTOT
        E2(JF) = EE(2)
        EE(3)  = log10(PP) + 17.
        E3(JF) = EE(3)
        EE(4)  = log10(TT) +  6.
        E4(JF) = EE(4)
        EE(5)  = log10(RHO)
        E5(JF) = EE (5)
        EE(6)  = ZL / G(2, MAXMEin)
        E6(JF) = EE(6)	
        EE(7)  = EPSN
        E7(JF) = EE(7)
        EE(8) = GRAVI
        EE(9) = EPNEU
        !     EE( 10 ) = DELTA
        EE(10) = GRSAD
        EE(11) = DRAD(JF)
        EE(12) = DAD
        EE(13) = CAP
        E13(JF) = EE(13)

        write(ioFISICA,109) EE , KECUIL(JF) , JF
     endif
  end do

  FG1 = GENPP  / G(2, MAXMEin)
  FG2 = GENCNO / G(2, MAXMEin)
  FG3 = GENHE  / G(2, MAXMEin)
  FG4 = GENGRA / G(2, MAXMEin)
  FG5 = GENNU  / G(2, MAXMEin)
  FG6 = GENCAR / G(2, MAXMEin)

  ! CALCOLA QUANTITA' AL CENTRO
  CV(1) = log10(G(3, 1)) + 17.
  CV(2) = log10(G(4, 1)) +  6.
  CV(3) = log10(ROC)
  CV(4) = XXX(1, 1)
  CV(5) = XXX(2, 1)
  CV(6) = XXX(3, 1)
  CV(7) = XXX(4, 1)
  CV(8) = XXX(5, 1)
  CV(9) = XXX(6, 1)

  ! inserisco overshooting
  if(kover == 1 .and. printover == 1) then
     XFE1 = XME*DEFAUFe
     dfe_fe = (XXX(21,1)-xfe1)/xfe1
     xzora = (1.0+dfe_fe)*XME
     CHT = CV(4)+CV(5)+CV(6)+xzora
     write(85,8453) NMD,AGE,CV(4),CV(5),CV(6),CV(7),CV(8),CV(9),CHT
  endif

  ! CALCOLA GRANDEZZE RELATIVE ALLA SHELL DI IDROGENO 
  if( XXX(1, 1) <= 0. ) then
     XHA = 0.01 * XXX(1, MAXMEin)
     XHB = 0.90 * XXX(1, MAXMEin)
     do K = 2,MAXMEin
        if( XXX(1, K) > 0. ) exit
     end do
     if( K < MAXMEin ) then
        CORHE = G(5, K) / EMTOT * 100.
        do LMA = 2,MAXMEin
           if( XXX(1, LMA) >= XHA ) exit
        end do
        if( LMA > MAXMEin ) LMA = MAXMEin
        do LMB = LMA,MAXMEin
           if( XXX(1, LMB) >= XHB ) exit
        end do
        if( LMB > MAXMEin ) LMB = MAXMEin
        
        FF(1) = G(5, NHMAX) / EMTOT * EMSOL
        FF(2) = AHMAX
        FF(3) = G(1, NHMAX) / RTOT
        FF(4) = G(2, NHMAX) / G(2, MAXMEin)
        FF(5) = log10(G(3, NHMAX)) + 17.
        FF(6) = log10(G(4, NHMAX)) +  6.
        FF(7) = log10(ROHMX)
        FF(8) = (G(5, LMB) - G(5, LMA)) / EMTOT * EMSOL
        KSHHY = 1
     endif

     ! CALCOLA GRANDEZZE RELATIVE ALLA SHELL DI ELIO 
     if( XXX(3, 1) <= 0. ) then
        XEA = 0.01
        XEB = 0.90
        do K = 2,MAXMEin
           if( XXX(3, K) > 0. ) exit
        end do
        if( K > MAXMEin ) K = MAXMEin
        CORCO = G(5, K) / EMTOT * 100.
        do LEA = 2,MAXMEin
           if( XXX(3, LEA) >= XEA ) exit
        end do
        if( LEA > MAXMEin ) LEA = MAXMEin
        do LEB = LEA,MAXMEin
           if( XXX(3, LEB) >= XEB ) exit
        end do
        if( LEB > MAXMEin ) LEB = MAXMEin
        HH(1) = G(5, NHEMAX) / EMTOT * EMSOL
        HH(2) = AHEMAX
        HH(3) = G(1, NHEMAX) / RTOT
        HH(4) = G(2, NHEMAX) / G( 2 , MAXMEin )
        HH(5) = log10(G(3, NHEMAX)) + 17.
        HH(6) = log10(G(4, NHEMAX)) +  6.
        HH(7) = log10(ROHEMX)
        HH(8) = (G(5, LEB) - G(5, LEA)) / EMTOT * EMSOL
        KSHHE = 1
     endif
  endif

  ! CHIMICA E VARIAZIONI TEMPORALI FISICA
  if( NABLA == 3 ) then
     ! chimica
     write(ioCHIMICA,'("#   MOD   AGE       MESH")')
     write(ioCHIMICA,'("# ", i5,1x,f14.10,1x,i5)') nmd, age, maxmein-1
     write(ioCHIMICA,111)
  endif

  do K = 2,MAXMEin
     do J = 1,MAXMV
        if(GG(5, J) >= G(5, K)) exit
     end do
     if( J > MAXMV ) J = MAXMV
     RAT = (G(5, K) - GG(5, J-1)) / (GG(5, J) - GG(5, J-1))
     do I = 1,4
        FISV = GG(I, J-1) + RAT * (GG(I, J) - GG(I, J-1))
        if( abs(FISV) >= 1.d-14 ) then
           DFIS(I) = (G(I, K) - FISV ) / FISV * 1.d2
        endif
     end do

     do I = 1,4
        if( K <= 2 ) then
           IVME(I) = K
           DFIX(I) = DFIS(I)
        endif
        if( abs(DFIS(I)) > abs(DFIX(I)) ) then
           DFIX(I) = DFIS(I)
           IVME(I) = K
        endif
     end do

     if( NABLA == 3 ) then
        JRONA = 0
        JRANA = K - ( K / KSD ) * KSD
        if( K == 2 .or. K == MAXMEin ) JRONA = 1
        if( JRANA /= 0 .and. JRONA == 0 ) cycle
        GM = (G(5, K-1) + G(5, K))/( 2. * EMTOT )
        KM = K - 1
        XI = XXX(1, KM)
        XHE3 = XXX(2, KM)
        XHE4 = XXX(3, KM)
        XC = XXX(4, KM)
        XN = XXX(5, KM)
        XO = XXX(6, KM)
        XO18 = XXX(7, KM)
        XNE20 = XXX(8, KM)
        XNE22 = XXX(9, KM)
        XMG24 = XXX(10, KM)
        XMG25 = XXX(11, KM)
        XMG26 = XXX(12, KM)
        XSI28 = XXX(13, KM)
        XNA23 = XXX(15, KM)
        XFE = XXX(21, KM)
        XLI6 = XXX(22, KM)
        XLI7 = XXX(23, KM)
        XBE = XXX(24, KM)
        XB = XXX(25, KM)
        XD = XXX(26, KM)

        if( XI <= 0. .and. XHE4 <= 0. ) then
           XI = XSERV(1, KM)
           XHE4 = XSERV(2, KM)
        endif
        ! scrittura CHIMICA.DAT
        write(ioCHIMICA,112) KM, Gm, XI,XHE3,XHE4,XC,XN,XO,XO18,XFE, &
             XLI6,XLI7,XBE,XB,XD,XSI28,XMG26,XNE22
     endif
  end do

  ! STAMPA FISICA VARIA
  if( NABLA == 3 ) call PLOTTA(MAXMEin, RTOT)
  write(2,113) CV
  if( ICOV == 1 .and. KCORE /= MAXMEin ) write(2, 114) BB
  if( NTMAX > 1 ) write(2, 115) AA
  if( KSHHE == 1 ) write(2, 116) HH
  if( N1 /= 0 ) then
     do LL = 1,N1
        write(2, 117 ) (CC(I, LL), I=1,11)
        write(2, 118 ) (DD(I, LL), I=1,6)
     end do
  endif
  if( KSHHY == 1 ) write(2, 119) FF
  if( IENV  == 1 ) write(2, 120) (CC(I, L1), I=1,11)
  write(2,121) GENPP,FG1,GENCNO,FG2,GENHE,FG3,GENGRA,FG4,GENNU,FG5, &
       GENCAR,FG6
  write(2,122) (IVME(LL), DFIX(LL), LL=1,4)

  write(*,200) A(2),A(3),A(4),CV(1),CV(2),CV(3),CORCO,CORHE,NMOD

  if( XXX(1,1)  >  0. ) then
     XCEN = XXX(1,1)
     VARI = DELCE
     write(*,201)BB(1),XCEN,VARI,FG1,FG2,FG3,FG4,EMSOL,MAXMEin
  else if( XXX(3,1)  >  0. ) then
     XCEN = XXX(3,1)
     VARI = BB(1) - GICO / EMTOT
     write(*,202) BB(1),XCEN,VARI,FG1,FG2,FG3,FG4,EMSOL,MAXMEin
  else
     XCEN = XXX(4,1)
     VARI = TMAX
     write(*,203) BB(1),XCEN,VARI,FG6,FG2,FG3,FG4,EMSOL,MAXMEin
  endif

  write(*,122) (IVME(LL), DFIX(LL), LL=1,4)
  dnorm = 2.8123d27
  GA_tot = 0.
  CL_TOT = 0.
  do k=1,8
     flussi(k) = ppneu(K)/dnorm
     GA_seg(k) = GAMATRIX(K)*flussi(K)
     GA_tot = GA_tot+GA_seg(K)
     CL_seg(K) = CLMATRIX(K)*flussi(K)
     CL_TOT = CL_TOT+CL_seg(K)
  end do

  ! MEMORIZZA VARIE QUANTITA' PER GRAFICI
  JRANA = NMOD - (NMOD / IFVA) * IFVA
  if( JRANA == 0 ) then
     if( L1 < 1 ) L1 = 1
     RPQ(1) = XXX(1, 1)
     RPQ(2) = XXX(3, 1)
     RPQ(3) = XXX(4, 1)
     RPQ(4) = XXX(5, 1)
     RPQ(5) = XXX(6, 1)
     RPQ(6) = FG1
     RPQ(7) = FG2
     RPQ(8) = FG3
     RPQ(9) = FG4
     lumin%Lpp = RPQ(6)
     lumin%Lcno = RPQ(7)
     lumin%L3a = RPQ(8)
     lumin%Lgra = RPQ(9)
     lumin%L = A(2)
     lumin%teff = A(3)

     SUMMA = 0.
     do JJ = 1,MELE
        ES(JJ) = XXX(JJ, 1)
        SUMMA = SUMMA + XXX(JJ, 1)
     end do

     if(ES(1) <= 0. .and. ES(3) <= 0.) then
        ES(1) = XSERV(1, 1)
        ES(3) = XSERV(2, 1)
        SUMMA = SUMMA + ES(1) + ES(2)
     endif
     write(10,300) NMOD,AGE,(A(J),J=1,4),A(6),CV(1),CV(2),CV(3),RPQ
     write(10,301) AA,FF,HH,BB(1),(BB(I),I=3,7),(BB(K),K=9,12), &
          CC(1,L1),(CC(K,L1),K=3,11)

     ! Scrittura OUT.DAT e BIGTAB.DAT 
     if(rpq(1) >= 1.e-29) then
        HorHe = rpq(1)
     else
        HorHe = rpq(2)
     endif

     ! questa parte sistema l'envelop convettivo
     if(kcore == MAxmein .or. (kcore == 0 .and. cc(1,l1) == 0.0) ) then
        ! struttura tutta convettiva
        envconv = a(6)
     else if( cc(1,l1) == 0.0 ) then
        envconv = 0.0
     else if (KCORE /= MAXMEin) then
        envconv = a(6) - cc(1,l1)*a(6)
     endif
     
     feh = log10(XXX(21,MAXMEin-1)/(XXX(1,MAXMEin-1)*55.847))

     ! [Fe/H] Asplund 2009
     write(ioOUT,307) NMOD,AGE,HorHe,a(2),a(3),a(6),rpq(9),rpq(8), &
          feh, feh+4.5, RAG, LOGG, sqrt(a(6)/rag**3),           &
          a(6)/(rag**2 * sqrt(10**a(3)/5777)), ALFA   ! 5777 = Teff_sun

     write(ioBIGTAB,304) NMOD,AGE,HorHe,a(2),a(3),CV(2),CV(3),aa(3),&
          aa(1)*a(6), &
          bb(1)*a(6),ff(1),envconv,rpq(6),rpq(7),rpq(8),rpq(9),CC(8,L1),&
          CC(5,L1),GICO, RAG, LOGG

     write(39,'(i6,3e14.5)') nmd, lumin%pn
  endif

  return

100 format(/,' AGE:',F15.11,' DT:',F8.4,' L:',F8.4,' TE:',F8.4,' R:',F8.4, &
       ' GR:',1P,e10.2,' M:',0P,F9.5,' E:',F6.3,' DM:',1P,e10.3,I2, &
       2I5,'  (ATMOSFERA ',a8)
101 format('#      ********** CORE CONVETTIVO **********')
102 format('#      ********** BORDO CONVETTIVO CENTRALE **********')
103 format('#      ********** BORDO INTERNO SHELL CONVETTIVA **********')
104 format('#      ********** BORDO ESTERNO SHELL CONVETTIVA **********')
105 format('#      ********** CORE DI CARBONIO - OSSIGENO **********')
106 format('#      ********** CORE DI ELIO **********')
107 format('#      ********** BORDO CONVEZIONE SUPERFICIALE **********')
108 format('#  M/MTOT      R/RTOT  log(P)  log(T)  log(RO)   L/LSUP&
       &     ENUC      EGRAV     ENEU      DELTA      RADIAT     ADIAB  &
       &    OPACITY   MESH')
109 format(F12.9,F10.7,F8.4,F7.4,F8.4,1P,e13.5,1P,3e10.3,1P,e9.2,1P,e12.4, &
       1P,e12.4,1P,e9.2,I2,I5)
111 format('#  MESH  M/Mtot     H         He3       He4       C12       N14 &
       &      O16       O18       Fe56      Li6       Li7       Be        B   &
       &      D         Si28      Mg26      Ne22')
112 format(1p,2x,I5,17e10.3)
113 format(1X,'VALORI CENTRALI:',3F8.4,1P,6e12.5)
114 format(1X,'CORE CONVETTIVO:',F10.7,F9.6,1P,e10.2,0P,3F7.3,1P,6e10.3)
115 format(1X,'TEMPERATURA MASSIMA NON CENTRALE:',F12.9,3F8.4)
116 format(1X,'HE-SHELL:',F10.7,1P,3e10.2,0P,3F7.3,1P,e10.3)
117 format(1X,'BORDO INT. SH. CONVET.:',2F10.7,1P,e10.2,0P,3F7.3,1P,5e10.3)
118 format(1X,'BORDO ESTERNO SHELL CONVETTIVA:',2F10.7,1P,e10.2,0P,3F7.3)
119 format(1X,' H-SHELL:',F10.7,1P,3e10.2,0P,3F7.3,1P,e10.3)
120 format(1X,'BORDO ENVELOPE CONVETTIVO:',2F10.7,1P,e10.2,0P,3F7.3,1P,6e10.3)
121 format(1X,1P,12e10.2)
122 format(1X,' DR(',I4,') =',F8.4, &
       ' DL(',I4,') =',F11.4, &
       ' DP(',I4,') =',F8.4, &
       ' DT(',I4,') =',F8.4)
200 format(' L:',F7.4,'  TE:',F7.4,'  R:',F7.4, &
       '  Pc:',F8.4,'  Tc:',F7.4,'  rho_c:',F7.3, &
       '  HE_SH:',F8.4,'  H_SH:',F8.4,'  NMOD:',I5) 
201 format(F7.3,'  Hc:',1P,e9.2,e9.2, &
       '  Lpp:',e9.2,'  Lcno:',e9.2,'  L3a:',e9.2, &
       '  Lgra:',e9.2,'  M:',0P,F9.5,'  NMESH:',I5)
202 format(F7.3,'  HEc:',1P,e9.2,e9.2, &
       '  Lpp:',e9.2,'  Lcno:',e9.2,'  L3a:',e9.2, &
       '  Lgra:',e9.2,'  M:',0P,F9.5,'  NMESH:',I5)
203 format(F7.3,'  Cc:',1P,e9.2,e9.2, &
       '  Lpp:',e9.2,'  Lcno:',e9.2,'  L3a:',e9.2, &
       '  Lgra:',e9.2,'  M:',0P,F9.5,'  NMESH:',I5)
300 format(I5,1P,e22.15,3E16.8,/,5E16.8,/,5E16.8,/,4E16.8)
301 format(1P,5E16.8)
302 format(2I5,1P,e22.15,3E16.8,/,2E16.8)
1981 format(i5,f12.8,12E14.6) 
123 format(1P,5e16.9)
124 format(I5)

304 format(I5,1X,F11.8,e12.5,6(1X,F8.5),2(1X,F8.5),3(1X,F10.3), &
       2(1X,E10.3),3(1X,F7.4),1X,F10.4,1X,F7.4)
307 format(I5,1x,F13.10,e12.5,1X,F13.10,1X,F13.10,1X,F8.5,2(1X,e10.3),1x,&
         9(f10.5,1x))
8453 format(I5,8e13.5)
end subroutine STAMPA


! questa routine calcola la media di CNO nel core di He
subroutine mix_cno(j, mx)
  use interfaccia
  use fisica
  use chimic
  use nummod
  use varie
  use mistura

  implicit none

  integer :: j, mx
  
  real :: dm , dmt, dmte, z, ze, xfe
  real,dimension(4) :: met
  real,dimension(5) :: mete
  integer :: i, k

  xfe = XME*DEFAUFe
  
  dmt = 0.
  met = 0.
  dmte = 0.
  mete = 0.

  ! calcolo media di C,N,O,Fe e z nel core di He
  do i=1,j
     dm = G(5,i+1) - G(5,i)
     do k=4,6
        met(k-3) = met(k-3) + xxx(k,i) * dm
     end do
     ! Fe = 21
     met(4) = met(4) + xxx(21,i) * dm
     dmt = dmt + dm
  end do
  
  met = met/dmt
  ! calcolo z  
  z = (met(4)/xfe)*XME
  
  ! calcolo media di C,N,O,Fe fuori dal core di He
  do i=j+1,mx-1
     dm = G(5,i+1) - G(5,i)
     do k=4,6
        mete(k-3) = mete(k-3) + xxx(k,i) * dm
     end do
     mete(4) = mete(4) + xxx(21,i) * dm
     mete(5) = mete(5) + xxx(2,i) * dm
     dmte = dmte + dm
  end do

  mete = mete/dmte
  ! calcolo z  
  ze = (mete(4)/xfe)*XME

  write(18,'(i5,1x,11(e13.6,1x))') nmd, met/z, z, mete/ze, ze
 
end subroutine mix_cno
