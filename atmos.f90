subroutine ATMOS(FLUMIN,TEFF,P,T,R,EMASS,ISCR,IRA) 
  use interfaccia
  use mistura
  use intero
  use fisica
  use preco
  use mesh
  use atmosfere
  use chimic
  use chim
  use numer
  use varie
  use nummod
  use costanti
  implicit none

  real :: FLUMIN, TEFF, P, T, R, EMASS
  integer :: ISCR, IRA

  real :: GIINI, TE 
  ! controlla save
  real,save :: C1(4), C2(4), C3(4), GA(5,300), GAV(5,300) 
  !   ATTENZIONE ALLA CHIMICA NEL CASO DI DREDGE-UP                       
  !   CLASSICAL VALUES : DTAU=.3 - PASSP=.3                                 
  real,parameter :: DTAU = .1, PASSP = .1, RATTO = .5

  ! disabilito il calcolo delle quantita' non volute in state
  real, parameter :: dum = 1.d99

  integer :: liala, jme, k, n_atm, num, nst, jjj, iter
  integer :: km, j, kk, maxa, lxa, lxo, ia
  real :: emassa, em, y, el, r7, rcm, rsol, a3, pmas, erre, elnuc
  real :: elnac, xfe, dfe_fe, tau_racc, pre, prad, ro, adiab, cp, pmu
  real :: cap, radiat, tau, tauv, gi, dtaup, dep, der, dem
  real :: ppp, tm, hr2, apm, arag, remma, emme, acca, passo, grad
  real :: dp, passu, ps, pmus, rs, ems, ts, grsad, rip, dr, dm
  real :: passt, passr, epsn, epsa, delta, hem, hp, ht, hro, hr
  real :: hgi, hcap, fh, hmas, hrag, zam, zem, zim, gradv, dt, t6, ecar
  character(len=15) :: caller="atmos         "
  logical :: vernazza=.false.
  

  EMASSA = EMASS / Msun  
  EM = 1.d3 * EMASS
  Y = 1.-XH-XME 
  if(Y < 0.) Y = 0. 

  EL = 1.d3 * (10.**FLUMIN) * Lsun 
  TE = 10.**TEFF 

  R = (1.d5 / sqrt(4.d0 * pigre * ssb) ) * (sqrt(EL)) / (TE * TE) 
  R7 = R / Rsun
  RCM = R * 1.d10 
  RSOL = R7 
  A3 = RCM 
  if(ISUB == 1) write(2,222) EMASSA,FLUMIN,TEFF,Y,XME,R7,BMAG,ALFA 
  PMAS = 0. 
  ERRE = 0.
  LIALA = 0 
  JME = 0 
  ELNUC = 0. 
  ELNAC = 0.
  
  XX(1:mele) = XXX(1:mele,MAXME-1) 
  

  ! PREPARO LE QUANTITA' CHE SERVONO PER IL CALCOLO DELLA NUOVA ATMOS
  X_sup = XXX(1,MAXME-1) 
  Y_sup = XXX(3,MAXME-1) 
  XFE = XME*DEFAUFe 
  ! variazione relativa di Fe(diffusione)
  DFE_FE = (XXX(21,MAXME-1)-XFE)/XFE 
  ! metallicita' riscalata col Ferro        
  Z_SUP = (1.+DFE_FE)*XME 

  if( .not. Z_SUP >= 0) then
     write(*,*) "160 - atmos: errore nel valore di Z_sup"
     write(66,*) "160 - atmos: errore nel valore di Z_sup"
     stop
  endif

  ! vale 1 se trovo valori nelle tabelle di Brott & Haus
  N_BH05 = 0 
  ! vale 1 se trovo valori nelle tabelle di Castelli & K
  N_CK03 = 0 

  GIINI = 1.d10 * Ggrav * EM / (R * R) 

  ! definisco il punto di raccordo atmosfera/sub-atmosfera        
  ! (tau = 10)                                      
  TAU_RACC = 10.0 

  ! cerca i valori nelle tabelle di Brott & Hauschildt (2005)        
  ! se non li trova, le cerca nelle tabelle Castelli & Kurucz (2003) 

  ! Brott & Haushildt (2005)
  if(do_bh_05 == 1 .or. nmd < 3) call ATM_BH05(GIINI,TE,TAU_RACC) 

  if(do_ck_03 == 1 .and. N_BH05 == 0) then 
     ! Castelli & Kurucz (2003)     
     call ATM_CK03(GIINI,TE,TAU_RACC) 
  end if

  N_ATM = 1

  ! controllo se ho trovato TE, Log G e [Z/X] nelle tabelle          
  ! Se N_ATM = 0 escludo atmosfere
  if(N_BH05 == 0 .and. N_CK03 == 0) N_ATM = 0 

  ! se uso i modelli di atmosfera allora specifico subito T, P alla base.
  if(N_ATM == 1) then 
     T = V_ATM(2) 
     ! pgas                                           
     PRE = V_ATM(3) 
     ! pradiativa                             
     PRAD = arad_3 * (T**4) 
     P = PRE + PRAD 
     call STATE(caller, PRE,T,RO,ADIAB,CP,PMU) 
     call KAPPA(RO,T,CAP,XX(3)) 
     ! gradiente radiativo
     !RADIAT = 2.632d12*CAP*EL*P/(R*R*GIINI*(T**4))
     RADIAT = cte_rad * CAP * EL * P / (R * R * GIINI * (T**4))
     ! per la stampa su print:                                          
     NUM = 1 
     PMAS = 0.
     TAU = TAU_RACC 
  else
     ! non dovrebbe più usare la vecchia T(tau) a meno che non ci siano
     ! dei punti non in tabella...                                      

     ! se vado avanti vuol dire che sono al primo modello
     ! CALCOLO CONDIZIONE AL CONTORNO: P0 PER TAU0
     TAU = 0.
     PRAD = arad_3 * (TE**4) 
     PRE = 0.1 * arad_3 * (.75*TE)**4 
     TAUV = TAU 
     GI = GIINI 
     NUM = 0 

     do
        NUM = NUM+1 
        if(NUM > 200) then
           write(*,*) '161 - atmos: troppe iterazioni KS'
           write(66,*) '161 - atmos: troppe iterazioni KS'
           stop
        endif
        ! calcolo la temperatura per il tau voluto:                         
        if(vernazza) then
           T = 3./4.*TE**4*(TAU+1.017-.3*exp(-2.54*TAU)-.291*exp(-30.*TAU))
        else
       ! KS originale                                                                           
           T = 3./4.*TE**4*(TAU+1.39-.815*exp(-2.54*TAU)-.025*exp(-30.*TAU))
        endif
        T = T**.25 
        if(isnan(T) .eqv. .true.) then
           write(*,*) '162 - atmos: T isnan'
           write(66,*) '162 - atmos: T isnsn'
           stop
        endif
        PRAD = 2.52d-15*(T**4) 
        P = PRE+PRAD 

        call STATE(caller,PRE,T,RO,dum,dum,dum) 
        call KAPPA(RO,T,CAP,XX(3)) 
        !RADIAT = 2.632d12*CAP*EL*P/(R*R*GIINI*(T**4))
        RADIAT = cte_rad * CAP * EL * P / (R * R * GIINI * (T**4)) 
        GI = GIINI*(1.-4.*RADIAT*PRAD/P) 
        TAU = PRE*CAP/GI 

        if(ISCR == 1) then 
           write(2,421)NUM,TAU,PRE,PRAD,T,GI,RADIAT,RO,CAP 
        endif

        if(abs((TAUV-TAU)/TAU)-1.d-4 <= 0 ) exit
        TAUV = TAU
     end do

     ! ==================================================
     !              INIZIO CALCOLO ATMOSFERA
     ! ==================================================
     if(iscr == 1) then
        write(2,423)PRE,TAU
     endif
     call STATE(caller,PRE,T,RO,dum,dum,dum)
     ERRE = P/(GI*RO)
     !PMAS = 12.566*A3*A3*RO*ERRE
     PMAS = 4.d0 * pigre * A3 * A3 * RO * ERRE
     if(TAU <= 1.d-8) TAU = 1.d-8
     NUM = 0

          tau_racc = 0.312  ! KS da MESA
     !     tau_racc = 2./3.  ! KS      
     if(vernazza) tau_racc = 2./3.
     
     do 
        NUM = NUM+1
         if(num > 1000) then
           write(*,*) '161 - atmos: troppe iterazioni KS'
           write(66,*) '161 - atmos: troppe iterazioni KS'
           stop
        endif
        DTAUP = DTAU*TAU
        
        if(TAU + DTAUP > tau_racc)  DTAUP = tau_racc - TAU
        
        TAU = TAU+DTAUP
        ! calcolo la temperatura per il tau voluto:
        ! Vernazza 1981
        if(vernazza) then
           T = 3./4.*TE**4*(TAU+1.017-.3*exp(-2.54*TAU)-.291*exp(-30.*TAU))
        else
       ! KS originale         
           T = 3./4.*TE**4*(TAU+1.39-.815*exp(-2.54*TAU)-.025*exp(-30.*TAU))
        endif
        T = T**.25

        PRAD = arad_3 * (T**4)
        DEP = GI*DTAUP/CAP
        PRE = PRE+DEP
        P = PRE+PRAD
        call STATE(caller,PRE,T,RO,ADIAB,dum,dum)
        DER = DEP/(GI*RO)
        !DEM = 12.566*A3*A3*RO*DER
        DEM = 4.d0 * pigre * A3 * A3 * RO * DER
        PMAS = PMAS+DEM
        ERRE = ERRE+DER
        call KAPPA(RO,T,CAP,XX(3))
        !RADIAT = 2.632d12*CAP*EL*P/(R*R*GIINI*T**4)
        RADIAT = cte_rad * CAP * EL * P / (R * R * GIINI * (T**4))
        if(iscr == 1) then
           ppp = pre+prad
           tm = T
           hr2 = R/(RSOL*Rsun)
           write(2,422) NUM,TAU,PRE,PRAD,T,RO,CAP,GI,RADIAT,ADIAB,PMAS*1d-33
        endif
        GI = GIINI*(1.-4.*RADIAT*PRAD/P)

        if(TAU >= TAU_RACC .or. RADIAT >= ADIAB) exit
     end do
     ! FINE VECCHIA ATMOS 
  endif

  ! ==================================================
  !           INIZIO CALCOLO SUBATMOSFERA
  ! ==================================================
  call STATE(caller,PRE,T,RO,ADIAB,CP,PMU)
  if(ISCR /= 0) then
     if(PMAS > 0.) APM = log10(PMAS)
     if(ERRE > 0.) ARAG = log10(ERRE)
     if(PMAS > 0. .and. ERRE > 0.) write(2,888) TAU,T,P,CAP,NUM,APM,ARAG
  endif
  T = T*1.d-6
  P = P*1.d-15
  REMMA = EM
  EMME = EM*FRAZ

56 continue
  NST = 1
  ACCA = 0.
  JJJ = 0
  ITER = 0
  PASSO = PASSP/10.
  GRAD = RADIAT
  GRADV = GRAD
  if(RADIAT > ADIAB) GRAD = ADIAB
  if(ISCR /= 0) then
     ! scrivo SUBATM.DAT
     write(23,'("#   MOD   AGE       MESH")')
     write(23,'("# ", i5,1x,f14.10,1x,i5)') nmd, 0., 0
     write(23,999)
  endif

  do while(JJJ /= 2)
     if(LIALA == 1) then
        JME = JME+1
        GA(1,JME) = R
        GA(2,JME) = EL/100.
        GA(3,JME) = P/100.
        GA(4,JME) = T
        GA(5,JME) = EM*1.d-3
     endif

     if(LIALA == 2) then
        JME = JME+1
        GAV(1,JME) = R
        GAV(2,JME) = EL/100.
        GAV(3,JME) = P/100.
        GAV(4,JME) = T
        GAV(5,JME) = EM*1.d-3
     endif

     DP = PASSO*P/2.
     PASSU = PASSO
     PS = P
     PMUS = PMU
     RS = R
     EMS = EM
     TS = T
     ITER = ITER+1
     do K=1,3
        call LOCAL(PS,TS,RS,EL,EMS,RO,GRAD,CP,PMU,ACCA, &
             CAP,RADIAT,ADIAB,GI,XX(3),GRSAD)
        !C1(K) = -DP*1.d5*RS*RS/(666.8*EMS*RO)
        !C2(K) = -DP*1884.57868*(RS**4)/EMS
        !C3(K) = DP*TS*GRAD/PS

        C1(K) = - 1.d-5 * DP * RS * RS / (Ggrav * EMS * RO)
        C2(K) = - (4.d-5 * pigre / Ggrav) * DP * (RS**4) / EMS
        C3(K) = DP*TS*GRAD/PS
        if(K-2 < 0) then 
           RIP = .5
        else
           RIP = 1
        endif
        PS = P+DP*RIP
        RS = R+C1(K)*RIP
        EMS = EM+C2(K)*RIP
        TS = T+C3(K)*RIP
     end do
     DR = (C1(1)+2.*C1(2)+C1(3))/4.
     DM = (C2(1)+2.*C2(2)+C2(3))/4.
     DT = (C3(1)+2.*C3(2)+C3(3))/4.
     PMAS = PMAS-DM*1.d30
     EM = EM+DM

     if(EM-EMME < 0) then
        JJJ = 1
        if(abs((EM-EMME)/EMME)-1.d-8 <= 0) then
           JJJ = 2
           T = T+DT
        else
           EM = EM-DM
           ITER = ITER-1
           PASSO  = PASSO*(EMME-EM)/DM
           cycle
        endif
     else if(EM-EMME == 0) then
        JJJ = 2
        T = T+DT
     else
        T = T+DT
     endif

     R = R+DR
     P = P+DP
     if(JJJ <= 0) then
        PASST = PASSP
        PASSR = PASSP
        ! CLASSICAL VALUES: PASSO*.03  AND PASSP*.04 
        if(abs(GRADV-GRAD) > 0. .and. ITER > 1) then
           PASST = PASSO*.03/(abs(GRADV-GRAD))
        endif
        PASSR = abs(PASSP*.04*R/DR)
        PASSO = PASSP
        if(PASSO > PASST) PASSO = PASST
        if(PASSO > PASSR) PASSO = PASSR
        if(PASSO < 1.d-3) PASSO = 1.d-3
        if(LIALA == 1 .or. LIALA == 2) PASSO = PASSO/3.
     endif
     EPSN = 0.
     T6 = T*1.d6
     call EPSI(RO,T6,1,0,0.0,EPSN,EPSA,ecar,DELTA)
     ELNUC = ELNUC-EPSN*DM
     if(ELNUC > 0.) ELNAC = log10(ELNUC/EL)

     if((ELNUC/EL) > 1.d-4 .and. ISUB /= 1 .and. ISCR == 0 .and. &
          IRA == 0 .and. KSF == 1)   then
        IRA = 1
        return
     endif
     GRADV = GRAD
     ERRE = ERRE-DR*1.d10

     if(JJJ /= 2 .or. ISCR < 1) then
        if(ISCR == 0) then
           if(JJJ-1 <= 0) then
              cycle
           else
              exit
           endif
        endif
        if(ISUB == 1) NST = 1
     endif
     HEM = EM/REMMA
     HP = log10(P)+15.
     HT = log10(T)+6.
     HRO = log10(RO)
     HR = R/(RSOL*Rsun)
     if(GI > 0.) HGI = log10(GI)
     HCAP = log10(CAP)
     FH = log10(CP)
     HMAS = log10(PMAS)
     HRAG = log10(ERRE)

     ! SUBATM.DAT
     write(23,1111)HEM,HP,HT,HRO,HR,PMU,HCAP,RADIAT,ADIAB,GRSAD, &
          ACCA,FH,HGI,HMAS,HRAG,ELNAC,XX(1),XX(3),ITER

     NST = NST+1
     if(NST >= 10) NST = 1
  enddo

  if(.not. (LIALA >= 1 .or. IRA < 3)) then
     LIALA = 1
     if(IRA == 4) LIALA = 2
     EMME = EMASS*(1.-(1.-FRAZ)/RATTO)*1.d3
     goto 56
  endif
  P = (log10(P))+15.
  T = (log10(T))+6.
  R = R/(RSOL*Rsun)
  HRO = log10(RO)
  if(ISUB==1 .and. IBAT==0) write(*,444) EMASSA,FLUMIN,TEFF,Y,XME,P,T,R,HRO

  if(LIALA-1 < 0) then
     return
  else if(LIALA-1 == 0) then
     JME = JME-2
     do K=1,JME
        KM = MAXME+K
        IN(KM) = 0
        do J=1,MELE
           XXX(J,KM) = XXX(J,MAXME)
        end do
        do j=1,5
           KK = JME-K+1
           G(J,KM) = GA(J,KK)
        end do
     end do

     MAXA = MAXME
     MAXME = MAXME+JME
     G(5,MAXME) = REMMA*FRAZ*1.d-3
     ZAM = G(1,MAXA+1)-G(1,MAXA)
     ZEM = G(3,MAXA+1)-G(3,MAXA)
     ZIM = G(4,MAXA+1)-G(4,MAXA)
     LXA = MAXA+1
     LXO = LXA+2
     do K=LXO,MAXME
        KK=K-2
        G(1,KK) = G(1,K)-ZAM*(G(5,K)/G(5,MAXME)-1.)/(G(5,LXA)/G(5,MAXME)-1.)
        G(3,KK) = G(3,K)-ZEM*(G(5,K)/G(5,MAXME)-1.)/(G(5,LXA)/G(5,MAXME)-1.)
        G(4,KK) = G(4,K)-ZIM*(G(5,K)/G(5,MAXME)-1.)/(G(5,LXA)/G(5,MAXME)-1.)
        G(5,KK) = G(5,K)
     end do

     MAXME = MAXME-2
     IA = MAXME-JME+2
     return
  else 
     JME = JME-2
     do K=1,JME
        KM = MAXMIV+K
        do J=1,MELE
           XXV(J,KM) = XXV(J,MAXMIV)
        end do
        do j=1,5
           KK = JME-K+1
           GG(J,KM) = GAV(J,KK)
        end do
     end do

     MAXA = MAXMIV
     MAXMIV = MAXMIV+JME
     GG(5,MAXMIV) = REMMA*FRAZ*1.d-3
     ZAM = GG(1,MAXA+1)-GG(1,MAXA)
     ZEM = GG(3,MAXA+1)-GG(3,MAXA)
     ZIM = GG(4,MAXA+1)-GG(4,MAXA)
     LXA = MAXA+1
     LXO = LXA+2
     do K=LXO,MAXMIV
        KK = K-2
        GG(1,KK) = GG(1,K)-ZAM*(GG(5,K)/GG(5,MAXMIV)-1.)/ &
             (GG(5,LXA)/GG(5,MAXMIV)-1.)
        GG(3,KK) = GG(3,K)-ZEM*(GG(5,K)/GG(5,MAXMIV)-1.)/ &
             (GG(5,LXA)/GG(5,MAXMIV)-1.)
        GG(4,KK) = GG(4,K)-ZIM*(GG(5,K)/GG(5,MAXMIV)-1.)/ &
             (GG(5,LXA)/GG(5,MAXMIV)-1.)
        GG(5,KK) = GG(5,K)
     end do

     MAXMIV = MAXMIV-2
     IA = MAXMIV-JME+2
     write(2,200)((GG(J,K),J=1,5),K=IA,MAXMIV)
     return
  endif

100 format(1X,'CALCOLO P0 DI TAU0 NON CONVERGE') 
888 format(/,1X,'FINE ATMOSFERA  TAU ',1P,E11.4,2X,'T ',E11.4,2X,'P ',&
       E11.4,2X,'OPAC ',E11.4,2X,'N ',I4,2X,2F7.3,/)                     
999 format('# MASSA    LOG_P LOG_T  LOG_RO  RAGGIO   MU    OPACITA G_RADIAT&
       &  G_AD  G_CONV   VELOC      CP   LOG_G  LOG_PMS  LOG_ER LOG_ELN  &
       &H    He   MESH')                     
1111 format(F9.6,F7.3,F6.3,1x,F7.3,1x,F7.4,1X,F6.3,F7.3,2x,1P,E10.3,0P,2F6.3,&
       2x,1P,E10.3,0P,F7.3,F6.3,1x,3(F7.3,1x),1x,2F5.2,1x,I4)          
222 format(///,1X,'MASSA=',F8.4,'  L=',F7.4,'  TE=',F7.4,'  Y=',      &
       1P,E10.3,'  Z=',E10.3,'  R=',E10.3,'  B=',E10.3,'  A=',0P,F4.1,3X &
       ,1P,2E10.3,/)
421 format(I5,1P,8E9.2) 
423 format(' CONDIZIONE AL CONTORNO P0 E TAU0:',1P,2E12.4)
422 format(I5,1P,10E9.2)
444 format(1X,2F8.4,F7.4,1P,E8.1,E8.1,0P,F8.4,3F7.4)
200 format(3X,1P,5D14.5)

end subroutine ATMOS
