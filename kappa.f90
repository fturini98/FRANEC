subroutine KAPPA(RHO,TEMP,CAP,Y) 
  use interfaccia
  use mistura
  use e
  use indicefile
  use chim
  use varie
  use numer
  use nummod

  ! questa versione legge la tabella Opacita.AF.OPAL.POT
  ! In combustione di H questa kappa e' simile alle vecchie,se non per il 
  ! fatto che non si creano piu' le tabelle di puro Y,C ed O.             
  ! Quando si inizia a bruciare l'He la kappa richiama la subroutine OPAC 
  ! (file xcotrin2006.f) che interpola le tabelle di Opal con C e O       
  ! variabili. La OPAC viene chiamata 10 volte per 10 valori di z e poi si
  ! interpola in metallicita` con la SPLINE. Poi la kappa chiama la      
  ! subroutine  POTEKHIN (file condint2006.f) che calcola l'opacita'      
  ! conduttiva. Infine si sommano gli inversi.                            

  implicit none

  real :: RHO, TEMP, CAP, Y

  character(len=1),dimension(47) :: INTEST, COMPO
  character(len=1),dimension(25) :: TESTO 
  ! blocco cache
  real,save,dimension(320) :: xABGDb
  real,dimension(80) :: xtABGD
  real,save :: cache21, cache4, cacheY, cache6, cache1

  real ::  Z1s
  real,save,dimension(21) :: ABBNUM
  real,save,dimension(26) :: xmetallog

  real,dimension(4) ::  TT, B
  real,dimension(80) :: CAPR
  real,dimension(20) ::  CAPS
  real,save,dimension(139) :: TK
  real,save,dimension(10,29,139,10) :: OP
  real,save,dimension(10,10) ::  YKZM 
  real,dimension(10) :: CAPZ, OPAVETT, y2tmp, YY

  integer,save :: IZ = 0
  real,parameter :: UM = 2.302585
  real,parameter,dimension(10) ::  YKZ = (/0.,.05,.1,.2,.3,.5,.65,.8,.9,1./) 
  integer :: JY1

  ! pesi atomici                                           
  real,parameter,dimension(21) :: PESIAT = (/1.00790,4.00260,12.01100, &
       14.00670,15.99940,20.17900,22.98977,24.30500,26.98154,28.08550, &
       30.97376,32.06000,35.45300,39.94800,39.09830,40.08000,47.90000, &
       51.99600,54.93800,55.84700,58.70000/)

  ! le cariche dei 21 ioni                          
  real,parameter,dimension(21) :: zet = (/1.0,2.0,6.0,7.0,8.0,10.0, &
       11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,22.0, &
       24.0,25.0,26.0,28.0/)
  ! per pesi e cariche, nell'ordine:
  ! H, He, C, N, O, Ne,
  ! Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Ti, 
  ! Cr, Mn, Fe, Ni

  real,dimension(21) :: logzet

  
  ! Abb. frac. in massa AGS (i 19 metalli sono norm. a 1)     
  real,parameter,dimension(21) :: FRACMASS09 = (/1.,1.,0.176993,0.051846, &
       0.429024,0.094032,0.002187,0.052975,0.004163,0.054575,0.000436, &
       0.023139,0.000614,0.005494,0.000229,0.004801,0.000234,0.001243, &
       0.000810,0.096689,0.005334/)
  real,parameter,dimension(21) :: FRACMASS05 = (/1.,1.,0.177050,0.050682, &
       0.439157,0.083833,0.002042,0.049455,0.003799,0.054575,0.000426, &
       0.026575,0.000674,0.003631,0.000281,0.004915,0.000230,0.001363, &
       0.000809,0.094517,0.005985/)
  real,parameter,dimension(21) :: FRACMASS09a3 = (/1.,1.,0.106859,0.031302, &
       0.516814,0.113274,0.001321,0.063815,0.002514,0.059939,0.000263, &
       0.027874,0.000371,0.006618,0.000138,0.005783,0.000282,0.000750, &
       0.000489,0.058375,0.003220/)
  real,parameter,dimension(21) :: FRACMASSGS98 = (/1.,1.,0.171060,0.050108, &  
       0.465234,0.104344,0.002114,0.039743,0.003425,0.042860,0.000375, &   
       0.029480,0.000482,0.004316,0.000222,0.003949,0.000216,0.001046, & 
       0.000580,0.075957,0.004490/)
  real,parameter,dimension(21) :: FRACMASSGN93 = (/1.,1.,0.173356,0.053190, &
       0.482519,0.098700,0.002000,0.037599,0.003239,0.040546,0.000355, & 
       0.021159,0.000456,0.005381,0.000210,0.003735,0.000204,0.000989, &
       0.000548,0.071849,0.003964/)

  ! lo inizializzo con le fracmass corrette la prima volta che giro
  real,dimension(21) :: FRACMASS  

  real :: CAPRAD, CAPCOND, ABBIDRO, pmedio
  real,save,dimension(10) :: xmetal 

  integer :: jump, m, i, j, k, l, ierror, isup, iinf, kt, ii, nz
  integer :: jyy, jota, mi, n, icache, ic, iop, iconta
  real :: xfe, dfe_fe
  real,save :: xzora, xzora_l
  real :: tmp1, r, t, rlg, deltacm, deltaom, t6
  real :: alfa1, beta1, gamma1, delta1, rap, xmetx, deltacz, deltaoz
  real :: zz, condrel, cond, zkappa, deltac, deltao, ro, akapparad

  integer :: cache_a, cache_b, cache_z, idx, offs, offs2, err
  integer :: start_i, end_i, n_i, izz
  integer,save :: old_start_i=-1, old_end_i=-1
  integer,save :: IIc, KTc 
  integer,save :: JYc(10) = (/1,1,1,1,1,1,1,1,1,1/)
  integer,save :: firsttime = 1
  real,save,dimension(21) :: FM_PA
  real,dimension(4,4) :: AA

  ! use_all_point = .true.  usa tutte le 10 xmetal per l'interpolazione
  ! spline quando ho esaurito idrogeno nel mesh
  ! use_all_point = .false. usa solo al max 3 xmetal prima di xzora e 
  ! 3 xmetal dopo tale valore.
  logical,parameter :: use_all_points = .true.

  if(firsttime == 1) then
     ! scelgo la mistura da utilizzare in base agli input
     if(misturaop == 1) then        ! AS09a0
        do i=1,21
           FRACMASS(i) = FRACMASS09(i)
        end do
     else if(misturaop == 2) then   ! AS05a0
        do i=1,21
           FRACMASS(i) = FRACMASS05(i)
        end do
     else if(misturaop == 3) then   ! AS09a3
        do i=1,21
           FRACMASS(i) = FRACMASS09a3(i)
        end do
     else if(misturaop == 4) then   ! GS98
        do i=1,21
           FRACMASS(i) = FRACMASSGS98(i)
        end do
     else if(misturaop == 5) then   ! GN93
        do i=1,21
           FRACMASS(i) = FRACMASSGN93(i)
        end do
     else
        write(*,*)'ATTENZIONE! Mistura non definita.' 
        write(66,*)'101 - kappa'
        write(66,*)'ATTENZIONE! Mistura non definita.' 
        stop 
     endif
     
     do i=3,21
        FM_PA(i) = FRACMASS(i)/PESIAT(i)
     end do
     do i=1,21
        logzet(i) = log10(zet(i))
     end do
     firsttime = 0
  endif

  jump = 0

  if(xx(21) == cache21 .and.  XX(4) == cache4 .and. Y == cacheY .and. &
       xx(6) == cache6 .and.  XX(1) == cache1 ) jump = 1

  if(jump == 0) then
     xfe = XME*DEFAUFe 
     ! variazione relativa di Fe(diffusione)    
     dfe_fe = (xx(21)-xfe)/xfe 
     if(iread == 4) dfe_fe = 0.0
     ! metallicita` riscalata col Ferro       
     xzora = (1.0+dfe_fe)*XME 
     xzora_l = log10(xzora)
     
     ! eventuale eccesso di C                
     deltaCm = XX(4)-DEFAUC*xzora 
     ! eventuale eccesso di O                
     deltaOm = XX(6)-DEFAUO*xzora 

     cache1 = XX(1)
     cache4 = XX(4) 
     cache6 = XX(6) 
     cache21 = XX(21) 
     cacheY = Y 

     pmedio = XX(1)/PESIAT(1)+Y/PESIAT(2) 
     pmedio = pmedio+deltaCm/PESIAT(3)+deltaOm/PESIAT(5) 

     do i=3,21 
        pmedio = pmedio+xzora*FM_PA(i) 
     end do
     pmedio = 1.0/pmedio 

     ! passo da abbondanze in massa ad abbondanze numeriche             
     ABBNUM(1) = XX(1)*pmedio/PESIAT(1) 
     ABBNUM(2) = Y*pmedio/PESIAT(2) 
     deltaC = deltaCm*pmedio/PESIAT(3) 
     deltaO = deltaOm*pmedio/PESIAT(5) 

     do i=3,21 
        ! abb rel numerica dell'iesimo elemento                   
        ABBNUM(i) = xzora*FM_PA(i)*pmedio
     end do

     ABBNUM(3) = ABBNUM(3)+deltaC 
     ABBNUM(5) = ABBNUM(5)+deltaO 
  endif

  !      leggo le tabelle: M=Z , J=H, K=temp, L=log(rho/T**3)             
  if(IZ == 0) then 
     open(UNIT=3, FILE='Opacita.AF.OPAL.POT',STATUS='OLD')

     if(misturaop == 1) then
        write(67,*) "================================"
        write(67,*) "Tabelle di opacita': Opacita.AF09.OPAL.POT06_alfa0p0"
        write(67,*) "Asplund 2009"
        write(67,*) "================================"
        write(67,*)
     else if(misturaop == 2) then
        write(67,*) "================================"
        write(67,*) "Tabelle di opacita': Opacita.AF05.OPAL.POT06_alfa0p0"
        write(67,*) "Asplund 2005"
        write(67,*) "================================"
        write(67,*)        
     else if(misturaop == 3) then
        write(67,*) "================================"
        write(67,*) "Tabelle di opacita': Opacita.AF09.OPAL.POT06_alfa0p3"
        write(67,*) "Asplund 2009 alpha enhanced 0.3"
        write(67,*) "================================"
        write(67,*)  
     else if(misturaop == 4) then
        write(67,*) "================================"
        write(67,*) "Tabelle di opacita': Opacita.AFGS98.OPAL.POT06_alfa0p0"
        write(67,*) "GS 1998"
        write(67,*) "================================"
        write(67,*) 
     endif

     ! lettura della prima riga       
     read(3,FMT='(47A1)',IOSTAT=IERROR) INTEST 
     if (IERROR /= 0) then 
        write(*,*)'ATTENZIONE! OPACITA NON PRESENTI' 
        write(66,*)'100 - kappa'
        write(66,*)'ATTENZIONE! OPACITA NON PRESENTI' 
        stop 
     endif

     ! Lettura delle tabelle
     do M=1,10 
        read(3,FMT='(D11.4)') xmetal(M) 
        read(3,FMT='(71A1)',IOSTAT=IERROR)  TESTO 

        !  Viene creata una tabella che ha come entrate i valori di Y per  
        ! l'M-simo valore di Z e il J-simo valore di X.                
        ! il valore corrisp a J=1 e' zero (X=1-Z). OSS: J e' letto al contrario
        ! rispetto all'ordine nel tabellone                                    
        ! Inoltre memorizzo nel vettore metallicita i 10 valori di Z         

        do J=10,1,-1 
           read(3,FMT='(47A1)') COMPO 
           read(3,FMT='(5X,F7.1)') TMP1 
           do K=1,139 
              read(3,FMT='(F6.3,29F7.3)') TK(K),(OP(M,L,K,J),L=1,29) 
           end do
        end do
     end do
     do M=1,10 
        if(M > 1) then
           xmetallog(M) = log10(xmetal(M)) 
        else
           xmetallog(M) = -20.   ! non posso fare log(0), imposto a -20.
        endif
        do J=2,10 
           YKZM(J,M) = YKZ(J)-xmetal(M) 
        end do
     end do
     close(3) 
     IZ = 1 
  endif

  !************************* CALCOLO OPACITA'**************************** 

  T6 = TEMP/1000000 
  R = RHO/(T6**3)
  T  = log10(TEMP) 
  RLG = log10(RHO) 
  RO = RLG - 3.0*T + 18.0 

  !**********************     BIFORCAZIONE    *************************** 
  !  Questo controllo mi dice che se sto ancora bruciando idrogeno devo   
  !  usare il vecchio metodo di interp, altrimenti uso le subroutine OPAC,
  !  e POTEKHIN e sommo gli inversi.                                      

  if(XX(1) > 1.d-30) then 
     !*********************   Se brucio idrogeno    ************************ 

     !******************** CERCA LE QUATTRO TEMPERATURE *********************
     if(TK(139) < T) then 
        K = 139 
     else if (TK(2) > T) then 
        K = 2 
     else 
        ISUP = 139 
        IINF = 2 
        do while (ISUP - IINF > 1)
           K = (ISUP + IINF)/2 
           if(TK(K) > T) then 
              ISUP = K 
           else 
              IINF = K 
           endif
        end do
        K = ISUP
     end if

     KT = K-2 
     if(KT < 1) KT = 1 
     if(KT > 136) KT = 136 

     
     !******************** CERCA LE QUATTRO DENSITA' *********************** 

     ! +8 perche' LogR comincia da -8          
     II = int((RO+8.)/0.5)+1 
     if(II < 2 ) II = 2 
     if(II > 27) II = 27 

     ! questa cache controlla se ci sono variazioni per II o KT
     cache_a = 1
     if(II /= IIc .or. KT /= KTc) then
        cache_a = 0
        IIc = II
        KTc = KT
     endif

     !***************** INTERPOLA CUBICAMENTE SU RO E T ******************** 
     !    lo fa per tutti gli Z  possibili, e  per gli H che gli servono     
     !********************************************************************** 
     if( use_all_points .eqv. .false.) then
        ! ricerca z nel vettore xmetal
        do izz=1,10
           if(xzora < xmetal(izz)) exit
        end do
        start_i = max(1,izz-3)
        end_i = min(izz+2,10)
        n_i = end_i - start_i + 1
!!$        if(old_start_i /= start_i .or. old_end_i /= end_i) then
!!$           old_start_i = start_i
!!$           old_end_i = end_i
!!$           cache_z = 0
!!$        endif
     else
        start_i = 1
        end_i = 10
        n_i = 10
     endif

     !  inizia ciclo sulle metallicita'                                      
     do NZ=start_i,end_i 

        ! cerco i valori di Y che intrappolano il valore corrente
        ! se i valori in cache non sono buoni
        if(.not. (YKZM(JYc(NZ),NZ) < Y .and. YKZM(JYc(NZ)+1,NZ) > Y)) then  
           if(YKZM(10,NZ) < Y) then 
              JY1 = 10 
           else if (YKZM(1,NZ) > Y) then 
              JY1 = 1 
           else 
              ISUP = 10 
              IINF = 1 
              do while(ISUP - IINF > 1)
                 JY1 = (ISUP + IINF)/2 
                 if(YKZM(JY1,NZ) >= Y) then 
                    ISUP = JY1 
                 else 
                    IINF = JY1 
                 endif
              end do
              JY1 = IINF 
           end if
           if (JY1 == 10) JY1 = 9 
           cache_b = 0
           JYc(NZ) = JY1
        else
           ! il valore in cache e' buono
           JY1 = JYc(NZ)
           cache_b = 1
        endif

        ! calcolo se ho i valori in cache (icache = 1)
        icache = cache_a*cache_b

        do JYY=1,2 
           ! interpolazione in logR: trova l'opacita'              
           ! per i 4 R corrispondenti alle 4 temperature di griglia ed ai 2 Y

           ! se non ho valori in cache devo fare il calcolo dei coeff.
           ! della cubica per RO
           if(icache /= 1) then 
              MI = II 
              JOTA = JY1 + JYY - 1
              do M=1,4
                 ! -8 perche' LogR inizia da -8 
                 B(m) = -8.+(MI-3+m)*.5
                 do N=1,4 
                    AA(N,m) = OP(NZ,MI-2+N,KT-1+M,JOTA)
                 end do
              end do
              idx = (JYY-1)*16+(NZ-1)*32
              call CUB2dr(AA,B, idx, xABGDb) 
           endif
           ! Valuto le 4 cubiche per i valori di JYY e NZ correnti.
           ! I valori vengono memorizzati nel vettore CAPR che contiene
           ! 20 blocchi = (2 * 10) di 4 valori ciascuno, corrispondenti 
           ! ai due valori di JYY e ai 10 valori di NZ
           offs = (JYY-1)*16+(NZ-1)*32
           offs2 = (JYY-1)*4 + (NZ-1)*8
           do m=1,4
              idx = (m-1)*4+offs
              CAPR(M+offs2) = ((xABGDb(idx+1)*RO+xABGDb(idx+2))*RO + &
                   xABGDb(idx+3))*RO+ xABGDb(idx+4)
           end do
        end do
     end do

     ! Interpolazione cubica in T
     ! Si usano gli 80 valori di CAPR calcolati nel ciclo precedente   
     ! e i valori in ascissa in TT
     ! Vengono fittate 20 cubiche
     ! Se use_all_points = .false. le cubiche sono meno (max 12)
     do K=1,4 
        TT(K) = TK(KT+K-1) 
     end do
     call CUB2dt(CAPR, TT, xtABGD, start_i, end_i) 

     ! Valuto i valori delle 20 cubiche per i 10 valori di NZ e i 
     ! due valori di JYY. I risultati sono 20 valori in CAPS
     ! Se use_all_points = .false. le cubiche sono meno (max 12)
     do NZ=start_i,end_i
        offs = (NZ-1)*8 + 1
        offs2 = (NZ-1)*2
        do JYY=1,2
           idx = (jyy-1)*4 + offs
           CAPS(JYY+offs2) = exp(UM*(((xtABGD(idx)*T+xtABGD(idx+1))*T+&
                xtABGD(idx+2))*T +xtABGD(idx+3)))
        end do
        ! Interpolazione lineare in Y delle 20 coppie di valori in CAPS
        ! per trovare l'opacita' a data metallicita'           
        JY1 = JYc(NZ)
        RAP = (CAPS(offs2+2)-CAPS(offs2+1))/(YKZM(JY1+1,NZ)-YKZM(JY1,NZ)) 
        CAPZ(NZ) = CAPS(offs2+1)+RAP*(Y-YKZM(JY1,NZ))
     end do

     ! Interpolo in Z con la spline   
     ! in x: i log delle metallicita'
     ! in y: i log delle CAPZ 
     do i=start_i,end_i 
        YY(i) = log10(CAPZ(i)) 
     end do

     ! costruisco e valuto la spline di Num.Rec.
     call makespline(xmetallog(start_i:end_i),yy(start_i:end_i),n_i,y2tmp)
     Z1s = xzora_l 
     call usesplint(xmetallog(start_i:end_i),yy(start_i:end_i),y2tmp, &
          n_i,z1s,CAP)

     CAP = exp(UM*CAP) 

!     if(T > 6.5) CAP = 1.025*CAP

!     if(T >= 6.5 .and. T <= 7.0) then
!        CAP = CAP / (-4459.304436 + 2630.259726*T -581.614061*T**2 +57.156195*T**3-2.106206*T**4)
!     endif

!!$     if(T >= 6.05 .and. T < 6.5) then
!!$        CAP = CAP / (1.262252e+05 -8.920206e+04*T + 2.393183e+04*T**2 -2.898511e+03*T**3 + 1.343693e+02*T**4  -1.028445e-05*T**10)
!!$     endif
!!$      if(T < 6.05) then
!!$         CAP = CAp*1.025
!!$         endif
     return
  else
     ! =====   Se non brucio piu' idrogeno  =====
     ! Chiamo la OPAC per tutti i valori di Z che ho e poi interpolo con la  
     ! spline                                                                

     ! per controllare se fa bene i cicli di lettura
     iconta = 0 
     ! Questa variabile viene messa a zero dalla OPAC
     ! quando la somma delle abbonanze in massa supera 1
     icntrl = 1 
     
     if( RO < 1. ) then
        err = 0
        if(ireadco /= 12345678 .or. use_all_points) then
           ! la prima volta faccio i giri completi
           start_i = 1
           end_i = 10
           n_i = end_i - start_i + 1
        else
           ! ricerca z nel vettore xmetal
           do izz=1,10
              if(xzora < xmetal(izz)) exit
           end do
           start_i = max(1,izz-3)
           end_i = min(izz+2,10)
           n_i = end_i - start_i + 1
        endif

        do M=start_i,end_i 
           ! indice serve alla OPAC e alla readco per aprire il file giusto 
           indice = M 
           xmetx = xmetal(indice) 
           ABBIDRO = XX(1) 
           deltaCz = XX(4)-DEFAUC*xmetx 
           deltaOz = XX(6)-DEFAUO*xmetx 
           if(deltaCz < 0) deltaCz = 0. 
           if(deltaOz < 0) deltaOz = 0. 

           call opac(xmetx,ABBIDRO,deltaCz,deltaOz,T6,R, err) 

           ! Il valore non e' in tabella
           ! Il primo giro continuo e leggo le tabelle;
           ! altrimenti esco e vado in akappetta
           if(ireadco == 12345678 .and. err == 1) exit

           if (icntrl == 1) then 
              OPAVETT(M) = OPACT 
              iconta = indice 
           else 
              OPAVETT(M) = OPAVETT(M-1) 
           endif
        end do
        !cosi la readco legge veramente da file solo la prima volta 
        ireadco = 12345678 

        if(err == 0) then
           ! costruisco e valuto la spline di Num.Rec.
           call makespline(xmetal(start_i:end_i),opavett(start_i:end_i), &
                n_i,y2tmp)
           call usesplint(xmetal(start_i:end_i),opavett(start_i:end_i), &
                y2tmp,n_i,xzora,CAPRAD)
        endif
     endif

     if(RO >= 1. .or. err == 1) then
        ! serve quando la OPAC va fuori dalle tabelle
        call akappetta(Y,XX(4),XX(6),T,RO,akapparad) 
        CAPRAD = exp(UM*akapparad) 
     else
        ! Per ora e' solo opacita' radiativa   
        CAPRAD = exp(UM*CAPRAD) 
     end if

     ! Chiamata della potekhin 
     cond = 0.
     do l=1,21 
        ZZ = logzet(l) 
        CONDREL = 0.
        ! oss: condrel e' una conducibilita' non  una opacita'      
        call potekhin(ZZ,RLG,T,CONDREL, l) 
        if(condrel == 99.) then 
           cond = 99. 
           exit
        else 
           cond = cond+(10.**condrel)*ABBNUM(l) 
        end if
     end do

     if(cond == 99.) then 
        ! questo perche' la zkappa si somma come inverso 
        zkappa = 99. 
     else 
        cond = log10(cond) 
        ZKAPPA = -3.51937915 + 3.0*T - RLG - cond 
     endif

     CAPCOND = exp(UM*ZKAPPA) 

     ! sommo gli inversi          
     CAP = (CAPRAD*CAPCOND)/(CAPRAD+CAPCOND) 

     return 
  endif
end subroutine KAPPA
