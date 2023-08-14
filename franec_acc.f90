!**************************************************************************
!*****                       FRANEC-2011 v3.3.0                      ******
!**************************************************************************

program STELEV

  use interfaccia
  use accrescimento
  use intero
  use mistura
  use fisica
  use strut
  use gratm
  use preco
  use nummod
  use atmosfere
  use chimic
  use tempe
  use chim
  use numer
  use varie
  use print1
  use fasievolut
  use overshoot

  implicit none

  real,parameter,dimension(5) :: TRAS = (/1.d-10,1.d-32,1.d-17,1.d-6,1.d-33/)
  real :: TEMPO = 1.d-9

  real ::  EN1 = 0., EN2 = 0., DMAS = 0.
  integer :: MCOR = 1, ITCC = 0
  integer :: INIZ = 1, NUM1 = 0
  integer :: IDIFF
  integer :: mesh

  integer :: lba, md, jsh, lpll, jupit, niter, niter1, ipunch, iora, nesci
  integer :: modello, ipp, nmod, maxme, ksa, ksb, ksd
  integer :: kse, ksg, ksh, ileg, k, j, iproi, maxmv, ip, maxmvv
  integer :: ines, numax, ichim, ierr, num2, kh1, aggiusta, primo
  integer :: kh2, nstp, i, imod, nabla1, nper, nper1, ntop, l, num0, dummy
  integer :: fase, fase1, nsmorza, npunti, test_iread4, logCNO, start_bump
  integer :: shellHoff, fase_start, fase_stop
  real :: scala, frat, prov, dam, error, yb, zb
  real :: ellogv, raz, emtov, emsol, xnemax, rat, provvv
  real :: ecnov, fg6v, tempo1, encno, div_pastem, max_pastem, TMAX
  real :: ETAFPR_print, mc, dm, mf, maxAGE, xHei
  real :: teff_in, teff_target  ! per ricerca RGB
  logical :: test_maxAGE

  ! Per gestire i MOD da generare con pepper per i modelli di HB
  real,dimension(60) :: MassHB
  integer :: indice_massa = 1
  character(len=15) :: massatxt

  ! Le misture gestite da fare uscire in scenario.log
  character(len=16) :: mix(32)

  ! da leggere da file e passare a optim
  real,dimension(4) :: URA, ULA, UPA, UTA, UMA 
  real,dimension(5) :: VMM

  type(luminosity) :: lumin

  ! per ciacio
  real :: temposole
  integer :: nsole
  ! se vale 1 vuol dire che si è esaurito 
  ! l'idrogeno centrale, siamo al TO, quindi si può levare la diffusione
  integer :: nfine_diff

  ! per la procedura di rilassamento del Modstart
  logical :: do_relax = .false.
  real :: He, Zeta, Alpha

  integer :: ipri = 0

  NABLA1 = -1
  fase = 0  ! indicatore della fase evolutiva raggiunta
  logCNO = 0  ! per monitorare il momento in cui loggare CNO per pepper
  start_bump = 0 ! per trovare lo start del bump
  shellHoff = 0
  ! valori che servono per la ricerca dell'inizio di RGB e bump
  teff_in = 0.0       
  teff_target = 0.0
  lumin%M_bordo_conv = -1.0
  
  ! apre i file di I/O e legge i valori dei parametri dinamici e mistura 
  call inout(IDIFF, div_pastem, max_pastem, URA, ULA, UPA, UTA, UMA, VMM, &
       mix, kover, par_OVER, maxAGE)

  nfine_diff = 0
  N_BH05 = 0
  N_CK03 = 0
  IBH05_LETTA = 0
  ICK03_LETTA = 0
  NMOD_PRE = 0

  nsole = 1199
  temposole = 4.57d9
  write(*,*) 'ATTENZIONE: si usa un Fe iniziale dato da:'
  write(*,'(5x,"XFE = Z*",F9.7)') DEFAUFE

  write(2,666)
  write(2,555)
  LBA = 4
  IBAT = 0
  MD = 12
  IPRALL = 0
  JSH = 0
  MAYA = 0
  SCALA = 1.
  LPLL = 0
  FRAT = 1.
  JUPIT = 1
  PROV = 0.

  ! IREAD:   1=MODANT   2=PRES   3=HOMOG   4=HB   5=WD
  read(4,*) IREAD
  read(4,*) NITER
  read(4,*) NPER
  read(4,*) ISUB
  read(4,*) IQUAT
  read(4,*) IHP
  read(4,*) IFAST
  read(4,*) IPUNCH
  read(4,*) IORA
  read(4,*) NESCI
  if(IREAD == 1) then   
     read(4,*) IDIFF    ! se sono in MODANT allora leggo da qui IDIFF
     if(IDIFF == 0 ) then
        write(*,*) "DIFFUSIONE NUOVA: OFF"
     endif
  else
     read(4,*) dummy    ! altrimenti scarto il valore
  endif

  NTOP = NITER-2
  NMD = 0
  MODELLO = 0
  IPP = 0
  NMOD = 0
  TEMPO = 1.d-9
  DAM = 0.
  ERROR = 0.
  MAIS = 0
  MAXME = 180

  ! KSA = IPRALL
  ! KSB=N PRINTA OGNI N MODELLI
  ! KSD=K SCRIVE OGNI K MESH
  ! KSE=1 NON USA LA OPTIM
  ! KSF=1 USA IL REZONING DELL'ATMOSFERA
  ! KSG=N GRAFICA TEMPORALE OGNI N MODELLI
  ! KSH=M GRAFICA STRUTTURA OGNI M MODELLI
  read(LBA,*) KSA
  read(LBA,*) KSB
  read(LBA,*) KSD
  read(LBA,*) KSE
  read(LBA,*) KSF
  read(LBA,*) KSG
  read(LBA,*) KSH

  if(IREAD == 1) then
     ! LETTURA Modstart PER RIPARTENZA
     read(LBA,*) MAXME
     read(LBA,*) MAXMV
     read(LBA,*) NMOD
     read(LBA,*) MCOR
     read(LBA,*) IPP
     read(LBA,*) IPROI
     read(LBA,*) JSH
     read(LBA,*) fase

     read(LBA,4809)(((U(J,K,L),L=1,4),K=1,4),J=1,3) &
          ,(((V(J,K,L),L=1,4),K=1,4),J=1,3),(EL(L),L=1,4),(TE(L),L=1,4), &
          ELLOG,ROCEN,ETAFPR,HECE,SCALA,EMTOT,TEMPO, &
          XH,XME,ALFA,BMAG,EMMA,HT1,HT1V,TEFF,TEFFV,PMA1,PMA2,EMAXH, &
          EMAHE,FRAZ,EMTOV,FG3,FG2
     read(LBA,4809) ((G(J,K),J=1,6),K=1,MAXME)
     read(LBA,4809,end=991) ((XXX(J,K),J=1,MELE),K=1,MAXME)
     read(LBA,4809) ((GG(J,K),J=1,6),K=1,MAXMV)
     read(LBA,4810) (IN(K),K=1,MAXME)
     read(LBA,4811) logCNO, start_bump, teff_in, teff_target, lumin%M_bordo_conv
     read(LBA,'(f9.5,i4)') Xhei, shellHoff

     write(67,'("Start dal modello: ", i5)') NMOD 

     do K=1,MAXMV
        do J=1,5
           GGV(J,K) = GG(J,K)
        end do
     end do
     do IP=1,MELE
        XX(IP) = XXX(IP,MAXME)
     end do
  else if(IREAD == 2) then
     ! FITTING PRESEQUENZA 
     call INNES(ITCC,IPP,MAXME,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax)

     if(ILEG == 1) goto 577   ! esce con errore di convergenza iniziale
     ELLOGV = ELLOG
     G(5,181) = G(5,180)
     do K=1,180
        do J=1,5
           G(J,K) = G(J,K)*TRAS(J)
           GG(J,K) = G(J,K)
           GGV(J,K) = G(J,K)
        end do
     end do
     TEFFV = TEFF
     EN1 = 0.
     do K=1,179
        EN1 = EN1+G(4,K)*(G(5,K+1)-G(5,K))
     end do
     HECE = XXX(3,1)
     ITCC = 1

     call INNES(ITCC,IPP,MAXME,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax)

     EN2 = 0.
     do K=1,179
        EN2 = EN2+1.d-6*G(4,K)*(GGV(5,K+1)-GGV(5,K))
     end do
     HT1 = abs(EN1-EN2)/(120.55*(10.**ELLOG+10.**ELLOGV))*2.5d8
     write(2,199) EN1,EN2,HT1
     IPP = 1
     call INNES(ITCC,IPP,MAXME,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax)

     RAZ = 10.**(ELLOGV-ELLOG)
     do K=1,180
        do J=1,5
           G(J,K) = G(J,K)*TRAS(J)
        end do
        GG(2,K) = G(2,K)*RAZ
        GGV(2,K) = G(2,K)*RAZ
        IN(K) = 1
     end do
     HT1V = HT1
     TEMPO = HT1
     XXV(3,1) = XXX(3,1)
     IPROI = 1
     MAXMV = MAXME
     EMTOV = EMTOT
     EMTOT_v = EMTOT
  else if (IREAD == 3 .or. IREAD == 4) then
     !  FITTING MS ED HB 
     call INNES(ITCC,IPP,MAXME,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax)
     write(17,*) 'Y=', YB,'Z=', ZB        

     if(ILEG == 1) goto 577    ! esce con errore di convergenza iniziale

     EMSOL = EMTOT/1.989
     do K=1,MAXME
        do J=1,5
           G(J,K) = G(J,K)*TRAS(J)
           GG(J,K) = G(J,K)
           GGV(J,K) = G(J,K)
           IN(K) = 1
        end do
     end do
     TEFFV = TEFF
     HECE = XXX(3,1)
     HT1 = 50000./(EMSOL**3*G(4,1)**2/225.)

     if(IREAD == 4) HT1 = HT1/500.
     if(IREAD == 3 .and. XXX(1,1) <= 1.d-5) HT1 = HT1/1000.

     HT1V = HT1
     IPROI = 0
     MAXMV = MAXME
     EMTOV = EMTOT
     EMTOT_v = EMTOT
  endif


  EMSOL = EMTOT/1.989
  EMMA = EMTOT-G(5,MAXME)
  HELI = 1.-XH-XME
  MAXMVV = MAXMV
  write(2,444) EMSOL,HELI,XME,ALFA,BMAG,IFAST
  primo = 1

  ! scrittura dei parametri M, He, z, alpha su file scenario.log
  ! Utile per ripartire con pepper in automatico
  open(64,file="scenario.log", status="unknown")
  write(64,'(f7.4,1x,f8.5,1x,e13.6,1x,f7.4,1x,e13.6,1x,a6)') &
       EMSOL, HELI, XME, ALFA, 100.*(1.-FRAZ), mix(misturaop)
  close(64)

  if(iread == 4) then
     ! leggo i punti massa con cui scandire HB
     ! da file si leggono le percentuali di accrescimento
     mc = emsol
     dm = (niter-11)*0.001*ETAFPR ! come da maslos.f90
     mf = mc + dm
     open(64,file="scansioneHB.in", status="old")
     read(64,*) npunti
     do i=1,npunti
        read(64,*) MassHB(i)                        ! percentuale ...
        MassHB(i) = mc + MassHB(i)*dm - 1.d-6       ! massa
        ! la sottrazione serve per non "sbordare" con l'ultimo valore
     end do
     close(64)
  endif

  !! ema: bisogna definire la massa totale precedente 
  !! serve per la epsig
  EMTOT_v = EMTOT

!! *** ema: controllo accrescimento *** !!
  if(acc_type == 0 .or. mdot == 0.d0)then
   acc_type = 0
   mdot = 0.0d0
  end if

  write(*,*)''
  open(64,file='Parametri_Accrescimento.dat')
   if(acc_type == 0)then
    write(64,'(" Accrescimento disabilitato!!!")')
   else

    write(64,'(" Riepilogo caratteristiche accrescimento:")')
    if(acc_type == 1)write(64,'(" Accrescimento da DISCO SOTTILE")')
    if(acc_type == 2)write(64,'(" Accrescimento SFERICO")')
    write(64,'(" Massa di partenza:",F8.4,", Massa finale:",F8.4)') &
                EMTOT/1.989,Mtot_acc
    write(64,'(" Rate di accrescimento (masse_solari/anno): ", &
               1PE11.4,0P)')mdot
   end if
  close(64)
!! *********************************** !!

  ! INIZIO ITERAZIONI 
  do
     if(primo == 1) then
        primo = 0
        goto 19  ! al primo giro salta la parte iniziale
     endif

     if(NMD >= 9999)then !! quando diventano troppo grandi
      NMD = 1
      MODELLO = 1 
     end if

     MAIS = MAIS+1
     INES = 0
     IPRALL = KSA
     NUMAX = 100
     if(NMD < 20)NUMAX = 1000
     if(MAIS > NUMAX) then
        ! SE LE ITERAZIONI SONO TROPPE SI FERMA
        write(2,4444)
        if(IBAT == 0) write(*,4444)
        write(*,*) 'Mesh interno:',NUM0,'con errore',ERROR
        write(*,*) 'Mesh atmosfera:',NUM2,'con errore',DAM
        write(66,*)'1 - main'
        write(66,*) 'Mesh interno:',NUM0,'con errore',ERROR
        write(66,*) 'Mesh atmosfera:',NUM2,'con errore',DAM
        stop
     endif
     if(MAIS <= 1) then
        call MIXING(MAXME,iread)

        if(XXX(1,1) > 1.d-5 .and. G(4,1) > 6.8) call EQLB(MAXME)
        do K=1,MAXME
           XNEMAX = xxx(4,k)/12.+xxx(5,k)/14.+xxx(6,k)/16.
           if(XXX(9,K)/22. > XNEMAX) XXX(9,K) = XNEMAX*22.
           do J=1,MELE
              if(XXX(J,K) < 1.d-40) XXX(J,K) = 0.
           end do
        end do

     !! questa cerca l'inviluppo convettivo e calcola cose utili
     !! alle subroutine di accrescimento !!
     !! se non c'e' accrescimento sistema anche la combustione
     !! degli elementi leggeri !!
       call conv(HT1,NMD,MAXME)   

 !! cerco il primo punto radiativo (dall'esterno) !!

       if(Lenv1 > 2 .and. G(4,1) < 1.d0)then
        write(*,'(" STRUTTURA NON COMPLETAMENTE CONVETTIVA !!",I4)')Lenv1
        write(*,'("     Sistemo la chimica !!")')
        do mesh=1,MAXME
         XME = Z_ini
         XXX(1,mesh) = X_ini
         XH = X_ini
         XXX(2,mesh) = XME*DEFAUHe3
         XXX(3,mesh)= Y_ini 
         XXX(4,mesh) = XME*DEFAUC 
         XXX(5,mesh) = XME*DEFAUN 
         XXX(6,mesh) = XME*DEFAUO
         XXX(7,mesh) = 0.
         XXX(8,mesh) = 0.
         XXX(9,mesh) = 0.
         XXX(10,mesh) = 0.
         XXX(11,mesh) = 0.
         XXX(12,mesh) = 0.
         XXX(13,mesh) = 0.
         XXX(14,mesh) = 0. 
         XXX(15,mesh) = 0.
         XXX(16,mesh) = 0.
         XXX(17,mesh) = 0.
         XXX(18,mesh) = 0.
         XXX(19,mesh) = 0. 
         XXX(20,mesh) = 0.
         XXX(21,mesh) = XME*DEFAUFe
         XXX(22,mesh) = XME*DEFAULi6 
         XXX(23,mesh) = XME*DEFAULi7 
         XXX(24,mesh) = XME*DEFAUBe 
         XXX(25,mesh) = XME*DEFAUB 
         XXX(26,mesh) = DEFAUD 
        end do

        do mesh=1,MAXMV
         XXV(1:MELE,mesh) = XXX(1:MELE,mesh)
        end do

!       stop
       end if

      !! accrescimento !!
        if(mdot > 0.0d0 .and. NMD >= NMD_acc) then
         call accresci(HT1,NMD,MAXME)
         call accresci_chim(HT1)
         EMSOL = EMTOT/1.989
        end if
         
     endif

     call QUATM(NABLA1,INIZ,MAXME,MAXMV,EMTOV)
     call HENYEY(MAXME,MAXMV,ERROR,NUM0,IERR)

     ! CONTROLLO CONVERGENZA
     aggiusta = 0
     if(MAIS > 1) then
        if(abs(ERROR) < 1.d-3 .and. abs(DAM) < 1.d-3) aggiusta = 1
        if(abs(ERROR) < 2.d-4 .or. abs(DAM) < 2.d-4) aggiusta = 1
        if(MAIS >= 7) then
           if(abs(DAM) <= 5.d-2 .or. abs(ERROR) <= 5.d-2) aggiusta = 1
        endif
        if(MAIS >= 10) then
           if(NUM1 == 2 .and. NUM2 > 500 .and. IERR == 3 .and. NUM0 > 500) &
                aggiusta = 1
        endif
     endif
     if(aggiusta == 0) then
        ! AGGIUSTO  CORREZIONI
        FRAT = 1.0
        if(XXX(1,1) <= 0. .and. XXX(3,1) <= 0. .and. abs(DAM) > 1.d1) then
           FRAT = 0.5
        endif

        if(NMD <= 3 .and. MAIS < 3) FRAT = 0.5
        if(MAIS <= 2 .and. MAYA /= 0) FRAT = 0.5
        if(MAIS > 12) FRAT = 0.5
        if(MAIS > 25) FRAT = 0.1

        call FITTA(MAXME,DAM,NUM1,NUM2,FRAT)
        cycle
     endif
19   continue

     ! CONTROLLA SE RAGGIO E/O PRESSIONE PRESENTANO INVERSIONI
     KH1 = 0
     KH2 = 0
     do K=1,Maxme-1
        if(G(1,K) >= G(1,K+1)) then
           write(2,4414) K
           ! meglio fermarsi...
           write(66,*) '3 - main'
           write(66,4414) k
           stop
           KH1 = K+1
           exit
        endif
     end do

     ! inversione pressione ...
     do K=1,Maxme-1
        if(G(3,K) <= G(3,K+1)) then
           write(2,4415) K
        endif
     end do

     MAYA = MAYA-1
     if(MAYA < 0) MAYA = 0
     MODELLO = MODELLO+1
     NSTP = MAIS
     NMD = NMOD+MODELLO

     !  STAMPA RISULTATI 
     if( ( NMD - ( NMD / KSB ) * KSB )  ==  0 ) NABLA1 = 3
     if( MODELLO  ==  1 ) NABLA1 = 3

     if (NMD  ==  nsole)  NABLA1 = 3

     if (nmd  ==  nsole .and. IREAD /= 4) then
        !     qui si stampano i risultati del modello nmd

        ! Se esaurisce l'idrogeno centrale
        ! nella ciacio viene messo nfine_diff = 1. 
        ! Da questo punto si può anche levare la diffusione

        !    IDIFF attiva o disattiva la diffusione
        if(IDIFF  ==  1 .and. nfine_diff /= 1) then   
           call risultati(nmd,tempo,maxme)
        endif
     endif

     call STAMPA(TEMPO,DMAS,NMD,NABLA1,KSD,MAXME,MAXMV,KSG,KSH, &
          lumin,TMAX, nsmorza, fase, logCNO, kover)

     if(ipri == 0 .and. (acc_type == 0 .or. mdot == 0.d0))then
       acc_type = 0
       mdot = 0.d0
       ipri = 1
       write(*,*)''
       write(*,'(5x,"NMD:",I4," ->  ARRIVATO ALLA MASSA VOLUTA!!")')NMD
       write(*,*)''
       write(*,'(5x,"   M(voluta) = ",F7.4,", M(iniziale) = ",F7.4, &
                 ", M(finale) = ",F7.4)')Mtot_acc,EMTOT_ini/1.989,EMTOT/1.989
       write(*,*)''
       write(*,'(5x,"        FERMO ACCRESCIMENTO!!")')
     end if

     ! PASSO TEMPORALE E MEMORIZZAZIONE
     if(MODELLO /= 1 .or. IREAD == 2) then
        call PASTEM(MAXME,MAXMV,FG3,PROVVV,FG2,ECNOV,FG6,FG6V,TMAX, &
             div_pastem, max_pastem, EMTOV)
      !! ema: al primo e secondo  modello metto un controllo sul 
      !! passo temporale, perche' ancora non ho attaccato l'accrescimento.
      !! HT1 deve essere almeno 1/10 del tempo scala dell'accrescimento !!
!        if(NMD <= NMD_acc .and. mdot > 0.0d0)then
        if(NMD <= 2 .and. mdot > 0.d0)then
         if(HT1 > 0.01*EMTOT/(1.989*mdot)) HT1 =  0.01*EMTOT/(1.989*mdot)
        endif
     endif

     ! controllo se devo scrivere un modello per pepper (test_iread4 = 1)
     test_iread4 = 0
     if(iread == 4) then
        if(emsol >= MassHB(indice_massa) .and. indice_massa <= npunti) &
             test_iread4 = 1
     endif

     ! controllo se devo scrivere un MOD, ossia:
     ! ogni nper modelli, se ipunch e' diverso da 1 e per ogni modello indicato
     ! da pepper
     if((mod(MODELLO,NPER) == 0 .and. IPUNCH /= 1) .or. test_iread4 == 1) then
        if(test_iread4 == 1 .and. modello >= 11) then
           ! pepper ha trovato un modello da scrivere in out
           ! scandendo il vettore MassHB

           ! Devo scegliere il nome corretto per l'out
           write(massatxt,'("MOD-M",f6.4,".DAT")')  emsol

           open(12, FILE=massatxt, status='unknown')
           md = 12
           indice_massa = indice_massa + 1
        else
           if( MD >= 14 ) MD = 12
           if(MD == 12) then
              open(UNIT=12,FILE='MOD1.DAT',status='unknown')
           else if(MD == 13) then
              open(UNIT=13,FILE='MOD2.DAT',status='unknown')
           end if
        endif

        ! scrivo quante iterazioni mancano. Se negativo, lascio invariata
        ! perche' indica la fase evolutiva
        if(niter > 0) then
           niter1 = niter - modello + 1
        else 
           niter1 = niter
        endif
        
        if(iread == 4) then
           ! se gira pepper lascio come fase bersaglio i pulsi
           niter1 = fase_pulsiT
        endif

        ! scrivo la fase evolutiva raggiunta, mettendo fase_postPepper
        ! se gira pepper
        fase1 = fase
        nper1 = nper
        if(iread == 4) then
           fase1 = fase_postPepper
           nper1 = 20
        endif

        ! scrivo 0 al posto di ETAFPR nel caso sia uno start da HB
        ! Cosi' posso rilanciare franec senza dover modificare 
        ! il file di uscita
        ETAFPR_print = ETAFPR
        if(IREAD == 4) ETAFPR_print = 0.0
        
        write(MD,*) 1
        write(MD,*) NITER1, "! NITER"
        write(MD,*) NPER1
        write(MD,*) ISUB 
        write(MD,*) IQUAT
        write(MD,*) IHP
        write(MD,*) IFAST
        write(MD,*) IPUNCH
        write(MD,*) IORA
        write(MD,*) NESCI
        if(IDIFF == 0 .or. nfine_diff == 1) then  ! flag per la diffusione
           write(MD,*) 0, "! IDIFF"
        else
           write(MD,*) 1, "! IDIFF"
        endif
        IMOD = NMD-1
        if(iread == 4) IMOD = 0    ! per ripartire con franec dopo pepper
        write(MD,*) KSA
        write(MD,*) KSB
        write(MD,*) KSD
        write(MD,*) KSE
        write(MD,*) KSF
        write(MD,*) KSG
        write(MD,*) KSH
        write(MD,*) MAXME, "! MAXME"
        write(MD,*) MAXMV
        write(MD,*) IMOD, "! modello di partenza"
        write(MD,*) MCOR
        write(MD,*) IPP
        write(MD,*) IPROI
        write(MD,*) JSH
        write(MD,*) fase1, "! fase evolutiva raggiunta"

        write(MD,4809) (((U(J,K,L),L=1,4),K=1,4),J=1,3) &
             ,(((V(J,K,L),L=1,4),K=1,4),J=1,3),(EL(J),J=1,4),(TE(J),J=1,4) &
             ,ELLOG,ROCEN,ETAFPR_print,HECE, &
             SCALA,EMTOT,TEMPO,XH,XME,ALFA,BMAG,EMMA,HT1,HT1V,TEFF,TEFFV, &
             PMA1,PMA2,EMAXH,EMAHE,FRAZ,EMTOV,FG3,FG2
        write(MD,4809) ((G(J,K),J=1,6),K=1,MAXME)
        write(MD,4809) ((XXX(J,K),J=1,MELE),K=1,MAXME)
        write(MD,4809) ((GG(J,K),J=1,6),K=1,MAXMV)
        write(MD,4810) (IN(K),K=1,MAXME)
        ! per stop_evolut
        write(MD,4811) logCNO, start_bump, teff_in, teff_target, &
             lumin%M_bordo_conv
        write(MD,'(f9.5,i4)') Xhei, shellHoff

        close(UNIT=MD)
        
        if(test_iread4 /= 1) then
           MD = MD + 1
        endif
     endif

     ! test per maxAGE
     test_maxAGE = (tempo*1e-9 > maxAGE .and. maxAGE > 0.)

     ! iterazioni terminate. USCITA.
     if(MODELLO == NITER .or. test_maxAGE) exit 
     ! controllo quale fase evolutiva e' stata raggiunta
     if(iread == 4) then ! partenza con pepper
        fase = fase_postPepper
     else
        call stop_evolut(fase, lumin, nmd, tempo, logCNO, start_bump, &
             teff_in, teff_target, XHei, shellHoff)
     endif
     ! idfase serve a state per il raccordo delle EOS
     idfase = fase 

     if(NITER < 0) then
        ! fase evolutiva bersaglio raggiunta. USCITA.
        if (fase == NITER .or. NMD > 15999 .or. test_maxAGE) exit 
        ! devo arrivare oltre, ma ho il flash...
        if(NITER < fase_flashHe .and. fase == fase_flashHe) then
           write(66,*) "-1 - SEQUENZA TERMINATA A CAUSA DEL FLASH..."
           write(*,*) "SEQUENZA TERMINATA A CAUSA DEL FLASH..."
           ! scrivo il file flash.log
           open(64,file="flash.log", status="unknown")
           write(64,'(f7.4,1x,f8.5,1x,e13.6,1x,f7.4,1x,f13.8,1x,a6)') &
                EMSOL, HELI, XME, ALFA, tempo*1e-9, mix(misturaop)
           close(64)
           stop
        endif
     endif

     ECNOV = FG2
     PROVVV = FG3
     FG6V = FG6
     do K=1,MAXMV
        do J=1,5
           GGV(J,K) = GG(J,K)
        end do
     end do
     TEFFVV = TEFFV
     MAXMVV = MAXMV
     do K=1,MAXME
        do J=1,6
           GG(J,K) = G(J,K)
        end do
     end do
     MAXMV = MAXME
     TEFFV = TEFF
     do K=1,MAXME
        do J=1,MELE
           XXV(J,K) = XXX(J,K)
        end do
     end do

     if(NMD >= 3 .and. KSE == 0 .and. MODELLO >= 2) then
        call OPTIM(MAXME,SCALA,FG3,ENCNO, URA,ULA,UPA,UTA,UMA,VMM,fase)
     endif

     ICHIM = 0
     do I=1,MELE
        XCE(I,1) = XXX(I,1)
     end do

     ! qui il valore del tempo con cui e' chiamata ciacio verra' usato SOLO
     ! per la scrittura in scrivi e in check.
     ! ATTENZIONE: il modello che si sta cercando di calcolare qui e'
     !             quello che verra' stampato come (nmd+1)-esimo

     ! Se esaurisce l'idrogeno centrale
     ! nella ciacio viene messo nfine_diff = 1. 
     ! Da questo punto si può anche levare la diffusione

     !      IDIFF attiva o disattiva la diffusione
     if(IDIFF  ==  1 .and. nfine_diff /= 1 .and. IREAD /= 4) then
        call ciacioLi(tempo,nmd,maxme,nsmorza, nsole, nfine_diff)
     endif

     if(iread /= 4) then
        ! metto la comp. al mesh maxme uguale a quella maxme-1
        do I=1,mele
           XXX(I,MAXME-1) = XXX(I,MAXME-2)
           XXX(I,MAXME) = XXX(I,MAXME-2)
        end do
     endif

 !! ema: mi salvo gli elementi prima di cambiare la chimica, ma dopo !!
 !! aver ottimizzato i mesh !!
     do K=1,MAXME
        do J=1,MELE
           XXV_a(J,K) = XXX(J,K)
        end do
     end do

     call EVOLUT(1,MAXME)

     if(do_relax) then
        call rilassamento(NMD, MAXME, He, Zeta, Alpha)
     endif

     do I=1,MELE
        XCE(I,2) = XXX(I,1)
     end do
     TEMPO = TEMPO+HT1
     EMTOV = EMTOT

     ! modificare fase_start e fase_stop per determinare le fasi in cui
     ! e' attiva la maslos
     fase_start = fase_exH
     fase_stop = fase_pulsiT

  !! ema: qui devo aggiungerci la subroutine di accrescimento !!
  !! prima mi salvo la massa totale del vecchio modello !!
     EMTOT_v = EMTOT

   !! dovrei dividere: se etafpr > 0 => accresci, altrimenti maslos
     call MASLOS(MAXME,DMAS, nmd,iread, fase, fase_start, fase_stop)

     EMSOL = EMTOT/1.989
     NABLA1 = -1
     MAIS = 0
  end do

  !  ESCE 
  if(IBAT == 0) write(*,9999)
  write(2,9999)
  write(66,*) "0 - SEQUENZA TERMINATA"
  stop

577 write(*,'(/,"Avvio impossibile con i valori in input",/)')
  write(66,*) "2 - Avvio impossibile con i valori in input" 
991 stop

111 format(I1,2I4,10I2)
222 format(11F7.4)
444 format(//,1X,'M=',F7.4,3X,'Y=',F7.5,3X,'Z=',1P,E8.1,3X,'ML=', &
       0P,F4.1,3X,'B MAG=',E8.1,1X,'IFAST=',I2,/)
666 format(1H1)
555 format(6X,'NON SOLO GALLEX...')
4444 format(1X,'STOP PER TROPPE ITERAZIONI')
7777 format(8I4,5I2)
9999 format(///,40X,'SEQUENZA TERMINATA',//)
952 format(3F5.3,e10.3)
951 format(9I2)
199 format(1X,1P,3e12.4)
4809 format(1P,4e20.13)
4810 format(80I1)
4811 format(2(i3,1x),3(f10.6,1x))  
4414 format(1X,'INVERSIONE DEL RAGGIO AL MESH',I5)
4415 format(1X,'INVERSIONE DELLA PRESSIONE AL MESH',I5)

end program STELEV
