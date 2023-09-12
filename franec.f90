!**************************************************************************
!*****                       FRANEC-2013 v3.4.0                      ******
!**************************************************************************

program STELEV

! Vengono importati i vari moduli che servono a definire l'equivalente di 
! variabili globali in C (Le loro definizioni si trovano in moduli.f90):
!     - fisica: modulo contenente la matrice G, che contiene raggio etc...
!     - chim: modulo contenente la variabile del passo temporale HT1
!     - chimic: modulo contenente la matrice delle abbondanze XXX(MELE,LIM)
!               dove MELE è il numero degli elementi mentre LIM è il numero
!               massimo di mesh.
  use interfaccia
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
  use sceltachim
  use mesh
  use costanti

  implicit none

  real,parameter,dimension(5) :: TRAS = (/1.d-10,1.d-32,1.d-17,1.d-6,1.d-33/)
  real :: TEMPO = 1.d-9

  real ::  EN1 = 0., EN2 = 0., DMAS = 0.
  integer :: MCOR = 1, ITCC = 0
  integer :: INIZ = 1, NUM1 = 0
  integer :: IDIFF

  integer :: lba, md, jsh, lpll, jupit, niter, niter1, ipunch, iora, nesci
  integer :: modello, ipp, nmod, ksa, ksb, ksd
  integer :: kse, ksg, ksh, ileg, k, j, iproi, maxmv, ip, maxmvv
  integer :: ines, numax, ichim, ierr, num2, kh1, aggiusta, primo
  integer :: kh2, nstp, i, imod, nabla1, nper, nper1, ntop, l, num0, dummy
  integer :: fase, fase1, nsmorza, npunti, test_iread4, logCNO, start_bump
  integer :: shellHoff, fase_start, fase_stop
  real :: scala, frat, prov, dam, error, yb, zb
  real :: ellogv, raz, emtov, emsol, xnemax, rat, provvv, oldeta
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
  real :: lastmass, svm

  ! per ciacio
  integer :: nsole
  ! se vale 1 vuol dire che si è esaurito 
  ! l'idrogeno centrale, siamo al TO, quindi si può levare la diffusione
  integer :: nfine_diff

  ! per la procedura di rilassamento del Modstart
  logical :: do_relax = .false., do_relax_subatm = .false.
  real :: He, Zeta, Alpha, firstHT1
  real :: diff, coeff, q
  
  ! gestione iterazioni in fitta per evitare strutture senza senso
  integer :: ierr_fitta, ntry

  integer :: relax_step, max_relax_step = 50

  !Abilita la presenza di Dark Matter stile WIMP//FRANCESCO
  integer :: on_off_DM=1 ! Se on_off_DM=1 c'é la DM, se on_off_DM=0, non c'è la DM
  if(on_off_DM==1) then
      write(*,*)"MATERIA OSCURA: on"
  else if(on_off_DM==0) then
      write(*,*)"MATERIA OSCURA: off"
  endif


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
       mix, kover, par_OVER, maxAGE)! La subrutine si trova in io.f90

  nfine_diff = 0
  N_BH05 = 0
  N_CK03 = 0
  IBH05_LETTA = 0
  ICK03_LETTA = 0
  NMOD_PRE = 0

  nsole = -2
  !!temposole = 4.57d9
  write(*,*) 'ATTENZIONE: si usa un Fe iniziale dato da:'
  write(*,'(5x,"XFE = Z*",F9.7)') DEFAUFE

  write(2,666)
  write(2,555)
  LBA = 4 ! Serve a dare la uinit ai read
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
!############## Leggo i Modstart.in
  ! IREAD:   1=MODANT   2=PRES   3=HOMOG   4=HB   5=WD
  read(4,*) IREAD !Tipo di avvio
  read(4,*) NITER !Numero di step temporali massimo (corrisponde ai cicli completi del mainloop)
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
  !#####################################
  !Qua sto leggendo sempre i Modstart.in
  !#####################################
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

     do IP=1,MELE
        XX(IP) = XXX(IP,MAXME)
     end do
   !#################################################################### 
   ! Calcolo del mdodello 0 da cui parte l'evoluzione per la presequenza 
   !####################################################################
   ! IREAD==2 -> Tale valore indica la presequenza
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
        EN2 = EN2+1.d-6*G(4,K)*tras(5)*(G(5,K+1)-G(5,K))
     end do
     HT1 = abs(EN1-EN2)/(120.55*(10.**ELLOG+10.**ELLOGV))*2.5d8
     firstHT1 = HT1
     write(2,199) EN1,EN2,HT1

     IPP = 1
     call INNES(ITCC,IPP,MAXME,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax)

     RAZ = 10.**(ELLOGV-ELLOG)
     do K=1,180
        do J=1,5
           G(J,K) = G(J,K)*TRAS(J)
        end do
        GG(2,K) = G(2,K)*RAZ
        IN(K) = 1
     end do
     HT1V = HT1
     TEMPO = HT1
     XXV(3,1) = XXX(3,1)
     IPROI = 1
     MAXMV = MAXME
     EMTOV = EMTOT
   !###########################################
   ! Qua fa il caso per la main sequence e l'HB
   !###########################################
  else if (IREAD == 3 .or. IREAD == 4) then
     !  FITTING MS ED HB 
     call INNES(ITCC,IPP,MAXME,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax)
     write(17,*) 'Y=', YB,'Z=', ZB        

     if(ILEG == 1) goto 577    ! esce con errore di convergenza iniziale

     EMSOL = EMTOT/Msun
     do K=1,MAXME
        do J=1,5
           G(J,K) = G(J,K)*TRAS(J)
           GG(J,K) = G(J,K)
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
  endif


  EMSOL = EMTOT/Msun
  EMMA = EMTOT-G(5,MAXME)
  MAXMVV = MAXMV
  write(2,444) EMSOL,yb,XME,ALFA,BMAG,IFAST
  primo = 1

  ! scrittura dei parametri M, He, z, alpha su file scenario.log
  ! Utile per ripartire con pepper in automatico
  open(64,file="scenario.log", status="unknown")
  write(64,'(f7.4,1x,f8.5,1x,e13.6,1x,f7.4,1x,e13.6,1x,a6)') &
       EMSOL, yb, XME, ALFA, 100.*(1.-FRAZ), mix(misturaop)
  close(64)

  if(iread == 4) then
     ! leggo i punti massa con cui scandire HB
     ! da file si leggono le percentuali di accrescimento
     mc = emsol
     dm = (niter-11)*0.001*ETAFPR ! come da maslos.f90
     mf = mc + dm
     write(*,*) "Aggiusto massa finale..."

     open(64,file="lastmass", status="old")
     read(64,*) lastmass
     close(64)

     open(64,file="scansioneHB.in", status="old")
     read(64,*) npunti
     do i=1,npunti
        read(64,*) MassHB(i)                        ! percentuale ...
        svm = MassHB(i)
        MassHB(i) = mc + MassHB(i)*dm        ! massa
        if(i == npunti .and. svm == 1.0)  MassHB(i) = lastmass
     end do
     close(64)
     
     write(*,*) "Output mass: "
     write(*,'(50(f8.4,1x))') MassHB(1:npunti)
     ! salvo etafpr per ripartire
     open(64,file="Modstart.in-preflash", status="old")
     do i=1,24
        read(64,*)
     end do
     read(64,*) oldETA
     close(64)
  endif
  !####################
  ! INIZIO ITERAZIONI 
  !####################
  !Ogni ciclo è un modello che avanza nel tempo
  mainloop: do
     if(primo == 1) then
        primo = 0
        goto 19  ! al primo giro salta la parte iniziale
     endif
     
     MAIS = MAIS+1 ! Contatore dei cicli per la convergenza della struttura
     INES = 0
     IPRALL = KSA
     ! massimo numero di iterazioni per arrivare a convergenza
     ! il controllo su nmd serve per partire in basso con Tc, altrimenti
     ! si blocca
     NUMAX = 100
     if(NMD < 40) NUMAX = 1000

     !#######################################################
     ! Sezione errore per troppe iterazioni senza convergenza
     !#######################################################
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
     !################################################
     !Mixa gli elementi nelle zone dove c'è convezione
     !################################################
     if(MAIS <= 1) then ! lo esegue solo al primo giro per ogni passo temporale
        call MIXING(MAXME,iread) 
        
        if(XXX(1,1) > 1.d-5 .and. G(4,1) > 6.8 .and. STDCHIM == 1) then
           call EQLB(MAXME)
        endif
        do K=1,MAXME
           XNEMAX = xxx(4,k)/12.+xxx(5,k)/14.+xxx(6,k)/16.
           if(XXX(9,K)/22. > XNEMAX) XXX(9,K) = XNEMAX*22.
           do J=1,MELE
              if(XXX(J,K) < 1.d-40) XXX(J,K) = 0.
           end do
        end do
     endif


     call QUATM(NABLA1,INIZ,MAXME,MAXMV,EMTOV)
     call HENYEY(MAXME,MAXMV,ERROR,NUM0,IERR)

     ! CONTROLLO CONVERGENZA 
     ! QUESTA SEZIONE E' QUELLA CHE USA IL MAIN LOOP SENZA AVANZARE DI PASSO TEMPORALE.
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
     if(aggiusta == 0) then ! Entro qui solo se la struttura non converge, e quindi provo a modificare i parametri fichè non funziona
        ! AGGIUSTO  CORREZIONI
        FRAT = 1.0
        if(XXX(1,1) <= 0. .and. XXX(3,1) <= 0. .and. abs(DAM) > 1.d1) then
           FRAT = 0.5
        endif

        if(NMD <= 3 .and. MAIS < 3) FRAT = 0.25
        if(MAIS <= 2 .and. MAYA /= 0) FRAT = 0.5
        if(MAIS > 12) FRAT = 0.5
        if(MAIS > 25) FRAT = 0.1

        ierr_fitta = 1
        ! il ciclo serve per controllare che las truttura calcolata abbia
        ! senso. Altrimenti la fitta ripristina G e riprova con frat piu'
        ! piccolo.
        ntry = 1
        do while( ierr_fitta /= 0) 
           call FITTA(MAXME,DAM,NUM1,NUM2,FRAT, ierr_fitta)
           if (ierr_fitta == 1) frat=frat*0.7
           ntry = ntry + 1
           if(ntry == 12) then
              write(*,*)  "4 - Troppe iterazini in fitta..."
              write(66,*) "4 - Troppe iterazini in fitta..."
              stop
           endif
        end do
        cycle ! cycle fa ripartire il loop dallo step successivo dopo che i parametri della stella sono stati cambiati perchè non convergeva
     endif
19   continue !Solo se chiami goto 19 viene eseguita questa riga, altrimenti viene ignorata


   !##################################
   ! Sezione comune per tutti i cicli
   !##################################

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

     if(do_relax_subatm .and. relax_step < max_relax_step .and. nmd > 3) then
        call rilassamento_subatm(nmd, fraz, ht1, emtot, maxme)
        relax_step = relax_step + 1
     endif

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
          lumin,TMAX, nsmorza, fase, logCNO, kover, nsole)

     ! PASSO TEMPORALE E MEMORIZZAZIONE
     if(MODELLO /= 1 .or. IREAD == 2) then
        call PASTEM(MAXME,MAXMV,FG3,PROVVV,FG2,ECNOV,FG6,FG6V,TMAX, &
             div_pastem, max_pastem)
        
        ! eta' sole
        if(Niter == fase_Sole .and. tempo < temposole .and. &
             tempo + HT1 > temposole) then
           HT1 = temposole - tempo
           nsole = nmd + 1 !! prossimo modello !!
        end if
     endif

     !Chiamo la funzione per calcolarmi la cattura della materia oscura.//FRANCESCO
     !Va chiamata dopo PASTEM o il passo temporale è per la vecchia struttura
     if(on_off_DM==1)then
      call Cattura_DM()
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
           ! se gira pepper lascio come fase bersaglio ex Hec
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
        if(IREAD == 4) ETAFPR_print = oldETA
        
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
             teff_in, teff_target, XHei, shellHoff, nsole)
     endif
     ! idfase serve a state per il raccordo delle EOS
     idfase = fase 

     if(NITER < 0) then
        ! fase evolutiva bersaglio raggiunta. USCITA.
        if (fase == NITER .or. NMD > 155999 .or. test_maxAGE) exit 
        ! devo arrivare oltre, ma ho il flash...
        if(NITER < fase_flashHe .and. fase == fase_flashHe) then
           write(66,*) "-1 - SEQUENZA TERMINATA A CAUSA DEL FLASH..."
           write(*,*) "SEQUENZA TERMINATA A CAUSA DEL FLASH..."
           ! scrivo il file flash.log
           open(64,file="flash.log", status="unknown")
           write(64,'(f7.4,1x,f8.5,1x,e13.6,1x,f7.4,1x,f13.8,1x,a6)') &
                EMSOL, yb, XME, ALFA, tempo*1e-9, mix(misturaop)
           close(64)
           stop
        endif
     endif

     ECNOV = FG2
     PROVVV = FG3
     FG6V = FG6
     TEFFVV = TEFFV
     MAXMVV = MAXMV
     do K=1,MAXME
        do J=1,6
           GG(J,K) = G(J,K)
        end do
     end do
     MAXMV = MAXME
     TEFFV = TEFF
     
     XXV(1:mele,1:maxme) = XXX(1:mele,1:maxme)
     if(NMD >= 3 .and. KSE == 0 .and. MODELLO >= 2) then
        call OPTIM(MAXME,SCALA,FG3,ENCNO, URA,ULA,UPA,UTA,UMA,VMM,fase)
     endif

     ICHIM = 0
     XCE(1:mele,1) = XXX(1:mele,1)
     ! qui il valore del tempo con cui e' chiamata ciacio verra' usato SOLO
     ! per la scrittura in scrivi e in check.
     ! ATTENZIONE: il modello che si sta cercando di calcolare qui e'
     !             quello che verra' stampato come (nmd+1)-esimo

     ! Se esaurisce l'idrogeno centrale
     ! nella ciacio viene messo nfine_diff = 1. 
     ! Da questo punto si può anche levare la diffusione

     call zone_convettive(MAXME-1, .true.)
          
     if(do_relax_subatm .and. relax_step < max_relax_step) then
        HT1 = firstHT1

        XCE(1:mele,2) = XXX(1:mele,1)

        TEMPO = TEMPO + HT1
        EMTOV = EMTOT
        EMSOL = EMTOT/Msun
        NABLA1 = -1
        MAIS = 0
        cycle mainloop
     endif

     !      IDIFF attiva o disattiva la diffusione
     if(IDIFF  ==  1 .and. nfine_diff /= 1 .and. IREAD /= 4) then 
        call ciacioLi(tempo,nmd,maxme,nsmorza, nsole, nfine_diff)
     endif

     
     if(iread /= 4) then
        ! metto la comp. al mesh maxme uguale a quella maxme-1
        XXX(1:mele, MAXME-1) = XXX(1:mele, MAXME-2)
        XXX(1:mele, MAXME) = XXX(1:mele, MAXME-2)
     endif

     call EVOLUT(1,MAXME)

     XCE(1:mele,2) = XXX(1:mele,1)

     TEMPO = TEMPO + HT1
     EMTOV = EMTOT

     ! modificare fase_start e fase_stop per determinare le fasi in cui
     ! e' attiva la maslos
     fase_start = fase_MS
     fase_stop = fase_pulsiT
     call MASLOS(MAXME,DMAS, nmd,iread, fase, fase_start, fase_stop)

     EMSOL = EMTOT/Msun
     NABLA1 = -1
     MAIS = 0
  end do mainloop

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
