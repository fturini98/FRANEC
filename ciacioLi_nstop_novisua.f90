! ***************** Si ferma solo se sbordano i metalli ********
! ***************** o X, Y. Non si ferma se sbordano Li, Be, B *****

subroutine ciacioLi(tempo,nmd,maxme,nsmorza, nsole, nfine_diff)
  use interfaccia
  use fisica
  use strut
  use chimic
  use tempe
  use chim

  implicit none

  ! questa e' la versione che fa i calcoli nei mesh e non nell'intermesh
  !     in data 26/08/95
  ! in data 13/09/95 e' stata adattata per il calcolo della diffusione dei
  ! metalli compreso il ferro (modifica in XXX(MELE,lim) )
  ! in data 6/9/2010 e' stata riscritta per gestire dinamicamente in 
  ! modo diretto gli elementi da diffondere

  !  M: numero di specie (INCLUSI GLI ELETTRONI) che si fanno diffondere.
  !  nmd+1: numero del modello `nuovo' che voglio calcolare 
  !  nabla: parametro che dice di stampare sul PRINT.DAT quando vale 3
  !  tempo: eta' del nuovo modello che verra' prodotto (in anni)
  !  tempo-deltat: eta' del modello `vecchio' (in anni) 
  !  NOTA: il tempo non mi serve per la diffusione. Uso solo deltat per sapere
  !       quanti passi (npassi) devo usare per far diffondere gli elementi
  !  maxme: numero di mesh della parte interna della nuova struttura che sara'
  !         prodotta

  ! NOTA: quando stampo qui,faccio la stampa prima di aver trovato la
  ! convergenza e il nuovo modello... (sono le vecchie quantita' diffuse...)

  real :: tempo
  integer :: nmd, maxme, nsmorza, nsole, nfine_diff

  integer,parameter :: M = 11
  ! vettore di selezione (1/0) degli elementi da diffondere
  integer,parameter,dimension(MELE) :: sel_ele = (/1,0,1,1,1, 1,0,0,0,0, &
       0,0,0,0,0, 0,0,0,0,0, 1,1,1,1,1, 0/)
  ! vettori Z e A per tutti gli elementi (all'ultimo posto ho gli elettroni)
  real,parameter,dimension(MELE+1) :: Zall = (/1.,2.,2., 6., 7., 8., 8.,10., &
       10.,12.,12.,12.,14.,0.,11., 6., 7., 8., 9.,10.,26.,3.,3.,4., 5.,1.,&
       -1.0/)
  real,parameter,dimension(MELE+1) :: Aall = (/1.,3.,4.,12.,14.,16.,18.,20., &
       22.,24.,25.,26.,28.,1.,23.,13.,15.,17.,19.,21.,56.,6.,7.,9.,11.,2., &
       0.000544617/)

  ! vettori riempiti dinamicamente al primo giro
  integer,dimension(M) :: ele
  real,dimension(M) :: Z, A
 
  real,parameter :: pigreco = 3.14159265
  real,dimension(lim) :: exx, pre, tem, grx, dm_dr, dm_dr_n
  real,dimension(m) :: X, AP, AT, C, XI
  real,dimension(m,m) :: AX, CL
  real :: ZXA, AC, NI, CZ, XIJ, NE, AO, LAMBDAD, LAMBDA
  real,dimension(m-1,lim) :: ELXXX, GRAX
  real,dimension(m,lim) :: emu, vel
  real,dimension(lim) :: GRAP, GRAT
  real :: TEMP1, TEMP2
  real,dimension(mele,lim) :: VELLO

  ! controllare save
  real,save :: dens(lim)
  integer, save :: icache = 0
  real,save,dimension(M,M) :: powZ
  ! ELLOG e' il Log (L/Lo);  TEFF e' il log(Te) (unita' cgs).
 
  real,parameter :: SUNM = 1.989d33, SUNR = 6.9599d10, TAU = 6.0d13
  real,dimension(lim) :: xsave
  real,dimension(lim) :: dm_dr_v, emutot
  real,dimension(lim) :: r, vl, dv_dr

  integer :: nscrivi_ema, k, i, mesh, j, npassi, nsmorzao, n, moda1, mstam1
  real :: tm, pr, roo, tmp, rho, t, sum1, sum2, tt, rrho, fatt
  real :: sup, deltadiff, deltao, rsmorza, y, tempomod, targetmass
  integer :: nscr, firstlog = 1
  integer,save :: firsttime = 1

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="ciacio        "

  if(firsttime == 1) then
     ! carico il vettore degli elementi da usare
     i = 1
     do k=1,mele
        if(sel_ele(k) == 1) then
           ele(i) = k
           z(i) = zall(k)
           a(i) = aall(k)
           i = i+1
        endif
     end do
     a(m) = aall(mele+1)
     z(m) = zall(mele+1)
     firsttime = 0
     
     if(i /= M) then ! controllo consistenza M e sel_ele
        write(66,*) '206 - ciacioLi'
        write(66,*) 'Il numero di elementi diffusi non concorda.'
        write(*,*) 'Il numero di elementi diffusi non concorda.'
        stop
     endif
     
     write(67,*) "@ ================================"
     write(67,*) "@             Diffusione          "
     write(67,*) "@ ================================"
     write(67,*) "@ Attiva su elementi:             "
     write(67,'(" @ ",26(i2,1x))') ele(1:m-1)
     write(67,*) "@ ================================"
     write(67,*)
  endif

  nfine_diff = 0
  nscrivi_ema = 0
  tempomod = tempo
   
  ! calcolo le grandezze da g(6,k)....
  !     G(1,k) e'      r/10^10          in cgs
  !     G(2,k) e' la   l/10^32          in cgs
  !     G(3,k) e' la  (P_totale)/10^17  in cgs
  !     G(4,k) e' la  T/10^6            in cgs     
  !     G(5,k) e' la  m/10^33           in cgs
  !     emas(k) e' la variabile di massa in Mo (non e' usata)
  !     emsol e' la massa totale in masse solari (non e' usata)
  !     ln(10)=2.302585
  !     RAG e' il raggio stellare in unita' solari
  !     RTOT e' il (raggio stellare)/10^10    in cgs 
  !     TEMP e' la temperatura efficace in k
  !     ELSU e' la (luminosita' superficiale)/10^32  in cgs
  !     TM  e' la temperatura in gradi kelvin
  !     PR  e' la pressione del gas in cgs
  !     ROO e' la densita' in cgs

  do k=1,maxme 
     TM = G(4,k) * 1.d6
     PR = G(3,k) * 1.d17 - 2.52d-15 * TM**4

     !  la STATE ha bisogno delle seguenti sei abbondanze per calcolare 
     ! la densita'.

     do i=1,6
        XX(i) = XXX(i, k)
     end do

     call STATE(caller, PR, TM, ROO, dum, dum, dum)

     !  preparo i vettori da usare per il calcolo dei gradienti:
     !  r in unita' solari
     r(k) = G(1,k)/6.9599
     pre(k) = log10(PR)
     tem(k) = log10(G(4,k)) + 6.
     dens(k) = log10(ROO)
  end do

  ! introduco le specie considerate delle XXX(mele,k) cioe':
  ! qelle che hanno 1 in sel_ele
  ! non metto il ferro in TBL (omogeneo):
  !          X(ferro)=56.*0.698*10**(-5.18)
  ! ma metto:
  !          X(6)=56. *0.698*10**(-4.49)   (vedi il soleciacio)

  ! calcolo la matrice dei logaritmi delle X(I) usate:

  ! Controllo sulla ciacio per vedere se succede del casino 
  ! Se qualche elemento sborda, o diventa negativo, fermo la ciacio 

  do k=1,maxme-1
     do i=1,m-1
        if (xxx(ele(i),k) > 1.) then
           write (*,*) 'elemento', ele(i) ,' > 1 al mesh', k, 'in ingresso'
           write(*,*) 'Mi FERMO la ciacio ha fatto casino...'
           write(66,*) '200 - ciacioLi'
           write(66,*) 'elemento', ele(i) ,' > 1 al mesh', k, 'in ingresso'
           stop
        else
           if (xxx(ele(i),k) <= 0.) then
              if(i==1) then ! idrogeno
                 ELXXX(i, k) = -40.
                 nfine_diff = 1
                 if(firstlog == 1) then
                    write(*,*) 'Diffusione OFF'
                    write(67,'(" # ",a32,1x,i5)') &
                         "Diffusione OFF al modello:", nmd
                    firstlog = 0
                 endif
              else if(ele(i) < 22) then
                 write (*,*) 'elemento ',ele(i),' < 0 al mesh',k,'in ingresso'
                 write(*,*)'Mi FERMO la ciacio ha fatto casino...'
                 write(66,*)'201 - ciacioLi'
                 write (66,*) 'elemento ',ele(i),' < 0 al mesh',k,'in ingresso'
              else
                 ELXXX(i, k) = -40.
              endif
           else
              ELXXX(i, k) = log(XXX(ele(i),k))
           endif
        end if
     end do
  end do


  !  CALCOLO I GRADIENTI ADIMENSIONALI DI  (ln T)  E (ln P):
  call deriva(maxme,r,pre,grap,0)
  call deriva(maxme,r,tem,grat,0)

  tmp = log(10.0)
  do k=1,maxme
     grap(k) = grap(k)*tmp
     grat(k) = grat(k)*tmp
  end do
  call rnmdna(grap,maxme)
  call rnmdna(grat,maxme)

  ! CALCOLO I GRADIENTI ADIMENSIONALI  delle abbondanze in massa (ln X):
  do i=1,m-1
     do k=1,maxme-1
        exx(k) = elxxx(i,k)
     end do
     call deriva(maxme-1,r,exx,grx,2)
     do k=1,maxme-1
        grax(i,k) = grx(k)
     end do
  end do

  ! dai gradienti di ABBONDANZA in massa calcolo i
  ! gradienti di CONCENTRAZIONE (che non usero'):

  !         do k=1,maxme
  !           do i=1,m-1
  !             GRADX(I)=grax(i,k)
  !           enddo
  !           TEMP=0.
  !           do I=1,M-1
  !             TEMP=TEMP+Z(I)*X(I)/A(I)*GRADX(I)
  !           enddo
  !           do J=1,M-1
  !             GRADC(J)=TEMP
  !           enddo
  !           TEMP=0.
  !           do I=1,M-1
  !             TEMP=TEMP+Z(I)*X(I)/A(I)
  !           enddo
  !           do J=1,M-1
  !             GRADC(J)=GRADX(J)-GRADC(J)/TEMP
  !             grac(J,k)=GRADC(J)
  !           enddo
  !         enddo

  !*********************************************************************
  do MESH =1, MAXME-1
     !     Nel seguito si calcolano i logaritmi coulombiani
     !     Si definiscono i valori di temperatura e densita
     !     per cui calcolare la densita elettronica (unita cgs)

     RHO = 10**dens(MESH)
     T = 10**tem(MESH)
 
     !   si definiscono i valori di X da usare tra quelli letti
     do i=1,m-1
        x(i) = xxx(ele(i),mesh)
        ! controllo per H e elem. legg.
        if(i == 1 .or. (ele(i) >= 22 .and. ele(i) <= 25)) then
           if(x(i) <= 0.) x(i) = 1.d-40
        endif
     end do

     xsave(mesh) = 0.
     do i=1,m-1
        xsave(mesh) = xsave(mesh) + X(i)
     end do

     ! calculate concentrations from mass fractions:
     ZXA = 0.
     do I=1,M-1
        ZXA = ZXA+Z(I)*X(I)/A(I)
     end do
     do I=1,M-1
        C(I) = X(I)/(A(I)*ZXA)
     end do
     C(M) = 1.
     ! calculate density of electrons (NE) from mass density (RHO):
     AC = 0.
     do I=1,M
        AC = AC+A(I)*C(I)
     end do
     NE = RHO/(1.6726d-24*AC)	
     ! calculate interionic distance (AO): 
     NI = 0.
     do I=1,M-1
        NI = NI+C(I)*NE
     end do
     AO = (0.23873/NI)**(1./3.)	
     ! calculate Debye length (LAMBDAD):
     CZ = 0.
     do I=1,M
        CZ = CZ+C(I)*Z(I)**2
     end do
     LAMBDAD = 6.9010*sqrt(T/(NE*CZ))
     ! calculate LAMBDA to use in Coulomb logarithm:
     LAMBDA = max(LAMBDAD, AO)
     ! calculate Coulomb logarithms:
     ! giada: metto in cache
     if(icache == 0) then 
        do J=1,M
           do I=1,M     
              powZ(i,j) = abs(Z(I)*Z(J))**1.2
           end do
        end do
        icache = 1
     endif
     
     tmp = (2.3939d3*T*LAMBDA)**1.2
     do J=1,M
        do I=1,M
           XIJ = tmp/powZ(i,j)
           CL(I,J) = 0.81245*log(1.+0.18769*XIJ)
        end do
     end do

     call DIFFUSION(M,A,Z,X,CL,AP,AT,AX)

     !  CALCOLO DELLE FUNZIONI DI DIFFUSIONE \XI(M):
     !******************
     !	GRADP=GRAP(MESH)
     !	GRADT=GRAT(MESH)
     !        do I =1, M-1
     !            GRADX(I)= GRAX(I,MESH)
     !            GRADC(I)= GRAC(I,MESH)
     !        enddo

     TEMP1 = 0.
     TEMP2 = 0.
     do I=1,M-1
        TEMP1 = TEMP1+Z(I)*X(I)/A(I)
        TEMP2 = TEMP2+Z(I)*X(I)/A(I)*GRAX(I,MESH)
     end do
     do I=1,M
        SUM1 = 0.
        SUM2 = 0.
        do J=1,M-1
           if (J /= 2) then
              SUM1 = SUM1+AX(I,J)*GRAX(J,MESH)
              SUM2 = SUM2+AX(I,J)
           endif
        end do
        XI(I) = AP(I)*GRAP(MESH)+AT(I)*GRAT(MESH)+SUM1-SUM2*TEMP2/TEMP1
     enddo

     !	CALCOLO DELLE VELOCITA DI DIFFUSIONE 
     !  To obtain the dimensionless diffusion velocity, simply multiply 
     !	\xi by T^{5/2}/\rho where T is in units of 10^7 K, and \rho is 
     !	in units of 100 g/cm^3.
     !	To obtain the velocity in cgs units, simply multiply the 
     !	dimensionless velocity by R_\sun/\tau_0 (see TBL).

     TT = T/1.0d7
     RRHO = RHO/1.0d2
     FATT = TT**2*sqrt(TT)/RRHO

     ! velocita adimensionale:
     do I=1,M
        if(MESH /= 1) then
           VEL(I,MESH) = XI(I)*FATT
        else
           VEL(I,MESH) = 0.
        endif
     enddo

     ! calcolo delle variabili mu totale e mu_i adimensionali:
     SUP = 4.*pigreco*(r(MESH)**2)
     emutot(MESH) = SUP * RRHO
     do j=1,m
        emu(j,MESH) = emutot(MESH) * X(j)
     end do
  end do

  ! ciclo di integrazioni con tali profili (per ogni elemento i)

  !     calcolo dt in unita' di 6*10^13 anni
  deltadiff = HT1/tau

  ! controllo sul valore di HT1 
  ! (massimo valore permesso = 10^-6 cioe' 60 My)
  deltao = 1.d-6
  npassi = 1
  if(deltadiff > deltao) then
     npassi = deltadiff/deltao+1
     deltadiff = deltadiff/npassi
     write(50,*) 'npassi diffusione=', npassi,'nmd=',nmd
  endif

  ! calcolo il valore di r da cui smorzare la velocita'
  nsmorzao = (maxme-nsmorza)/4
  rsmorza = r(maxme-nsmorzao)

  if(nsmorzao == maxme/4.) then
     nsmorzao = maxme-50
     rsmorza = r(maxme-50)
     nscrivi_ema = 1
  endif

  !! new matt
  targetmass = 0.99*g(5,maxme)
  do MESH=MAXME, 1, -1
     if(G(5,mesh) < targetmass) exit
  end do
     
!!$  if(nsmorzao < maxme-50) then
!!$     nsmorzao = maxme-50
!!$     rsmorza = r(maxme-50)
!!$     nscrivi_ema = 1
!!$  endif
  if(nsmorzao < mesh) then
     nsmorzao = mesh
     rsmorza = r(mesh)
     nscrivi_ema = 1
  endif
  
  write(*,'("  DIFF:  DT = ",F7.4," MYR ", &
       &" RSMORZA = ",F7.4,1x," NSMORZAO = ",I7,1x, &
       &" NSMORZA = ",I7)') deltadiff*6d7,rsmorza,nsmorzao,nsmorza
  if(nscrivi_ema == 1) then
     write(*,'("   Manca envelope esterno convettivo")')
     write(*,'("      Vedi ciacio")') 
  endif

  do i=1,M-1
     nscr = i
     do k=1,maxme-1
        vl(k) = vel(i,k)
        VELLO(i,k) = vel(i,k)
        dm_dr(k) = emu(i,k)
     enddo

     ! modifico i profili di velocita per r>=rsmorza in modo da portarli a zero
     ! all'ultimo mesh considerato (maxme-1)
     do j=1,maxme-1
        if (r(j) >= rsmorza) then
           ! aggiungo un controllo su j
           if(j < (maxme-1)) then 
              y = (r(j)-r(maxme-1) ) * 0.1 / (r(maxme-1) - rsmorza  )
           else
              y = 0. ! perche' tanto se j=maxme-1 il numeratore fa zero,
              ! ma potrebbe essere zero anche in denominatore
              ! allora fa casino.
           end if
           vl(j) = vl(j)*f(y)
           VELLO(i,j) = vl(j)
        endif
        if((nmd+1) == nsole) write(18,*) r(j),vl(j)
     enddo

     call cambia(dm_dr,dm_dr_v)
     n = 1

     ! controllo eventuali perdite numeriche di massa
     moda1 = 70
     mstam1 = mod(nmd,moda1)
     if(mstam1 == 0. .or. nmd+1 == nsole) then
        call check(nmd-1,maxme, nscr, tempomod, dm_dr_v, emutot, r, vl, &
             dv_dr, m)
     endif

     do while(n-1 < npassi)
        call integra(dm_dr_v,dm_dr_n,deltadiff,maxme, r, vl, dv_dr)
        ! qui ho i nuovi valori di dm; posso farci i  check e scrivere
        call cambia(dm_dr_n,dm_dr_v)
        n = n+1
     end do
     
     call decidi(nmd,maxme,xsave, nsole, nscr, tempomod, dm_dr_v,emutot,&
          r, vl, dv_dr, m, ele)
  end do

  return
end subroutine ciacioLi


real function f(y)
  !       questa serve per correggere l'andamento del campo di
  !       velocita' fra r=0.6 ed r(maxme-1) ;
  !       la funzione f vale zero per y=0, vale 1 per y=-.1 e
  !       la sua derivata e' nulla in questo punto

  implicit none
  real :: y
  f = -100.0*y*y-20.0*y
  return
end function f


subroutine deriva(max,x,y,dy,nsmo)
  use interfaccia
  use parametri

  implicit none

  integer :: max, nsmo
  real,dimension(max) :: x, y, dy

  real,dimension(lim) :: xa, ya, dya, erin, upsilon

  integer :: k, nuno, ki, num 
  real :: delmax, dx, ypn, yp1, delta

  ! controllo dimensioni
  if (lim < max) then
     write(*,*) 'aumenta lim nella subroutine deriva'
  endif

  ! immagazzinamento dati da derivare
  do k=1,max
     erin(k) = x(k)
     upsilon(k) = y(k)
  enddo

  ! calcolo minimo intervallo
  delta = x(1)
  do k=1,max-1
     dx = x(k+1)-x(k)
     if (dx < delta) then
        delta = dx
     endif
  enddo

  ! calcolo massimo intervallo
  delmax = x(1)
  do k=1,max-1
     dx = x(k+1)-x(k)
     if (dx > delmax) then
        delmax = dx
     endif
  enddo

  ! calcolo intervallo medio
  delta = (delta+delmax)/2.

  ! calcolo numero di punti necessari
2 num = (x(max)-x(1))/delta
  if (num > lim) then
     delta = 2.*delta
     write (*,*)'Delta raddoppiato.'
     goto 2
  endif

  ! calcolo intervallo corrispondente
  delta = (x(max)-x(1))/num

  !	in tal modo ho suddiviso l'intervallo tra r(1) ed r(max) in num
  !	intervalli, equispaziati di delta, e individuati da (num+1) punti
  !	di cui il primo e l'ultimo coincidenti esattamente con r(1) ed r(max)

  ! dimensionamento vettori
  num = num+1

  ! da EVITARE lo smooth dei dati

  ! interpolo con spline cubica
  yp1 = y(2)/x(2)
  ypn = (y(max)-y(max-1))/(x(max)-x(max-1))

  call SPLININT(x,y,max,yp1,ypn,num,xa,ya,0, erin, upsilon, delta)

  ! smooth dei dati "splinati"
  if (nsmo == 3) then
     call smooth(ya,num)
     call smooth(ya,num)
  else if ((nsmo == 2).or.(nsmo == 0)) then
     call smooth(ya,num)
     call smooth(ya,num)
  endif

  ! calcolo della derivata a 5 punti
  do k=1,num
     if (k == 1) then
        dya(k) = (-25.*ya(k)+48.*ya(k+1)-36.*ya(k+2)+ &
             16.*ya(k+3)-3.*ya(k+4))/(12.*delta)
     else if (k == 2) then
        ki = 1
        dya(k) = (-3.*ya(ki)-10.*ya(ki+1)+18.*ya(ki+2) &
             -6.*ya(ki+3)+ya(ki+4))/(12.*delta)
     else if (k == num) then
        ki = num-4
        dya(k) = (3.*ya(ki)-16.*ya(ki+1)+36.*ya(ki+2)- &
             48.*ya(ki+3)+25.*ya(ki+4))/(12.*delta)
     else if (k == (num-1)) then
        ki=num-4
        dya(k) = (-ya(ki)+6.*ya(ki+1)-18.*ya(ki+2)+ &
             10.*ya(ki+3)+3.*ya(ki+4))/(12.*delta)
     else
        ki=k-2
        dya(k) = (ya(ki)-8.*ya(ki+1)+8.*ya(ki+3)-ya(ki+4))/(12.*delta)
     endif
  end do

  ! smooth delle derivate
  if (nsmo == 1 .or. nsmo == 2) then
     call smooth(dya,num)
     call smooth(dya,num)
  else if (nsmo == 4) then
     call rnmdna(dya,num)
  endif

  if (nsmo == 1 .or. nsmo == 0 .or. nsmo == 3) then
     call rnmdna(dya,num)
  endif
  ! REinterpolo con spline cubica 
  yp1 = dya(2)/xa(2)
  ypn = (dya(num)-dya(num-1))/(xa(num)-xa(num-1))
  call SPLININT(xa,dya,num,yp1,ypn,max,x,dy,1, erin, upsilon, delta)

  ! smooth delle derivate reinterpolate
  if (nsmo == 1 .or. nsmo == 2) then
     call smooth(dy,max)
     call smooth(dy,max)
  endif
  return
end subroutine deriva


subroutine smooth(y,npti)
  implicit none
  integer :: npti
  real,dimension(npti) ::  y

  integer :: imax,k
  real :: yi, ynew

  imax = npti-1
  yi = y(1)

  do  k=2,imax
     ynew = (yi+2.*y(k)+y(k+1))/4.
     yi = y(k)
     y(k) = ynew
  end do
  return
end subroutine smooth


!!$subroutine smoothb(y,npti)
!!$  implicit none
!!$  integer :: npti
!!$  real,dimension(npti) ::  y
!!$
!!$  integer :: imax,k
!!$  real :: yi, ynew
!!$
!!$  imax = npti-1
!!$  yi = y(1)
!!$
!!$  do  k=1,imax
!!$     ynew = (yi+2.*y(k)+y(k+1))/4.
!!$     yi = y(k)
!!$     y(k) = ynew
!!$  end do
!!$  return
!!$end subroutine smoothb


!!$subroutine smoothc(y,npti)
!!$  implicit none
!!$  integer :: npti
!!$  real,dimension(npti) ::  y
!!$
!!$  integer :: imax,k
!!$  real :: yi, ynew
!!$
!!$  imax = npti-1
!!$  yi = y(1)
!!$
!!$  do  k=1,imax
!!$     ynew = (yi+y(k)+y(k+1))/3.
!!$     yi = y(k)
!!$     y(k) = ynew
!!$  end do
!!$  return
!!$end subroutine smoothc


!!$subroutine smdia(y,npti)
!!$  implicit none
!!$  integer :: npti
!!$  real,dimension(npti) ::  y
!!$
!!$  real,dimension(3000) :: yy
!!$  integer :: k
!!$
!!$  if (npti > 3000) then
!!$     write(*,*) 'controlla smdia'
!!$  endif
!!$
!!$  do k=1,npti
!!$     if (k > 2 .and. k < (npti-1)) then
!!$        yy(k) = (y(k+1)+y(k-1)+y(k+2)+y(k-2)+y(k))/5.
!!$     else if (k == 1) then
!!$        yy(k) = y(k)
!!$     else if (k == 2 .or. k == (npti-1)) then
!!$        yy(k) = (y(k)+y(k+1)+y(k-1))/3.
!!$     else if (k == npti) then
!!$        yy(k) = y(k)
!!$     endif
!!$  enddo
!!$  do k=1,npti
!!$     y(k) = yy(k)
!!$  end do
!!$  return
!!$end subroutine smdia



subroutine SPLININT(x,y,n,yp1,ypn,num,xa,ya,niscri, erin, upsilon, delta)
  use parametri

  implicit none

  integer :: n, num, niscri
  real :: x(n), y(n)
  real :: xa(num),ya(num), yp1,ypn
  real,dimension(lim) :: erin, upsilon
  real :: delta
  
  real,dimension(lim) :: u, ydue
  real :: sig, p, qn, un, xint
  real :: a, b, rh
  integer :: k, j, klo, khi, i

  !*********************************************************
  !               CALCOLO DELLE DERIVATE SECONDE
  !*********************************************************

  !  Calcola la derivata 2 in i=1 usando il valore dato in i=1 della
  !  derivata prima (se la derivata prima in i=1 e' >=10^30
  !  la der. 2 in 1 e' posta a zero )

  if (yp1 > 0.99d30) then
     ydue(1) = 0.0
     u(1) = 0.0
  else
     ydue(1) = -0.5
     u(1) = (3./(x(2)-x(1))*(y(2)-y(1))/(x(2)-x(1))-yp1) 
  endif

  do i=2,n-1
     sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
     p = sig*ydue(i-1)+2.
     ydue(i) = (sig-1.)/p
     u(i) = (6.* ((y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1)) &
          / (x(i)-x(i-1))) / (x(i+1)-x(i-1)) -sig*u(i-1)) /p
  end do

  !  Calcola la derivata 2 in i=n usando il valore dato in i=n della
  !  derivata prima (se la derivata prima in  e' >=10^30
  !  la der. 2 in n e' posta a zero )

  if (ypn > 99d30) then
     qn = 0.
     un = 0.
  else
     qn = 0.5
     un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif

  ydue(n) = (un-qn*u(n-1))/(qn*ydue(n-1)+1.)
  do k=n-1,1,-1
     ydue(k) = ydue(k)*ydue(k+1)+u(k)
  end do

  !*********************************************************
  !         CALCOLO DELL'INTERPOLAZIONE tra j-1 e j+1
  !*********************************************************

  do j=1,num
     ! scelta del punto da interpolare
     if (niscri == 0) then
        if (j == 1) then
           xa(j) = x(j)
           ya(j) = y(j)
           xint = x(j)
           cycle
        else if (j /= 1 .and. j /= num) then
           xint = xint+delta
           xa(j) = xint
        else if (j == num) then
           xa(num) = x(n)
           ya(num) = y(n)
           cycle
        endif
     else if(niscri == 1) then
        xa(j) = erin(j)
        xint = erin(j) 

        if (j == 1) then
           ya(j) = y(j)
           cycle
        else if (j == num) then
           ya(num) = y(n)
           cycle
        endif
     endif

     ! ricerca della posizione nella tabella
     klo = 1
     khi = n
     do while (khi-klo > 1)
        k = (khi+klo)/2
        if (x(k) > xint) then
           khi = k
        else
           klo = k
        endif
     end do

     ! controllo sui punti forniti
     rh = x(khi)-x(klo)
     if(rh == 0.) then
        write(*,*) 'i dati in x devono esser distinti'
        write(66,*)'210 - splinint'
        write(66,*) 'i dati in x devono esser distinti'
        stop
     endif

     ! calcolo dell'interpolazione cubica
     a = (x(khi)-xint)/rh
     b = (xint-x(klo))/rh
     ya(j) = a*y(klo)+b*y(khi)+ &
          ((a**3-a)*ydue(klo)+(b**3-b)*ydue(khi))*(rh**2)/6.
  end do

  return
end subroutine SPLININT


subroutine rnmdna(y,max)
  use parametri

  implicit none
  integer :: max
  real :: y(max)

  real :: b(7), yy(lim)
  integer :: imed, kernel, mesh, k, j, i
  real :: ba

  ! controllo dimensioni
  if (max > lim) then
     write(*,*) 'controlla rnmdna'
     write(66,*)'220 - rnmdna'
     write(*,*) 'max > MAX MESH'
     stop
  endif

  ! determina kernel
  kernel = 9
  imed = kernel/2+1

  do mesh=2,max-1
     if (mesh == 2 .or. mesh == max-1) then
        kernel = 3
        imed = kernel/2+1
     else if (mesh == 3 .or. mesh == max-2) then
        kernel = 5
        imed = kernel/2+1
     else if (mesh == 4 .or. mesh == max-3) then
        kernel = 7
        imed = kernel/2+1
     endif
     ! seleziona schiera
     do i=1,kernel
        b(i) = y(mesh-imed+i)
     end do

     ! ordina la schiera
     do j=1,kernel-1
        do k=j+1,kernel
           if (b(k) < b(j)) then
              ba = b(j)
              b(j) = b(k)
              b(k) = ba
           endif
        end do
     end do

     ! seleziona mediana
     yy(mesh) = b(imed)
  end do

  ! riscrivi mediane
  do k=2,max-1
     y(k) = yy(k)
  end do

  return
end subroutine rnmdna


subroutine cambia(X,Y)
  use parametri

  implicit none

  real,dimension(lim) :: X,Y
  integer :: i

  do i=1,lim
     y(i) = X(i)
  end do
  return
end subroutine cambia


subroutine integra(dm_dr_v,dm_dr_n,deltadiff,maxme, r, vl, dv_dr)
  use parametri

  implicit none

  real,dimension(lim) :: dm_dr_v, dm_dr_n, r, vl, dv_dr
  real :: deltadiff
  integer :: maxme

  integer :: mai, i

  mai = maxme-1

  i = 2
  dm_dr_n(i) = dm_dr_v(i)-deltadiff*(dm_dr_v(i)+dm_dr_v(i+1))* &
       (vl(i)+vl(i+1))/(2.*(3.*r(i)+r(i+1)))

  do i=3,mai-1
     dm_dr_n(i) = dm_dr_v(i)+deltadiff*(((dm_dr_v(i)+dm_dr_v(i-1))* &
          (vl(i)+vl(i-1)))-((dm_dr_v(i)+dm_dr_v(i+1))*(vl(i)+vl(i+1)))) &
          /(2.*(r(i+1)-r(i-1)))
  end do

  i = mai
  dm_dr_n(i) = dm_dr_v(i)+(dm_dr_v(i-1)+dm_dr_v(i))* &
       (vl(i-1)+vl(i))*deltadiff/(2.*(r(i)-r(i-1)))

  return
end subroutine integra


subroutine decidi(nmd,maxme,xsave, nsole, nscr, tempomod, dm_dr_v, emutot, &
     r, vl, dv_dr, m, ele)
  use interfaccia
  use chimic

  implicit none

  integer :: nmd, maxme, nsole, nscr, m
  real :: tempomod
  real,dimension(lim) :: xsave, dm_dr_v, emutot, r, vl, dv_dr
  integer,dimension(m) :: ele

  ! si ferma l'integrazione dopo n passi e si chiama questa subroutine
  real,dimension(m-1,lim) ::  xnew

  integer :: mai, k, i

  ! calcolo delle nuove abbondanze chimiche dalle dm_dr_v:
  mai = maxme-1

  do k=2,mai
     xnew(nscr,k) = dm_dr_v(k)/emutot(k)
     !******   modifica fatta da Scilla ********
     if (xnew(nscr,k) > 1.) then
        write(*,*) 'ATTENZIONE AL MODELLO', nmd 
        write(*,*) 'ELEMENTO n.',nscr,'> 1 !!!'
        !****** si ferma quando fa casino *****************
        write(*,*)'Mi FERMO la ciacio ha fatto casino...'
        write(66,*)'230 - decidi'
        write(66,*) 'ATTENZIONE AL MODELLO', nmd 
        write(66,*) 'ELEMENTO n.',nscr,'> 1 !!!'
        stop   
     endif
  end do
  xnew(nscr,1) = xnew(nscr,2)

  ! scrivo le nuove abbondanze nel modulo CHIMIC:
  do k=1,mai
     xxx(ele(nscr),k) = xnew(nscr,k)
     if (nscr == 1 .and. xnew(nscr,k) < 1.d-39) then
        xxx(nscr,k) = 0.
     endif
  end do
  
  ! calcolo l'elio come differenza (per la conservazione della massa)
  ! alla fine:
  if(nscr == m-1) then
     do k=1,maxme-1
        xxx(3,k) = xsave(k)-xxx(1,k)-xxx(4,k)-xxx(5,k) &
             -xxx(6,k)-xxx(21,k)-xxx(22,k)-xxx(23,k)-xxx(24,k)-xxx(25,k)
     end do
  endif
  !   Scilla: ATTENZIONE!!!!!  se l'elio e' < 0 ...
  do k=1,maxme-1
     if(xxx(3,k) < 0. ) then
        write(*,*) 'ATTENZIONE Y < 0' 
        write(*,*) 'al mesh',k
        !****** ... lo faccio fermare ************
        write(*,*)'Mi FERMO la ciacio ha fatto casino...'
        write(66,*)'231 - decidi'
        write(66,*) 'ATTENZIONE Y < 0' 
        write(66,*) 'al mesh',k
        stop
     endif
  end do

  if(nmd+1 == nsole) then
     do i=1,maxme-1
        write(nscr+69,101) nmd+1,r(i),dm_dr_v(i),xxx(ele(nscr),i)
     end do

     call check(nmd,maxme, nscr, tempomod, dm_dr_v, emutot, r, vl, dv_dr, m)
  endif
101 format(1x,i6,3(1x,1pd15.6))

  return
end subroutine decidi


subroutine check(nmd,maxme, nscr, tempomod, dm_dr_v, emutot, r, vl, dv_dr, m)
  use parametri

  implicit none

  integer :: nmd, maxme, nscr, m
  real :: tempomod
  real,dimension(lim) :: dm_dr_v, emutot, r, vl, dv_dr

  ! qui il calcolo della massa totale include solo quella interna
  !       devo dividere per il fattore M_o/(R_o^3 * 100) per avere le masse
  !       in unita' solari
  real,parameter :: fatt = 19.89/6.9599**3
  real, parameter :: TAU = 6.d13

  integer :: i, mai
  real :: emlim, dra, sum, sum_pro, sumtot, t9

  !   emlim = ~massa `limite' entro cui ho combustione 
  ! (corrispondente a 0.3 Mo)

  emlim = 0.3*fatt
  i = 2
  dra = (3.*r(i)+r(i+1))/2.
  sum = dm_dr_v(i)*dra
  sum_pro = sum
  sumtot = emutot(i)*dra

  mai = maxme-1

  ! da calcolare mescon che qui tengo fisso ................
  do i=3,mai-1
     sum = sum+dm_dr_v(i)*(r(i+1)-r(i-1))/2.
     sumtot = sumtot+emutot(i)*(r(i+1)-r(i-1))/2.
     if(sumtot <= emlim) sum_pro = sum_pro+dm_dr_v(i)*(r(i+1)-r(i-1))/2.
  end do

  dra = (r(mai)-r(mai-1))/2.

  sum = (sum+dm_dr_v(mai)*dra)/fatt
  sumtot = (sumtot+emutot(mai)*dra)/fatt
  sum_pro = sum_pro/fatt

  t9 = tempomod/tau*6.d4
  write(50,102) nscr,nmd+1,t9,sum,sum_pro,sumtot
  write(6,102) nscr,nmd+1,t9,sum,sum_pro,sumtot

102 format(2(2x,i6),4(1x,1pd14.6))
  return
end subroutine check


!!$subroutine risulta(nmd,tempo,maxme)
!!$  use fisica
!!$  use strut
!!$  use chimic
!!$
!!$  implicit none
!!$
!!$  integer :: nmd, maxme
!!$  real :: tempo
!!$
!!$  real,parameter :: pigreco = 3.14159265
!!$  integer,parameter :: m = 11
!!$  real,dimension(mele) :: xtot, dmtot, dmpro
!!$  integer :: menv, l1, i, max, l
!!$  real :: t9, tuttam, ematoc, emlim 
!!$
!!$  !  calcolo la massa xtot(i) dell'elemento i (in Mo) contenuta nella regione
!!$  !  convettiva (vedi MIXING)
!!$
!!$  t9 = tempo/1.d9
!!$
!!$  tuttam = G(5,2)
!!$  MENV = 1
!!$  L1 = maxme
!!$  if(G(6,1) > 0.) then
!!$     L1 = 1
!!$     do I = 1 , mele
!!$        XTOT(I) = XXX(I,1)*G(5,2)
!!$     end do
!!$  else
!!$     do I = 1 , mele
!!$        XTOT(I) = 0.
!!$     end do
!!$  endif
!!$  MAX = MAXME-1
!!$  do L = 2,MAX
!!$     tuttam = tuttam+(G(5,L+1)-G(5,L))
!!$     if( G(6,L) >= 0. .and. G(6,L-1) < 0. ) then
!!$        do I = 1,mele
!!$           XTOT(I) = XXX(I,L)*(G(5,L+1)-G(5,L))
!!$        end do
!!$        L1 = L
!!$     elseif( G(6,L) >= 0. .and. G(6,L-1) >= 0. ) then
!!$        do I = 1,mele
!!$           XTOT(I) = XTOT(I)+XXX(I,L)*(G(5,L+1)-G(5,L))
!!$        end do
!!$        MENV =  L
!!$     endif
!!$  end do
!!$
!!$  if(MENV == MAX) then
!!$     do I = 1,mele
!!$        XTOT(I) = XTOT(I)+XXX(I,MAX)*(EMTOT-G(5,MAXME))
!!$     end do
!!$  endif
!!$
!!$  ! massa convettiva in masse solari:
!!$  do i=1,mele
!!$     xtot(i) = xtot(i)/1.989
!!$  end do
!!$
!!$  ! massa totale convettiva in masse solari:
!!$  if(MENV == MAX) then
!!$     ematoc = (EMTOT-G(5,L1))/1.989
!!$  else
!!$     ematoc = G(5,MENV)-G(5,L1)/1.989
!!$  endif
!!$
!!$  !  calcolo la massa totale dmtot(i) e la massa entro 0.3Mo dmpro(i)
!!$  !  dell'elemento i (in Mo)
!!$
!!$  emlim = 0.3*1.989
!!$
!!$  do i=1,mele
!!$     DMTOT(i) = XXX(i,1)* G(5,2)
!!$     dmpro(i) = DMTOT(i)
!!$
!!$     MAX = MAXME-1
!!$     do  L = 2 , MAX
!!$        DMTOT(i) = DMTOT(i)+xxx(i,L)*(G(5,L+1)-G(5,L))
!!$        if(G(5,L) < emlim) dmpro(i) = dmpro(i)+xxx(i,L)*(G(5,L+1)-G(5,L))
!!$     end do
!!$
!!$     DMTOT(i) = (DMTOT(i)+xxx(i,max)*(EMTOT-G(5,MAXME)))/1.989
!!$     dmpro(i) = dmpro(i)/1.989
!!$  end do
!!$
!!$  write(68,*) 'massa totale convettiva =',  ematoc
!!$  write(68,*) 'L1=',L1,'MENV=',MENV
!!$  write(68,*) 'r base zona convettiva (in Ro)=', G(1,L1)/6.9599
!!$  write(68,*) 'elem.  ','nmd  ','  t9  ', '  massa convettiva  ' &
!!$       ,'  massa totale  ', '  massa entro 0.3Mo'
!!$
!!$  write(68,*) 'H' ,nmd,t9,xtot(1), dmtot(1), dmpro(1)
!!$  write(68,*) 'He',nmd,t9, xtot(3), dmtot(3), dmpro(3)
!!$  write(68,*) 'C' ,nmd,t9, xtot(4), dmtot(4), dmpro(4)
!!$  write(68,*) 'N' ,nmd,t9, xtot(5), dmtot(5), dmpro(5)
!!$  write(68,*) 'O' ,nmd,t9, xtot(6), dmtot(6), dmpro(6)
!!$  write(68,*) 'tuttam=',tuttam
!!$
!!$  return
!!$end subroutine risulta


subroutine risultati(nmd,tempo,maxme)
  use interfaccia
  use fisica
  use strut
  use chimic
  use tempe
  use chim
  
  implicit none
  integer :: nmd, maxme
  real :: tempo

  integer,parameter :: m = 11
  real,parameter :: pigreco = 3.14159265
  real,parameter :: fatt = 19.89/6.9599**3
  real,parameter :: TAU = 6.d13
  real,dimension(lim) :: tuttamu, eragg
  real, dimension(m) :: X, sum1, sum_pro
  real,dimension(m,lim) :: emu
  integer :: k, i, kk, mai
  real :: t9, tm1, pr1, roo1, rrho1, sup, dra
  real :: sumtot, emlim

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="risultati     "

  t9 = tempo/1.d9
  !                       TM  e' la temperatura in gradi kelvin
  !                       PR  e' la pressione del gas in cgs
  !                       ROO e' la densita' in cgs
  do k=1,maxme-1
     !   si definiscono i valori di X da usare tra quelli letti
     X(1) = XXX(1,k)
     X(2) = XXX(3,k)
     X(3) = XXX(4,k)
     X(4) = XXX(5,k)
     X(5) = XXX(6,k)
     X(6) = XXX(21,k)
     X(7) = XXX(22,k)
     X(8) = XXX(23,k)        
     X(9) = XXX(24,k)
     X(10) = XXX(25,k)

     !  la STATE ha bisogno delle seguenti sei abbondanze 
     ! per calcolare la densita'.
     !  Sono passate tramite modulo CHIM
     do  i = 1,6
        XX(i) = XXX(i, k)
     end do

     kk = k
     TM1 = G(4,kk) * 1.d6
     PR1 = G(3,kk) * 1.d17 - 2.52d-15 * TM1**4

     call STATE(caller, PR1, TM1, ROO1, dum, dum, dum)

     RRHO1 = ROO1/1.0d2

     !  r in unita' solari 
     eragg(k) = G(1,k)/ 6.9599

     ! calcolo della  variabile mu totale e delle variabili mu_i adimensionali:
     SUP = 4.*pigreco* eragg(k)**2
     tuttamu(k) = SUP * RRHO1
     do i=1,m-1
        emu(i,k) = tuttamu(k) * X(i)
     end do
  end do

167 format(1x,i6,13(1x,1pd12.5))

  ! calcolo delle masse per controllo check
  ! qui la massa totale include solo quella interna

  !       devo dividere per il fattore M_o/(R_o^3 * 100) per avere le masse
  !       in unita' solari
  !       emlim = ~massa `limite' entro cui ho combustione 
  !       (corrispondente a 0.3 Mo)


  emlim = 0.3*fatt

  k = 2
  dra = (3.* eragg(k)+eragg(k+1))/2.
  do i=1,m-1
     sum1(i) = emu(i,k)*dra
     sum_pro(i) = sum1(i)
  end do
  sumtot = tuttamu(k)*dra

  mai = maxme-1
  do k=3,mai-1
     dra = (eragg(k+1)-eragg(k-1))/2.
     sumtot = sumtot+tuttamu(k)*dra
     do i=1,m-1
        sum1(i) = sum1(i)+emu(i,k)*dra
        if(sumtot <= emlim) then
           sum_pro(i) = sum_pro(i)+emu(i,k)*dra
        endif
     end do
  end do

  dra = (eragg(mai)-eragg(mai-1))/2.
  do i=1,m-1
     sum1(i) = sum1(i)+emu(i,mai)*dra
  end do
  sumtot = sumtot+tuttamu(mai)*dra

  sumtot = sumtot/fatt
  do i=1,m-1
     sum1(i) = sum1(i)/fatt
     sum_pro(i) = sum_pro(i)/fatt
     write(69,169) nmd,t9,sum1(i),sum_pro(i),sumtot
  enddo

169 format(2x,i6,4(1x,1pd14.6))

  return
end subroutine risultati
