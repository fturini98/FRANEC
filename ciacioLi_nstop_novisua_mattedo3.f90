! ***************** Si ferma solo se sbordano i metalli ********
! ***************** o X, Y. Non si ferma se sbordano Li, Be, B *****

subroutine ciacioLi(tempo,nmd,maxme,nsmorza, nsole, nfine_diff)
  use interfaccia
  use fisica
  use strut
  use chimic
  use tempe
  use chim
  use costanti

  implicit none

  ! questa e' la versione che risolve il sistema di abbondanze in
  ! modo implicito

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
  ! convergenza e il nuovo modello (sono le vecchie quantita' diffuse).

  integer,intent(in) ::  maxme

  real :: tempo
  integer :: nmd, nsmorza, nsole, nfine_diff, ierr

  integer,parameter :: M = 11  ! numero elementi che diffondono + 1 (elettroni)
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
  integer,save,dimension(M) :: ele
  real,save,dimension(M) :: Z, A

  real,dimension(maxme) :: dens, pre, rho, tem, grx
  real,dimension(m-1,maxme) :: ELXXX, GRAX
  real,dimension(m,maxme) :: emu, vel, vgt, vel_correct
  real,dimension(maxme) :: GRAP, GRAT
  real,dimension(m,m,maxme) :: sigma

  real, dimension((m-1)*(maxme-1)) :: eleX

  ! ELLOG e' il Log (L/Lo);  TEFF e' il log(Te) (unita' cgs).

  real,parameter :: SUNM = 1.989d33, SUNR = 6.9599d10, TAU = 6.0d13
  real,dimension(maxme) :: xsave
  real,dimension(maxme) :: emutot
  real,dimension(maxme) :: r

  integer :: nscrivi_ema, k, i, n_passi_diff, m1,m2
  real :: tm, pr, roo
  real :: deltadiff, deltadiff_0, deltao, tempomod, tempo_ad
  integer,save :: firsttime = 1

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="ciacio        "

  real,dimension(mele,maxme) :: XXXsave, XXXor
  logical :: do_grad_chim_smooth = .true.
  logical, parameter :: rinormalizza_abund = .true.
  logical, parameter :: smooth_abund = .true.

  real :: mmin, mmax, mtarget
  integer :: dim, repeat
  
  real, dimension(:,:),allocatable :: a_matrix, diag, al
  integer :: nsup
  real, save :: tgt = 1
  integer :: indx

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

  m1 = 2*(m-1)-1
  m2 = m1
  dim = (m-1)*(maxme-1)
  allocate(diag(dim,1+m1+m2), a_matrix(dim,1+m1+m2), al(dim,1+m1+m2))

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
     XX(1:6) = XXX(1:6, k)

     call STATE(caller, PR, TM, ROO, dum, dum, dum)

     !  preparo i vettori da usare per il calcolo dei gradienti:
     !  r in unita' solari
     r(k) = G(1,k)/6.9599
     pre(k) = log(PR)
     tem(k) = log(TM)
     dens(k) = log(ROO)
     rho(k) = ROO/1.d2
  end do

  if(nsmorza > 0 .and. (1-G(5,nsmorza)/G(5,maxme)) > 1.0d-7 ) then  ! 1.d-7
     mmin = G(5,nsmorza)
     mmax = G(5,maxme)
     mtarget = mmin  + (mmax-mmin)/4.0
  else
     mtarget = (1.0 - 1.0d-7)*G(5,maxme)
  endif


  do k=maxme,1,-1 
     if(G(5,k) < mtarget) exit
  end do
  k = k+1
  nsup = min(k, maxme)
  write(*,*) "nsup, maxme =", nsup, maxme
  write(*,*) 1.-G(5,nsup)/G(5,maxme)

  call init_elements(maxme, M, ele, nfine_diff, XXXsave)
  XXXor(:,:) = XXXsave(:,:)

  ! abbondanze logaritmiche
  !  CALCOLO I GRADIENTI ADIMENSIONALI DI  (ln T)  E (ln P):
  call derivami(maxme,r,pre,grap)
  call derivami(maxme,r,tem,grat)

  !  calcolo dt in unita' di 6*10^13 anni (in sec)
  deltadiff_0 = HT1*sec_anno

  ! controllo sul valore di HT1 (massimo valore permesso = 60 My)
  deltao = 60d6*sec_anno

  if (deltadiff_0 > deltao) then
     n_passi_diff = int(deltadiff_0/deltao)+1
  else
     n_passi_diff = 1
  endif

  ierr = 0
  tempo_ad = 0
  repeat = 1

  ! ciclo di integrazione delle abbondanze diffuse
  integrazione: do while(tempo_ad < deltadiff_0)  
     deltadiff = deltadiff_0/n_passi_diff

     do k = 1,maxme-1
        do i = 1,M-1
           ELXXX(i, k) = log(XXXsave(ele(i),k))
           eleX((m-1)*(k-1)+i) = XXXsave(ele(i),k)
        end do
     end do

     ! CALCOLO I GRADIENTI ADIMENSIONALI  delle abbondanze in massa (ln X):
     call grad_chim(maxme, m, elxxx, r, grx, grax)
     ! smooth del gradiente della chimica per evitare bruschi salti al
     ! bordo delle zone convettive
     if (do_grad_chim_smooth) then
        do k=1,m-1
           call smooth(grax(k,:), maxme-1)
           call smooth(grax(k,:), maxme-1)
        end do
     endif

     call get_diffusion_v(maxme, M, Z, A, dens, tem, r, grax, grap, grat, &
          ele, XXXsave, vel, emu, emutot, xsave, vgt, sigma)
     call insert_turb(m, maxme, dens, tem, vel, grax, vel_correct)
     do k=2,maxme-1
        if (vel_correct(8,k)>vel(8,k)) then
           do i=1,m-1
              vel(i,k)=vel_correct(i,k)
           enddo
        endif
     enddo

     ! riempio la matrice a bande per la soluzione del sistema: 
     ! a_matrix * Abb_nuove = Abb_vecchie
     call get_matrix_4(m, maxme, dim, r, 1.d2*rho, vel, sigma, deltadiff, &
          a_matrix)
!!$     call get_matrix_4_S(nmd, m, maxme, dim, r, 1.d2*rho, vgt, sigma, eleX, deltadiff, a_matrix)

     INDX = 0

     ! risolvo il sistema ottenendo le nuove abbondanze in eleX
     call solve_system(maxme, a_matrix, m, dim, m1, m2, eleX)

     ! rinormalizza le abbondanze 
     ! scaricando la differenza sull'elemento piu' abbondante del mesh
     if(rinormalizza_abund) call rinormalizza_abbondanze(m, maxme, eleX, xsave)

     do i = 1,m-1
        do k = 1, maxme-1  
           if ((eleX((m-1)*(k-1)+i) >1.) .or. (eleX((m-1)*(k-1)+i)<0.) .or. &
                (isnan(eleX((m-1)*(k-1)+i)) .eqv. .true.)) then
              n_passi_diff = n_passi_diff*2
              XXXsave(:,1:maxme) = XXXor(:,1:maxme)
              write(*,*) 'problema: ricomincio iterazione'
              write(*,'(i3,2x,i3,2x,e11.5,2x,e20.14)') i,k,eleX((m-1)*(k-1)+i), G(5,k)/G(5,maxme)
              tempo_ad = 0
              repeat = repeat + 1 
              if(repeat < 10) then
                 cycle integrazione
              else
                 stop "troppe iterazioni"
              endif
           endif
        enddo
     enddo
     
     ! aggiorno le abbondanze
     do i = 1,m-1
        do k = 1,maxme-1
           XXXsave(ele(i),k) = eleX((m-1)*(k-1)+i)
        enddo
     enddo
     
     ! aumento il tempo e riparto con il prossimo giro
     tempo_ad = tempo_ad + deltadiff

  enddo integrazione

  ! carico le nuove abbondanze nella matrice della chimica
  do k = 1,maxme-1
     do i = 1,m-1      
        XXX(ele(i),k) = eleX((m-1)*(k-1)+i)
     enddo
  enddo

  if(nmd == 500 ) then
     do k = 1,maxme-1
        write(136,'(2(i5,1x),18(e15.8,1x))') nmd, k, G(5,k), (grax(i,k), i=1,m-1)
     end do
  endif

  ! smoothing delle abbondanze, per ovviare a problemi numerici del solutore
  if(smooth_abund) then
     ! smooth ogni 200 Myr
     if(tempomod*1d-9*5. - tgt > 1) then
        tgt = tgt + 1 
        do k = 2,maxme-1
           do i = 1,m-1
              XXXsave(ele(i),k) = 1./8.*XXX(ele(i),k-1) + 3./4.*XXX(ele(i),k) &
                   +1./8.*XXX(ele(i),k+1)
           end do
        end do
        do i = 1,m-1
           XXXsave(ele(i),1) =  7./8.*XXX(ele(i),1) +1./8.*XXX(ele(i),2)
           XXXsave(ele(i),maxme) =  7./8.*XXX(ele(i),maxme) + &
                1./8.*XXX(ele(i),maxme-1)
        end do
        
        XXX(1:mele,1:maxme) = XXXsave(1:mele,1:maxme)
     endif
  endif
  
  do k = 1,maxme-1
     do i = 1,m-1    
        if(nmd == 500 .and. i == 1) then
           write(135,'(2(i5,1x),18(e15.8,1x))') nmd, k, G(5,k), XXXor(1,k), eleX((m-1)*(k-1)+i), XXX(1,k), grap(k), grat(k), vel(1,k), rho(k)
        endif
     end do
  end do


  deallocate(a_matrix, al, diag)

  return

end subroutine ciacioLi




subroutine get_new_abundances(it, m, maxme, vel, emu, r, nsole, emutot, &
     xsave, ele, tempomod, nsmorza, deltadiff, ierr, nsup)
  use nummod

  implicit none

  integer :: m, maxme, nsole, nsmorza, it, ierr, nsup
  real,dimension(m,maxme) :: emu, vel
  real,dimension(maxme) :: r, emutot, xsave
  integer,dimension(M) :: ele
  real :: tempomod, deltadiff

  integer :: i, k, nsmorzao, nscrivi_ema
  real,dimension(maxme) :: dm_dr,  dm_dr_v,  dm_dr_n, dv_dr, vl
  integer,parameter :: moda1 = 70
  real :: rsmorza

  nscrivi_ema = 0
  ! calcolo il valore di r da cui smorzare la velocita'
  nsmorzao = (maxme-nsmorza)/4
  rsmorza = r(maxme-nsmorzao)

  if(nsmorzao == maxme/4.) then
     nsmorzao = maxme-50
     rsmorza = r(maxme-50)
     nscrivi_ema = 1
  endif

  if(it == 1) then
     write(*,'("  DIFF:  DT = ",F7.4," MYR ", &
          &" RSMORZA = ",F7.4,1x," NSMORZAO = ",I7,1x, &
          &" NSMORZA = ",I7)') deltadiff*6d7,rsmorza,nsmorzao,nsmorza
     if(nscrivi_ema == 1) then
        write(*,'("   Manca envelope esterno convettivo")')
        write(*,'("      Vedi ciacio")') 
     endif
  endif

  do i=1,M-1
     do k=1,maxme-1
        vl(k) = vel(i,k)
        dm_dr(k) = emu(i,k)
     enddo

     call cambia(dm_dr, dm_dr_v, maxme)

     ! controllo eventuali perdite numeriche di massa
     if(mod(nmd, moda1) == 0. .or. nmd+1 == nsole) then
        call check(nmd-1, maxme, i, tempomod, dm_dr_v, emutot, r, vl, &
             dv_dr, m)
     endif

     call integra2(dm_dr_v, dm_dr_n, deltadiff, maxme, r, vl, dv_dr)
     ! qui ho i nuovi valori di dm; posso farci i  check e scrivere
     call cambia(dm_dr_n, dm_dr, maxme)

     call decidi(nmd, maxme, xsave, nsole, i, tempomod, dm_dr, emutot,&
          r, vl, dv_dr, m, ele, ierr, nsup)
  end do

end subroutine get_new_abundances


real function f(y)
  !       questa serve per correggere l'andamento del campo di
  !       velocita' fra r=0.6 ed r(maxme-1) ;
  !       la funzione f vale zero per y=0, vale 1 per y=-.1 e
  !       la sua derivata e' nulla in questo punto

  implicit none
  real :: y
  f = -100.0*y*y-20.0*y
!  f = exp(-150*y-7.5)/(1+exp(-150*y-7.5))
  return
end function f


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




subroutine cambia(X,Y, maxme)

  implicit none

  integer :: maxme
  real,dimension(maxme) :: X,Y
  integer :: i

  do i=1,maxme
     y(i) = X(i)
  end do

  return
end subroutine cambia


subroutine integra(dm_dr_v,dm_dr_n,deltadiff,maxme, r, vl, dv_dr)

  implicit none

  real,dimension(maxme) :: dm_dr_v, dm_dr_n, r, vl, dv_dr
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


subroutine integra2(dm_dr_v,dm_dr_n,deltadiff,maxme, r, vl, dv_dr)
  use parametri
  use costanti

  implicit none

  real,dimension(maxme) :: dm_dr_v, dm_dr_n, r, vl, ef, def_dr, dv_dr
  real :: deltadiff
  integer :: maxme

  integer,parameter :: num = 4
  integer :: i, numnum
  real, dimension(num+1) :: roxv, droxv_dr
  real, parameter :: cost_4pi = 4.0*pigre

  numnum = num + 1
  roxv(1) = 0.
  do i=2,numnum
     roxv(i) = ( dm_dr_v(i)*vl(i) ) / ( cost_4pi*r(i)**2 )
  end do
  call derivami(numnum,r,roxv,droxv_dr)

  dm_dr_n(1) = dm_dr_v(1) - 3.*deltadiff*droxv_dr(1)
  do i=2,num
     dm_dr_n(i) = dm_dr_v(i) - 3.*deltadiff*cost_4pi*r(i)**2*droxv_dr(i)
  end do

  ef(1) = 0.
  do i=2,maxme-1
     ef(i) = dm_dr_v(i) * vl(i)
  end do
  call derivami(maxme-1,r,ef,def_dr)

  do i=numnum,maxme-1
     dm_dr_n(i) = dm_dr_v(i) - deltadiff*def_dr(i)
  end do

end subroutine integra2


subroutine decidi(nmd,maxme,xsave, nsole, nscr, tempomod, dm_dr, emutot, &
     r, vl, dv_dr, m, ele, ierr, nsup)
  use interfaccia
  use parametri
  use chimic
  use fisica
  use strut

  implicit none

  integer :: nmd, maxme, nsole, nscr, m, ierr, nsup
  real :: tempomod
  real,dimension(maxme) :: xsave, dm_dr, emutot, r, vl, dv_dr
  integer,dimension(m) :: ele

  ! si ferma l'integrazione dopo n passi e si chiama questa subroutine
  real,dimension(m-1,maxme) ::  xnew

  integer :: j, k, i

  ierr = 0

  ! calcolo delle nuove abbondanze chimiche dalle dm_dr:
  if(nsup > maxme-1) nsup = maxme-1
!  do k=2,maxme-1
  do k=2,nsup
     xnew(nscr,k) = dm_dr(k)/emutot(k)
  end do
  xnew(nscr,1) = xnew(nscr,2)

  

  ! scrivo le nuove abbondanze nel modulo CHIMIC:
!  do k=1,maxme-1
   do k=1,nsup
     xxx(ele(nscr),k) = xnew(nscr,k)
     if (nscr == 1 .and. xnew(nscr,k) < 1.d-39) then
        xxx(nscr,k) = 0.
     endif
  end do

  ! alla fine:  (sull'elemento m-1)
  if(nscr == m-1) then
     ! calcolo l'elio come differenza (per la conservazione della massa)
     do k=1,maxme-1
        xxx(3,k) = xsave(k)-xxx(1,k)-xxx(4,k)-xxx(5,k) &
             -xxx(6,k)-xxx(21,k)-xxx(22,k)-xxx(23,k)-xxx(24,k)-xxx(25,k)
     end do
     
     do j=1,mele
        !       do k=1,maxme-1
        do k=1,nsup
           if (xxx(j,k) > 1.) then
              write(*,*) 'ELEMENTO n.',j,'> 1 !!!', xxx(j,k)
              write(*,*) 'al mesh',k
              write(*,*) 'coordinata in massa: ', G(5,k)/emtot
              !****** si ferma quando fa casino *****************
              write(*,*)'Mi FERMO la ciacio ha fatto casino...'
              write(66,*)'230 - decidi'
              write(66,*) 'ELEMENTO n.',j,'> 1 !!!'
              write(66,*) 'al mesh',k
              write(66,*) 'coordinata in massa: ', G(5,k)/emtot
              ierr = 1
              return
!!$              stop   
           endif
           
           ! controllo < 0 per H e He4
!!$           if((j == 1 .or. j == 3) .and. xxx(j,k) < 0. ) then
           if (xxx(j,k) < 0. ) then
              write(*,*) 'ELEMENTO n.',j,'< 0 !!!', xxx(j,k) 
              write(*,*) 'al mesh',k
              write(*,*) 'coordinata in massa: ', G(5,k)/emtot
              !****** ... lo faccio fermare ************
              write(*,*)'Mi FERMO la ciacio ha fatto casino...'
              write(66,*)'231 - decidi'
              write(66,*) 'ELEMENTO n.',j,'< 0 !!!' 
              write(66,*) 'al mesh',k
              write(66,*) 'coordinata in massa: ', G(5,k)/emtot
              ierr = 1
              return
!!$              stop
           endif
        end do
     end do

  endif

  if(nmd+1 == nsole) then
     do i=1,maxme-1
        write(nscr+69,101) nmd+1,r(i),dm_dr(i),xxx(ele(nscr),i)
     end do

     call check(nmd,maxme, nscr, tempomod, dm_dr, emutot, r, vl, dv_dr, m)
  endif
101 format(1x,i6,3(1x,1pd15.6))

  return
end subroutine decidi


subroutine check(nmd,maxme, nscr, tempomod, dm_dr_v, emutot, r, vl, dv_dr, m)
  use parametri

  implicit none

  integer :: nmd, maxme, nscr, m
  real :: tempomod
  real,dimension(maxme) :: dm_dr_v, emutot, r, vl, dv_dr

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
  real,dimension(maxme) :: tuttamu, eragg
  real, dimension(m) :: X, sum1, sum_pro
  real,dimension(m,maxme) :: emu
  integer :: k, i, kk, mai
  real :: t9, tm1, pr1, roo1, rrho1, sup, dra
  real :: sumtot, emlim

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="ciacio        "


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


subroutine derivami(max, x, y, dy)

  implicit none

  integer, intent(in) :: max
  real, dimension(max), intent(in) :: x, y
  real, dimension(max), intent(out) :: dy

  ! Variabili locali
  integer :: i

  if ( max < 3 ) stop 'IMPOSSIBILE DERIVARE CON derivami: POCHI PUNTI!'

  dy(1) = ((y(2)-y(1))*(x(3)-x(1))*(x(3)-x(1)) - &
       (y(3)-y(1))*(x(2)-x(1))*(x(2)-x(1))) /  &
       ((x(2)-x(1))*(x(3)-x(1))*(x(3)-x(2)))

  do i=2,max-1
     dy(i) = ((y(i)-y(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i)) + &
          (y(i+1)-y(i))*(x(i)-x(i-1))*(x(i)-x(i-1))) / &
          ((x(i+1)-x(i))*(x(i)-x(i-1))*(x(i+1)-x(i-1)))
  end do

  dy(max) = ((y(max)-y(max-1))*(x(max)-x(max-2))*(x(max)-x(max-2)) - &
       (y(max)-y(max-2))*(x(max)-x(max-1))*(x(max)-x(max-1))) /      &
       ((x(max)-x(max-1))*(x(max)-x(max-2))*(x(max-1)-x(max-2)))

end subroutine derivami


subroutine init_elements(maxme, M, ele, nfine_diff, XXXsave)
  use parametri
  use chimic
  use nummod
  
  implicit none
  integer :: maxme, M, nfine_diff

  integer,dimension(M) :: ele

  integer :: k, i
  real,dimension(mele,maxme) :: XXXsave
  integer,save :: firstlog

  ! introduco le specie considerate delle XXX(mele,k) cioe':
  ! qelle che hanno 1 in sel_ele
  ! non metto il ferro in TBL (omogeneo):
  !          X(ferro)=56.*0.698*10**(-5.18)
  ! ma metto:
  !          X(6)=56. *0.698*10**(-4.49)   (vedi il soleciacio)

  ! calcolo la matrice dei logaritmi delle X(I) usate:

  ! Controllo sulla ciacio per vedere se succede del casino 
  ! Se qualche elemento sborda, o diventa negativo, fermo la ciacio 

  ! salvataggio abbondanze originarie in ingresso
  XXXsave(:,1:maxme) = XXX(:,1:maxme)

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
                 XXXsave(ele(i), k) = 1d-40
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
                 XXXsave(ele(i), k) = 1d-40
              endif
           endif
        end if
     end do
  end do

end subroutine init_elements


subroutine grad_chim(maxme, m, elxxx, r, grx, grax)

  implicit none

  integer :: maxme, m 
  real,dimension(m-1,maxme) :: ELXXX, GRAX

  integer :: i, k
  real,dimension(maxme) :: exx
  real,dimension(maxme) :: r, grx

  do i=1,m-1
     do k=1,maxme-1
        exx(k) = elxxx(i,k)
     end do

     call derivami(maxme-1,r,exx,grx)
!     call deriva(maxme-1,r,exx,grx,2)
     do k=1,maxme-1
        grax(i,k) = grx(k)
     end do

  end do
  return
end subroutine grad_chim


subroutine get_diffusion_v(maxme, M, Z, A, dens, tem, r, grax, grap, grat, &
     ele, XXXsave, vel, emu, emutot, xsave, vgt, S)
  use parametri
  use costanti
  
  implicit none

  integer :: maxme, M, k
  real,dimension(maxme) :: dens, tem, r
  real,dimension(M) :: Z, A
  real,dimension(m-1,maxme) :: GRAX
  real,dimension(maxme) :: GRAP, GRAT, emutot
  integer,dimension(M) :: ele
  real,dimension(mele,maxme) :: XXXsave
  real,dimension(m,maxme) :: emu, vel, vgt
  real,dimension(m,m,maxme) :: S
  
  integer :: mesh, i, j
  real,dimension(maxme) :: xsave
  real,dimension(m) :: X, AP, AT, C, XI, XInchim
  real,dimension(m, m) :: AX, CL
  real :: ZXA, AC, NI, CZ, XIJ, NE, AO, LAMBDAD, LAMBDA, ax_sum
  real :: tmp, TEMP1, TEMP2, SUM1, SUM2
  real :: TT, RRho, fatt, sup, RHO, T, dv
  real,save,dimension(mele,mele) :: powZ
  integer :: im
  
  real, parameter ::  SUNR = 6.9599d10, TAU = 6.0d13
  
  rewind(126)
  rewind(127)

  !*********************************************************************
  do MESH=1, MAXME-1

     RHO = exp(dens(MESH))
     T = exp(tem(MESH))
     
     !     Nel seguito si calcolano i logaritmi coulombiani
     !     Si definiscono i valori di temperatura e densita
     !     per cui calcolare la densita elettronica (unita cgs)

     !   si definiscono i valori di X da usare tra quelli letti
     do i=1,m-1
        x(i) = xxxsave(ele(i),mesh)
        ! controllo per H e elem. legg.
        if(i == 1 .or. (ele(i) >= 22 .and. ele(i) <= 25)) then
           if(x(i) <= 0.) x(i) = 1.d-40
        endif
     end do
     
     xsave(mesh) = sum(X(1:(m-1)))

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
     do J=1,M
        do I=1,M     
           powZ(i,j) = abs(Z(I)*Z(J))**1.2
        end do
     end do
     
     tmp = (2.3939d3*T*LAMBDA)**1.2
     do J=1,M
        do I=1,M
           XIJ = tmp/powZ(i,j)
           CL(I,J) = 0.81245*log(1.+0.18769*XIJ)
        end do
     end do
     
     call DIFFUSION(M,A,Z,X,CL,AP,AT,AX)

     !  CALCOLO DELLE FUNZIONI DI DIFFUSIONE \XI(M):
 
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
        XInchim(I) = AP(I)*GRAP(MESH)+AT(I)*GRAT(MESH)
        XI(I) = XInchim(I)+SUM1-SUM2*TEMP2/TEMP1
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
     ! passo in cgs
     FATT = FATT * SUNR/ (TAU * sec_anno)

     do I=1,M
        ! attivo anche la prima
        if(MESH /= -1) then
           VEL(I,MESH) = XI(I)*FATT
           vgt(i,mesh) = XInchim(i)*FATT
        else
           VEL(I,MESH) = 0.
           vgt(i,mesh) = 0.
        endif
     enddo
     
     ! final fixup for velocity of most abundant to give exact local mass conservation
     im = maxloc(X(1:m-1),dim=1)
     dv = -dot_product(vel(1:m-1,mesh), x(1:m-1))/x(im)
     vgt(im, mesh) = vgt(im, mesh) + dv
     vel(im, mesh) = vel(im, mesh) + dv
     
!!$     write(126,'(2(i5,1x),20(e14.6,1x))') mesh, m, (vgt(i, mesh), i=1,m-1) 
!!$     write(127,'(2(i5,1x),20(e14.6,1x))') mesh, m, (vel(i, mesh), i=1,m-1)

     do I=1,M-1
        ax_sum = 0.
        do k = 1,m-1
           ax_sum = ax_sum+AX(i,k)
        enddo
        do J=1,M-1
           S(i,j,mesh) = FATT*SUNR*(-AX(i,j)+(Z(J)*X(J)/(A(J)*TEMP1))*ax_sum)
        end do
     end do
     
     ! calcolo delle variabili mu totale e mu_i adimensionali:
!     if (MESH /= 1) then
        SUP = 4*pigre*r(MESH)**2
        emutot(MESH) = SUP * RRHO
 !    else
 !       emutot(MESH) = RRHO
 !    endif
     do j=1,m
        emu(j,MESH) = emutot(MESH) * X(j)
     end do
  end do
  return

end subroutine get_diffusion_v



subroutine deriva(max,x,y,dy,nsmo)
  use interfaccia
  use parametri
  use nummod

  implicit none

  integer :: max, nsmo
  real,dimension(max) :: x, y, dy

  real,dimension(max) :: xa, ya, dya, erin, upsilon

  integer :: k, ki, num 
  real :: delmax, dx, ypn, yp1, delta

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
  if (num > max) then
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

  if(nmd == 69) then
     write(56,'(i5,1x,12(e13.6,1x))') num, yp1, ypn, delta
  endif

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

  real :: b(7), yy(max)
  integer :: imed, kernel, mesh, k, j, i
  real :: ba

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

                   

                   
subroutine get_matrix_4(m, maxme, dim, r, rho, vgt, sigma, dt, MA)
  use interfaccia
  use fisica
  use strut
  use chimic
  use tempe
  use chim
  use nummod
  use costanti
  
  implicit none

  integer :: m, maxme, dim
  real :: dt
  real, dimension(maxme) :: r, rho
  real, dimension(m,maxme) :: vgt
  real,dimension(m,m,maxme) :: sigma

  integer :: i, j, k, m1
  real, dimension(dim, 4*(m-1)-1) :: MA
  real, dimension(m-1, maxme-1) :: GT
  real, dimension(maxme-1) :: cell_dm, r2_mid

 
  real, parameter :: SUNM = 1.989d33, SUNR = 6.9599d10
  logical, parameter :: armonic_vel = .false.
 
  real, dimension(m-1,m-1,maxme-1) :: e00, ep1, em1

  !real, dimension(:,:,:),allocatable :: e00, ep1, em1 
  !allocate(e00(m-1,m-1,maxme-1), ep1(m-1,m-1, maxme-1), em1(m-1,m-1,maxme-1))

  do k = 1, maxme-1
     ! grandezze in cgs
     r2_mid(k) = (((SUNR*r(k))**3 + (SUNR*r(k+1))**3)/2.)**(2./3.)
!     rmid(k) = (((SUNR*r(k))**3 + (SUNR*r(k+1))**3)/2.)**(1./3.)
!     cell_dm(k) = (G(5,k+1) - G(5, k))*1.d33
     cell_dm(k) = 4./3.*pigre*rho(k)*SUNR**3*(r(k+1)**3 - r(k)**3)

!     if(k>1) cell_dm(k) = cell_dm(k) + 0.5*(G(5,k) - G(5, k-1))*1d33
!     if(k<maxme-1) cell_dm(k) =  cell_dm(k)+ 0.5*(G(5,k+2) - G(5, k+1))*1d33
  end do
  
  do k = 1, maxme-2
     do i = 1, m-1
        if(armonic_vel) then
           ! media armonica vel.
           GT(i,k) = -dt*pigre*2.*(rho(k)+ rho(k+1))*r2_mid(k) &
                * (vgt(i,k)*vgt(i,k+1))/(vgt(i,k) + vgt(i,k+1))
        else
           ! media aritmetica vel.
           GT(i,k) = -dt*pigre/2.*(rho(k)+ rho(k+1))* r2_mid(k) &
                * (vgt(i,k) + vgt(i,k+1))
        endif
     end do
  end do
  GT(:, maxme-1) = GT(:, maxme-2)
 
  e00 = 0.0
  ep1 = 0.0
  em1 = 0.0

  do k=1, maxme-1
     do i=1,m-1
        e00(i,i,k) = 1.0
     end do
     
     if(k > 1) then
        do i=1,m-1
           e00(i,i,k) = e00(i,i,k) + GT(i, k-1)/cell_dm(k)
           em1(i,i,k) = GT(i, k-1)/cell_dm(k)
        end do
     endif
     
     if(k < maxme-1) then
        do i=1,m-1
           e00(i,i,k) = e00(i,i,k) - GT(i, k)/cell_dm(k)
           ep1(i,i,k) = - GT(i, k)/cell_dm(k)
        end do
     endif
     
  end do

  ! riempio la matrice a banda
  m1 = 2*(m-1)-1
  MA = 0.0
  i = 1
  do k=1,maxme-1
     do j = 1, m-1
        MA(i,m1+1) = e00(j,j,k)
        MA(i,m1+1 + (m-1)) = ep1(j,j,k)
        MA(i,m1+1 - (m-1)) = em1(j,j,k)
        i = i +1
     end do
  end do
  
  if(nmd == 500) then
     do k=1,maxme
        write(137,'(i5,1x,20(e15.7,1x))') k, G(5,k), GT(1,k), cell_dm(k), r2_mid(k), (rho(k)+ rho(k+1))/2., (vgt(i,k) + vgt(i,k+1))/2.
     end do
  endif

  return
end subroutine get_matrix_4






subroutine bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
  implicit none

  integer :: m1,m2,mp,mpl,n,np,indx(n)
  real :: d,a(np,mp),al(np,mpl)
  real,parameter :: myTINY=1.d-20
  
!!$  Given an n × n band diagonal matrix A with m1 subdiagonal rows and m2 superdiagonal
!!$  rows, compactly stored in the array a(1:n,1:m1+m2+1) as described in the comment for
!!$  routine banmul, this routine constructs an LU decomposition of a rowwise permutation
!!$  of A. The upper triangular matrix replaces a, while the lower triangular matrix is returned
!!$  in al(1:n,1:m1). indx(1:n) is an output vector which records the row permutation
!!$  eﬀected by the partial pivoting; d is output as ±1 depending on whether the number of
!!$  row interchanges was even or odd, respectively. This routine is used in combination with
!!$  banbks to solve band-diagonal sets of equations.
  integer :: i,j,k,l,mm
  real :: dum
  
  mm = m1+m2+1
  if(mm > mp .or. m1 > mpl .or. n > np) then
     write(*,*) 'bad args in bandec'
     stop
  endif
  l = m1
  do i=1,m1
!!$     Rearrange the storage a bit.
     do j=m1+2-i,mm
        a(i,j-l) = a(i,j)
     enddo
     l=l-1
     do j=mm-l,mm
        a(i,j) = 0.
     enddo
  enddo
  d = 1.
  l = m1
  do k=1,n
!!$     For each row...
     dum = a(k,1)
     i = k
     if(l < n) l = l+1
     do j=k+1,l
!!$        Find the pivot element.
        if(abs(a(j,1)) > abs(dum)) then
           dum = a(j,1)
           i = j
        endif
     enddo
     indx(k) = i
     if(dum == 0.) a(k,1) = myTINY
!!$     Matrix is algorithmically singular, but proceed anyway with TINY pivot (desirable in some applications).
     if(i /= k) then
!!$        Interchange rows.
        d = -d
        do j=1,mm
           dum = a(k,j)
           a(k,j) = a(i,j)
           a(i,j) = dum
        enddo
     endif
     do i=k+1,l
!!$        Do the elimination.
           dum = a(i,1)/a(k,1)
           al(k,i-k) = dum
           do j=2,mm
              a(i,j-1) = a(i,j)-dum*a(k,j)
           enddo
           a(i,mm) = 0.
        enddo
     enddo
     
  return
end subroutine bandec


subroutine banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
  implicit none

  integer :: m1,m2,mp,mpl,n,np,indx(n)
  real :: a(np,mp),al(np,mpl),b(n)
  !!$Given the arrays a, al, and indx as returned from bandec, and given a right-hand side
  !!$vector b(1:n), solves the band diagonal linear equations A · x = b. The solution vector x
  !!$overwrites b(1:n). The other input arrays are not modified, and can be left in place for
  !!$successive calls with different right-hand sides.
  integer :: i,k,l,mm
  real :: dum

  mm = m1+m2+1
  if(mm > mp .or. m1 > mpl .or. n > np) then
     write(*,*) 'bad args in banbks'
     stop
  endif

  l = m1
  do k=1,n
!!$     Forward substitution, unscrambling the permuted rows as we go.
     i = indx(k)
     if(i /= k) then
        dum = b(k)
        b(k) = b(i)
        b(i) = dum
     endif
     if(l < n) l = l+1
     do i=k+1,l
        b(i) = b(i)-al(k,i-k)*b(k)
     end do
  enddo
  l = 1
  do i=n,1,-1
!!$     Backsubstitution.
     dum = b(i)
     do k=2,l
        dum = dum-a(i,k)*b(k+i-1)
     end do
     b(i) = dum/a(i,1)
     if(l < mm) l = l+1
  end do
  
  return
end subroutine banbks




subroutine banmul(a,n,m1,m2,np,mp,x,b)
  integer ::m1,m2,mp,n,np
  real a(np,mp),b(n),x(n)
!!$Matrix multiply b = A x, where A is band diagonal with m1 rows below the
!!$diagonal and m2 rows above. The input vector x and output vector b are stored
!!$as x(1:n) and b(1:n), respectively. The array a(1:n,1:m1+m2+1) stores A as
!!$follows: The diagonal elements are in a(1:n,m1+1). Subdiagonal elements are in
!!$a(j:n,1:m1) (with j > 1 appropriate to the number of elements on each
!!$subdiagonal). Superdiagonal elements are in a(1:j,m1+2:m1+m2+1) with j < n
!!$appropriate to the number of elements on each superdiagonal.
  integer i,j,k
  
  do i=1,n
     b(i)=0.
     k=i-m1-1
     do j = max(1,1-k),min(m1+m2+1,n-k)
        b(i) = b(i)+a(i,j)*x(j+k)
     enddo
  end do
  return
end subroutine banmul


subroutine solve_system(maxme, a_matrix, m, dim, m1, m2, eleX)

  implicit none

  integer :: m, m1, m2, dim, maxme
  real ::  a_matrix(dim,1+m1+m2), eleX((m-1)*(maxme-1))

  real :: meleX((m-1)*(maxme-1))
  real :: d
  integer, dimension((m-1)*(maxme-1)) :: INDX
  integer :: i, k

  logical,parameter :: do_matrix_mult = .true.
  real, dimension(:,:),allocatable :: diag, al

  allocate(diag(dim,1+m1+m2), al(dim,1+m1+m2))

  INDX = 0

  if (do_matrix_mult) then
     diag = a_matrix
     diag(:,1+m1) =  diag(:,1+m1) - 1.0
     call banmul(diag, dim, m1, m2, dim, 1+m1+m2, eleX, meleX)
     meleX = -meleX
  endif

  call bandec(a_matrix,dim,m1,m2,dim,1+m1+m2,&
       al, m1, indx, d)

  if (do_matrix_mult) then
     call banbks(a_matrix,dim,m1,m2,dim,1+m1+m2, &
          al, m1, indx, meleX)
     eleX = eleX + meleX
  else
     call banbks(a_matrix,dim,m1,m2,dim,1+m1+m2, &
          al, m1, indx, eleX)
  endif


end subroutine solve_system


subroutine solve_system_L(maxme, a_matrix, m, dim, m1, m2, eleX)

  implicit none

  integer :: m, m1, m2, dim, maxme
  real ::  a_matrix(dim,1+m1+m2), eleX((m-1)*(maxme-1))

  integer :: ldab
  integer :: n, info, k, i, j
  integer,dimension((m-1)*(maxme-1)) :: ipiv
  real,dimension(2*m1+m2+1, (m-1)*(maxme-1)) :: a_matrix_aug

  n = (m-1)*(maxme-1)
  ldab = 2*m1+m2+1

  k = m1 + m2 + 1
  
  do i=1,n
     do j=max(i-m1,1),min(i+m2,n)
        a_matrix_aug(k+i-j,j) = a_matrix(i,j-i+m1+1)
     end do
  end do 

  ! LAPACK
  call DGBSV(n, m1, m2, 1, a_matrix_aug, ldab, ipiv, eleX, n, INFO)

end subroutine solve_system_L


subroutine rinormalizza_abbondanze(m, maxme, eleX, xsave)
  implicit none

  integer :: m, maxme
  real,dimension(maxme) :: xsave
  real, dimension((m-1)*(maxme-1)) :: eleX

  integer :: k, i, im1, im2, im
  real :: somma

  do k=1, maxme-1
     somma = 0.0
     do i = 1, m-1
        somma = somma + eleX((m-1)*(k-1)+i)
     end do
     im1 = (m-1)*(k-1) + 1
     im2 = (m-1)*(k-1) + m-1
     im = maxloc(eleX(im1:im2),dim=1)
     eleX((m-1)*(k-1) + im) = eleX((m-1)*(k-1) + im) + xsave(k) - somma
  end do

end subroutine rinormalizza_abbondanze



subroutine insert_turb(n_el, n_point, dens, tem, vel_diff, grax, vel_correct)
  
  use fisica
  use nummod

  implicit none

  integer :: n_el, n_point
  real, dimension(n_point) :: dens, tem
  real, dimension(n_el,n_point) :: vel_diff, vel_correct
  real, dimension(n_el-1,n_point-1) :: grax

  integer, parameter :: n=3
  real, parameter :: T_0=7.d5, f=400.

  integer :: i, k
  real :: RHO_0, D_He_approx, RHO
  
  if (exp(tem(n_point)) .gt. T_0) then
     write(*,*) "attenzione! Il punto esterno e` piu` caldo di T_0!"
     stop
  else
     do k=1,n_point
        if (exp(tem(k)) .lt. T_0) exit
     enddo
     RHO_0=exp(dens(k-1))+((exp(dens(k))-exp(dens(k-1)))/(exp(tem(k))-exp(tem(k-1))))*(T_0-exp(tem(k-1)))
  endif
  D_He_approx=3.3d-15*(T_0**2.5)/(4.*RHO_0*log(1.+(1.125d-16)*(T_0**3)/RHO_0))
  do k=1,n_point-1 
     if (k==1) then
        do i=1,n_el
           vel_correct(i,k)=0.
        enddo
     else
        RHO=exp(dens(k))     
        do i=1,n_el
           write(155,*) nmd, i,k, f*(grax(i,k)/6.9599d10)*D_He_approx*((RHO_0/RHO)**n)/vel_diff(i,k)
           vel_correct(i,k) = vel_diff(i,k)-f*(grax(i,k)/6.9599d10)*D_He_approx*((RHO_0/RHO)**n)
        enddo
     endif
  enddo
  return
  
end subroutine insert_turb
