subroutine stop_evolut(fase, lumin, nmd, age, logCNO, start_bump, teff_in, &
     teff_target, XHei, shellHoff, nsole)
  use interfaccia
  use chimic
  use varie
  use strut
  use fasievolut ! qui sono definite le fasi evolutive

  implicit none

  integer :: fase, nmd, logCNO, start_bump, shellHoff, nsole
  type(luminosity) :: lumin
  real :: age, teff_in, teff_target, XHei

  real,save :: old_L = 0, age_old
  real,save,dimension(10) :: dL, Lum
  real :: d_age
  integer :: i
  logical :: test

  ! metto da parte Teff iniziale (al modello 20)
  if(nmd == 20) teff_in = lumin%teff

  if((fase == fase_exH .or. fase == fase_RGB .or. fase == fase_bump1) .and. logCNO >= 0) then
     ! controllo per bloccare C.N,O per pepper
     if( lumin%L3a > 1.0d-8 ) logCNO = 1
     if( lumin%L3a > 1.0d-6 ) logCNO = -1
  endif

  if((fase == fase_exH .or. fase == fase_RGB) .and. start_bump >= 0) then
     if(lumin%nuclear > 10. .and. start_bump == 0 ) then
        start_bump = 1
        fase = fase_bump1
        write(68,'(i5,1x,i5)') fase_bump1, nmd
         write(67,'(" # ",a32,1x,i5)') "Bump RGB al modello:", nmd 
     endif
     if(lumin%nuclear < 1. .and. start_bump == 1) then
        start_bump = -1
        write(68,'(i5,1x,i5)') fase_bump2, nmd
     endif
  endif

  if( fase == 0 ) then          ! controllo per MS
     if(XXX(1,1) < 0.99*XH ) then
        fase = fase_MS

        write(67,'(" # ",a32,1x,i5)') "Inizio MS al modello:", nmd 
        write(68,'(i5,1x,i5)') fase, nmd
     endif
  else if(fase == fase_MS .or. fase == fase_sole) then      
     ! controllo esaustione H centrale
     if(XXX(1,1) < 1.0d-30 ) then
        fase = fase_exH
        ! questo valore lo uso per la ricerca di RGB
        teff_target = teff_in + 0.90*(lumin%teff - teff_in)
        write(67,'(" # ",a32,1x,i5)') "Esaurimento Hc al modello:", nmd 
        write(68,'(i5,1x,i5)') fase, nmd
        XHei = XXX(3,1)
     endif

  else if(fase == fase_exH .and. lumin%teff < teff_target ) then
     ! Cerco l'inizio di RGB:
     ! Trovo il punto in cui la L corrente e' superiore a 
     ! quella precedente, e cosi' via fino a 4 passi indietro
     do i=10,2,-1
        dL(i) = dl(i-1)
     end do
     dL(1) = lumin%L - old_L
     old_L = lumin%L

     ! test L
     test = (dL(1) > dL(2) .and. dL(2) > dL(3) .and. dL(3) > dL(4) .and. &
          dL(4) > dL(5))
     ! test Lgra
     if(test) then
        fase = fase_RGB 
        write(68,'(i5,1x,i5)') fase, nmd
        write(67,'(" # ",a32,1x,i5)') "Inizio RGB al modello:", nmd 
        write(67,*)
     endif
  endif

  ! Trovato il sole 
  if(nmd == nsole) then
     fase = fase_Sole
  end if

  if(fase == fase_RGB .or. fase == fase_exH .or. fase == fase_bump1) then   
     ! controllo flash He
     if(lumin%L3a/(lumin%Lpp+lumin%Lcno) > 4.0d1 .and. lumin%L3a > 4.0d1) then
        fase = fase_flashHe
        write(68,'(i5,1x,i5)') fase, nmd
        write(67,*) "==== STOP SIMULAZIONE ====" 
        write(67,'(" # ",a32,1x,i5)') "He flash al modello:", nmd 
        ! per aggiustare l'ultima massa in modo esatto alla ripartenza
        open(64, file="lastmass", status="unknown")
        write(64,*) emtot/1.989
        close(64)
       endif
       
       ! controlo innesco He per stelle che non fanno il flash
       if(XXX(3,1) < 0.99*XHei) then
          fase = fase_innHe
          write(67,'(" # ",a32,1x,i5)') "Innesco He al modello:", nmd 
          write(68,'(i5,1x,i5)') fase, nmd
       endif
    endif

  ! per stelle che flashano si arriva al massimo al flash 
  ! poi riparte con fase fase_postPepper
  if(fase == fase_postPepper .or. fase == fase_innHe) then 
     ! controllo esaustione He centrale
     if(XXX(3,1) < 1.0d-30 ) then
        fase = fase_exHe
        write(68,'(i5,1x,i5)') fase, nmd
        write(67,'(" # ",a32,1x,i5)') "Esaurimento He_c al modello:", nmd 
     endif
  else if(fase == fase_exHe) then
     ! cerco i pulsi

     ! questo va bene per masse piccole (sicuramente <= 1.10)
     if( lumin%L3a > 5.0 ) then 
        fase = fase_pulsiT
        write(68,'(i5,1x,i5)') fase, nmd
        write(67,'(" # ",a32,1x,i5)') "Pulsi termici al modello:", nmd 
        write(67,'(a50)') "Condizione su L3a > 5.0"
     end if

     ! questo va bene per masse che innescano il carbonio (probabilmente >= 6)
     if(lumin%pn(1) > 0 ) then
        fase = fase_pulsiT
        write(68,'(i5,1x,i5)') fase, nmd
        write(67,'(" # ",a32,1x,i5)') "Pulsi termici al modello:", nmd 
        write(67,'(a50)') "Condizione su NUCL C > 0"
     end if

     ! questo va bene per le masse che pulsano quando si accende la
     ! shell H
     if(lumin%pn(3) > lumin%pn(2) .and. shellHoff == 1 ) then
        shellHoff = 0
        fase = fase_pulsiT
        write(68,'(i5,1x,i5)') fase, nmd
        write(67,'(" # ",a32,1x,i5)') "Pulsi termici al modello:", nmd 
        write(67,'(a50)') "Condizione su accensione shell H"
     end if
     
     ! la shell H e' spenta
     if(lumin%pn(3) < 1e-1) then
        shellHoff = 1
     endif

  endif
  
end subroutine stop_evolut
