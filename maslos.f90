subroutine MASLOS(MAXME,DM,nmd,iread, fase, fase_start, fase_stop) 
  use interfaccia
  use fisica
  use strut
  use tempe
  use chim
  use varie
  use costanti
  implicit none

  integer :: MAXME, nmd, iread, fase, fase_start, fase_stop
  real :: DM

  integer :: LA, L
  real :: FINMAS, RAT

  
  if(ETAFPR <= 0.) return 
  if(iread == 4) then
     if(nmd <= 10) return   ! primi 10 modelli...
     DM = 0.001*ETAFPR*Msun
  else
     if( fase <= fase_start .and. fase >= fase_stop ) then
        DM = -5.3d-5*HT1*ETAFPR*(10.**(ELLOG*1.5 - 2.*TEFF))/EMTOT 
     else
        return
     endif
  endif
  EMTOT = EMTOT+DM 
  FINMAS = EMTOT*FRAZ 
  RAT = (FINMAS-G(5,MAXME-30))/(G(5,MAXME)-G(5,MAXME-30)) 
  LA = MAXME-29 
  do L=LA,MAXME 
     G(5,L) = G(5,MAXME-30)+RAT*(G(5,L)-G(5,MAXME-30)) 
  end do
  return 
end subroutine MASLOS
