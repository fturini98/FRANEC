! La routine check_convettivo agisce sulle zone convettive
! cumulando i valori di yf e Srho, pesandoli con
! la massa del mesh (yfs e Ssum)
! Al termine della zona convettiva media tali valori.
! Queste quantita' possono essere usate per evolvere il blocco convettivo


subroutine check_convettivo(convettivo, jf, Ssum, yfs, Srho, yf, jfstart)
  use parametri
  use fisica
  use overshoot

  implicit none

  integer :: convettivo, jf, jfstart
  real,dimension(nsezi+2) :: SSum, Srho
  real,dimension(MELE) :: yfs, yf

  integer, save :: jfstop

  if (convettivo == 1) then
     ! inizio zona convettiva
     jfstart = jf
     SSum = 0
     yfs = 0

     ! gestisco il core convettivo
     if(jfstart > 1) then
        SSum = Srho*(G(5,JF+1)-G(5,JF))
        yfs = yf*(G(5,JF+1)-G(5,JF))
     else
        SSum = Srho*G(5,JF+1)
        yfs = yf*G(5,JF+1)
     endif
  else if (convettivo == 2) then
     ! interno zona convettiva
     Ssum = SSum + Srho*(G(5,JF+1)-G(5,JF))
     yfs = yfs + yf*(G(5,JF+1)-G(5,JF))
  else if (convettivo == 3) then
     ! termine zona convettiva
     jfstop = jf
     Ssum = SSum + Srho*(G(5,JF+1)-G(5,JF))
     yfs = yfs + yf*(G(5,JF+1)-G(5,JF))

     ! gestisco il core convettivo
     if(jfstart > 1) then
        SSum = Ssum/(G(5,jfstop+1) - G(5,jfstart))
        yfs = yfs/(G(5,jfstop+1) - G(5,jfstart))
     else
        SSum = Ssum/G(5,jfstop+1)
        yfs = yfs/G(5,jfstop+1)
     endif
  endif

end subroutine check_convettivo


