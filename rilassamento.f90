subroutine rilassamento(MAXME, M, Y, fit, nstep)
  use chimic
  use varie
  use mistura
  use costanti

  implicit none
  
  integer :: nmd, maxme
  real :: M, Y

  integer,save :: firsttime = 1
  integer :: i
  real,save :: s_M, s_HELI, s_XME, s_XMEl, s_ALFA

  ! numero max di modelli in cui rilassare
  integer :: nstep 
  real, dimension(8) :: fit

  if(firsttime == 1) then
     s_M = (fit(1)*Msun - M)/nstep
     s_XMEl = (log10(fit(2)) - log10(XME))/nstep
     s_HELI = (fit(3) - Y)/nstep
     s_ALFA = (fit(4) - ALFA)/nstep

     firsttime = 0
  endif

  ! massa
  M = M + s_M

  ! Alpha
  ALFA = ALFA + s_ALFA

  ! He
  do i=1,maxme
     xxx(3,i) = xxx(3,i) + s_HELI
     xxx(1,i) = xxx(1,i) - s_HELI
  end do
  XH = XH - s_HELI

  ! z
  s_XME = 10.**(log10(XME) + s_XMEl) - XME
  do i=1,maxme
     xxx(1,i) = xxx(1,i) - s_XME
     xxx(2, i) = xxx(2,i) + s_XME*DEFAUHe3
     xxx(4,i) = xxx(4,i) + s_XME*DEFAUC
     xxx(5,i) = xxx(5,i) + s_XME*DEFAUN
     xxx(6,i) = xxx(6,i) + s_XME*DEFAUO
     xxx(21,i) = xxx(21,i) + s_XME*DEFAUFe
     xxx(22,i) = xxx(22,i) + s_XME*DEFAULi6 
     xxx(23,i) = xxx(23,i) + s_XME*DEFAULi7 
     xxx(24,i) = xxx(24,i) + s_XME*DEFAUBe 
     xxx(25,i) = xxx(25,i) + s_XME*DEFAUB
  end do
  XME = XME + s_XME
  XH = XH - s_XME

end subroutine rilassamento


subroutine rilassamento_subatm(nmd, fraz, ht1, emtot, maxme)
  use fisica

  implicit none

  integer :: nmd, maxme
  real :: fraz, ht1, emtot

  real :: mult, del
  real :: TARGET_FRAZ

  TARGET_FRAZ = 1e-5

  ! soppressione progressiva della subatmosfera
  ! la scelta dei valori soglia e' puramente empirica e non ottimizzata...
  if(nmd > 3 .and. nmd < 50 .and. fraz < 1.0) then
     mult = 0.95
     if(nmd > 40) then
        mult = 0.5
     else if(nmd > 30) then
        mult = 0.7
     else if(nmd > 20) then
        mult = 0.8 
     else if(nmd > 12) then
        mult = 0.9
     endif

     ! riporto la frazione eliminata nell'ultimo mesh
     if(1 - fraz > TARGET_FRAZ) then
        del = (1 - fraz)*mult
        fraz = 1 - del
        G(5, maxme) = EMTOT*fraz
     endif
  endif

end subroutine rilassamento_subatm
