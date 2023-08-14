subroutine FITTA(M,DAM,NUM1,NUM2,FRAT)
  use fisica
  use second
  use atm
  use tempe
  use chim
  use numer
 
  use nummod
  
  implicit none

  integer :: M, NUM1, NUM2
  real :: DAM, FRAT

  real,dimension(4) :: DA, DAV 
  integer :: k, l, i, j, it
  real :: ddr, ddp, ddt, den, tao, dao, cao, dtef, del1, dav3, dav4

  DDR = log(RB/G(1,M))
  DDP = log(PB/G(3,M))
  DDT = log(TB/G(4,M))
  DEN = 1.-BET(2,M)*DPL-GAM(2,M)*DTL 
  TAO = (BET(2,M)*DPT+GAM(2,M)*DTT)/DEN 
  DAO = (ALF(2,M)+BET(2,M)*DDP+GAM(2,M)*DDT)/DEN 
  CAO = BET(1,M)*(DPL*TAO+DPT)+GAM(1,M)*(DTL*TAO+DTT)-DRL*TAO-DRT 
  DTEF = -ALF(1,M)-DAO*(BET(1,M)*DPL+GAM(1,M)*DTL-DRL)+DDR-BET(1,M)*  &
       DDP-GAM(1,M)*DDT
  DTEF = DTEF/CAO 
  DEL1 = TAO*DTEF+DAO 
  TEFF = TEFF+DTEF*FRAT/2.30258509 
  
  DA(3) = DDP+DPL*DEL1+DPT*DTEF 
  DA(4) = DDT+DTL*DEL1+DTT*DTEF 
  DAM = 0. 
  do K=1,M
     L = M-K+1
     DA(1) = ALF(1,L)+BET(1,L)*DA(3)+GAM(1,L)*DA(4)
     DA(2) = ALF(2,L)+BET(2,L)*DA(3)+GAM(2,L)*DA(4)
     DAV3 = DA(3)
     DAV4 = DA(4)

     do i=1,4
        if(abs(DA(I))-abs(DAM) > 0) then 
           DAM = DA(I)*1.000001
           NUM1 = I
           NUM2 = L
        endif
     end do

     do J=1,4
        DA(J) = DA(J)*G(J,L)
        G(J,L) = G(J,L)+DA(J)*FRAT
        if(L == M) cycle
        DG(J,L) = DG(J,L)+DAV(J)*FRAT-DA(J)*FRAT
     end do

     if(IPRALL == 1) then
        write(2,112) (G(IT,L),IT=1,4)
        write(2,112) (DA(IT),IT=1,4)
        write(2,112) (DG(IT,L),IT=1,4)
     endif
     do I=1,4
        DAV(I) = DA(I)
     end do
     DA(3) = ALF(3,L)+BET(3,L)*DAV3+GAM(3,L)*DAV4
     DA(4) = ALF(4,L)+BET(4,L)*DAV3+GAM(4,L)*DAV4
  end do

  if(G(2,M) <= 0.0) then
     write(*,*) 'fitta : G(2,M) negativo'
     write(66,*) '170 - fitta'
     write(66,*) 'G(2,M) negativo'
     stop
  endif
  ELLOG = log10(G(2,M)/38.26)
  
  return 
112 format(1X,1P,4E14.5) 
end subroutine FITTA
