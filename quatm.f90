subroutine QUATM(NABLA,INI,MAXMEin,MAXMV,EMTOV) 
  use interfaccia
  use fisica
  use gratm
  use strut
  use mesh
  use preco
  use atm
  use tempe
  use numer
  use varie
  use costanti
  implicit none

  integer :: NABLA, INI, MAXMEin, MAXMV
  real :: EMTOV

  real,save :: SERV(4),TT1(2),EL1(2),C1(4),C2(4),AT(4,2,2),PTR(3,4,4)
  real,save,dimension(3,4) :: alfac, betac, gammac, deltac
  real,parameter :: RAPPO = 0.10, RAPMA = 1.20, DRATT = 0.0005, DRAPP = 0.005
  real, parameter :: RATTO = .5
  integer, parameter :: LALLA = 1 
  real,dimension(4) :: tmp1, tmp2, A,B

  integer :: num , icont, ira, nana, incas, jl, jt, i, it, j, k, l, ite, ilu
  real :: emtut, pd, td, rd, pm11, pm22, vartmp, ddratt, ddrapp
  real :: alfa1, beta1, gamma1, delta1, ratt, ratl, aluce, ema

  NUM = 0 
  ICONT = 0 
  IRA = 0 
  NANA = 0 
  INCAS = 0 
  MAXME = MAXMEin 
  MAXMIV = MAXMV 

  !*****SE DEVE CALCOLARE UNA SOLA ATM. O DEVE PRINTARE *********
  if(ISUB == 1 .or. NABLA == 3) then
     EMTUT = EMTOT 
     call ATMOS(ELLOG,TEFF,PD,TD,RD,EMTUT,1,IRA) 
     PB = PD 
     TB = TD 
     RB = RD 
     return 
  endif

  !*************** SINGOLA ATMOSFERA *************************            
  EL(1) = ELLOG 
  EL(2) = EL(1)+0.0001 
  TE(1) = TEFF 
  TE(2) = TE(1)+0.0001 
  PMA1 = EMTOT 
  PMA2 = PMA1*0.99 
  
  do JL=1,2 
     do JT=1,2 
        if(JL == 2 .and. JT == 2) cycle
        call ATMOS(EL(JL),TE(JT),U(1,JL,JT),U(2,JL,JT),U(3,JL,JT),PMA1,0,IRA)
        if(IRA == 1) goto 38 
        U(3,JL,JT) = log10(U(3,JL,JT)*sqrt(1.d1*Lsun*(10.**EL(JL))/ &
             (7.125d-16*(10.**(4.*TE(JT)))))) 
     end do
  end do
  PM11 = PMA1/Msun
  PB = 10.**(U(1,1,1)-17) 
  TB = 10.**(U(2,1,1)-6) 
  RB = 10.**U(3,1,1) 
  DPL = (U(1,2,1)-U(1,1,1))/0.0001 
  DPT = (U(1,1,2)-U(1,1,1))/0.0001 
  DTL = (U(2,2,1)-U(2,1,1))/0.0001 
  DTT = (U(2,1,2)-U(2,1,1))/0.0001 
  DRL = (U(3,2,1)-U(3,1,1))/0.0001 
  DRT = (U(3,1,2)-U(3,1,1))/0.0001 
  return 

  !*********************************************************************  
39 continue 
  PMA1 = EMTOT 
  PMA2 = PMA1*0.99 
  TE(1) = 3.560 
  TE(2) = 3.580 
  TE(3) = 3.600 
  TE(4) = 3.620 
  EL(1) = 4.51 
  EL(2) = 4.53 
  EL(3) = 4.55 
  EL(4) = 4.57 

  ! giada
  do i=1,4
     tmp1(i) = 1.d1*Lsun * 10.**EL(i)
     tmp2(i) = 7.125d-16 * 10.**(4.*TE(i))
  end do

  do JL=1,4 
     do JT=1,4 
        call ATMOS(EL(JL),TE(JT),U(1,JL,JT),U(2,JL,JT),U(3,JL,JT),PMA1,0,IRA)
        if(IRA == 1) goto 38 
        ! giada
        vartmp = sqrt(tmp1(JL)/tmp2(JT))
        U(3,JL,JT) = log10(U(3,JL,JT)*vartmp) 
        PM11 = PMA1/Msun
        if(ETAFPR <= 1.D-10) cycle
        call ATMOS(EL(JL),TE(JT),V(1,JL,JT),V(2,JL,JT),V(3,JL,JT),PMA2,0,IRA)
        if(IRA == 1) goto 38 
        V(3,JL,JT) = log10(V(3,JL,JT)*vartmp) 
        PM22 = PMA2/Msun 
        write(2,333) EL(JL),TE(JT),(V(IT,JL,JT),IT=1,3),PM22 
        write(*,333) EL(JL),TE(JT),(V(IT,JL,JT),IT=1,3),PM22 
     end do
  end do
  INI = 2 
  ICONT = 1 

  do J=1,3 
     do K=1,4 
        do L=1,4 
           PTR(J,K,L) = U(J,K,L) 
        end do
     end do
  end do

  !****************CALCOLO GRIGLIATINO************************************
  DDRATT = 2.0*DRATT 
  DDRAPP = 2.0*DRAPP 
  TT1(1) = TEFF+DRATT 
  TT1(2) = TEFF-DRATT 
  EL1(1) = ELLOG+DRAPP 
  EL1(2) = ELLOG-DRAPP 
  do ITE=1,2 
     do ILU=1,2 
        do I=1,3 
           do J=1,4 
              if(ITE == 1 .and. ILU == 1) then
                 do K=1,4 
                    A(K) = PTR(I,J,K) 
                    B(K) = TE(K) 
                 end do
                 call CUB(A,B,ALFA1,BETA1,GAMMA1,DELTA1) 
                 alfac(i,j) = alfa1
                 betac(i,j) = beta1
                 gammac(i,j) = gamma1
                 deltac(i,j) = delta1
              else
                 alfa1 = alfac(i,j)
                 beta1 = betac(i,j)
                 gamma1 = gammac(i,j)
                 delta1 = deltac(i,j)
              endif
              SERV(J) = ALFA1*TT1(ITE)*TT1(ITE)*TT1(ITE)+BETA1*TT1(ITE)* & 
                   TT1(ITE)+GAMMA1*TT1(ITE)+DELTA1    
           end do
           do J=1,4 
              A(J) = SERV(J) 
              B(J) = EL(J) 
           end do
           call CUB(A,B,ALFA1,BETA1,GAMMA1,DELTA1) 
           AT(I,ITE,ILU) = ALFA1*EL1(ILU)*EL1(ILU)*EL1(ILU)+BETA1*EL1(ILU)* &
                EL1(ILU)+GAMMA1*EL1(ILU)+DELTA1     
        end do
     end do
  end do
  RATT = (TEFF-TT1(2))/DDRATT 
  RATL = (ELLOG-EL1(2))/DDRAPP 
  do K=1,4 
     C1(K) = AT(K,2,1)+RATT*(AT(K,1,1)-AT(K,2,1)) 
     C2(K) = AT(K,2,2)+RATT*(AT(K,1,2)-AT(K,2,2)) 
  end do
  PB = 10.**(C2(1)+RATL*(C1(1)-C2(1))-17.) 
  TB = 10.**(C2(2)+RATL*(C1(2)-C2(2))-6.) 
  RB = 10.**(C2(3)+RATL*(C1(3)-C2(3))) 
  DPL = (C1(1)-C2(1))/DDRAPP 
  DTL = (C1(2)-C2(2))/DDRAPP 
  DRL = (C1(3)-C2(3))/DDRAPP 
  do K=1,3 
     C1(K) = AT(K,1,2)+RATL*(AT(K,1,1)-AT(K,1,2)) 
     C2(K) = AT(K,2,2)+RATL*(AT(K,2,1)-AT(K,2,2))
  end do
  DPT = (C1(1)-C2(1))/DDRATT 
  DTT = (C1(2)-C2(2))/DDRATT 
  DRT = (C1(3)-C2(3))/DDRATT 
  if(INCAS == 1) then
     EMTUT = EMTOT 
     IRA = 3 
     call ATMOS(ELLOG,TEFF,PD,TD,RD,EMTUT,0,IRA) 
     IRA = 4 
     ALUCE = log10(GG(2,MAXMV)/38.2) 
     call ATMOS(ALUCE,TEFFV,PD,TD,RD,EMTOV,0,IRA) 
     MAXMEin = MAXME 
     MAXMV = MAXMIV 
     return 
  endif
  return 

  !*********RIDUCE ZONA SUBATMOSFERICA************************************
38 continue 

  EMA = PMA1 
  if(ETAFPR > 1.d-10) EMA = PMA2 
  IRA = 2 
  FRAZ = 1.-(1.-FRAZ)*RATTO 
  EMMA = EMTOT*(1.-FRAZ) 
  NANA = NANA+1 
  if(NANA >= 30) then
     write(66,*)'120 - quatm'
     write(66,*)'NANA>= 30'
     stop 'NANA>= 30' 
  endif
  write(2,100) FRAZ 
  write(*,100) FRAZ 

  INCAS = 1 
  MAYA = 5 
  goto 39 

100 format(/,3X,'ATTENZIONE-REZONING ATMOSFERA-FRAZ=',1P,D12.5,/) 
333 format(6F8.4,1P,3E10.3) 
end subroutine QUATM
