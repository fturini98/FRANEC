subroutine EPSI (RO,T,JF,NABLA,HT,EPS,EPSA,ECAR,DELTA) 
  use interfaccia
  use mistura
  use fisica
  use chequ
  use nummod
  use equil
  use neutri
  use chim
  use varie
  use chimic
  use overshoot
  use mesh
  use sceltachim
  use tempi_scala
  use accrescimento
  use zone_conv

  implicit none

  real :: RO, T, HT, EPS, EPSA, ECAR, DELTA
  integer :: JF, NABLA

  real,dimension(35) :: E
  integer,parameter :: IDEBUG = 0

  real,parameter,dimension(MELE) ::  Z = (/2.,1.,2., 6., 7., 8., 8.,10., &
       10.,12.,12.,12.,14.,0.,11., 6., 7., 8., 9.,10.,26.,3.,3.,4., 5.,1./)
  real,parameter,dimension(MELE) :: AT = (/3.,1.,4.,12.,14.,16.,18.,20., &
       22.,24.,25.,26.,28.,1.,23.,13.,15.,17.,19.,21.,56.,6.,7.,9.,11.,2./)

  ! controllare Q13 e Q32
  real,dimension(35) :: Q = (/ &
       1.943,   7.162,  7.551,  2.216,  7.297, &
       4.415,   4.966,  4.014,  0.600,  4.730, &
       1.192,   0.586,  7.993,  9.667,  8.114, &
       1.675,   2.431,  9.316,  6.739,  8.794, &
       -0.478, 10.615,  2.377, 11.693,  2.242, &
       4.621,   9.986,  7.274,  6.671, 12.859, &
       1.0,    11.004,  14.031, 3.616,  5.493 /)

  real,parameter :: AN = 6.022d23
  integer,parameter :: ISHORT = 1
  
  integer :: i, ir, k
  real :: epsp, epsc, xhe3
  real :: dmass, rpp, zeta, shep, schep, xhi
  real :: t6, t612, t9, t92, t93, t912, t913, t923, t932
  real :: roan

  integer :: ncache

  real,dimension(nsezi) :: S

  real,dimension(MELE) :: YV, YF 

  type(pow_of_T) :: powT


  EPSP = 0. 
  EPSC = 0. 
  EPSA = 0. 
  ECAR = 0. 
  EPS = 0. 
  DELTA = 1. 
  if( JF == 1 ) then 
     PPNEU = 0.
  endif
  
  if(T < 1.d5) then 
     return 
  end if
  if(T < 1.d6 .and. NMD >= 40) then 
     return 
  end if
  if(XX(1) <= 0. .and. XX(3) > 0. .and. T <= 4.d7) then 
     return 
  end if
  if(XX(1) <= 0. .and. XX(3) <=  0. .and. T <= 4.d8) then 
     return 
  end if
  if(XX(4) <= 0. .and. T > 3.d9) then 
     write(*,150) JF
     write(2,150) JF
     write(66,*) '50 - epsi'
     write(66,150) JF
     stop
  end if
  ! ******** IL RESTO DEL PROGRAMMA INVERTE IDROGENO E ELIO3 ************ 
  XHE3 = XX(2) 
  XX(2) = XX(1) 
  XX(1) = XHE3 
  !************************ INIZIALIZATION RAP-NEW ***********************

  ! controllo di non avere abbondanze di elementi leggeri negative
  do I=22,MELE 
     if(xx(i) < 0.) then
        xx(i) = 0.0
     endif
  end do

  ! ** CONTROLLO FASE EVOLUTIVA:                                        **
  ! ** IR=1 --> COMBUST. IDROGENO     IR=2 --> COMBUST. IDROGENO + ELIO **
  ! ** IR=3 --> COMBUST. ELIO         IR=4 --> COMBUST. CARBONIO        **
  if(XX(2) > 0. .and. T <= 5.d8) then 
     if(T < 1.d8) then 
        IR = 1 
     else 
        IR = 2 
     endif
  endif
  if(XX(2) <= 0. .and. T <= 5.d8) then 
     IR = 3 
  endif
  if(T > 5.d8) then 
     IR = 4 
  endif

  E(1:31) = 0. 

  do I=1,MELE 
     YV(I) = XX(I)/AT(I) 
     YF(I) = YV(I) 
  end do

  ! calcolo potenze di temperatura utili in cross
  T6 = T*1.d-6   
  T612 = sqrt(T6)
  T9 = T*1.d-9
  T92 = T9*T9
  T93 = T92*T9
  T912 = sqrt(T9)
  T913 = T9**(1./3.)
  T923 = T913*T913
  T932 = T9*T912

  powt%T9 = T9
  powt%T92 = T92
  powt%T93 = T93
  powt%T912 = T912
  powt%T913 = T913
  powt%T923 = T923
  powt%T932 = T932
  powt%T6 = T6
  powt%T612 = T612

  call CROSS(RO,T, YV, IR, ISHORT, S, powt, nabla)
  
  ncache = 1
  ROAN = RO*AN
  ! ************ CALCOLO COEFFICIENTE ENERGETICA NUCLEARE  ************** 
  if(IR <=  2) then 
     E(1) = Q(29)*S(29)*ROAN*YV(2)**2/2. 
     E(2) = Q(30)*S(30)*ROAN*YV(1)**2/2. 
     Q(31) = 1.586+(S(33)*11.407+17.397)/(1.+S(33))

     E(3) = Q(31)*S(31)*ROAN*YV(1)*YV(3)

     ! energia da deuterio S(34)
     E(35) = Q(35)*S(34)*ROAN*YV(2)*YV(26)

     ! ********* se c'e' deuterio devo de-commentare questa *****
     EPSP = 1.6022d-6*(E(1)+E(2)+E(3)+E(35))

     !******************** CALCOLO NEUTRINI PROTONE PROTONE ******
     DMASS = (G(5,JF+1)-G(5,JF))*1.d33
     PPNEU(1) = PPNEU(1)+S(29)*ROAN*YV(2)**2/2.*DMASS
     PPNEU(2) = PPNEU(2)+S(31)*ROAN*YV(1)*YV(3)/(S(33)+1.0)*DMASS
     PPNEU(3) = PPNEU(3)+S(31)*ROAN*YV(1)*YV(3)*S(33)/(S(33)+1.0)*DMASS
     E(4) = Q(32)*S(1)*ROAN*YV(2)*YV(4)
     E(5) = Q(33)*S(5)*ROAN*YV(2)*YV(5)
     E(6) = Q(34)*S(9)*ROAN*YV(2)*YV(6)
     EPSC = 1.6022d-6*(E(4)+E(5)+E(6))
     !****************** CALCOLO NEUTRINI pep ******************************
     RPP = S(29)*ROAN*YV(2)**2/2.
     ! giada
     PPNEU(7) = PPNEU(7)+1.102d-4*RO*(1.+XX(2))/2./T612*(1.+0.02 &
          *T6)*RPP*DMASS
     !******** CALCOLO SCHERMAGGIO ELETTRONICO PER Hep *******************
     ZETA = sqrt(sum(YV(1:MELE)*Z(1:MELE)*(Z(1:MELE) + 1.)))
     SCHep = exp(0.188*2*ZETA*sqrt(RO)/(T6*T612))
     !*************** CALCOLO NEUTRINI Hep ********************************
     SHep = SCHep*1.41d-13/T923*exp(-6.141/T913)
     PPNEU(8) = PPNEU(8)+SHep*ROAN*YV(2)*YV(1)*DMASS
     !************************* CALCOLO NEUTRINI CNO ***********************
     PPNEU(4) = PPNEU(4)+S(1)*ROAN*YV(2)*YV(4)*DMASS
     PPNEU(5) = PPNEU(5)+S(5)*ROAN*YV(2)*YV(5)*DMASS
     PPNEU(6) = PPNEU(6)+S(9)*ROAN*YV(2)*YV(6)*DMASS
  endif
  if(IR > 1) then
     E(7) = Q(28)*S(28)*RO**2*AN*YV(3)**3/6.
     E(8) = Q(2)*S(2)*ROAN*YV(3)*YV(4)
     E(9) = Q(6)*S(6)*ROAN*YV(3)*YV(5)
     E(10) = Q(10)*S(10)*ROAN*YV(3)*YV(6)
     E(11) = Q(14)*S(14)*ROAN*YV(3)*YV(7)
     E(12) = Q(18)*S(18)*ROAN*YV(3)*YV(8)
     E(13) = Q(21)*S(21)*ROAN*YV(3)*YV(9)
     E(14) = Q(22)*S(22)*ROAN*YV(3)*YV(9)
     E(15) = Q(27)*S(27)*ROAN*YV(3)*YV(10)
     do K=7,15
        EPSA = EPSA+1.6022d-6*E(K)
     end do
  endif
  if(IR == 4) then
     E(16) = Q(3)*S(3)*ROAN*YV(2)*YV(16)
     E(17) = Q(4)*S(4)*ROAN*YV(3)*YV(16)
     E(18) = Q(7)*S(7)*ROAN*YV(2)*YV(17)
     E(19) = Q(8)*S(8)*ROAN*YV(3)*YV(17)
     E(20) = Q(11)*S(11)*ROAN*YV(2)*YV(18)
     E(21) = Q(12)*S(12)*ROAN*YV(3)*YV(18)
     E(22) = Q(13)*S(13)*ROAN*YV(2)*YV(7)
     E(23) = Q(15)*S(15)*ROAN*YV(2)*YV(19)
     E(24) = Q(16)*S(16)*ROAN*YV(3)*YV(19)
     E(25) = Q(17)*S(17)*ROAN*YV(2)*YV(8)
     E(26) = Q(19)*S(19)*ROAN*YV(2)*YV(20)
     E(27) = Q(20)*S(20)*ROAN*YV(2)*YV(9)
     E(28) = Q(23)*S(23)*ROAN*YV(2)*YV(15)
     E(29) = Q(24)*S(24)*ROAN*YV(2)*YV(15)
     E(30) = Q(25)*S(25)*ROAN*YV(4)**2/2.
     E(31) = Q(26)*S(25)*ROAN*YV(4)**2/2.
     E(32) = Q(1)*S(1)*ROAN*YV(2)*YV(4)
     E(33) = Q(5)*S(5)*ROAN*YV(2)*YV(5)
     E(34) = Q(9)*S(9)*ROAN*YV(2)*YV(6)
     do K=16,34
        ECAR = ECAR+1.6022d-6*E(K)
     end do
  endif
  if(EPSP > 0. .and. EPSC > 0.) then
     DELTA = EPSP/(EPSP+EPSC)
  else
     DELTA = 1.
  endif
  EPS = ECAR+EPSA+EPSC+EPSP
  if( NABLA == 0 ) then
     XHI = XX(2)
     XX(2) = XX(1)
     XX(1) = XHI
     return
  endif
  
150 format(1X,'SUPIN RICADDE E PIU NON PARVE FORA ALLO MESCIO',I4)

end subroutine EPSI
