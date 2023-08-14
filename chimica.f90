! queste routine risolvono il sistema di equazioni differenziali del
! bruciamento degli elementi (Y).
! il modulo netw_chimica definisce le derivate delle abbondanze (ydot) e il 
! jacobiano del sistema (PD).

module mod_f
  use parametri
  real :: S(nsezi+2)
end module mod_f

module netw_chimica
contains

subroutine FEXc(NEQ,T,Y,YDOT)
    use mod_f

    implicit none
    
    integer :: NEQ
    real :: T, Y(NEQ), YDOT(NEQ)
    intent(IN)  :: NEQ, T, Y
    intent(OUT) :: YDOT

    ydot = 0.0

!! He3 :
    ydot(1) = (Y(2)**2*S(29)/2.-Y(1)**2*S(30)-Y(3)*Y(1)*S(31) + &
               Y(2)*y(14)*S(35) + Y(2)*y(10)*S(36)) 
!! H :
    ydot(2) = (Y(1)**2*S(30)-3.*Y(2)**2*S(29)/2.-Y(3)*Y(1)*S(31) - &
               2.*Y(2)*(Y(4)*S(1)+Y(5)*S(5)+Y(6)*S(9)) - Y(2)*y(14)*S(35) - &
               Y(2)*y(10)*S(36) - Y(2)*y(11)*S(37) - Y(2)*y(12)*S(38) - &
               Y(2)*y(13)*S(39))
!! He4 :
    ydot(3) = (-Y(3)*(Y(4)*S(2)+Y(5)*S(6)+Y(6)*S(10)+Y(7)*S(14) + &
                Y(8)*S(18)+Y(9)*(S(21)+S(22))+Y(3)**2*S(28)/2.) + &
                Y(1)**2*S(30)/2.+Y(3)*Y(1)*S(31) + Y(5)*Y(2)*S(nsezi+1)+ &
                Y(6)*Y(2)*S(9) + Y(2)*y(10)*S(36) + 2.0*Y(2)*y(11)*S(37) + &
                2.0*Y(2)*y(12)*S(38) + 3.0*Y(2)*y(13)*S(39))
!! C12 :  
    ydot(4) = (-Y(4)*Y(3)*S(2)+Y(3)**3*S(28)/6.+ Y(2)*(-Y(4)*S(1)+&
         Y(5)*S(nsezi+1))) 
!! N14 :    
    ydot(5) = (-Y(5)*Y(3)*S(6) + Y(2)*(-Y(5)*S(5)+Y(6)*S(9)+Y(4)*S(1)))
!! O16 :
    ydot(6) = (-Y(6)*Y(3)*S(10)+Y(3)*Y(4)*S(2) + Y(2)*(-Y(6)*S(9)+&
         Y(5)*S(nsezi+2))) 
!! O18 :
    ydot(7) = (-Y(7)*Y(3)*S(14)+Y(3)*Y(5)*S(6)) 
!! Ne20 :
    ydot(8) = (-Y(8)*Y(3)*S(18)+Y(3)*Y(6)*S(10)) 
!! Ne22 :
    ydot(9) = (-Y(9)*Y(3)*(S(21)+S(22))+Y(3)*Y(7)*S(14))
!! Li6 :
    ydot(10) = -Y(2)*y(10)*S(36)
!! Li7 :
    ydot(11) = -Y(2)*y(11)*S(37)
!! Be9 :
    ydot(12) = -Y(2)*y(12)*S(38)
!! B11 :
    ydot(13) = -Y(2)*y(13)*S(39)
!! D :
    ydot(14) = -Y(2)*y(14)*S(35)!! + Y(2)*y(12)*S(38)

    return
  end subroutine FEXc

  subroutine JEXc(NEQ,T,Y,ML,MU,PD,NRPD)
    use mod_f

    implicit none

    integer :: NEQ, ML, MU, NRPD
    real :: PD(NRPD,NEQ), T, Y(NEQ)

    PD = 0.0

!! He3 :
    PD(1,1) = -(2.*Y(1)*S(30)+Y(3)*S(31)) 
    PD(1,2) = Y(2)*S(29) + y(14)*S(35) + y(10)*S(36)
    PD(1,3) = -Y(1)*S(31) 
    PD(1,10) = Y(2)*S(36)
    PD(1,14) = Y(2)*S(35)

!! H :
    PD(2,1) = (2.*Y(1)*S(30)-Y(3)*S(31)) 
    PD(2,2) = -(3.*Y(2)*S(29)+2.*(Y(4)*S(1)+Y(5)*S(5)+Y(6)*S(9)) + &
                y(14)*S(35) + y(10)*S(36) + y(11)*S(37) + y(12)*S(38) + &
                y(13)*S(39))
    PD(2,3) = -Y(1)*S(31) 
    PD(2,4) = -2.*Y(2)*S(1) 
    PD(2,5) = -2.*Y(2)*S(5) 
    PD(2,6) = -2.*Y(2)*S(9)
    PD(2,10) = -Y(2)*S(36)
    PD(2,11) = -Y(2)*S(37)
    PD(2,12) = -Y(2)*S(38)
    PD(2,13) = -Y(2)*S(39)
    PD(2,14) = -Y(2)*S(35) 

!! He4 :
    PD(3,1) = (Y(1)*S(30)+Y(3)*S(31)) 
    PD(3,2) = (Y(5)*S(nsezi+1)+Y(6)*S(9) + y(10)*S(36) + 2.0*y(11)*S(37) + &
               2.0*y(12)*S(38) + 3.0*y(13)*S(39))
    PD(3,3) = -(Y(4)*S(2)+Y(5)*S(6)+Y(6)*S(10)+Y(7)*S(14)+ & 
          Y(8)*S(18)+Y(9)*(S(21)+S(22))+Y(3)**2*S(28)*3./2.) 
    !-y(22)*S(27)
    PD(3,4) = -Y(3)*S(2) 
    PD(3,5) = (Y(2)*S(nsezi+1)-Y(3)*S(6)) 
    PD(3,6) = -(Y(3)*S(10)-Y(2)*S(9)) 
    PD(3,7) = -Y(3)*S(14) 
    PD(3,10) = Y(2)*S(36)
    PD(3,11) = 2.0*Y(2)*S(37)
    PD(3,12) = 2.0*Y(2)*S(38)
    PD(3,13) = 3.0*Y(2)*S(39)

!! C12 :
    PD(4,2) = (Y(5)*S(nsezi+1)-Y(4)*S(1)) 
    PD(4,3) = -(Y(4)*S(2)-Y(3)**2*S(28)/2.) 
    PD(4,4) = -(Y(2)*S(1)+Y(3)*S(2)) 
    PD(4,5) = Y(2)*S(nsezi+1) 

!! N14 :
    PD(5,2) = (Y(4)*S(1)-Y(5)*S(5)+Y(6)*S(6)) 
    PD(5,3) = -Y(5)*S(6) 
    PD(5,4) = Y(2)*S(1) 
    PD(5,5) = -(Y(2)*S(5)+Y(3)*S(6)) 
    PD(5,6) = Y(2)*S(9) 
    
!! O16 :
    PD(6,2) = (Y(5)*S(nsezi+2)-Y(6)*S(9)) 
    PD(6,3) = -(Y(6)*S(10)-Y(4)*S(2)) 
    PD(6,4) = S(2)*Y(3)
    PD(6,5) = Y(2)*S(nsezi+2) 
    PD(6,6) = -(Y(2)*S(9)+Y(3)*S(10)) 

!! O18 :
    PD(7,3) = -(Y(7)*S(14)-Y(5)*S(6)) 
    PD(7,5) = Y(3)*S(6) 
    PD(7,7) = -Y(3)*S(14) 

!! Ne20 :
    PD(8,3) = (Y(6)*S(10)-Y(8)*S(18))
    PD(8,6) = Y(3)*S(10) 
    PD(8,8) = -Y(3)*S(18)

!! Ne22 :
    PD(9,3) = (Y(7)*S(14)-Y(9)*(S(21)+S(22)))
    PD(9,7) = Y(3)*S(14) 
    PD(9,9) = (-Y(3)*(S(21)+S(22)))

!! Li6 :
    PD(10,2) = -y(10)*S(36)
    PD(10,10) = -Y(2)*S(36)

!! Li7 :
    PD(11,2) = -y(11)*S(37)
    PD(11,11) = -Y(2)*S(37)

!! Be9 :
    PD(12,2) = -y(12)*S(38)
    PD(12,12) = -Y(2)*S(38)

!! B11 :
    PD(13,2) = -y(13)*S(39)
    PD(13,13) = -Y(2)*S(39)

!! D :
    PD(14,2) = -y(14)*S(35)
    PD(14,14) = -Y(2)*S(35)

    return
  end subroutine JEXc

  ! questo routine attuamlmente non e' in uso.
  subroutine struct_sparse(ia,nia, ja, nja)
    implicit none
    integer :: ia(nia), ja(nja), nia, nja

    ia(1) = 1
    ia(2) = ia(1) + 3
    ia(3) = ia(2) + 11
    ia(4) = ia(3) + 9
    ia(5) = ia(4) + 5
    ia(6) = ia(5) + 6
    ia(7) = ia(6) + 5
    ia(8) = ia(7) + 3
    ia(9) = ia(8) + 1
    ia(10) = ia(9) + 1
    ia(11) = ia(10) + 4
    ia(12) = ia(11) + 3
    ia(13) = ia(12) + 3
    ia(14) = ia(13) + 3
    ia(15) = ia(14) + 3

    ja(1) = 1
    ja(2) = 2
    ja(3) = 3

    ja(4) = 1
    ja(5) = 2
    ja(6) = 3
    ja(7) = 4
    ja(8) = 5
    ja(9) = 6
    ja(10) = 10
    ja(11) = 11
    ja(12) = 12   
    ja(13) = 13
    ja(14) = 14

    ja(15) = 1
    ja(16) = 2
    ja(17) = 3   
    ja(18) = 4
    ja(19) = 5
    ja(20) = 6
    ja(21) = 7
    ja(22) = 8
    ja(23) = 9

    ja(24) = 2
    ja(25) = 3
    ja(26) = 4
    ja(27) = 5
    ja(28) = 6
    
    ja(29) = 2
    ja(30) = 3
    ja(31) = 4
    ja(32) = 5
    ja(33) = 6
    ja(34) = 7

    ja(35) = 2
    ja(36) = 3
    ja(37) = 5
    ja(38) = 6
    ja(39) = 8

    ja(40) = 3
    ja(41) = 7
    ja(42) = 9

    ja(43) = 8

    ja(44) = 9

    ja(45) = 1
    ja(46) = 2
    ja(47) = 3
    ja(48) = 10
    
    ja(49) = 2
    ja(50) = 3
    ja(51) = 11

    ja(52) = 2
    ja(53) = 3
    ja(54) = 12
    
    ja(55) = 2
    ja(56) = 3
    ja(57) = 13

    ja(58) = 1
    ja(59) = 2
    ja(60) = 14

  end subroutine struct_sparse

end module netw_chimica



subroutine chimica(PAS, N, Y, Si)
  use parametri
  use mod_f
  use netw_chimica
  use DVODE_F90_M
  
  implicit none
  
  integer :: n, j
  real :: roi, pas
  real,dimension(n) :: y 
  real,dimension(nsezi+2) :: Si

  real :: RTOL, T, TOUT, RSTATS(22)
  integer :: ISTATS(31), IOUT, IERROR, I, ISTATE, ITASK, myn
  real :: ATOL(n), locY(14)
  type (VODE_OPTS),save :: OPTIONS

  integer, save :: firsttime = 1
  
  integer :: ia(15), ja(60)
  
  S = Si

  T = 0.0
  TOUT = PAS
  
  RTOL = 1.d-8
  ATOL = 1.d-14
      
  ITASK = 1
  ISTATE = 1
  
  ! compatto la chimica
  locY(1:9) = y(1:9)
  locY(10:14) = y(22:26)
  ! numero di equazioni compattate
  myn = 14

  if(firsttime == 1) then
     OPTIONS = SET_NORMAL_OPTS(DENSE_J=.true.,ABSERR_VECTOR=ATOL,      &
          RELERR=RTOL,USER_SUPPLIED_JACOBIAN=.true.)
  
     firsttime = 0
  endif

  call DVODE_F90(FEXc,myN,locY,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEXc)
  
  if (ISTATE<0) then
     write (*,*) "errore  -  istate: ",ISTATE
     stop
  end if

  ! ripristino i valori nel vettore della chimica
  y(1:9) = locY(1:9)
  y(22:26) = locY(10:14)
  
end subroutine chimica
