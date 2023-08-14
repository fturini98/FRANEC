subroutine IDRELI (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S) 
  use fisica

  implicit none

  real :: RO, T, PAS
  integer :: JF, INI, IFI, N
  real,dimension(MELE) :: YV, YF 
  integer,dimension(MELE) :: NEQUI
  real,dimension(nsezi) :: S

  real,parameter :: CR = 0.5
  integer :: ni, ih, kk, j, ib, ig, l, lk, jl, ks    
  real :: the3, beta, gamma, ef, dd, rapp1, rapp2 

  real,dimension(MELE,MELE) :: DXX
  real,dimension(MELE) :: Y, DX, DB, B
  real,dimension(MELE*MELE) :: A 

  ! Efficienza relativa del ramo CN rispetto a CN+NO
  RAPP1 = S(7)/(S(7)+S(32))
  RAPP2 = 1.-RAPP1

  !********************* TEMPO SCALA EQUILIBRIO HE3 **********************
  if( G(6,JF) < 0.0 ) then 
     THE3 = 1./((YV(1)*S(30)+YV(3)*S(31))*RO) 
     if(THE3 < PAS) NEQUI(1) = 0 
  endif

  ! azzero la matrice jacobiana del sistema da risolvere
  do KK=ini,ifi 
     do IH=ini,ifi
        DXX(IH,KK) = 0.
     end do
  end do

  do NI=0,20
     do IH=INI,IFI 
        Y(IH) = .5*(YF(IH)+YV(IH)) 
     end do
     ! ******* CALCOLO DERIVATE DY(K)/DT IN T+PAS/2 SE RADIATIVO   **********
     ! ******* - - - - - - - - - - - -   IN    T    SE CONVETTIVO  **********
     DX(1) = RO*(Y(2)**2*S(29)/2.-Y(1)**2*S(30)-Y(3)*Y(1)*S(31)) 
     DX(2) = RO*(Y(1)**2*S(30)-3.*Y(2)**2*S(29)/2.-Y(3)*Y(1)*S(31)       &
          -2.*Y(2)*(Y(4)*S(1)+Y(5)*S(5)+Y(6)*S(9)))
     DX(3) = RO*(-Y(3)*(Y(4)*S(2)+Y(5)*S(6)+Y(6)*S(10)+Y(7)*S(14)+        &
          Y(3)**2*S(28)*RO/2.)+Y(1)**2*S(30)/2.+Y(3)*Y(1)*S(31)+RAPP1* &
          Y(5)*Y(2)*S(5)+Y(6)*Y(2)*S(9))
     DX(4) = RO*(-Y(4)*Y(3)*S(2)+Y(3)**3*S(28)*RO/6.+Y(2)*(-Y(4)*S(1)+ &
          RAPP1*Y(5)*S(5)))
     DX(5) = RO*(-Y(5)*Y(3)*S(6)+Y(2)*(-Y(5)*S(5)+Y(6)*S(9)+Y(4)*S(1))) 
     DX(6) = RO*(-Y(6)*Y(3)*S(10)+Y(3)*Y(4)*S(2)+Y(2)*(-Y(6)*S(9)+RAPP2  &
          *Y(5)*S(5)))      
     DX(7) = RO*(-Y(7)*Y(3)*S(14)+Y(3)*Y(5)*S(6)) 
     DX(8) = RO*Y(3)*Y(6)*S(10) 
     DX(9) = RO*Y(3)*Y(7)*S(14) 

     !************* SE CONVETTIVO INTEGRA LINEARMENTE  ******************   
!!$     if( G(6,JF) >= 0. ) then 
!!$        do KK=INI,IFI 
!!$           YF(KK) = Y(KK)+DX(KK)*PAS 
!!$        end do
!!$        return 
!!$     endif
     !********************************************************************** 
     !************************* RAHPSON - NEWTON *************************** 
     !********************************************************************** 
     !************* CALCOLO TERMINI NOTI E CONTROLLO CONVERGENZA *********** 
     BETA = 0.0 
     do J=INI,IFI 
        DB(J) = -(YF(J)-YV(J))/PAS*NEQUI(J)+DX(J) 
        if(.not. (YV(J) < 1.d-9 .or. NI == 0 .or. NEQUI(J) == 0)) then
           EF = abs(1.0-(YF(J)-YV(J))/(PAS*DX(J))) 
           if( EF > BETA ) then 
              IB = J 
              BETA = EF 
           endif
        endif
     end do
     if( NI /= 0 ) then  
        if( NI > 3 ) then 
           write(2,444) IB,BETA,IG,GAMMA,NI,JF 
        endif
        if( (BETA < 1.d-4 .or. GAMMA < 1.d-4) ) return 
     endif

     ! ************************  DERIVATE -DB(J)/DYF  *********************
     DXX(1,1) = 1./PAS+RO*(2.*Y(1)*S(30)+Y(3)*S(31))*CR 
     DXX(1,2) = -RO*Y(2)*S(29)*CR 
     DXX(1,3) = RO*Y(1)*S(31)*CR 
     DXX(2,1) = RO*(-2.*Y(1)*S(30)+Y(3)*S(31))*CR 
     DXX(2,2) = 1./PAS+RO*(3.*Y(2)*S(29)+2.*(Y(4)*S(1)+Y(5)*S(5)+Y(6)*   &
          S(9)))*CR
     DXX(2,3) = RO*Y(1)*S(31)*CR 
     DXX(2,4) = RO*2.*Y(2)*S(1)*CR 
     DXX(2,5) = RO*2.*Y(2)*S(5)*CR 
     DXX(2,6) = RO*2.*Y(2)*S(9)*CR 
     DXX(3,1) = RO*(-Y(1)*S(30)-Y(3)*S(31))*CR 
     DXX(3,2) = RO*(-RAPP1*Y(5)*S(5)-Y(6)*S(9))*CR 
     DXX(4,2) = RO*(-RAPP1*Y(5)*S(5)+Y(4)*S(1))*CR 
     DXX(4,5) = -RO*RAPP1*Y(2)*S(5)*CR 
     DXX(5,2) = RO*(-Y(4)*S(1)+Y(5)*S(5)-Y(6)*S(6))*CR 
     DXX(5,4) = -RO*Y(2)*S(1)*CR 
     DXX(5,6) = -RO*Y(2)*S(9)*CR 
     DXX(6,2) = RO*(-RAPP2*Y(5)*S(5)+Y(6)*S(9))*CR 
     DXX(6,5) = -RO*RAPP2*Y(2)*S(5)*CR 
     DXX(3,3) = 1./PAS+RO*(Y(4)*S(2)+Y(5)*S(6)+Y(6)*S(10)+Y(7)*S(14)+ &
          Y(3)**2*S(28)*3.*RO/2.)*CR
     DXX(3,4) = RO*Y(3)*S(2)*CR 
     DXX(3,5) = RO*(-RAPP1*Y(2)*S(5)+Y(3)*S(6))*CR 
     DXX(3,6) = RO*(Y(3)*S(10)-Y(2)*S(9))*CR 
     DXX(3,7) = RO*Y(3)*S(14)*CR 
     DXX(4,3) = RO*(Y(4)*S(2)-Y(3)**2*S(28)*RO/2.)*CR 
     DXX(4,4) = 1./PAS+RO*(Y(2)*S(1)+Y(3)*S(2))*CR 
     DXX(5,3) = RO*Y(5)*S(6)*CR 
     DXX(5,5) = 1./PAS+RO*(Y(2)*S(5)+Y(3)*S(6))*CR 
     DXX(6,3) = RO*(Y(6)*S(10)-Y(4)*S(2))*CR 
     DXX(6,4) = -RO*S(2)*Y(3)*CR 
     DXX(6,6) = 1./PAS+RO*(Y(2)*S(9)+Y(3)*S(10))*CR 
     DXX(7,3) = RO*(Y(7)*S(14)-Y(5)*S(6))*CR 
     DXX(7,5) = -RO*Y(3)*S(6)*CR 
     DXX(7,7) = 1./PAS+RO*Y(3)*S(14)*CR 
     DXX(8,3) = -RO*Y(6)*S(10)*CR 
     DXX(8,6) = -RO*Y(3)*S(10)*CR 
     DXX(8,8) = 1./PAS 
     DXX(9,3) = -RO*Y(7)*S(14)*CR 
     DXX(9,7) = -RO*Y(3)*S(14)*CR 
     DXX(9,9) = 1./PAS 
     ! *************  CALCOLO DELLE DIFFERENZE FINITE  *******************
     LK = 0 
     do J=INI,IFI 
        KK = 0 
        LK = LK+1 
        B(LK) = DB(J) 
        do L=INI,IFI 
           KK = KK+1 
           JL = (J-INI)*N+KK 
           A(JL) = DXX(L,J) 
        end do
     end do
     call SIMQ(A,B,N,KS) 
     if( KS == 1 ) then
        write(*,130) 
        write(2,130) 
        write(66,*)'71 - idreli'
        write(66,130) 
        stop
     endif
     GAMMA = 0.0 
     KK = 0 
     do J=INI,IFI 
        KK = KK+1 
        YF(J) = YF(J)+B(KK) 
        !*********** CALCOLA VARIAZIONI PERCENTUALI ELEMENTI IMPORTANTI **
        if( YV(J) > 1.d-9 ) then 
           DD = abs(B(KK)/YF(J)) 
           if( DD > GAMMA ) then 
              IG = J 
              GAMMA = DD 
           endif
        endif
     end do

     ! questo patch serve per impedire H, He negative in uscita
     if(yf(2) < 0.) yf(2) = 0.
     if(yf(3) < 0.) yf(3) = 0.

  end do

  ! se sono qui vuol dire che sono uscito dal ciclo per troppe iterazioni
  write(*,250) 
  write(2,250) 
  write(2,113) T,RO,JF 
  write(66,*)'70 - idreli'
  write(66,250) 
  write(66,113) T,RO,JF 
  stop 

113 format(1X,'T =',1P,E10.3,2X,'RO =',1P,E10.3,2X,'MESH =',I4) 
130 format(1X,'ATTENZIONE QUALCOSA NON VA NELLA EPSI. IL DET = 0 NEL  &
         &RAPHSON/NEWTON - vedi idreli')
250 format(1X,'MORTO: LA EPSI ITERA TROPPO NEL RAPHSON/NEWTON ') 
444 format(1X,'ITERAZIONI EPSI: DF =',I3,E10.3,' DX =',I3,E10.3,      &
         ' N. ITER.=',I3,' N. MESH=',I4)                     

end subroutine IDRELI
