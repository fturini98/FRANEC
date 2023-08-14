subroutine ELIO (RO,T,PAS,JF,INI,IFI,N,ISHORT, YV, YF, NEQUI, S) 
  use fisica

  implicit none

  real :: RO, T, PAS
  integer :: JF, INI, IFI, N, ISHORT
  real,dimension(MELE) :: YV, YF 
  integer,dimension(MELE) :: NEQUI
  real,dimension(nsezi) :: S

  real,parameter :: CR = 0.5
  integer :: ni, ih, kk, j, ib, ig, l, lk, ks, jl
  real :: beta, gamma, ef, dd

  real,dimension(MELE,MELE) :: DXX
  real,dimension(MELE) :: Y, DX, DB, B
  real,dimension(MELE*MELE) :: A 

  logical,parameter :: linConv = .false.
      
  ! azzero la matrice jacobiana del sistema da risolvere
  do KK=ini,ifi 
     do IH=ini,ifi
        DXX(IH,KK) = 0.
     end do
  end do

  do NI=0,20
     do IH=1,mele 
        Y(IH) = .5*(YF(IH)+YV(IH)) 
     end do
     ! ******* CALCOLO DERIVATE DY(K)/DT IN T+PAS/2 SE RADIATIVO   **********
     ! ******* - - - - - - - - - - - -   IN    T    SE CONVETTIVO  **********
     DX(3) = RO*(-Y(3)*(Y(4)*S(2)+Y(5)*S(6)+Y(6)*S(10)+Y(7)*S(14)+       &
          Y(8)*S(18)+Y(9)*(S(21)+S(22))+Y(10)*S(27)+Y(3)**2*S(28)*RO/2.)    &
          +Y(1)**2*S(30)/2.+Y(3)*Y(1)*S(31))
     DX(4) = RO*(-Y(4)*Y(3)*S(2)+Y(3)**3*S(28)*RO/6.) 
     DX(5) = -RO*Y(5)*Y(3)*S(6) 
     DX(6) = RO*(-Y(6)*Y(3)*S(10)+Y(3)*Y(4)*S(2)) 
     DX(7) = RO*(-Y(7)*Y(3)*S(14)+Y(3)*Y(5)*S(6)) 
     DX(8) = RO*(-Y(8)*Y(3)*S(18)+Y(3)*Y(6)*S(10)) 
     DX(9) = RO*(-Y(9)*Y(3)*(S(21)+S(22))+Y(3)*Y(7)*S(14)) 
     if(ISHORT /= 1) then 
        DX(10) = RO*(Y(8)*Y(3)*S(18)-Y(3)*Y(10)*S(27)) 
        DX(11) = RO*Y(9)*Y(3)*S(21) 
        DX(12) = RO*Y(9)*Y(3)*S(22) 
        DX(13) = RO*Y(3)*Y(10)*S(27) 
        DX(14) = RO*Y(9)*Y(3)*S(21) 
     endif
     !************* SE CONVETTIVO INTEGRA LINEARMENTE   ******************   

     if( G(6,JF) >= 0. .and. linConv ) then 
        do KK=INI,IFI 
           YF(KK) = Y(KK)+DX(KK)*PAS 
        end do
        return 
     endif
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
        if( BETA < 1.d-4 .or. GAMMA < 1.d-4) return 
     endif

     ! **********************  DERIVATE -DB(J)/DYF  *************************
     DXX(3,3) = 1./PAS+RO*(Y(4)*S(2)+Y(5)*S(6)+Y(6)*S(10)+Y(7)*S(14)+ & 
          Y(8)*S(18)+Y(9)*(S(21)+S(22))+Y(10)*S(27)+Y(3)**2*S(28)*3.*RO/2.)*CR
     DXX(3,4) = RO*Y(3)*S(2)*CR 
     DXX(3,5) = RO*Y(3)*S(6)*CR 
     DXX(3,6) = RO*Y(3)*S(10)*CR 
     DXX(3,7) = RO*Y(3)*S(14)*CR 
     DXX(4,3) = RO*(Y(4)*S(2)-Y(3)**2*S(28)*RO/2.)*CR 
     DXX(4,4) = 1./PAS+RO*Y(3)*S(2)*CR 
     DXX(5,3) = RO*Y(5)*S(6)*CR 
     DXX(5,5) = 1./PAS+RO*Y(3)*S(6)*CR 
     DXX(6,3) = RO*(Y(6)*S(10)-Y(4)*S(2))*CR 
     DXX(6,4) = -RO*S(2)*Y(3)*CR 
     DXX(6,6) = 1./PAS+RO*Y(3)*S(10)*CR 
     DXX(7,3) = RO*(Y(7)*S(14)-Y(5)*S(6))*CR 
     DXX(7,5) = -RO*Y(3)*S(6)*CR 
     DXX(7,7) = 1./PAS+RO*Y(3)*S(14)*CR 
     DXX(8,3) = RO*(Y(8)*S(18)-Y(6)*S(10))*CR 
     DXX(8,6) = -RO*Y(3)*S(10)*CR 
     DXX(8,8) = 1./PAS+RO*Y(3)*S(18)*CR 
     DXX(9,3) = RO*(Y(9)*(S(21)+S(22))-Y(7)*S(14))*CR 
     DXX(9,7) = -RO*Y(3)*S(14)*CR 
     DXX(9,9) = 1./PAS+RO*Y(3)*(S(21)+S(22))*CR 
     if(ISHORT /= 1) then 
        DXX(3,8) = RO*Y(3)*S(18)*CR 
        DXX(3,9) = RO*Y(3)*(S(21)+S(22))*CR 
        DXX(3,10) = RO*Y(3)*S(27)*CR 
        DXX(10,3) = RO*(-Y(8)*S(18)+Y(10)*S(27))*CR 
        DXX(10,8) = -RO*Y(3)*S(18)*CR 
        DXX(10,10) = 1./PAS+RO*Y(3)*S(27)*CR 
        DXX(11,3) = -RO*Y(9)*S(21)*CR 
        DXX(11,9) = -RO*Y(3)*S(21)*CR 
        DXX(11,11) = 1./PAS 
        DXX(12,3) = -RO*Y(9)*S(22)*CR 
        DXX(12,9) = -RO*Y(3)*S(22)*CR 
        DXX(12,12) = 1./PAS 
        DXX(13,3) = -RO*Y(10)*S(27)*CR 
        DXX(13,10) = -RO*Y(3)*S(27)*CR 
        DXX(13,13) = 1./PAS 
        DXX(14,3) = -RO*Y(9)*S(21)*CR 
        DXX(14,9) = -RO*Y(3)*S(21)*CR 
        DXX(14,14) = 1./PAS 
     endif
     ! *************  CALCOLO DELLE DIFFERENZE FINITE  ********************* 
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
        write(66,*)'41 - elio'
        write(66,130) 
        stop 
     endif
     GAMMA = 0.0 
     KK = 0 
     do J=INI,IFI 
        KK = KK+1 
        YF(J) = YF(J)+B(KK) 
        !*********** CALCOLA VARIAZIONI PERCENTUALI ELEMENTI IMPORTANTI *****
        if( YV(J) > 1.d-9 ) then 
           DD = abs(B(KK)/YF(J)) 
           if( DD > GAMMA ) then 
              IG = J 
              GAMMA = DD 
           endif
        endif
     end do
  end do

  ! se sono qui vuol dire che sono uscito dal ciclo per troppe iterazioni
  write(*,250) 
  write(2,250) 
  write(2,113) T,RO,JF 
  write(66,*)'40 - elio'
  write(66,250)
  write(66,113) T,RO,JF 
  stop 

113 format(1X,'T =',1P,E10.3,2X,'RO =',1P,E10.3,2X,'MESH =',I4) 
130 format(1X,'ATTENZIONE QUALCOSA NON VA NELLA EPSI. IL DET = 0 NEL  &
         &RAHPSON/NEWTON - vedi elio')                                     
250 format(1X,'MORTO: LA EPSI ITERA TROPPO NEL RAPHSON/NEWTON ') 
444 format(1X,'ITERAZIONI EPSI: DF =',I3,E10.3,' DX =',I3,E10.3,      &
         ' N. ITER.=',I3,' N. MESH=',I4) 
end subroutine ELIO
