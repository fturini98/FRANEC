subroutine CARBON (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S) 
  use interfaccia
  use fisica

  implicit none

  real :: RO, T, PAS
  integer :: JF, INI, IFI, N
  real,dimension(MELE) :: YV, YF 
  integer,dimension(MELE) :: NEQUI
  real,dimension(nsezi) :: S

  real,parameter :: CR = 0.5
  integer :: ni, ih, j, ib, ig, lk, kk, l, jl, ks
  real :: beta, gamma, ef, dd

  real,dimension(MELE,MELE) :: DXX
  real,dimension(MELE) :: Y, DX, DB, B
  real,dimension(MELE*MELE) :: A 

  NI = 0 
  
  ! azzero la matrice jacobiana del sistema da risolvere
  do KK=ini,ifi 
     do IH=ini,ifi
        DXX(IH,KK) = 0.
     end do
  end do
  
  do
     do IH=1,mele
        Y(IH) = .5*(YF(IH)+YV(IH)) 
     end do
     ! ******* CALCOLO DERIVATE DY(K)/DT IN T+PAS/2 SE RADIATIVO   **********
     ! ******* - - - - - - - - - - - -   IN    T    SE CONVETTIVO  **********
     DX(2) = RO*(-YF(2)*(YF(4)*S(1)+YF(16)*S(3)+YF(5)*S(5)+YF(17)*S(7)        &
          +YF(6)*S(9)+YF(18)*S(11)+YF(7)*S(13)+YF(19)*S(15)+YF(8)*S(17)+      &
          YF(20)*S(19)+YF(9)*S(20)+YF(15)*(S(23)+S(24)))+YF(4)**2*S(25)/2.+   &
          YF(3)*YF(19)*S(16))
     DX(3) = RO*(-YF(3)*(YF(4)*S(2)+YF(16)*S(4)+YF(5)*S(6)+YF(17)*S(8)+       &
          YF(6)*S(10)+YF(18)*S(12)+YF(7)*S(14)+YF(19)*S(16)+YF(8)*S(18)+YF(9) &
          *(S(21)+S(22))+YF(10)*S(27))+YF(2)*(YF(18)*S(11)+YF(7)*S(13)+       &
          YF(19)*S(15)+YF(15)*S(23)+YF(17)*S(7))+YF(4)**2*S(40)/2.)
     DX(4) = RO*(-Y(4)*(Y(2)*S(1)+Y(3)*S(2)+Y(4)*(S(25)+S(40)))+ &
          Y(2)*Y(17)*S(7)) 
     DX(5) = RO*(-Y(5)*(Y(2)*S(5)+Y(3)*S(6))+Y(2)*Y(16)*S(3)+Y(2)*Y(18)*S(11))
     DX(6) = RO*(-Y(6)*(Y(2)*S(9)+Y(3)*S(10))+Y(2)*Y(19)*S(15)+               &
          Y(3)*Y(4)*S(2)+Y(3)*Y(16)*S(4))
     DX(7) = RO*(-Y(7)*(Y(2)*S(13)+Y(3)*S(14))+Y(3)*Y(5)*S(6)) 
     DX(8) = RO*(-Y(8)*(Y(2)*S(17)+Y(3)*S(18))+Y(3)*Y(6)*S(10)+Y(3)*Y(18)     &
          *S(12)+Y(2)*Y(15)*S(23)+Y(4)**2*S(40)/2.) 
     DX(9) = RO*(-Y(9)*(Y(2)*S(20)+Y(3)*(S(21)+S(22)))+Y(3)*Y(7)*S(14)+       &
          Y(3)*Y(19)*S(16)+Y(2)*Y(20)*S(19))                               
     DX(10) = RO*(Y(8)*Y(3)*S(18)+Y(15)*Y(2)*S(24)-Y(3)*Y(10)*S(27)) 
     DX(11) = RO*Y(9)*Y(3)*S(21) 
     DX(12) = RO*Y(9)*Y(3)*S(22) 
     DX(13) = RO*Y(3)*Y(10)*S(27) 
     DX(14) = RO*(Y(16)*Y(3)*S(4)+Y(18)*Y(3)*S(12)+Y(9)*Y(3)*S(21)) 
     DX(15) = RO*(-Y(15)*Y(2)*(S(23)+S(24))+Y(2)*Y(9)*S(20)+Y(4)**2*S(25)/2.)
     DX(16) = RO*(-Y(16)*(Y(2)*S(3)+Y(3)*S(4))+Y(2)*Y(4)*S(1)) 
     DX(17) = RO*(-Y(17)*(Y(2)*S(7)+Y(3)*S(8))+Y(2)*Y(5)*S(5)+Y(2)*Y(7)*S(13))
     DX(18) = RO*(-Y(18)*(Y(2)*S(11)+Y(3)*S(12))+Y(2)*Y(6)*S(9)) 
     DX(19) = RO*(-Y(19)*(Y(2)*S(15)+Y(3)*S(16))+Y(3)*Y(17)*S(8)) 
     DX(20) = RO*(-Y(20)*Y(2)*S(19)+Y(2)*Y(8)*S(17)) 
     !************* SE CONVETTIVO INTEGRA LINEARMENTE  *********&*********   
     !     IF( G(6,JF) >=  0. ) THEN                                        
     !       DO KK=INI,IFI                                                 
     !       YF(KK)=Y(KK)+DX(KK)*PAS                                         
     !     END DO
     !       RETURN                                                          
     !     ENDIF                                                             
     !********************************************************************** 
     !************************* RAHPSON - NEWTON *************************** 
     !********************************************************************** 
     !************* CALCOLO TERMINI NOTI E CONTROLLO CONVERGENZA *********** 
     BETA = 0.0 
     do J=INI,IFI 
        DB(J) = -(YF(J)-YV(J))/PAS*NEQUI(J)+DX(J) 
        if(YV(J) < 1.d-9 .or. NI == 0 .or. NEQUI(J) == 0) cycle
        EF = abs(1.0-(YF(J)-YV(J))/(PAS*DX(J))) 
        if( EF > BETA ) then 
           IB = J 
           BETA = EF 
        endif
     end do
     if( NI /= 0 ) then
        if( NI > 5 ) then 
           write(2,444) IB,BETA,IG,GAMMA,NI,JF 
        endif
        if( BETA < 1.d-4 .or. GAMMA < 1.d-4) return 
        if( NI > 20 ) then
           write(*,250) 
           write(2,250) 
           write(2,113)T,RO,JF 
           write(66,*) '30 - carbon'
           write(66,250) 
           write(66,113)T,RO,JF 
           stop 
        endif
     endif
     NI = NI+1 
     ! ************************  DERIVATE DB(J)/DY  *************************
     DXX(2,2) = RO*(YF(4)*S(1)+YF(16)*S(3)+YF(5)*S(5)+YF(17)*S(7)+YF(6)*   &
          S(9)+YF(18)*S(11)+YF(7)*S(13)+YF(19)*S(15)+YF(8)*S(17)+          &
          YF(20)*S(19)+YF(9)*S(20)+YF(15)*(S(23)+S(24))) 
     DXX(2,3) = -RO*YF(19)*S(16) 
     DXX(2,4) = -RO*(-YF(2)*S(1)+YF(4)*S(25)) 
     DXX(2,5) = RO*YF(2)*S(5) 
     DXX(2,6) = RO*YF(2)*S(9) 
     DXX(2,7) = RO*YF(2)*S(13) 
     DXX(2,8) = RO*YF(2)*S(17) 
     DXX(2,9) = RO*YF(2)*S(20) 
     DXX(2,15) = RO*YF(2)*(S(23)+S(24)) 
     DXX(2,16) = RO*YF(2)*S(3) 
     DXX(2,17) = RO*YF(2)*S(7) 
     DXX(2,18) = RO*YF(2)*S(11) 
     DXX(2,19) = -RO*(-YF(2)*S(15)+YF(3)*S(16)) 
     DXX(2,20) = RO*YF(2)*S(19) 
     DXX(3,2) = -RO*(YF(17)*S(7)+YF(18)*S(11)+YF(7)*S(13)+YF(19)*S(15)+    &
          YF(15)*S(23))
     DXX(3,3) = RO*(YF(4)*S(2)+YF(16)*S(4)+YF(5)*S(6)+YF(17)*S(8)+YF(6)*   &
          S(10)+YF(18)*S(12)+YF(7)*S(14)+YF(19)*S(16)                      &
          +YF(8)*S(18)+YF(9)*(S(21)+S(22))+YF(10)*S(27))
     DXX(3,4) = RO*(-YF(3)*S(2)+YF(4)*S(40)) 
     DXX(3,5) = RO*YF(3)*S(6) 
     DXX(3,6) = RO*YF(3)*S(10) 
     DXX(3,7) = -RO*(-YF(3)*S(14)+YF(2)*S(13)) 
     DXX(3,8) = RO*YF(3)*S(18) 
     DXX(3,9) = RO*YF(3)*(S(21)+S(22)) 
     DXX(3,10) = RO*YF(3)*S(27) 
     DXX(3,16) = RO*YF(3)*S(4) 
     DXX(3,17) = -RO*(-YF(3)*S(8)+YF(2)*S(7)) 
     DXX(3,18) = -RO*(-YF(3)*S(12)+YF(2)*S(11)) 
     DXX(3,19) = -RO*(-YF(3)*S(16)+YF(2)*S(15)) 
     DXX(3,15) = -RO*YF(2)*S(23) 
     DXX(4,2) = -RO*(-Y(4)*S(1)+Y(17)*S(7))*CR 
     DXX(4,3) = RO*(Y(4)*S(2)-Y(3)**2*S(28)*RO/2.)*CR 
     DXX(4,4) = 1./PAS+RO*(Y(2)*S(1)+Y(3)*S(2)+2.*Y(4)*(S(25)+S(40)))*CR 
     DXX(4,17) = -RO*Y(2)*S(7)*CR 
     DXX(5,2) = -RO*(-Y(5)*S(5)+Y(16)*S(3)+Y(18)*S(11))*CR 
     DXX(5,3) = RO*Y(5)*S(6)*CR 
     DXX(5,5) = 1./PAS+RO*(Y(2)*S(5)+Y(3)*S(6))*CR 
     DXX(5,16) = -RO*Y(2)*S(3)*CR 
     DXX(5,18) = -RO*Y(2)*S(11)*CR 
     DXX(6,2) = -RO*(-Y(6)*S(9)+Y(19)*S(15))*CR 
     DXX(6,3) = -RO*(-Y(6)*S(10)+Y(4)*S(2)+Y(16)*S(4))*CR 
     DXX(6,4) = -RO*S(2)*Y(3)*CR 
     DXX(6,6) = 1./PAS+RO*(Y(2)*S(9)+Y(3)*S(10))*CR 
     DXX(6,16) = -RO*S(4)*Y(3)*CR 
     DXX(6,19) = -RO*S(15)*Y(2)*CR 
     DXX(7,2) = RO*Y(7)*S(13)*CR 
     DXX(7,3) = -RO*(-Y(7)*S(14)+Y(5)*S(6))*CR 
     DXX(7,5) = -RO*Y(3)*S(6)*CR 
     DXX(7,7) = 1./PAS+RO*(Y(3)*S(14)+Y(2)*S(13))*CR 
     DXX(8,2) = -RO*(-Y(8)*S(17)+Y(15)*S(23))*CR 
     DXX(8,3) = -RO*(-Y(8)*S(18)+Y(6)*S(10)+Y(18)*S(12))*CR 
     DXX(8,4) = -RO*Y(4)*S(40)*CR 
     DXX(8,6) = -RO*Y(3)*S(10)*CR 
     DXX(8,8) = 1./PAS+RO*(Y(2)*S(17)+Y(3)*S(18))*CR 
     DXX(8,15) = -RO*Y(2)*S(23)*CR 
     DXX(8,18) = -RO*Y(3)*S(12)*CR 
     DXX(9,2) = -RO*(-Y(9)*S(20)+Y(20)*S(19))*CR 
     DXX(9,3) = -RO*(-Y(9)*(S(21)+S(22))+Y(7)*S(14)+Y(19)*S(16))*CR 
     DXX(9,7) = -RO*Y(3)*S(14)*CR 
     DXX(9,9) = 1./PAS+RO*(Y(2)*S(20)+Y(3)*(S(21)+S(22)))*CR 
     DXX(9,19) = -RO*Y(3)*S(16)*CR 
     DXX(9,20) = -RO*Y(2)*S(19)*CR 
     DXX(10,2) = -RO*Y(15)*S(24)*CR 
     DXX(10,3) = RO*(-Y(8)*S(18)+Y(10)*S(27))*CR 
     DXX(10,8) = -RO*Y(3)*S(18)*CR 
     DXX(10,10) = 1./PAS+RO*Y(3)*S(27)*CR 
     DXX(10,15) = -RO*Y(2)*S(24)*CR 
     DXX(11,3) = -RO*Y(9)*S(21)*CR 
     DXX(11,9) = -RO*Y(3)*S(21)*CR 
     DXX(11,11) = 1./PAS 
     DXX(12,3) = -RO*Y(9)*S(22)*CR 
     DXX(12,9) = -RO*Y(3)*S(22)*CR 
     DXX(12,12) = 1./PAS 
     DXX(13,3) = -RO*Y(10)*S(27)*CR 
     DXX(13,10) = -RO*Y(3)*S(27)*CR 
     DXX(13,13) = 1./PAS 
     DXX(14,3) = -RO*(Y(16)*S(4)+Y(18)*S(12)+Y(9)*S(21))*CR 
     DXX(14,9) = -RO*Y(3)*S(21)*CR 
     DXX(14,14) = 1./PAS 
     DXX(14,16) = -RO*Y(3)*S(4)*CR 
     DXX(14,18) = -RO*Y(3)*S(12)*CR 
     DXX(16,2) = -RO*(-Y(16)*S(3)+Y(4)*S(1))*CR 
     DXX(16,3) = RO*Y(16)*S(4)*CR 
     DXX(16,4) = -RO*Y(2)*S(1)*CR 
     DXX(16,16) = 1./PAS+RO*(Y(2)*S(3)+Y(3)*S(4))*CR 
     DXX(15,2) = -RO*(-Y(15)*(S(23)+S(24))+Y(9)*S(20))*CR 
     DXX(15,4) = -RO*Y(4)*S(25)*CR 
     DXX(15,9) = -RO*Y(2)*S(20)*CR 
     DXX(15,15) = 1./PAS+RO*Y(2)*(S(23)+S(24))*CR 
     DXX(17,2) = -RO*(-Y(17)*S(7)+Y(5)*S(5)+Y(7)*S(13))*CR 
     DXX(17,3) = RO*Y(17)*S(8)*CR 
     DXX(17,5) = -RO*Y(2)*S(5)*CR 
     DXX(17,7) = -RO*Y(2)*S(13)*CR 
     DXX(17,17) = 1./PAS+RO*(Y(2)*S(7)+Y(3)*S(8))*CR 
     DXX(18,2) = -RO*(-Y(18)*S(11)+Y(6)*S(9))*CR 
     DXX(18,3) = RO*Y(18)*S(12)*CR 
     DXX(18,6) = -RO*Y(2)*S(9)*CR 
     DXX(18,18) = 1./PAS+RO*(Y(2)*S(11)+Y(3)*S(12))*CR 
     DXX(19,2) = RO*Y(19)*S(15)*CR 
     DXX(19,3) = -RO*(-Y(19)*S(16)+Y(17)*S(8))*CR 
     DXX(19,17) = -RO*Y(3)*S(8)*CR 
     DXX(19,19) = 1./PAS+RO*(Y(2)*S(15)+Y(3)*S(16))*CR 
     DXX(20,2) = -RO*(-Y(20)*S(19)+Y(8)*S(17))*CR 
     DXX(20,8) = -RO*Y(2)*S(17)*CR 
     DXX(20,20) = 1./PAS+RO*Y(2)*S(19)*CR 
     ! *************  CALCOLO DELLE DIFFERENZE FINITE  ********************* 
     LK = 0 
     do J=INI,IFI 
        KK =0 
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
        write(66,*)'31 - carbon'
        write(66,130) 
        stop 
     endif
     GAMMA = 0.0 
     KK=0 
     do J=INI,IFI 
        KK = KK+1 
        YF(J) = YF(J)+B(KK) 
        !*********** CALCOLA VARIAZIONI PERCENTUALI ELEMENTI IMPORTANTI ****
        if( YV(J) > 1.d-9 ) then 
           DD = abs(B(KK)/YF(J)) 
           if( DD > GAMMA ) then 
              IG = J 
              GAMMA = DD 
           endif
        endif
     end do
  end do

113 format(1X,'T =',1P,E10.3,2X,'RO =',1P,E10.3,2X,'MESH =',I4) 
130 format(1X,'ATTENZIONE QUALCOSA NON VA NELLA EPSI. IL DET = 0 NEL  &
       &RAHPSON/NEWTON')                                                  
250 format(1X,'MORTO: LA EPSI ITERA TROPPO NEL RAPHSON/NEWTON ') 
444 format(1X,'ITERAZIONI EPSI: DF =',I3,E10.3,' DX =',I3,E10.3,      &
       ' N. ITER.=',I3,' N. MESH=',I4)        
end subroutine CARBON
