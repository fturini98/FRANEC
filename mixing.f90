subroutine MIXING(MAXME, iread)
  use interfaccia
  use fisica
  use strut
  use nuconv
  use chimic
  use chim
  use overshoot
  use nummod
  use costanti
  use zone_conv

  implicit none

  integer :: MAXME, iread

  real,save,dimension(mele) :: XTOT 
  real,dimension(mele,lim) :: XSERV 
  real,save,dimension(2) :: GRAB
  real,save :: YECON 

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="mixing        "

  !*************IBLK=1 BLOCCA I PULSI: IBLK=0 NO *************************
  !* SE IBLK=1 ALLORA ISMOTH=1 AMMORBIDISCE I BP, ISMOTH = 0 BLOCCA I BP
  !* SE IVOV=1 ATTIVA UNDERSHOOTING ENVELOPE CONVETTIVO ******************
  integer :: LBCON = 2 , IBLK = 1, ISMOTH = 0
  integer ::  IVOV = 0

  integer :: k, igro, ibor, mcore, l1, lbit, kk, jj, l2, j, kl, ln
  integer :: max1, la, i, iconv, ia, menv, l, n, esci
  real :: ddm, zl, pp, tt, em, xh, xhe, tm, pr, rho, dad, csp
  real :: cap, drad, xxmax, grohe, gica, aconv, dmc
  integer :: meshstart, meshstop, jf
  logical, parameter :: test_version = .false.

  DDM = 0.
  K  = 0
  IGRO = 0
  IBOR = 0
  MCORE = 0
  GICO = 0.
  GSEMI = 0.
  L1 = MAXME
  LBIT = 1

  if(XXX(1,1) > 0.) then
     if(test_version) then
        do l=1,n_conv
           meshstart =  b_conv(1,l)
           meshstop =  b_conv(2,l)
           if(meshstop == maxme-1) meshstop = maxme

           do jf=meshstart,meshstop
              if(jf == meshstart) then
                 xtot = 0
                 dmc = 0
              endif

              if(meshstart > 1) then
                 do I = 1,MELE
                    XTOT(I) = XTOT(I) + XXX(I,jf)*(G(5,jf+1)-G(5,jf))
                 end do
                 DMC = DMC+G(5,jf+1)-G(5,jf)
              else
                 do I = 1,MELE
                    XTOT(I) = XXX(I,1)*G(5,2)
                 end do
                 DMC = G(5,2)
              endif
           end do

           do I = 1,MELE
              XXX(I,meshstart:meshstop) = XTOT(I)/DMC
           end do

        end do
        return
        ! torno via qui...
     endif
     ! se ho H al centro
     MENV = 1

     ! overshooting
     if(KOVER == 1 .and. NMD > 1) then
        if(G(6,1) > 0.0) then
           do LN = 1,L3
              G(6,LN) = 1.0
           end do
        endif
     endif

     if(G(6,1) > 0.) then
        L1 = 1
        do I = 1,MELE
           XTOT(I) = XXX(I,1)*G(5,2)
        end do
        DMC = G(5,2)
     else
        do I = 1,MELE
           XTOT(I) = 0.
        end do
        DMC = 0.
     endif
     MAX1 = MAXME-1
     do L = 2 , MAX1
        if(G(6,L) >= 0. .and. G(6,L-1) < 0.) then
           do I = 1,MELE
              XTOT(I) = XXX(I,L)*(G(5,L+1)-G(5,L))
           end do
           DMC = G(5,L+1)-G(5,L)
           L1 = L
        elseif(G(6,L) >= 0. .and. G(6,L-1) >= 0.) then
           do I = 1,MELE
              XTOT(I) = XTOT(I)+XXX(I,L)*(G(5,L+1)-G(5,L))
           end do
           DMC = DMC+G(5,L+1)-G(5,L)
           do N = L1,L
              do I = 1,MELE
                 XXX(I,N) = XTOT(I)/DMC
              end do
           end do
           if(L1 == 1) then
              GICO = G(5,L+1)
              MENV = L
           endif
        endif
     end do
     if(MENV == MAX1) then
        do I = 1,MELE
           XTOT(I) = XTOT(I)+XXX(I,MAX1)*(EMTOT-G(5,MAXME))
        end do
        DMC = DMC+EMTOT-G(5,MAXME)
        do N = L1,MAX1
           do I = 1,MELE
              XXX(I,N) = XTOT(I)/DMC
           end do
        end do

        ! metto la comp. al mesh maxme uguale a quella a maxme-1
        do I = 1,MELE
           XXX(I,MAXME) = XXX(I,MAXME-1)
        end do

     endif
     return
  endif

  if(G(6,1) < 0.) goto 25

  do KK=1,MAXME
     do JJ=1,MELE
        XSERV(JJ,KK) = XXX(JJ,KK)
     end do
  end do
  LBCON = 64
  LBIT = 64
  L1 = 1
  L2 = LBCON
  if( L2  <  2 ) L2 = 2

10 continue
  if( LBIT  <  1 ) LBIT = 1
  if( LBCON  <  2 ) LBCON = 2
  !************************ AZZERAMENTO VETTORI **************************
  do J = 1,MELE
     XTOT(J) = 0.
  end do
  DMC = 0.
  !************************** MIXING DA L1 A L2 **************************
  do K = L1,L2
     do J = 1  , MELE
        XTOT(J) = XTOT(J) + XXX(J,K) * (G(5,K+1) - G(5,K))
     end do
     DMC = DMC + G(5,K+1) - G(5,K)
  end do

  do
     do K = L1,L2
        do J = 1,MELE
           XXX(J,K) = XTOT(J)/DMC
        end do
     end do

     !******************** TEST: RAD > AD DA L1 A L2 ************************
     do K = L1,L2
        ZL = (G(2, K) + G(2, K+1))/2.
        PP = (G(3, K) + G(3, K+1))/2.
        TT = (G(4, K) + G(4, K+1))/2.
        EM = (G(5, K) + G(5, K+1))/2.
        do J = 1,MELE
           XX(J) = XXX(J, K)
        end do
        XH  = XX(1)
        XHE = XX(3)
        TM = TT * 1.d6
        PR = PP * 1.d17 - arad_3 * TM**4
        call STATE(caller,PR, TM, RHO, DAD, CSP, dum)
        call KAPPA(RHO, TM, CAP, XHE)
        !!DRAD = 39.46 * CAP * PP * ZL / ( EM * TT**4 )
        DRAD = cte_grad * CAP * PP * ZL / ( EM * TT**4 )
        esci = 0
        if( DRAD  < DAD ) then
           esci = 1
           exit
        endif
     end do
     if(esci == 1) exit
     !***************** AGGIUNGE LBIT MESH AL CONVETTIVO ********************
     if( L2 >= MAXME-1 ) return
     L2 = L2+LBIT

     if(LBIT > 1) goto 10
     do  J = 1 , MELE
        XTOT(J) = XTOT(J) + XXX(J,L2) * ( G(5,L2+1) - G(5,L2) )
     end do
     DMC = DMC + G(5,L2+1) - G(5,L2)
  end do

  !**************** MESH K RADIATIVO *************************************
  if(L1 == 1 .and. LBIT > 1) then
     L2 = L2-LBIT/2
     LBIT = LBIT/2
     do KK = 1,MAXME
        do JJ = 1,MELE
           XXX(JJ,KK) = XSERV(JJ,KK)
        end do
     end do
     goto 10
  end if
  if(L1 == 1) then
     GICO = G(5,K+1)
     MCORE = K+1
  endif
  XXMAX = XXV(3,1)
  !***********************************************************************
  if(XCE(1,1) <= 0. .and. L1 == 1 .and. XXX(3,1) > XXMAX .and. IBLK ==  1 &
       .and. XXMAX > 0. .and. XXMAX < 0.15) then
     IBOR = 1
  else
     IBOR = 0
  endif
  if(IBOR == 1) then
     do J = 1 , MELE
        XTOT(J) = 0.
     end do
     DMC = 0.
     do KK = 1,MAXME
        do JJ = 1,MELE
           XXX(JJ,KK) = XSERV(JJ,KK)
        end do
     end do
     do KL = 1 , L2
        do J = 1  , MELE
           XTOT(J) = XTOT(J) + XXX(J,KL) * (G(5,KL+1) - G(5,KL))
        end do
        DMC = DMC + G(5,KL+1) - G(5,KL)
        GROHE = XTOT(3)/DMC
        if(GROHE > XXV(3,1)) then
           if(ISMOTH == 1) then
              !  QUESTO BLOCCHETTO SERVE PER SMOOTTARE IL CONVETTIVO ********
              IGRO = KL
              GICA = GICO
              GICO = G(5,KL+1)
           else
              ! QUESTO BLOCCHETTO SERVE PER BLOCCARE IL CONVETTIVO *********
              IGRO = KL-1
              GICA = GICO
              GICO = G(5,KL)
           endif
           write(2,200) L2,K,GICA,KL,GICO
           goto 32
        endif
     end do
     write(66,*) '110 - mixing'
     stop
32   continue
     !********** QUESTO BLOCCHETTO SERVE PER BLOCCARE IL CONVETTIVO *********
     if(IBLK == 1 .and. ISMOTH == 0) then
        do J=1,MELE
           XTOT(J) = XTOT(J) - XXX(J,IGRO+1)*(G(5,IGRO+2)-G(5,IGRO+1))
        end do
        DMC = DMC - (G(5,IGRO+2)-G(5,IGRO+1))
     endif

     do KL = 1,IGRO
        do J = 1, MELE
           XXX(J,KL) = XTOT(J)/DMC
        end do
     end do
  endif
  !************************ AZZERAMENTO VETTORI **************************
25 continue
  do J = 1,MELE
     XTOT(J) = 0.
  end do
  DMC = 0.
  !******** CONTROLLA PRESENZA DI SHELL CONVETTIVE ***********************
  LBIT = 1
  MAX1 = MAXME - 1
  if(IBOR == 1) goto 43
  LA = K + 1
  I = 0
  do K = LA , MAX1
     if( G(6,1) < 0. .and. G(6,K) < 0.) cycle
     if( G(6,1) >= 0. .and. G(6,K) < 0. .and. K > L2) cycle
     if( G(6,1) < 0. .and. G(6,K) >= 0. .and. G(4,K) >= 80.) then
        L1 = K
        L2 = L1+1
        goto 10
     end if
     if( G(6,K) >= 0. .and. XXX(1,K) > 0.) cycle
     if(G(6,1) < 0. .and. G(6,K) >= 0. .and. G(4,K) < 80.) cycle
     ZL = ( G(2 , K) + G(2, K+1) )/2.
     PP = ( G(3 , K) + G(3, K+1) )/2.
     TT = ( G(4 , K) + G(4, K+1) )/2.
     EM = ( G(5 , K) + G(5, K+1) )/2.
     do J = 1,MELE
        XX(J) = XXX(J, K)
     end do
     XH  = XX( 1 )
     XHE = XX( 2 ) + XX( 3 )
     TM = TT * 1.d6
     PR = PP * 1.d17 - arad_3 * TM**4
     call STATE(caller,PR, TM, RHO, DAD, CSP, dum)
     call KAPPA(RHO, TM, CAP, XHE)
     !!DRAD = 39.46 * CAP * PP * ZL / ( EM * TT**4 )
     DRAD = cte_grad * CAP * PP * ZL / ( EM * TT**4 )
     if( DRAD < DAD ) cycle
     if(XXX(1,K) > 0.) cycle
     L1 = K
     if(G(6,1) < 0.) then
        L2 = L1+1
     else
        L2 = L2+1
        L2 = max(L2,L1+1)
     end if
     goto 10
  end do
  !************CONTROLLA PRESENZA ENVELOPE CONVETTIVO ********************
43 continue
  ICONV = 0
  do J=1,MAX1
     IA = MAX1-J+1
     if(G(6,IA) > 0.) then
        ICONV = ICONV+1
        L1 = IA
     else
        exit
     endif
  end do

  if(ICONV == 0) return
58 continue
  L2 = MAXME - 1
  YECON = G(5,L1)
  do K = L1,L2
     do J = 1  ,  MELE
        XTOT(J) = XTOT(J) + XXX(J,K) * ( G(5,K+1) - G(5,K) )
     end do
     DMC = DMC + G(5,K+1) - G(5,K)
  end do
  do K = L1,L2
     do J = 1, MELE
        XXX(J,K) = XTOT(J)/DMC
     end do
  end do

  !**** metto la composizione chimica a maxme uguale a quella a maxme-1
  if(iread /= 4) then
     do J=1,MELE
        XXX(J,MAXME) = XXX(J,MAXME-1)
     end do
  endif

  if(abs(DDM) > 1.d-10) return

  !********** CONTROLLO OVERSHOOTING ENVELOPE CONVETTIVO *****************
  if(IVOV == 1) then
     do K=1,2
        ZL = G(2 , L1)
        PP = G(3 , L1)
        TT = G(4 , L1)
        EM = G(5 , L1)
        do J = 1,MELE
           XX(J) = XXX(J,MAXME-1)
        end do
        XH  = XX(1)
        XHE = XX(2) + XX(3)
        TM = TT * 1.d6
        PR = PP * 1.d17 - arad_3 * TM**4
        call STATE(caller,PR, TM, RHO, DAD, CSP, dum)
        call KAPPA(RHO, TM, CAP, XHE)
        !!DRAD = 39.46 * CAP * PP * ZL / ( EM * TT**4 )
        DRAD = cte_grad * CAP * PP * ZL / ( EM * TT**4 )
        GRAB(K) = DRAD
        L1 = L1+1
     end do
     DDM = (DAD-GRAB(1))/(GRAB(2)-GRAB(1))*(G(5,L1-1)-G(5,L1-2))
     ACONV = G(5,L1-2)+DDM

     write(2,*) DDM,ACONV
     do K=1,MAXME
        if(G(5,K) > ACONV) then
           L1 = K
           if(DDM < -1.d-10) goto 58
        endif
     end do
     write(66,*)'111 - mixing'
     stop
  endif
  return

33 continue


200 format(1X,'************* ELIO CENTRALE SUPERATO *************',/, &
       ' MIX. FINO A:',I4,' BORDO CC:',I4,'M-CC:',1P,E12.5,' STOP A:',I4 &
       ,' M-CC:',1P,E12.4)                
end subroutine MIXING
  
