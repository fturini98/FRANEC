subroutine MIXING(MAXME, iread)
  use interfaccia
  use fisica
  use strut
  use nuconv
  use chimic
  use chim
  use overshoot
  use nummod
  use zone_conv
  use costanti
  implicit none

  integer :: MAXME, iread

  real,dimension(mele) :: XTOT 
  real,dimension(mele,lim) :: XSERV 
  real,save,dimension(2) :: GRAB
  real,save :: YECON 

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99

  !*************IBLK=1 BLOCCA I PULSI: IBLK=0 NO *************************
  !* SE IBLK=1 ALLORA ISMOTH=1 AMMORBIDISCE I BP, ISMOTH = 0 BLOCCA I BP
  !* SE IVOV=1 ATTIVA UNDERSHOOTING ENVELOPE CONVETTIVO ******************
  integer :: LBCON = 2 , IBLK = 1, ISMOTH = 0
  integer ::  IVOV = 0

  integer :: k, igro, ibor, mcore, l1, lbit, kk, jj, l2, j, kl, ln, l2old
  integer :: max1, la, i, iconv, ia, menv, l, n, esci, stp, mycore, ext
  real :: ddm, zl, pp, tt, em, xh, xhe, tm, pr, rho, dad, csp
  real :: cap, drad, xxmax, grohe, gica, aconv, dmc

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
        XTOT = 0.
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

  write(65,*) "START", nmd, n_conv
  write(65, *) b_conv(1:2, 1:n_conv)
  write(65, *) g(5,b_conv(2,1))/EMTOT, g(5,b_conv(2,1)+1)/EMTOT
  write(65, *) g(6,b_conv(2,1)), g(6,b_conv(2,1)+1)
  
  xserv(1:mele,1:maxme) = xxx(1:mele,1:maxme)
!!$  l1 = 1
!!$  l2 = b_conv(2,1)
!!$  call do_mix_core(l1, l2, maxme)
!!$  
!!$  b_conv(2,1) = l2
!!$  do i=1,n_conv
!!$     XTOT = 0.
!!$     DMC = 0.
!!$     do K = b_conv(1,i), b_conv(2,i)
!!$        do J = 1  , MELE
!!$           XTOT(J) = XTOT(J) + Xserv(J,K) * (G(5,K+1) - G(5,K))
!!$        end do
!!$        DMC = DMC + G(5,K+1) - G(5,K)
!!$     end do
!!$  
!!$     do K = b_conv(1,i), b_conv(2,i)
!!$        do J = 1  , MELE
!!$           XXX(j,k) = XTOT(j)/dmc
!!$        end do
!!$     end do
!!$     
!!$  end do
!!$  GICO = G(5,b_conv(2,1)+1)
!!$  return

!  LBCON = 64
  LBIT = 64
  LBCON = b_conv(2,1)
  L1 = 1
  L2 = LBCON
  l2old = l2
  if( L2  <  2 ) L2 = 2

!!$  call do_mix_core(l1, l2, maxme)
!!$  call mixa(l1,l2)
!!$  mycore = l2
!!$  lbit = 1
!!$  k = l2
!!$  la = 2
!!$  goto 26
  
10 continue
  if( LBIT  <  1 ) LBIT = 1
  if( LBCON  <  2 ) LBCON = 2
  !************************ AZZERAMENTO VETTORI **************************
  XTOT = 0.
  DMC = 0.
  !************************** MIXING DA L1 A L2 **************************
  do K = L1,L2
     do J = 1  , MELE
        XTOT(J) = XTOT(J) + XXX(J,K) * (G(5,K+1) - G(5,K))
     end do
     DMC = DMC + G(5,K+1) - G(5,K)
  end do
  write(65,*) "mix ", l1, l2

  do
     do K = L1,l2
        do J = 1,MELE
           XXX(J,K) = XTOT(J)/DMC
        end do
     end do

     !******************** TEST: RAD > AD DA L1 A L2 ************************
     do K = L2,L1,-1
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
        call STATE(PR, TM, RHO, DAD, CSP, dum)
        call KAPPA(RHO, TM, CAP, XHE)
        !!DRAD = 39.46 * CAP * PP * ZL / ( EM * TT**4 )
        DRAD = cte_grad * CAP * PP * ZL / ( EM * TT**4 )
        write(65,*) "last ", l2, k, drad - dad
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
  write(65,*) "exit ", l1,l2, lbit

  
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
  endif
26 continue
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
              write(65,*) "reset ", igro, gico/EMTOT
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
  if(l1 == 1) write(65,*) "gico", gico/EMTOT, l1,l2
  write(65,*) "shell", k, la
  XTOT = 0.
  DMC = 0.
  !******** CONTROLLA PRESENZA DI SHELL CONVETTIVE ***********************
  LBIT = 1
  MAX1 = MAXME - 1
  if(IBOR == 1) goto 43
  LA = K + 1
  I = 0
  do K = LA , MAX1
     write(65,'("shell 2 ",2i5,2e10.3)') k, l2, G(6,1), g(6,k)
     if( G(6,1) < 0. .and. G(6,K) < 0.) cycle
     if( G(6,1) >= 0. .and. G(6,K) < 0. .and. K > L2) cycle
     if( G(6,1) < 0. .and. G(6,K) >= 0. .and. G(4,K) >= 80.) then
        write(65,*) "g4"
        L1 = K
        L2 = L1+1
        goto 10
     end if
     write(65,'("shell 2b ",i5,e10.3)') k, xxx(1,k)
     if( G(6,K) >= 0. .and. XXX(1,K) > 0.) cycle
     if(G(6,1) < 0. .and. G(6,K) >= 0. .and. G(4,K) < 80.) cycle
     ZL = ( G(2 , K) + G(2, K+1) )/2.
     PP = ( G(3 , K) + G(3, K+1) )/2.
     TT = ( G(4 , K) + G(4, K+1) )/2.
     EM = ( G(5 , K) + G(5, K+1) )/2.
     do J = 1,MELE
        XX(J) = XXX(J, K)
     end do
     write(65,*) "shell 3", k
     XH  = XX( 1 )
     XHE = XX( 2 ) + XX( 3 )
     TM = TT * 1.d6
     PR = PP * 1.d17 - arad_3 * TM**4
     call STATE(PR, TM, RHO, DAD, CSP, dum)
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
  write(65,*) "envelop"
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
  write(65,*) "env2", l1,l2
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
        call STATE(PR, TM, RHO, DAD, CSP, dum)
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
  

subroutine do_mix_core(l1, l2, maxme)
  use interfaccia
  use fisica
  use strut
  use nuconv
  use chimic
  use chim
  use overshoot
  use nummod
  use zone_conv

  implicit none

  integer :: l1, l2, maxme

  real,dimension(mele,lim) :: XSERV 

  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99

  integer :: k, j
  integer :: esci
  real :: zl, pp, tt, em, xh, xhe, tm, pr, rho, dad, csp
  real :: cap, drad
  logical:: convettivo

  convettivo = .true.
  Xserv(1:mele,1:maxme) = XXX(1:mele,1:maxme)

  do while(convettivo) 
     !************************** MIXING DA L1 A L2 **************************
     
     call mixa(l1, l2)
!!$     XTOT = 0.
!!$     DMC = 0.
!!$
!!$     do K = L1,L2
!!$        do J = 1, MELE
!!$           XTOT(J) = XTOT(J) + XXX(J,K) * (G(5,K+1) - G(5,K))
!!$        end do
!!$        DMC = DMC + G(5,K+1) - G(5,K)
!!$     end do
!!$     write(65,*) "mix ", l1, l2
!!$
!!$     do K = L1,L2
!!$        do J = 1  , MELE
!!$           XXX(j,k) = XTOT(j)/dmc
!!$        end do
!!$     end do
     !******************** TEST: RAD > AD DA L1 A L2 ************************
     do K = L2,L1,-1
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
        call STATE(PR, TM, RHO, DAD, CSP, dum)
        call KAPPA(RHO, TM, CAP, XHE)
        !!DRAD = 39.46 * CAP * PP * ZL / ( EM * TT**4 )
        DRAD = cte_grad * CAP * PP * ZL / ( EM * TT**4 )
        write(65,*) "last ", l2, k, drad - dad
        esci = 0
        if( DRAD  < DAD ) then
           convettivo = .false.
           exit
        endif
     end do

     !***************** AGGIUNGE 1 MESH AL CONVETTIVO ********************
     XXX(1:mele,l1:l2) = Xserv(1:mele,l1:l2)
     if(convettivo) L2 = L2+1
  end do

  XXX(1:mele,l1:l2) = Xserv(1:mele,l1:l2)
  
end subroutine do_mix_core



subroutine mixa(i1, i2)
  use interfaccia
  use fisica
  use strut
  use chimic
  use chim
    
  implicit none
  
  integer :: i1, i2

  integer :: k, j
  real,dimension(mele) :: XTOT 
  real :: dmc

  XTOT = 0.
  DMC = 0.
  do K = i1, i2
     do J = 1, MELE
        XTOT(J) = XTOT(J) + XXX(J,K) * (G(5,K+1) - G(5,K))
     end do
     DMC = DMC + G(5,K+1) - G(5,K)
  end do
  
  do K = i1, i2
     do J = 1, MELE
        XXX(j,k) = XTOT(j)/dmc
     end do
  end do

end subroutine mixa
