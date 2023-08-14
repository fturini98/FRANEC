! *****************************************************************
! ATTENZIONE: cambiare le referenze in ref_cross se cambiano le cross!!!!
! *****************************************************************

subroutine CROSS (RO,T, XX, IR,ISHORT, S, powt, nabla)
  use interfaccia
  use parametri
  use sceltachim

  implicit none

  real :: RO, T
  real,dimension(mele) :: xx
  integer :: IR, ISHORT, nabla
  real,dimension(nsezi) :: S
  type(pow_of_T) :: powt

  integer :: ncache, Z1, Z2, i
  real,dimension(32) :: schermo
  real :: T94, T914, T943, T953, tmpvar, sk1, seb1, sk2, seb2
  real :: skc12, skn14, sko, skneon, g20a, gt9, t9a, t9a13, t9a56, t9a23
  real :: ft9a, fpt9a, skmg, skc, skca, skop, skoa, skfa, sknap, skcc
  real :: skfp, sknep, sbep, sbee, skbe7, skli, skb, cbe1, cbe2
  real :: T9, T92, T93, T912, T913, T923, T932, T954, T915, T6, T612
  integer, save :: firsttime = 1

  TMPVAR = 0.1

  do i=1,nsezi
     S(i) = 0.
  end do

  T9 = powt%T9
  T92 = powt%T92
  T93 = powt%T93
  T912 = powt%T912
  T913 = powt%T913
  T923 = powt%T923
  T932 = powt%T932
  T6 = powt%T6
  T612 = powt%T612

  T94 = T93*T9
  T914 = sqrt(T912)
  T943 = T923*T923
  T953 = T943*T913

  ncache = 0
  
  if(firsttime == 1) then
     ! Stampo le referenze delle cross
     ! Ricordarsi di modificarle se si cambia una sezione d'urto !
     call ref_cross
     firsttime = 0
  endif

  ! CALCOLO SEZIONI D'URTO
  if(IR <= 2 .or. stdchim == 0) then
     !     ******************  PPI: P + P => D + P => HE3 + G  **************
     Z1 = 1
     Z2 = 1
     SCHERMO(1) = SK(RO,T,XX,Z1,Z2,ncache)
     ncache = 1
     ! GIADA:  rate NACRE
     S(29) = SCHERMO(1)*(4.08d-15/T923*exp(-3.381/T913) * &
          (1.0 + 3.82*T9 + 1.51*T92 +0.144*T93 -1.14d-2*T94))
    

     !     *****************   PPI: HE3 + HE3 => HE4 + 2P     ***************
     Z1 = 2
     Z2 = 2
     ! GIADA:  rate NACRE
     SCHERMO(2) = SK(RO,T,XX,Z1,Z2,ncache)
     S(30) = SCHERMO(2)* (5.59d10/T923 * exp(-12.277/T913) * &
          (1.0 - 0.135*T9 +2.54d-2*T92 - 1.29d-3*T93 ))

     !     ****************** PPII : HE3 + HE4 => BE7 + G ECC. **************
     Z1 = 2
     Z2 = 2
     SCHERMO(3) = SCHERMO(2)
     ! GIADA: rate di Cyburt e Davids 2008, Phys.Rev.C 78, 064614
     S(31) = SCHERMO(3)*(exp(15.609867-12.82707707/T913-2./3.*log(T9)) * &
          ((1-0.020478*T923+0.211995*T943)/(1+0.255059*T923+0.338573*T943)))


     !     ******************   N15 + P => O16 + G  *************************
     Z1 = 7
     Z2 = 1
     SCHERMO(4) = SK(RO,T,XX,Z1,Z2,ncache)
     ! GIADA:  rate NACRE : controlla Q = 12.127
     if(T9 <= 3.5) then
        S(32) = SCHERMO(4) * (1.08d9/T923*exp(-15.254/T913 -(T9/0.34)**2 ) &
             * (1. +6.15*T9 +16.4*T92 ) + 9.23d3/T932 * exp(-3.597/T9)+ &
             3.27d6/T932 * exp(-11.024/T9) )
     else
        S(32) = SCHERMO(4) * (3.54d4*T9**0.095*exp(-2.306/T9))
     end if

     ! ****************  CALCOLO Q DEL BRANCHING DEL BE7 ******************* 
     Z1 = 4 
     Z2 = 1 
     skbe7 = SK(RO,T,XX,Z1,Z2,ncache)
     ! GIADA: rate NACRE    
     SBEP = skbe7*( 2.61d5/T923*exp(-10.264/T913) *   &
          (1.0 -5.11d-2*T9 + 4.68d-2*T92- 6.60d-3*T93 +3.12d-4*T9**4) +  &
          2.05d3/T932*exp(-7.345/T9) )
     ! rate Adelberger 1998 
     ! Rev.Mod.Phys. vol.70, n.4, 1265
     SBEE = 5.60d-9/T612 * (1.0+0.004*(T6-16.0))

     S(33) = SBEP*XX(2)/(SBEE*(1.+XX(2))/2.)

     ! ============== D, Li, Be, B ========================
     ! ------------------- D + p -> He3 + g
     ! serve sia non schermato S(34) che schermato S(35)
     ! rate NACRE
     if(T9 < 0.11) then
        S(34) = 1.81d3/T923*exp(-3.721/T913)*(1.+14.3*T9-90.5*T92+ &
             395.*T93)
     else
        S(34) = 2.58d3/T923*exp(-3.721/T913)*(1.+3.96*T9+0.116*T92)
     endif
     S(35) = SCHERMO(1)*S(34)

     if(nabla == 1) then
        !********  6Li+p ----> 3He + alfa
        Z1 = 3
        Z2 = 1
        skli = SK(RO,T,XX,Z1,Z2,ncache)
        ! rate NACRE
        S(36) = skli*(3.54d10/T923*exp(-8.415/T913)*(1.-0.137*T9+ &
             2.41d-2*T92-1.28d-3*T93))

        !*********  7Li +p ----> 4He + alfa
        ! rate NACRE
        S(37) = skli*(7.20d8/T923*exp(-8.473/T913 -(T9/6.5)**2)*(1.+1.05*T9 &
             -0.653*T92+0.185*T93-2.12d-2*T9**4+9.30d-4*T9**5)+ &
             9.85d6*T9**0.576*exp(-10.415/T9))

        ! ------------------- Be9 + p -> Li6 + alpha, con la reazione
        ! ------------------- Be9 + p -> Be8 + D -> 2alpha + D
        !*********** 9Be + p ----> 6Li + alfa
        ! rate NACRE
        CBE1 = 2.11d11/T923*exp(-10.361/T913-(T9/0.4)**2)*(1.-0.189*T9 &
             + 35.2*T92)+5.24d8/T932*exp(-3.446/T9)+4.65d8/T9**0.293* &
             exp(-4.396/T9)
        !*********** 9Be + p ----> 8Be+d ----> 2alfa+d  
        ! rate NACRE
        CBE2 = 2.18d11/T923*exp(-10.361/T913 -(T9/0.42)**2)*(1.-0.427*T9 &
             +34.055*T92)+6.24d8/T932*exp(-3.446/T9)+3.53d8/T9**0.205 &
             *exp(-3.889/T9)
        S(38) = skbe7*(CBE1+CBE2)

        !******** 11B + p ---> alfa + 2alfa
        Z1 = 1
        Z2 = 5
        skb = SK(RO,T,XX,Z1,Z2,ncache)
        ! rate NACRE
        S(39) = skb*(2.68d12/T923*exp(-12.097/T913)*(1.+1.62*T9-1.31*T92+ &
             0.260*T93)+2.12d6/T932*exp(-1.724/T9))
     endif
  endif                     ! IF(IR <= 2)

  if(IR /= 3 .or. stdchim == 0) then
     !     *******************  C12 + P => N13 + G  *************************
     Z1 = 6
     Z2 = 1
     SCHERMO(5) = SK(RO,T,XX,Z1,Z2,ncache)
     ncache = 1
     ! GIADA:  rate NACRE
     S(1) = SCHERMO(5) * (2.0d7/T923 * exp(-13.692/T913 -       &
          (T9/0.46)**2 ) * ( 1.0+9.89*T9-59.8*T92 +266.*T93)+   &
          1.0d5/T932*exp(-4.913/T9) + 4.24d5/T932 * exp(-21.62/T9) )

     !     *******************  N14 + P => O15 + G  **************************
     Z1 = 7
     Z2 = 1
     if(IR <= 2) then
        SCHERMO(6) = SCHERMO(4) 
     else
        SCHERMO(6) = SK(RO,T,XX,Z1,Z2,ncache)
     end if
     ! GIADA: rate Imbriani et al. 2005, Eur.Phys. J.A 25, 455
     S(5) = SCHERMO(6)*(3.12d7/T923*exp(-15.193/T913-(T9/0.486)**2)*        &
          (0.782-1.50*T9+17.97*T92-3.32*T93) +2.11d3/T932*exp(-2.998/T9)+   &
          8.42d2*(T9**0.0682)*exp(-4.891/T9))

     !     ****************    N15 + P => A + C12   *************************
     Z1 = 7
     Z2 = 1
     SCHERMO(7) = SCHERMO(6)
     !   GIADA:  rate NACRE
     if(T9 <= 2.5)then
        S(7) = SCHERMO(7) * (1.12d12/T923*exp(-15.253/T913 -  &
             (T9/0.28)**2) * (1.0+4.95*T9+143.0*T92) + &
             1.01d8/T932*exp(-3.643/T9) + 1.19d9/T932 * exp(-7.406/T9) )
     else
        S(7) = SCHERMO(7) * ( 4.17d7*T9**0.917*exp(-3.292/T9))
     end if
     ! jina
!!$     S(7) = exp(2.747640e+01 -1.525300e+01/T913 + &
!!$          1.593180e+00*T913 + &
!!$          2.447900e+00*T9 -2.197080e+00*T9**(5/3)-6.666670e-01*log(T9))
!!$     S(7) = S(7) + exp(-4.873470e+00 -2.021170e+00/T9 + 3.084970e+01*T913-8.504330e+00*T9-1.544260e+00*T9**(5/3)-1.500000e+00*log(T9))
!!$     S(7) = S(7) + exp(2.089720e+01-7.406000e+00/T913-1.500000e+00*log(T9))
!!$     S(7) = S(7) + exp(-6.575220e+00-1.163800e+00/T9+2.271050e+01*T913-2.907070e+00*T9+2.057540e-01*T9**(5/3)-1.500000e+00*log(T9))
!!$     S(7) = S(7) * SCHERMO(7)


     ! *****************   O16 + P => F17 + G    *************************
     Z1 = 8
     Z2 = 1
     SCHERMO(8) = SK(RO,T,XX,Z1,Z2,ncache)
     !   GIADA:  rate NACRE
     S(9) = SCHERMO(8)*7.37d7*exp(-16.696/T913)*T9**(-0.82)
  endif

  if(IR == 2 .or. IR == 3 .or. stdchim == 0) then

     !     ****************    HE4 + 2A => C12 + G *************************
     Z1 = 2
     Z2 = 2
     if(IR <= 2) then
        SK1 = SCHERMO(2)
     else
        SK1 = SK(RO,T,XX,Z1,Z2,0)
        ncache = 1
     endif
     ! GIADA: rate NACRE
     SEB1 = SK1*(2.43d9/T923*exp(-13.490/T913-(T9/0.15)**2)* &
          (1.0+74.5*T9)+6.09d5/T932*exp(-1.054/T9)) 

     Z1 = 2
     Z2 = 4
     SK2 = SK(RO,T,XX,Z1,Z2,ncache)
     SEB2 = SK2*(2.76d7/T923*exp(-23.570/T913-(T9/0.4)**2)* &
          (1.0+5.47*T9+326*T92)+130.7/T932*exp(-3.338/T9)+   &
          2.51d4/T932*exp(-20.307/T9))

     if(T9 <= 0.03)then
        S(28) = SEB1*SEB2*3.07d-16*(1.0-29.1*T9+1308*T92)
     else
        S(28) = SEB1*SEB2*3.44d-16*(1.0+0.0158*T9**(-0.65))   
     end if
  endif

  if(IR > 1 .or. stdchim == 0) then
     !     ********************  C12 + A => O16 + G  ************************
     Z1 = 6
     Z2 = 2
     SKC12 = SK(RO,T,XX,Z1,Z2,ncache)
     ncache = 1
     ! GIADA: rate Hammer et al. 2005, Nucl.Phys.A 758, 363
     S(2) = SKC12*(1.51d8/(T92*(1+0.0666/T923)**2)*            &
          exp(-32.12/T913-(T9/1.03)**2)+(1.11d9/((T92)*        &
          (1.0+0.735/T923)**2)*exp(-32.12/T913)+16200.0/T923*  &
          (1.0+2.19d6*T913)*exp(-38.814/T913)))

     !     ****************    N14 + A => F18 + G   ***********************
     Z1 = 7
     Z2 = 2
     SKN14 = SK(RO,T,XX,Z1,Z2,ncache)
     ! GIADA: rate NACRE
     if(T9 <= 2)then
        S(6) = SKN14*(7.93d11/T923*exp(-36.035/T913    &
             -(T9/0.07)**2)+1.85d-10/T932*exp(-2.750/T9)+ &
             2.62/T932*exp(-5.045/T9)+2.93d3*T9**0.344*exp(-10.561/T9))
     else
        S(6) = SKN14*(1.52d2*T9**1.567*exp(-6.315/T9))
     end if

     !     *****************   O16 + A => NE20 + G   *************************
     Z1 = 8
     Z2 = 2
     SKO = SK(RO,T,XX,Z1,Z2,ncache)
     ! GIADA: rate NACRE
     S(10) = SKO*(2.68d10/T923*exp(-39.76/T913-(T9/1.6)**2)+   &
          51.1/T932*exp(-10.32/T9)+616.1/T932*exp(-12.2/T9)+   &
          0.41*T9**2.966*exp(-11.9/T9))

     ! ****************    O18 + A => NE22 + G     *************************
     Z1 = 8
     Z2 = 2
     ! GIADA: rate NACRE
     if(T9 <= 6) then
        S(14) = SKO*(1.95d-13/T932*exp(-2.069/T9)+1.56d-2/T932*  &
             exp(-4.462/T9)+10.1/T932*exp(-6.391/T9)+            &
             44.1/T932*exp(-7.389/T9)+3.44d5/T912*exp(-22.103/T9))
     else
        S(14) = SKO*(3.31d5*T9**(-0.221)*exp(-24.99/T9))
     end if

     if(IR == 2 .or. (IR == 3 .and. ISHORT == 1)) return

     !     *****************   NE20 + A => MG24 + G   *************************
     Z1 = 10
     Z2 = 2
     SKNEON = SK(RO,T,XX,Z1,Z2,ncache)
     ! rate NACRE
     if(T9 <= 1) then
        S(18) = SKNEON*(8.72*T9**(-0.532)*exp(-8.995/T9))
     else
        S(18) = SKNEON*(3.72d2*T9**(2.229)*exp(-12.681/T9))
     endif

     !     *****************   NE22 + A => MG25 + N    ************************
     Z1 = 10
     Z2 = 2
     ! rate NACRE
     if(T9 <= 2) then
        S(21) = SKNEON*(7.4*exp(-7.79/T9)+1.3d-4*T9**0.83*exp(-5.52/T9)+ &
             9.41d3*T9**2.78*exp(-11.7/T9)+8.59d6*T9**0.892*exp(-24.4/T9))
     else
        S(21) = SKNEON*(1.51d5*T9**2.879*exp(-16.717/T9))
     endif

     !     *******************     NE22 + A => MG26 + G  ********************
     Z1 = 10
     Z2 = 2
     ! rate NACRE
     if(T9 <= 1.25) then
        S(22) = SKNEON*(3.55d-9/T932*exp(-3.927/T9) + 7.07d-1*T9**(-1.064)* &
             exp(-7.759/T9) + 1.27d-3*T9**(-2.556)*exp(-6.555/T9))
     else
        S(22) = SKNEON*(1.76*T9**3.322*exp(-12.412/T9))
     endif

     !     ******************* MG24 + A => SI28 + G *************************
     Z1 = 12
     Z2 = 2
     SKMG = SK(RO,T,XX,Z1,Z2,ncache)
     ! rate CF 1988
     GT9 = 1.+5.*exp(-15.882/T9)
     S(27) = SKMG*((4.78d1/T932*exp(-13.506/T9)+2.38d3/T932   &
          *exp(-15.218/T9)+2.47d2*T932*exp(-15.147/T9)+       &
          TMPVAR*1.72d-9/T932*exp(-5.028/T9)+1.25d-3/T932*    &
          exp(-7.929/T9)+2.43d1/T9*exp(-11.523/T9))/GT9)
  endif

  if(IR == 4) then
     !     ********************  C13 + P => N14 + G  *************************
     Z1 = 6
     Z2 = 1
     SKC = SCHERMO(5)
     !  GIADA: rate NACRE
     S(3) = SKC*(9.57d7/T923*(1+3.56*T9)*exp(-13.720/T913-T92)+ &
          1.50d6/T932*exp(-5.930/T9)+6.83d5*T9**(-0.864)*exp(-12.057/T9))

     !     *******************  C13 + A => O16 + N ***************************
     Z1 = 6
     Z2 = 2
     SKCA = SK(RO,T,XX,Z1,Z2,ncache)
     ncache = 1
     !  GIADA: rate NACRE
     if(T9 <= 4) then
        S(4) = SKCA*(3.78d14/T92*exp(-32.333/T913-(T9/0.71)**2)*   &
             (1.+46.8*T9-292*T92+738*T93)+2.30d7*T9**0.45*exp(-13.03/T9))
     else
        S(4) = SKCA*(7.59d6*T9**1.078*exp(-12.056/T9))
     end if

     !     *****************   N15 + A => F19 + G    *************************
     Z1 = 7
     Z2 = 2
     !  GIADA: rate NACRE
     S(8) = SKN14*(1.1d11/T923*exp(-36.214/T913-(T9/0.6)**2)+  &
          1.65d-4/T932*exp(-4.224/T9)+2.66/T932*exp(-6.22/T9)+ &
          1.56d2/T932*exp(-7.764/T9)+3.92d4*T9**(-0.333)*exp(-14.522/T9))

     !     ****************    O17 + P => N14 + A      ***********************
     Z1 = 8
     Z2 = 1
     SKOP = SCHERMO(8)
     ! GIADA: rate NACRE
     if(T9 <= 6) then
        S(11) = SKOP*(9.20d8/T923*exp(-16.715/T913-(T9/0.06)**2)*   &
             (1-80.31*T9+2211*T92)+9.13d-4/T932*exp(-0.7667/T9)+    &
             9.68/T932*exp(-2.083/T9)+8.13d6/T932*exp(-5.685/T9)+   &
             1.85d6*T9**1.591*exp(-4.848/T9))
     else
        S(11) = SKOP*(8.73d6*T9**0.950*exp(-7.508/T9))
     end if

     !     ****************    O17 + A => NE20 + N     ***********************
     Z1 = 8
     Z2 = 2
     SKOA = SKO
     ! GIADA: rate NACRE
     S(12) = SKOA*(4.38d17/T923*exp(-39.918/T913-(T9/1.1)**2)+  &
          1.73d3/T932*exp(-8.55/T9)+7.5d5*T9**1.83*exp(-13.8/T9))

     !     ****************    O18 + P => N15 + A     *************************
     Z1 = 8
     Z2 = 1
     ! GIADA: rate NACRE : controllare Q = 3.981
     S(13) = SKOP*(5.58d11/T923*exp(-16.732/T913-(T9/0.51)**2) &
          *(1+3.2*T9+21.8*T92)+9.91d-14/T932*exp(-0.232/T9)+    &
          2.58d4/T932*exp(-1.665/T9)+3.24d8*T9**(-0.378)*exp(-6.395/T9))

     !     ****************    F19 + P => O16 + A      ***********************
     Z1 = 9
     Z2 = 1
     SKFP = SK(RO,T,XX,Z1,Z2,ncache)
     ! GIADA: rate NACRE
     S(15) = SKFP*(2.62d11/T923*exp(-18.116*T913-(T9/0.185)**2)*        &
          (1.+6.26d-2*T9+0.285*T92+4.94d-3*T93+11.5*T9**4+7.4d4*T9**5)+ &
          3.8d6/T932*exp(-3.752/T9)+3.27d7*T9**(-0.193)*exp(-6.587/T9)+ &
          7.3d8*T9**(-0.201)*exp(-16.249/T9))

     !     *****************   F19 + A => NE22 + P    ************************
     Z1 = 9
     Z2 = 2
     SKFA = SK(RO,T,XX,Z1,Z2,ncache)
     S(16) = SKFA*(4.50d18/T923*exp(-43.467/T913-(T9/0.637)**2)+ &
          7.98d4*T932*exp(-12.760/T9))

     !     *****************   NE20 + P => NA21 + G   *************************
     Z1 = 10
     Z2 = 1
     SKNEP = SK(RO,T,XX,Z1,Z2,ncache)
     ! rate NACRE
     S(17) = SKNEP*(2.35d7/T9**(-1.84)*exp(-19.451/T913)*(1.0+10.80*T9) + &
          18.0/T932*exp(-4.247/T9) + 9.83/T932*exp(-4.619/T9) + &
          6.74d4*T9**(-0.641)*exp(-11.922/T9))

     !     *****************   NE21 + P => NA22 + G   *************************
     Z1 = 10
     Z2 = 1
     ! rate NACRE
     if(T9 <= 2.) then
        S(19) = SKNEP*(4.68d8/T923*exp(-19.465/T913 - (T9/0.2)**2) + &
             8.18d-4/T932*exp(-1.085/T9)+6.11/T932*exp(-1.399/T9) + &
             1.34d4/T932*exp(-3.009/T9)+1.26d5*T9**(-0.128)*exp(-4.962/T9))
     else
        S(19) = SKNEP*(3.04d4*T9**(0.42)*exp(-2.650/T9))
     endif

     ! *****************   NE22 + P => NA23 + G   *************************
     Z1 = 10
     Z2 = 1
     ! rate NACRE
     if(T9 <= 2.) then
        S(20) = SKNEP*(1.11d-9/T932*exp(-0.422/T9)+6.83d-5/T932* &
             exp(-0.810/T9)+9.76d-3/T932*exp(-1.187/T9)+ &
             1.06d-1/T932*exp(-1.775/T9)+8.51d4*T9**(0.725)*exp(-4.315/T9))
     else
        S(20) = SKNEP*(6.30d4*T9**(0.816)*exp(-3.910/T9))
     endif

     !     *******************     NA23 + P => NE20 + A  ********************
     Z1 = 11
     Z2 = 1
     SKNAP = SK(RO,T,XX,Z1,Z2,ncache)
     ! rate NACRE
     if(T9 <= 5.) then
        S(23) = SKNAP*(8.39d9/T923*exp(-20.770/T913 - (T9/0.1)**2)* &
             (1.0+45.2*T9)+3.09d-13/T932*exp(-0.420/T9)+8.12d-3/T932* &
             exp(-1.601/T9)+4.37/T932*exp(-1.934/T9)+7.50d3*T9**(-1.48)* &
             exp(-3.150/T9)+1.05d6*T9**(1.456)*exp(-4.482/T9))
     else
        S(23) = SKNAP*(3.96d6*T9**(1.291)*exp(-9.277/T9))
     endif

     !     *******************     NA23 + P => MG24 + G  *********************
     Z1 = 11
     Z2 = 1    
     ! rate NACRE
     if(T9 <= 5.) then
        S(24) = SKNAP*(9.55d7/T923*exp(-20.770/T913 - (T9/0.3)**2)* &
             (1.0-10.80*T9+61.08*T92)+8.20d-2/T932*exp(-1.601/T9)+ &
             85.2/T932*exp(-2.808/T9)+1.70d4/T932*exp(-3.458/T9)+ &
             5.94d4*exp(-5.734/T9))
     else
        S(24) = SKNAP*(5.60d3*T9**(1.112)*exp(-2.337/T9))
     endif

     !     *******************  C12 + C12 ----- NE20  NA23  ******************
     Z1 = 6
     Z2 = 6
     SKCC = SK(RO,T,XX,Z1,Z2,ncache)
     !     *******************  C12 + C12 => P + NA23  ******************
     ! rate JINA CF88
     S(25) = SKCC*exp(6.096490d1 -8.416500d1/T913 -1.4191*T913 &
          -1.146190d-1*T9 -7.030700d-2*T953 -6.666670d-1*log(T9))

     !     *******************  C12 + C12 => A + NE20  ******************
     ! rate JINA CF88
     S(40) = SKCC*exp(6.128630d1 -8.416500d1/T913 -1.56627*T913	&
          -7.360840d-2*T9 -7.279700d-2*T953 -6.666670d-1*log(T9))
  endif

  return
end subroutine CROSS


! Questa routine stampa su runlog le referenze delle cross usate
! *****************************************************************
! ATTENZIONE: cambiare qui le referenze se cambiano le cross!!!!
! *****************************************************************
subroutine ref_cross

  implicit none

  ! Numero di reazioni.
  ! Aggiornare se se ne introducono altre.
  integer,parameter :: sezi = 39

  character(len=55),dimension(sezi) :: lab, ref
  integer :: i

  lab(1) = "P + P => D + P => He3 + G"
  ref(1) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(2) = "He3 + He3 => He4 + 2P"
  ref(2) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(3) = "He3 + He4 => Be7 + G ecc."
  ref(3) = "Cyburt, Davids Phys.Rev.C 78, 064614, 2008"

  lab(4) = "N15 + P => O16 + G"
  ref(4) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(5) = "Branching Be7"
  ref(5) = "NACRE 1999 + Adelberger Rev.Mod.Phys. 70, 4, 1265, 1998"

  lab(6) = "C12 + P => N13 + G"
  ref(6) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(7) = "N14 + P => O15 + G"
  ref(7) = "Imbriani et al. Eur.Phys. J.A 25, 455, 2005"

  lab(8) = "N15 + P => A + C12"
  ref(8) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(9) = "O16 + P => F17 + G"
  ref(9) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(10) = "He4 + 2A => C12 + G"
  ref(10) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(11) = "C12 + A => O16 + G"
  ref(11) = "Hammer et al. Nucl.Phys.A 758, 363, 2005"

  lab(12) = "N14 + A => F18 + G"
  ref(12) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(13) = "O16 + A => Ne20 + G"
  ref(13) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(14) = "O18 + A => Ne22 + G"
  ref(14) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(15) = "Ne20 + A => Mg24 + G"
  ref(15) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(16) = "Ne22 + A => Mg25 + N"
  ref(16) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(17) = "Ne22 + A => Mg26 + G"
  ref(17) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(18) = "Mg24 + A => Si28 + G"
  ref(18) = "Caughlan, Fowler. Atomic Data & Nuc.Data Tab., 40, 1988"

  lab(19) = "C13 + P => N14 + G"
  ref(19) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(20) = "C13 + A => O16 + N"
  ref(20) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(21) = "N15 + A => F19 + G"
  ref(21) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(22) = "O17 + P => N14 + A" 
  ref(22) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(23) = "O17 + A => Ne20 + N"
  ref(23) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(24) = "O18 + P => N15 + A" 
  ref(24) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(25) = "F19 + P => O16 + A"
  ref(25) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(26) = "F19 + A => Ne22 + P"
  ref(26) = "Caughlan, Fowler. Atomic Data & Nuc.Data Tab., 40, 1988"

  lab(27) = "Ne20 + P => Na21 + G" 
  ref(27) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(28) = "Ne21 + P => Na22 + G"
  ref(28) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(29) = "Ne22 + P => Na23 + G"
  ref(29) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(30) = "Na23 + P => Ne20 + A"
  ref(30) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(31) = "Na23 + P => Mg24 + G"
  ref(31) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(32) = "C12 + C12 => P + NA23"
  ref(32) = "JINA: Caughlan, Fowler 1988"

  lab(33) = "C12 + C12 => A + NE20"
  ref(33) = "JINA: Caughlan, Fowler 1988"

  lab(34) = "D + P => He3 + G"
  ref(34) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(35) = "Li6 + P => He3 + A"
  ref(35) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(36) = "Li7 + P => He4 + A"
  ref(36) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(37) = "Be9 + P => Li6 + A"
  ref(37) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(38) = "Be9 + P => Be8 + D => 2A + D"
  ref(38) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  lab(39) = "B11 + P => A + 2A"
  ref(39) = "NACRE: Angulo et al. Nucl.Phys.A 656, 3, 1999"

  write(67,*) "@ ================================"
  write(67,*) "@            Cross                "
  write(67,*) "@ ================================"
  do i=1,sezi
     write(67,'(" @ ",a55)') lab(i) 
     write(67,'(" @ ",4x,a55)') ref(i) 
  end do
  write(67,*) "@ ================================"
  write(67,*)

end subroutine ref_cross
