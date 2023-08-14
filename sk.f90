! Questa funzione calcola gli schermi per le reazioni nucleari

! X e' l'abbondanza in numero passata da epsi

function SK(RO,T,X,Z1,Z2, cache) 
  use parametri

  ! cache = 1 serve a eliminare il calcolo di quantita' ripetute
  ! che la funzione ha appena calcolato al giro precedente
  implicit none 

  real :: sk, RO,T,X(MELE)
  integer :: Z1, Z2, cache 

  real :: SUM,F(MELE),PMU,FER,FL12,PASS1,PASS2,PASS3,SCR1,SCR2,EF
  real :: A1,A2,G12,CHK 
  real :: SKtmp, ZMt,TV1,TV2, THE, tmp

  ! questi blocchi data precalcolano delle quantita' funzioni di Z1 e Z2
  ! velocizzando il calcolo della routine
  ! CZ = Z**(1./3.) con Z = (1,..., 10)
  !real,parameter,dimension(10) :: CZ = (/1.000000,1.259921,1.442250, &
  !     1.587401,1.709976,1.817121,1.912931,2.000000,2.080084,2.154435/)
  real,parameter,dimension(20) :: CZ = (/1.000000,1.259921,1.442250, &
       1.587401,1.709976,1.817121,1.912931,2.000000,2.080084,2.154435, &
       2.223980091,2.289428485,2.351334688,2.410142264,2.466212074, &
       2.5198421,2.571281591,2.620741394,2.668401649,2.714417617/)
  ! PZ1 = Z**1.86  con Z = (1,..., 20)
  real,parameter,dimension(20) :: PZ1 = (/1.000000,3.630077,7.716947, &
       13.177456,19.956492,28.013110,37.314906,47.835176,59.551275, &
       72.443596,86.494879,101.689735,118.014300,135.455969,154.003199, &
       173.645354,194.372578,216.175693,239.046119,262.975804/)
  ! PZ2 = Z**(5./3.)   con Z = (1,..., 20)
  real,parameter,dimension(20) :: PZ2 = (/1.000000,3.174802,6.240251, &
       10.079368,14.620089,19.811563,25.615140,32.000000,38.940738, &
       46.415888,54.406962,62.897793,71.874073,81.323000,91.233030, &
       101.593667,112.395313,123.629138,135.286980,147.361260/)    
  ! PZ3 = Z**(4./3.)   con Z = (1,..., 20)
  real,parameter,dimension(20) :: PZ3 = (/1.000000,2.519842,4.326749, &
       6.349604,8.549880,10.902724,13.390518,16.000000,18.720754, &
       21.544347,24.463781,27.473142,30.567351,33.741992,36.993181, &
       40.317474,43.711787,47.173345,50.699631,54.288352/) 
  ! PZ4 = Z**(2./3.)   con Z = (1,..., 20)
  real,parameter,dimension(20) :: PZ4 = (/1.000000,1.587401,2.080084, &
       2.519842,2.924018,3.301927,3.659306,4.000000,4.326749,4.641589, &
       4.946087,5.241483,5.528775,5.808786,6.082202,6.349604,6.611489, &
       6.868285,7.120367,7.368063/)
  ! PZ5 = Z**(1./3.)   con Z = (1,..., 20)
  real,parameter,dimension(20) :: PZ5 = (/1.000000,1.259921,1.442250, &
       1.587401,1.709976,1.817121,1.912931,2.000000,2.080084,2.154435, &
       2.223980,2.289428,2.351335,2.410142,2.466212,2.519842,2.571282, &
       2.620741,2.668402,2.714418/)
  ! Z**1.58
  real,parameter,dimension(MELE) :: Z158 = (/2.989698,1.000000,2.989698, &
       16.962075,21.639883,26.722813,26.722813,38.018940,38.018940,50.711490, &
       50.711490,50.711490,64.696725,0.000000,44.197779,16.962075, &
       21.639883,26.722813,32.188680,38.018940,172.051273,5.673507,5.673507, &
       8.938297,12.716647,1.000000/)
  integer :: K
  real, save :: Zms,ZTs,Atmp,FL,ZT,ZM,ZAV,T12,Told, Flt, Fls, T13
  real :: MU12

  real,parameter,dimension(MELE) :: Z = (/2.,1.,2.,6.,7.,8.,8.,10.,10., &
       12.,12.,12.,14.,0.,11.,6.,7.,8.,9.,10.,26.,3.,3.,4.,5.,1./) 

  real,parameter :: OVER = 200., PI = 3.141592654, H = 1.660531d-24, &
       E = 4.80325d-10 
  real,parameter :: HTAG = 1.055d-27, KB = 1.38062d-16, INFSTR = 0.2
  real,parameter :: PI4H = 7.567682d24, E2 = 2.307121d-19, V = 66.61983 
  real,parameter :: E4 = 5.322808d-38, V2 = 1.536665d-70 
  real,parameter :: coeff = 5351.69

  integer,save :: firsttime = 1

  if(firsttime == 1) then
     ! *****************************************************************
     ! ATTENZIONE: cambiare qui le referenze se cambia l'algoritmo!!!!
     ! *****************************************************************
     write(67,*) "@ ================================"
     write(67,*) "@            Screening            "
     write(67,*) "@ ================================"
     write(67,*) "@ WEAK-INTERMEDIATE: Graboske et al., AP.J., 181, 457, 1973"
     write(67,*) "@ INTERMEDIATE-STRONG: Dewitt et al., AP.J., 181, 439, 1973"
     write(67,*) "@ STRONG: Itoh et al., AP.J., 218, 477, 1977"
     write(67,*) "@ STRONG: Itoh et al., AP.J., 234, 1079, 1979"
     write(67,*)
     firsttime = 0
  endif

  !     WEAK - INTERMEDIATE AND INTERMEDIATE-STRONG SCREENING:            
  !     GRABOSKE, DE WITT, GROSSMAN AND COOPER, AP. J., 181, 457-474, 1973
  !     DEWITT, GRABOSKE AND COOPER, AP. J., 181, 439-456, 1973           
  !     STRONG SCREENING:                                                 
  !     ITOH,TOTSUJI AND ICHIMARU, AP. J., 218, 477-483, 1977             
  !     ITOH, TOTSUJI, ICHIMARU AND DE WITT, AP. J., 234, 1079-1084 , 1979
  
  ! eseguo il blocco seguente solo se non ho una cache
  if(cache /= 1) then
     sum = dot_product(Z, X)
     Atmp = 1.442250/(PI4H*RO*SUM)**(1./3.) 
     T13 = T**(1./3.)
  endif

  A1   = CZ(Z1)*Atmp 
  A2   = CZ(Z2)*Atmp 
  G12  = (Z1*Z2*E2)/(0.5*(A1+A2)*KB*T) 
  T12 = Z1*Z2*coeff/(PZ5(Z1+Z2)*T13)
  CHK  = 3.*G12/T12 

  ! STRONG (ITHO)                         
  if(CHK >= INFSTR) then 
     EF = 1.25*G12-0.095*T12*CHK**2
     SK = exp(EF) 
     return 
  endif

  if(cache /= 1) then
     ! CALCOLO PESO MOLECOLARE
     SUM = 0.
     do K=1,MELE 
        ! HE3 E H SPOSTATI !!!                  
        F(K) = X(K) 
        SUM = SUM+X(K) 
     end do
     PMU = 1./SUM 
     ! CALCOLO DI TETA 
     FER = 5.4885d7*RO*(1.0+X(2))/(T*sqrt(T)) 
     if(FER <= 3.) then
        THE = 1./(1.+0.39716*FER-0.00929*FER**2) 
     else
        THE = 1.1447/((FER**(2./3.))*(1.+0.28077/FER))
     endif
     ! Teta  
     ! CALCOLO ZETA MEDIO E ZETA TILDE
     do K=1,MELE 
        ! ni/nI                                 
        F(K) = F(K)*PMU 
     end do
     ZM = 0.
     ZT = 0. 
     do K=1,MELE 
        ! z Medio   
        tmp = Z(K)*F(K)
        ZM = ZM+tmp 
        ZT = ZT+tmp*(Z(K)+THE) 
     end do
     ZMs = ZM**(0.28) 
     ! z Tilde                               
     ZT = sqrt(ZT) 
     ZTs = ZT**(0.58)
     ! CALCOLO LAMBDA ZERO
     FL = 1.88d8*sqrt(RO/(PMU*T**3)) 
     FLs = FL**(0.86)
     FLt = FL**(2./3.) 
     ! CALCOLO  < z**(3b-1) > 
     ZAV = 0. 
     do K=1,MELE 
        ! < z**(3b-1) >                    
        ZAV = ZAV+F(K)*Z158(K) 
     end do
  endif

  ! Lambda 12                             
  FL12 = FL*Z1*Z2*ZT 

  !  CALCOLO SCREENING 
  if(FL12 <= 0.1) then 
     ! WEAK                                  
     SK = exp(FL12) 
  else if(FL12 <= 5.) then 
     PASS1 = PZ1(Z1+Z2)-PZ1(Z1)-PZ1(Z2) 
     TV1 = ZTs*ZMs 
     TV2 = FLs
     SCR1 = 0.38*ZAV/TV1 * PASS1 * TV2 
     SCR1 = exp(SCR1) 
     if(FL12 <= 2.) then 
        ! INTERMEDATE                           
        SK = SCR1 
     else 
        PASS1 = PZ2(Z1+Z2)-PZ2(Z1)-PZ2(Z2) 
        PASS2 = PZ3(Z1+Z2)-PZ3(Z1)-PZ3(Z2) 
        PASS3 = ( PZ4(Z1+Z2)-PZ4(Z1)-PZ4(Z2) )/FLt 
        ZMt = ZM**(1./3.) 
        SCR2 = 0.624*ZMt * (PASS1 + 0.316 * ZMt * PASS2 + &
             0.737/ZM * PASS3) * FLt                           
        if(SCR2 > OVER) SCR2 = OVER 
        SCR2 = exp(SCR2) 
        ! INTERMEDIATE-STRONG              
        SK = min(SCR1,SCR2) 
     endif
  else 
     PASS1 = PZ2(Z1+Z2)-PZ2(Z1)-PZ2(Z2) 
     PASS2 = PZ3(Z1+Z2)-PZ3(Z1)-PZ3(Z2) 
     PASS3 = ( PZ4(Z1+Z2)-PZ4(Z1)-PZ4(Z2) )/FLt 
     ZMt = ZM**(1./3.) 
     SCR2 = 0.624*ZMt * (PASS1 + 0.316 * ZMt * PASS2 + &
          0.737/ZM*PASS3) * FLt                                
     if(SCR2 > OVER) SCR2 = OVER 
     ! STRONG (GRABOSKE)                     
     SK = exp(SCR2) 
  endif

  return 
end function SK
