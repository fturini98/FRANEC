!     questa e' la neutr di Raffelt et al.
! *****************************************************************
! ATTENZIONE: cambiare le referenze in write se cambia l'algoritmo!!!!
! *****************************************************************
subroutine NEUTR(RO,T,QT) 
  use interfaccia
  use chim
  
  implicit none

  real :: RO, T, QT
 
  ! PRODUZIONE DI NEUTRINI
  real, save, dimension(15) :: TT
  real, save, dimension(51) :: RR
  real, save, dimension(15,51) :: Q
  real, save :: QQQ1, QTRASV, EFREQ, ROME 

  real,parameter,dimension(3) :: A0 = (/6.002d19,4.886d10,2.320d-7/)
  real,parameter,dimension(3) :: A1 = (/2.084d20,7.580d10,8.449d-8/)
  real,parameter,dimension(3) :: A2 = (/1.872d21,6.023d10,1.787d-8/)
  real,parameter,dimension(3) :: B1 = (/9.383d-1,6.290d-3,2.581d-2/)
  real,parameter,dimension(3) :: B2 = (/-4.141d-1,7.483d-3,1.734d-2/)
  real,parameter,dimension(3) :: B3 = (/5.829d-2,3.061d-4,6.990d-4/)
  real,parameter,dimension(3) :: C1 = (/5.5924,1.5654,0.56457/)
  real,dimension(3) :: F1 

  ! NN NUMERO DI NEUTRINI OLTRE QUELLI ELETTRONICI
  integer :: NN = 2, IP = 0

  integer :: i, j, k, kt, kr
  real :: cvpa, cvma, vme, vla, vla2, vla3, vla4, rome13, csi, chk, chk2
  real :: g1, q1, q2, q3, q4, q5, qpa, xt, xrm, dt, dr, qq1, qq2, freq, freq2
  real :: xgam, xgam2, ft, fl, tmp1, tmp2, xics, xy, dxx, var, fxy
  real :: q1t, qappr, am, zm, s1, s2, t8, efer, bb, ba2, vlb, d2, eXGAM
  real :: vld2, bb1, bb2, bb3, fb, gb, hh1, hh2, hh3, hh4, hh5, sXGAM, sVLA

  ! PRODUZIONE NEUTRINI DA:    PLASMA (Q1) - PHOTO (Q2) - COPPIE (Q3)     
  !  Q3 DA MUNAKATA,H., KOHYAMA,Y., ITOH,N. AP.J. 296,197 - 1985          
  !                + ERRATA DEL 1986                                      
  !  Q2 DA ITOH  ET AL. AP.J 339, 354 1989                                
  !  Q1 da haft et al. ApJ 425,222 1994                                   
  
  if(IP /= 1) then 
     ! LEGGE LA TABELLA TABEL.in DEI FOTONEUTRINI
     read(25,102) (TT(I),I=1,15) 
     read(25,103) (RR(I),I=1,51) 
     do J=1,15 
        read(25,104) (Q(J,K),K=1,51) 
     end do
     IP = 1 

     ! *****************************************************************
     ! ATTENZIONE: cambiare qui le referenze se cambia l'algoritmo!!!!
     ! *****************************************************************
     write(67,*) "@ ================================"
     write(67,*) "@            Neutrini             "
     write(67,*) "@ ================================"
     write(67,*) "@ PLASMA (Q1): Haft et al. ApJ 425, 222, 1994"
     write(67,*) "@ PHOTO (Q2): Itoh et al. AP.J 339, 354, 1989"
     write(67,*) "@ COPPIE (Q3): Munakata et al. AP.J. 296, 197, 1985"
     write(67,*) "@              + Errata 1986"
     write(67,*) "@ BREMSSTRAHLUNG: Dicus et al. AP.J. 210, 481, 1976"
     write(67,*) "@ ================================"
     write(67,*)
     
  endif
  QT = 0. 
  QQQ1 = 0. 

  if(T < 1.d8 .and. XX(1) > 1.d-5) return 

  CVPA = 1.122356+NN*0.254356 
  CVMA = 0.622356-NN*0.245644 
  VME = 2./(1.+XX(1)) 
 
  VLA = T/5.93097d9 
  VLA2 = VLA**2
  VLA3 = VLA**3 
  VLA4 = VLA**4 
  ROME = RO/VME 
  ROME13 = ROME**(1./3.)
  CSI = (ROME13/1.d3)/VLA 
  do K = 1,3 
     F1(K) = 0. 
     CHK = -C1(K)*CSI 
     if(abs(CHK) <= 70.) then
        F1(K) = (A0(K)+A1(K)*CSI+A2(K)*CSI**2)*exp(CHK) 
        F1(K) = F1(K)/(CSI**3+B1(K)/VLA+B2(K)/VLA2+B3(K)/VLA3) 
     endif
  end do
  G1 = 1.-13.04*VLA2+133.5*VLA4+1534.*VLA4*VLA2+918.6*VLA4*VLA4 

  ! COPPIE  
  Q3 = 0. 
  CHK2 = abs(2./VLA) 
  if(CHK2 <= 70.) then
     sVLA = sqrt(VLA)
     QPA = 1./(10.748*VLA2+.3967*sVLA+1.005)*(1.+ROME/(7.692E7*VLA3+  &
          9.715E6*sVLA))**(-.3)                                          
     Q3 = .5*CVPA*(1+CVMA/CVPA*QPA)*G1*exp(-2./VLA)*F1(1) 
  endif

  ! FOTONEUTRINI  
  ! INTERPOLA TABELLE ITOH                      
  Q2 = 0. 
  XT = log10(T) 
  if(XT >= 7.) then 
     XRM = log10(ROME) 
     if(XRM >= 0.) then 
        ! ricerca indice comprensiva di gestione di bordo superiore
        ! i dati in TT iniziano da 7 e sono spaziati di 0.1
        KT = min(int((XT-7.)/0.1) + 2, 15)

        ! ricerca indice comprensiva di gestione di bordo superiore
        KR = min(int(XRM/0.2) + 2, 51)

        DT = (XT-TT(KT-1))/(TT(KT)-TT(KT-1)) 
        QQ1 = Q(KT-1,KR-1)+(Q(KT,KR-1)-Q(KT-1,KR-1))*DT 
        QQ2 = Q(KT-1,KR)+(Q(KT,KR)-Q(KT-1,KR))*DT 
        DR = (XRM-RR(KR-1))/(RR(KR)-RR(KR-1)) 
        Q2 = QQ1+(QQ2-QQ1)*DR 
     endif
  endif
  
  ! PLASMA 
  ! CALCOLO LA FREQUENZA DI PLASMA --> E= FREQ/2                          
  FREQ2 = 8.223d-4*ROME/sqrt(1.+ 1.054d-4*ROME13**2) 
  FREQ = sqrt(FREQ2)
  EFREQ = FREQ/2. 

  ! CALCOLO LE PERDITE ENERGETICHE STANDARD PER PLASMA NEUTRINI          
  ! FORMULE DI HAFT et al. 1994                                       
  Q1 = 0. 
  XGAM2 = (1.1095d11*ROME)/( T**2*sqrt(1.+0.0001012627*ROME13**2))
  XGAM = sqrt(XGAM2)
  sXGAM = sqrt(XGAM)
  FT = 2.4+0.6*sXGAM+0.51*XGAM+1.25*sXGAM*XGAM
  FL = (8.6*XGAM2+1.35*sXGAM*XGAM**3)/(225.-17.*XGAM+XGAM2) 
  
  tmp1 = log10(2.*ROME)

  tmp2 = 3.*XT
  XICS = (1./6.)*(17.5 + tmp1 - tmp2) 
  XY = (1./6.)*(-24.5 + tmp1 + tmp2) 
  DXX = abs(XICS) 
  VAR = XY-1.6+1.25*XICS 
  if(DXX > 0.7 .or. XY < 0.) then 
     FXY = 1. 
  else 
     FXY = 1.05+(0.39-1.25*XICS-0.35*sin(4.5*XICS)-  &
          0.3*exp(-(4.5*XICS+0.9)**2))*             &
          exp(-(min(0.,VAR)/(0.57-0.25*XICS))**2) 
  endif

  ! calcolo anche solo l'emissione trasversa per poi poter contare i      
  ! neutrini      
  eXGAM = exp(-XGAM)
  Q1T = 0.9325*3.00d21*VLA**9*XGAM**6*eXGAM*FT*FXY 
  QAPPR = 3.00d21*VLA**9*XGAM**6*eXGAM*(FT+FL)*FXY 

  Q1 = 0.9325*QAPPR 

  ! RICOMBINAZIONE 
  ZM = XX(1)+2.*XX(2)+2.*XX(3)+ 6.*XX(4)+ 7.*XX(5)+ 8.*XX(6) 
  AM = XX(1)+3.*XX(2)+4.*XX(3)+12.*XX(4)+14.*XX(5)+16.*XX(6) 
  S1 = 1.85d-4*(ZM**6)*RO**2*VLA2/AM**2 
  S2 = -1.57d5*ZM*ZM/T-2.428d-5*(RO**(2./3.))/VLA 
  Q4 = S1*exp(S2) 

  ! NEUTRINI DI BREMSSTRAHLUNG - DICUS ET AL. (1976) AP.J. 210,481 
  T8 = T/1.d8 
  EFER = sqrt(1.018d-4 * ROME13**2 + 1. )
  BB = 1./sqrt(1.-1./(EFER*EFER)) 
  BA2 = BB*BB 
  VLB = log(abs((BB+1.)/(BB-1.))) 
  D2 = BB/215. 
  VLD2 = log((2.+D2)/D2) 
  BB1 = -2.0/3.0+BA2+.50*BB*(1.-BA2)*VLB 
  BB2 = -4.+2.*(1.+D2)*VLD2 
  BB3 = VLD2-2./(2.+D2) 
  FB = 0.14*(BB1*BB2-BB3*(BA2-1.)*(2./3.+.5*BA2-BB/4.*(BA2+1.)*VLB)) 
  GB = 0.14*BB3*(BA2-1.)*(2./3.-5./2.*BA2-BB/4.*(3.-5.*BA2)*VLB) 
  Q5 = 2.*ZM*ZM/AM*T8**6*(.5*CVPA*FB-.5*CVMA*GB)*RO 
  QT = -(Q1+Q2+Q3+Q4+Q5)/RO  ! 16*   mult
  QQQ1 = -Q1/RO 
  QTRASV = -Q1T/RO 
  HH1 = Q1 
  HH2 = Q2 
  HH3 = Q3 
  HH4 = Q4 
  HH5 = Q5 

102  format(15F4.1) 
103  format(15F5.1) 
104  format (8E11.4) 

  return 
end subroutine NEUTR

