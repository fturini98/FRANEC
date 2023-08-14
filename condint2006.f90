! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! Subroutine di Potekhin: legge la tabella con il logaritmo Coulombiano 
! efficace per gli elettroni (condall06.d) per vari Z, RHO e T.         
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

subroutine potekhin(Zion,RLG,TLG,CK, index) 
  use interfaccia
  
  implicit none

  real :: Zion, RLG, TLG, CK
  integer :: index

  real :: DRK, DTK
  real,parameter :: TLGmin=3., TLGmax=9.
  real,parameter :: RLGmin=-6., RLGmax=9.75

  if (TLG < TLGmin .or. TLG > TLGmax) then 
     print*,'choose Lg T between',TLGmin,'  and',TLGmax 
     write(66,*)'301 - potekhin'
     write(66,*)'choose Lg T between',TLGmin,'  and',TLGmax
     stop 
  endif

  if (RLG < RLGmin .or. RLG > RLGmax) then 
     print*,'choose Lg R between',RLGmin,'  and',RLGmax 
     CK = 99.000 
  else 
     call CONINTER(Zion,TLG,RLG,CK, index) 
  endif

  return 
END subroutine potekhin


!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!  This subroutine interpolates the electron thermal conductivity       
!               from the data file "condall06.d"                        
!  ---------------------------------------------------  Version 23.05.99
! Modificata in modo da accettare il logaritmo di Zion e un indice che 
! evita il calcolo delle hunt per gli elementi 2...21
subroutine CONINTER(Zion,TLG,RLG,CK, index) 
  ! Input: Zion - ion charge, TLG - lg(T[K]), RLG - lg(rho[g/cc])         
  ! Output: CK - Log_{10} thermal conductivity (kappa) [CGS units]        
  !**      (it is also possible to obtain all second derivatives)      ***
  implicit none

  real :: Zion, TLG, RLG, CK
  integer :: index

  integer,parameter :: MAXT=19, MAXR=64, MAXZ=15
!!! NB: These parameters must be consistent with the table "condall06.d"
  real,save,dimension(MAXT) :: AT
  real,save,dimension(MAXR) :: AR
  real,save,dimension(MAXZ) :: AZ
  real,save,dimension(MAXT,MAXR,MAXZ) :: AKAP 
  integer,save ::  KRUN = -1

  real,save :: CKTMZ0, CKT0Z0,  CKT1Z0, CKTPZ0, RLGs
  integer,save :: ITs, IRs, IZs
  real,save :: CKTMZ1, CKT0Z1, CKT1Z1, CKTPZ1 

  integer,save :: iz, it, ir, firsttime = 1
  integer :: itm, itp, irm, irp, isalta, icache 
  real :: z, zlg, xr, xz0, xz1, cktm, drktm, dr2ktm, ckt0, ckt1, drkt0
  real :: dr2kt0, drkt1, dr2kt1, cktp, drktp, dt2k, xt, drtk, drt2k, dr2ktp

  real,dimension(4) :: a
  real,dimension(8) :: r
  real,dimension(32) :: b

  real,save,dimension(21) :: izsv

  ! Reading                                 
  if (KRUN /= 12345) then 
     open(1,file='condall06.d',status='OLD') 
     print*,'Reading thermal conductivity data...' 
     ! skip the first line                            
     read(1,'(A)') 
     do IZ=1,MAXZ 
        read(1,*) Z,(AT(IT),IT=1,MAXT) 
        AZ(IZ) = log10(Z) 
        do IR=1,MAXR 
           read(1,*) AR(IR),(AKAP(IT,IR,IZ),IT=1,MAXT) 
        enddo
     enddo
     close(1) 
     KRUN = 12345 
     IZ = MAXZ/2+1 
     IT = MAXT/2+1 
     IR = MAXR/2+1 
     print*,'done.' 
  endif
  ! matt non e' piu' necessario perche' la kappa passa il logaritmo di Zion
  !ZLG = log10(Zion) 
  ZLG = Zion

  ! matt: gestisco la ricerca solo la prima volta
  ! altrimenti prendo la cache
  if(firsttime == 1) then
     call HUNT(AZ,MAXZ,ZLG,IZ) 
     if (IZ == 0 .or. IZ == MAXZ) then
        write(66,*)'310 - cointer'
        write(66,*)'Z out of range'
        stop'CONINTER: Z out of range' 
     endif
     izsv(index) = iz
     if(index == 21) firsttime = 0
  else
     iz = izsv(index)
  endif
  
  if(index == 1) then
     call HUNT(AT,MAXT,TLG,IT) 
     if (IT == 0 .or. IT == MAXT) then
        write(66,*)'311 - cointer'
        write(66,*)'T out of range'
        stop'CONINTER: T out of range' 
     endif
     call HUNT(AR,MAXR,RLG,IR) 
     if (IR == 0 .or. IR == MAXR) then
        write(66,*)'312 - cointer'
        write(66,*)'rho out of range'
        stop'CONINTER: rho out of range' 
     endif
  endif

  ITM = max(1,IT-1) 
  ITP = min(MAXT,IT+2) 
  IRM = max(1,IR-1) 
  IRP = min(MAXR,IR+2) 

  ! giada
  isalta = 0
  if(IT  ==  ITs .and. IR  ==  IRs .and. IZ  ==  IZs) then 
     if( RLG == RLGs ) then
        isalta = 1
     endif
  endif

  if( isalta == 0 ) then
     ITs = IT 
     IRs = IR 
     IZs = IZ 
     RLGs = RLG 

     a(1) = AR(IRM)
     a(2) = AR(IR)
     a(3) = AR(IR+1)
     a(4) = AR(IRP)
     b(1) = AKAP(ITM,IRM,IZ)
     b(2) = AKAP(ITM,IR,IZ)
     b(3) = AKAP(ITM,IR+1,IZ)
     b(4) = AKAP(ITM,IRP,IZ)
     b(5) = AKAP(IT,IRM,IZ)
     b(6) = AKAP(IT,IR,IZ)
     b(7) = AKAP(IT,IR+1,IZ)
     b(8) = AKAP(IT,IRP,IZ)
     b(9) = AKAP(IT+1,IRM,IZ)
     b(10) = AKAP(IT+1,IR,IZ)
     b(11) = AKAP(IT+1,IR+1,IZ)
     b(12) = AKAP(IT+1,IRP,IZ)
     b(13) = AKAP(ITP,IRM,IZ)
     b(14) = AKAP(ITP,IR,IZ)
     b(15) = AKAP(ITP,IR+1,IZ)
     b(16) = AKAP(ITP,IRP,IZ)
     b(17) = AKAP(ITM,IRM,IZ+1)
     b(18) = AKAP(ITM,IR,IZ+1)
     b(19) = AKAP(ITM,IR+1,IZ+1)
     b(20) = AKAP(ITM,IRP,IZ+1)
     b(21) = AKAP(IT,IRM,IZ+1)
     b(22) = AKAP(IT,IR,IZ+1)
     b(23) = AKAP(IT,IR+1,IZ+1)
     b(24) = AKAP(IT,IRP,IZ+1)
     b(25) = AKAP(IT+1,IRM,IZ+1)
     b(26) = AKAP(IT+1,IR,IZ+1)
     b(27) = AKAP(IT+1,IR+1,IZ+1)
     b(28) = AKAP(IT+1,IRP,IZ+1)
     b(29) = AKAP(ITP,IRM,IZ+1)
     b(30) = AKAP(ITP,IR,IZ+1)
     b(31) = AKAP(ITP,IR+1,IZ+1)
     b(32) = AKAP(ITP,IRP,IZ+1)

     ! Cubic interpolation in RLG: 
     call CINTERP3b(a,RLG,b, IR,MAXR, r)

     CKTMZ0 = r(1)
     CKT0Z0 = r(2)
     CKT1Z0 = r(3)
     CKTPZ0 = r(4)
     CKTMZ1 = r(5)
     CKT0Z1 = r(6)
     CKT1Z1 = r(7)
     CKTPZ1 = r(8)
  endif

  ! Linear interpolation in ZLG:                                          
  XZ1 = (ZLG-AZ(IZ))/(AZ(IZ+1)-AZ(IZ)) 
  XZ0 = 1.-XZ1 
  CKTM = XZ0*CKTMZ0+XZ1*CKTMZ1 
  CKT0 = XZ0*CKT0Z0+XZ1*CKT0Z1 
  CKT1 = XZ0*CKT1Z0+XZ1*CKT1Z1 
  CKTP = XZ0*CKTPZ0+XZ1*CKTPZ1 

  ! Cubic interpolation in TLG:                                           
  call CINTERP3(AT(ITM),AT(IT),AT(IT+1),AT(ITP),TLG,IT,MAXT,        &
       CKTM,CKT0,CKT1,CKTP,CK)
  ! input: values of lg kappa               
  ! lg kappa, d lg k / d lg T, d2 lg k / d2 lg T 
  return 
END subroutine CONINTER


subroutine CINTERP3b(a,Z,b, N0,MXNV,r)
  ! Questa routine svolge lo stesso lavoro di cinterp3, ma in modo
  ! vettoriale. 
  implicit none

  real :: Z
  real,dimension(4) :: a
  real,dimension(8) :: r
  real,dimension(32) :: b
  integer :: N0, MXNV

  real :: X,H,H2,XH2,HM,HP,HM2,HP2,HTMP1,HTMP2, XH
  real :: vv, v01, v11, c2, c3
  integer :: i, idx

  if (N0 <= 0 .or. N0 >= MXNV) then
     write(66,*)'320 - CINTERP3'
     write(66,*)'N0 out of range'
     stop'CINTERP: N0 out of range' 
  endif

  X = Z-a(2)
  ! basic interval                        
  H = a(3)-a(2) 
  H2 = H**2
  XH = X/H 
  XH2 = XH*XH 
  ! left adjoint interval                 
  HM = a(2)-a(1)
  ! right adjoint interval                
  HP = a(4)-a(3) 
  HM2 = HM*HM 
  HP2 = HP*HP 
  HTMP1 = (1./H+1./HM) 
  HTMP2 = (1./H+1./HP) 

  if (N0 > 1 .and. N0 < MXNV-1) then 
     do i = 1,8
        idx = (i-1)*4
        VV = b(3+idx)-b(2+idx) 
        ! left derivative                
        V01 = (VV/H2+(b(2+idx)-b(1+idx))/HM2)/HTMP1 
        ! right derivative               
        V11 = (VV/H2+(b(4+idx)-b(3+idx))/HP2)/HTMP2 
        ! Cubic interpolation        
        C2 = 3.*VV - H*(V11 + 2.*V01) 
        C3 = H*(V01+V11) - 2.*VV 
        r(i) = b(2+idx) + V01*X + C2*XH2 + C3*XH2*XH 
     end do
  else
     ! Quadratic interpolation
     if (N0 == 1) then 
        do i = 1,8
           idx = (i-1)*4
           V11 = (VV/H2+(b(4+idx)-b(3+idx))/HP2)/HTMP2
           C2 = -VV + V11*H 
           r(i) = b(3+idx) - V11*(H-X) + C2*(1.-XH)**2 
        end do
     else 
        do i = 1,8
           idx = (i-1)*4
           V01 = (VV/H2+(b(2+idx)-b(1+idx))/HM2)/HTMP1 
           C2 = VV - V01*H 
           r(i) = b(2+idx) + V01*X+C2*XH2 
        end do
     endif
  end if

  return 
END subroutine CINTERP3b


subroutine CINTERP3(ZM,Z0,Z1,ZP,Z,N0,MXNV,VM,V0,V1,VP,VF)
  ! Given 4 values of Z and 4 values of V, find VF corresponding to 5th Z 
  !                                                       Version 23.05.99
  !   Output: VF - interpolated value of function                         
  implicit none

  real :: ZM,Z0,Z1,ZP,Z,VM,V0,V1,VP,VF
  integer :: N0, MXNV

  real :: X,H,H2,HM,HP,HM2,HP2,HTMP1,HTMP2, XH
  real :: vv, v01, v11, c2, c3

  X = Z-Z0 
  ! basic interval                        
  H = Z1-Z0 
  H2 = H**2
  XH = X/H 

  ! left adjoint interval                 
  HM = Z0-ZM 
  ! right adjoint interval                
  HP = ZP-Z1 
  HTMP1 = (1./H+1./HM) 
  HTMP2 = (1./H+1./HP) 
  HM2 = HM*HM 
  HP2 = HP*HP 
 
  VV = V1-V0 

  if (N0 > 1 .and. N0 < MXNV-1) then 
     ! left derivative                
     V01 = (VV/H2+(V0-VM)/HM2)/HTMP1 
     ! right derivative               
     V11 = (VV/H2+(VP-V1)/HP2)/HTMP2 
     ! Cubic interpolation        
     C2 = 3.*VV - H*(V11 + 2.*V01) 
     C3 = H*(V01+V11) - 2.*VV 
     VF = V0 + V01*X + C2*XH**2 + C3*XH**3 
     return
  else
     ! Quadratic interpolation                     
     if (N0 == 1) then 
        V11 = (VV/H2+(VP-V1)/HP2)/HTMP2 
        C2 = -VV + V11*H 
        VF = V1 - V11*(H-X) + C2*(1.-XH)**2 
     else 
        V01 = (VV/H2+(V0-VM)/HM2)/HTMP1
        C2 = VV - V01*H 
        VF = V0 + V01*X+C2*XH**2 
     endif
  endif
  
  return 
END subroutine CINTERP3

subroutine HUNT(XX,N,X,JLO) 
  !   W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling              
  !   Numerical Receipes(Cambridge Univ., 1986)                           
  !     Given an array XX of length N, and given a value X,               
  !     returns a value JLO such that X is between XX(JLO) and XX(JLO+1). 
  !     XX must be monotonic, either increasing or decreasing.            
  !     JLO=0 or JLO=N is returned to indicate that X is out of range.    
  !     JLO on input is taken as the initial guess for JLO on output.     
  implicit none
  
  real,dimension(*) :: XX
  integer :: N,JLO
  real :: X

  logical :: ASCND 
  integer :: jhi, inc, jm

  ! true if ascending order, false otherwise   
  ASCND = XX(N) > XX(1) 
  ! Input guess not useful.          
  if (JLO <= 0 .or. JLO > N) then 
     JLO = 0 
     JHI = N+1 
     ! go immediately to bisection                           
     goto 3 
  endif
  ! set the hunting increment                                 
  INC = 1 
  ! Hunt up:                       
  if (X >= XX(JLO) .eqv. ASCND) then 
1    JHI = JLO+INC 
     ! Done hunting, since off end of table       
     if (JHI > N) then 
        JHI = N+1 
        ! Not done hunting         
     else if (X >= XX(JHI) .eqv. ASCND) then 
        JLO = JHI 
        INC = INC+INC 
        goto 1 
     endif
     ! Hunt down:                                                 
  else 
     JHI = JLO 
2    JLO = JHI-INC 
     ! Done hunting, since off end of table       
     if (JLO < 1) then 
        JLO=0 
        ! Not done hunting         
     elseif (X < XX(JLO).eqv.ASCND) then 
        JHI = JLO 
        ! so double the increment                        
        INC = INC+INC 
        ! and try again                                  
        goto 2 
        ! Done hunting, value bracketed                           
     endif
  endif
  !   Hunt is done, so begin the final bisection phase:                   
3 if (JHI-JLO == 1) return 
  JM = (JHI+JLO)/2. 
  if (X >= XX(JM) .eqv. ASCND) then 
     JLO = JM 
  else 
     JHI = JM 
  endif
  goto 3 
END subroutine HUNT
