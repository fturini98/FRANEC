!*************************************************************          
! This routine was written by Anne A. Thoul, at the Institute           
! for Advanced Study, Princeton, NJ 08540.                              
! See Thoul et al., Ap.J. 421, p. 828 (1994)                            
! The subroutines LUBKSB and LUDCMP are from Numerical Recipes.         
!*************************************************************          
! This routine inverses the burgers equations.                          
!                                                                       
! The system contains N equations with N unknowns.                      
! The equations are: the M momentum equations,                          
!                    the M energy equations,                            
!                    two constraints: the current neutrality            
!                                     the zero fluid velocity.          
! The unknowns are: the M diffusion velocities,                         
!                   the M heat fluxes,                                  
!                   the electric field E                                
!                   the gravitational force g.                          
!                                                                       
!**************************************************                     
subroutine DIFFUSION(M,A,Z,X,CL,AP,AT,AX) 
  use interfaccia

  ! The parameter M is the number of species considered.                  
  !                                                                       
  ! Fluid 1 is the hydrogen                                               
  ! Fluid 2 is the helium                                                 
  ! Fluids 3 to M-1 are the heavy elements                                
  ! Fluid M is the electrons                                              
  !                                                                       
  ! The vectors A,Z and X contain the atomic mass numbers,                
  ! the charges (ionization), and the mass fractions, of the elements.    
  ! NOTE: Since M is the electron fluid, its mass and charge must be      
  !      A(M)=m_e/m_u                                                     
  !      Z(M)=-1.                                                         
  !                                                                       
  ! The array CL contains the values of the Coulomb Logarithms.           
  ! The vector AP, AT, and array AX contains the results for the diffusion
  ! coefficients.                                                         

  implicit none     
  integer,parameter:: MMAX=20, NMAX=42
  integer :: M
  real,dimension(M) :: A,Z,X,AP,AT
  real,dimension(M,M) :: AX, CL 

  integer :: I,J,L,N

  ! giada                                                                 
  real,save :: XX(MMAX,MMAX),Y(MMAX,MMAX),YY(MMAX,MMAX),VTMP1(MMAX,MMAX)
  real,save :: Ao(MMAX)
  integer :: icache                      

  real :: TEMP, KO, D, VARTMP,CC, AC  

  integer,dimension(2*M+2) :: INDX 
  real,dimension(2*M+2) :: gammatmp, ALPHA, NU, GA
  real,dimension(M) :: C
  real,dimension(M,M) :: K
  real,dimension(2*M+2, 2*M+2) :: GAMMA, DELTA

  logical,parameter :: use_lapack = .false.
  integer :: ierr

  ! The vector C contains the concentrations                              
  ! CC is the total concentration: CC=sum(C_s)                            
  ! AC is proportional to the mass density: AC=sum(A_s C_s)               
  ! The arrays XX,Y,YY and K are various parameters which appear in       
  ! Burgers equations.                                                    
  ! The vectors and arrays ALPHA, NU, GAMMA, DELTA, and GA represent      
  ! the "right- and left-hand-sides" of Burgers equations, and later      
  ! the diffusion coefficients.                                           

  ! Initialize parameters:                                                

  KO = 2. 
  N = 2*M+2 
  indx = 0
  
  C = 0. 
  CC = 0. 
  AC = 0. 

  ! Calculate concentrations from mass fractions:                         
  TEMP=0. 
  do I=1,M-1 
     TEMP = TEMP+Z(I)*X(I)/A(I) 
  end do
  do I=1,M-1 
     C(I) = X(I)/A(I)/TEMP 
  end do
  C(M) = 1. 

  ! Calculate CC and AC:                                                  

  do I=1,M 
     CC = CC+C(I) 
     AC = AC+A(I)*C(I) 
  enddo

  ! Calculate the mass fraction of electrons:                             
  X(M) = A(M)/AC 

  ! Calculate the coefficients of the burgers equations   

  ! giada
  icache = 1
  do i=1,M
     if(A(i) /= Ao(i)) then
        icache = 0
        do j=1,M
           Ao(i) = A(i)  ! metto in cache i valori
        end do
        exit
     endif
  end do
                
  do I=1,M 
     do J=1,M 
        ! giada 
        if(icache /= 1) then 
           XX(I,J) = A(J)/(A(I)+A(J)) 
           Y(I,J) = A(I)/(A(I)+A(J)) 
           YY(I,J) = 3.0*Y(I,J)+1.3*XX(I,J)*A(J)/A(I) 
           VTMP1(I,J) = sqrt(A(I)*A(J)/(A(I)+A(J)))*Z(I)**2*Z(J)**2 
        endif
        K(I,J) = 1.*CL(I,J)*VTMP1(I,J)*C(I)*C(J) 
     end do
  end do
  ! Write the burgers equations and the two constraints as                
  ! alpha_s dp + nu_s dT + sum_t(not 2 or M) gamma_st dC_t                
  !                     = sum_t delta_st w_t                              

  do I=1,M 
     ALPHA(I) = C(I)/CC 
     NU(I) = 0. 
     do J=1,M 
        GAMMA(I,J) = 0. 
     enddo
     do J=1,M 
        if (J /= 2 .and. J /= M) then 
           ! giada    
           if(i == 1) then 
              gammatmp(j) = -C(J)/CC+C(2)/CC*Z(J)*C(J)/Z(2)/C(2) 
              GAMMA(I,J) = gammatmp(j) 
           else 
              GAMMA(I,J) = gammatmp(J) 
           endif
           if (J == I) then 
              GAMMA(I,J) = GAMMA(I,J)+1. 
           endif
           if (I == 2) then 
              GAMMA(I,J) = GAMMA(I,J)-Z(J)*C(J)/Z(2)/C(2) 
           endif
           GAMMA(I,J) = GAMMA(I,J)*C(I)/CC 
        endif
     end do

     do J=M+1,N 
        GAMMA(I,J) = 0. 
     end do
  end do

  do I=M+1,N-2 
     ALPHA(I) = 0. 
     NU(I) = 2.5*C(I-M)/CC 
     do J=1,N 
        GAMMA(I,J) = 0. 
     end do
  end do

  ALPHA(N-1) = 0. 
  NU(N-1) = 0. 
  do J=1,N 
     GAMMA(N-1,J) = 0. 
  enddo

  ALPHA(N) = 0. 
  NU(N)=0. 
  do J=1,N 
     GAMMA(N,J) = 0. 
  end do

  do I=1,N 
     do J=1,N 
        DELTA(I,J) = 0. 
     end do
  end do

  do I=1,M 
     do J=1,M 
        if (J == I) then 
           do L=1,M 
              if(L /= I) then 
                 DELTA(I,J) = DELTA(I,J)-K(I,L) 
              endif
           end do
        else 
           DELTA(I,J) = K(I,J) 
        endif
     end do

     do J=M+1,N-2 
        if(J-M == I) then 
           do L=1,M 
              if (L /= I) then 
                 DELTA(I,J) = DELTA(I,J)+0.6*XX(I,L)*K(I,L) 
              endif
           end do
        else 
           DELTA(I,J) = -0.6*Y(I,J-M)*K(I,J-M) 
        endif
     end do

     DELTA(I,N-1) = C(I)*Z(I) 

     DELTA(I,N) = -C(I)*A(I) 
  end do

  do I=M+1,N-2 
     do J=1,M 
        if (J == I-M) then 
           do L=1,M 
              if (L/=I-M) then 
                 DELTA(I,J) = DELTA(I,J)+1.5*XX(I-M,L)*K(I-M,L) 
              endif
           end do
        else 
           DELTA(I,J) = -1.5*XX(I-M,J)*K(I-M,J) 
        endif
     end do

     do J=M+1,N-2 
        if (J-M == I-M) then 
           do L=1,M 
              if (L /= I-M) then 
                 DELTA(I,J) = DELTA(I,J)-Y(I-M,L)*K(I-M,L)*           &
                      (1.6*XX(I-M,L)+YY(I-M,L))                     
              endif
           end do
           DELTA(I,J) = DELTA(I,J)-0.8*K(I-M,I-M) 
        else 
           DELTA(I,J) = 2.7*K(I-M,J-M)*XX(I-M,J-M)*Y(I-M,J-M) 
        endif
     end do

     DELTA(I,N-1) = 0. 
     DELTA(I,N) = 0. 
  end do

  do J=1,M 
     DELTA(N-1,J) = C(J)*Z(J) 
  enddo
  do J=M+1,N 
     DELTA(N-1,J) = 0. 
  end do

  do J=1,M 
     DELTA(N,J) = C(J)*A(J) 
  end do
  do J=M+1,N 
     DELTA(N,J) = 0. 
  enddo

  ! Inverse the system for each possible right-hand-side, i.e.,           
  ! if alpha is the r.h.s., we obtain the coefficient A_p                 
  ! if nu    ---------------------------------------- A_T                 
  ! if gamma(i,j) ----------------------------------- A_Cj                
  !                                                                       
  ! If I=1, we obtain the hydrogen diffusion velocity                     
  ! If I=2, ------------- helium   ------------------                     
  ! If I=3,M-1, --------- heavy element -------------                     
  ! If I=M, ------------- electrons -----------------                     
  ! For I=M,2M, we get the heat fluxes                                    
  ! For I=N-1, we get the electric field                                  
  ! For I=N, we get the gravitational force g                             

  if(use_lapack) then
     call dgetrf(n, n, delta, n, indx, ierr)
     if (ierr /= 0) then
        stop 'ierr'
     endif
     
     call dgetrs( 'n', n, 1, delta, n, indx, alpha, n, ierr )
     if (ierr /= 0) stop 'ierr2'
     
     call dgetrs( 'n', n, 1, delta, n, indx, nu, n, ierr )
     if (ierr /= 0) stop 'ierr3'
     
      call dgetrs( 'n', n, n, delta, n, indx, gamma, n, ierr )
  else
     call LUDCMP(DELTA,N,N,INDX,D) 
     
     call LUBKSB(DELTA,N,N,INDX,ALPHA) 
     call LUBKSB(DELTA,N,N,INDX,NU)

     do J=1,N 
        do I=1,N 
           GA(I) = GAMMA(I,J) 
        enddo
        call LUBKSB(DELTA,N,N,INDX,GA) 
        do I=1,N 
           GAMMA(I,J) = GA(I) 
        end do
     end do
  endif
  
  ! The results for the coefficients must be multiplied by p/K_0:         
  VARTMP = KO*AC*CC 
  do I=1,M 
     ALPHA(I) = ALPHA(I)*VARTMP 
     NU(I) = NU(I)*VARTMP 
     do J=1,M 
        GAMMA(I,J) = GAMMA(I,J)*VARTMP 
     end do
  end do

  do I=1,M 
     AP(I) = ALPHA(I) 
     AT(I) = NU(I) 
     do J=1,M 
        AX(I,J) = GAMMA(I,J) 
     end do
  end do

  return 
end subroutine DIFFUSION


!********************************************************************   
subroutine lubksb(a,n,np, indx,b)
  implicit none
  integer :: N, np
  real :: A(Np,Np),B(N) 
  integer,dimension(N) :: INDX 

  real :: SUM
  integer :: I,k,J,L  

  k = 0
  do i=1,n
     l = indx(i)
     sum = b(l)
     b(l) = b(i)
     if (k > 0) then
        sum = sum - dot_product(a(i,k:i-1),b(k:i-1))
     else if (sum /= 0.0) then
        k = i
     endif
     b(i) = sum
  enddo

  do i=n,1,-1
     b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
  end do
  return
end subroutine lubksb



!********************************************************               
subroutine LUDCMP(A,N,NP,INDX,D) 
  implicit none                                                    

  real :: D                                                           
  integer :: N,NP 
  real,dimension(NP,NP) :: A
  integer,dimension(N) :: INDX

  integer,parameter :: NMAX = 100
  real,parameter :: TINY = 1.0d-20

  real :: AAMAX,DUM,SUM                                               
  integer I,IMAX,J,K 
  real,dimension(NMAX) :: VV 

  D = 1. 
  do I = 1,N 
     AAMAX = 0. 
     do J = 1,N 
        if (abs(A(I,J)) > AAMAX) AAMAX = abs(A(I,J)) 
     end do
     if (AAMAX == 0.) then 
        write(*,*) 'Singular matrix.' 
        write(66,*)'400 - LUDCMP'
        write(66,*) 'Singular matrix.' 
        stop 
     endif
     VV(I) = 1./AAMAX 
  end do
  do J = 1,N 
     if (J > 1) then 
        do I = 1,J - 1 
           SUM = A(I,J) 
           if (I > 1) then 
              do K = 1,I - 1 
                 SUM = SUM - A(I,K)*A(K,J) 
              end do
              A(I,J) = SUM 
           end if
        end do
     end if

     AAMAX = 0. 
     do I = J,N 
        SUM = A(I,J) 
        if (J > 1) then 
           do K = 1,J - 1 
              SUM = SUM - A(I,K)*A(K,J) 
           end do
           A(I,J) = SUM 
        end if

        DUM = VV(I)*abs(SUM) 
        if (DUM >= AAMAX) then 
           IMAX = I 
           AAMAX = DUM 
        end if
     end do
     if (J /= IMAX) then 
        do K = 1,N 
           DUM = A(IMAX,K) 
           A(IMAX,K) = A(J,K) 
           A(J,K) = DUM 
        end do
        D = -D 
        VV(IMAX) = VV(J) 
     end if

     INDX(J) = IMAX 
     if (J /= N) then 
        if (A(J,J) == 0.) A(J,J) = TINY 
        DUM = 1./A(J,J) 
        do I = J + 1,N 
           A(I,J) = A(I,J)*DUM 
        end do
     end if
  end do
  if (A(N,N) == 0.) A(N,N) = TINY 
  return 

end subroutine LUDCMP


