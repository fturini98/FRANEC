subroutine SUPERA(RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCO,RADIAT)
  use interfaccia
  use parametri
  use varie
  use costanti
  implicit none

  real :: RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCO,RADIAT

  integer :: L, stato
  real :: GUESS1, B, AML, CV, E, FAC1, FAC2
  real :: HP, PMU1, P1, PRE1, Q,  STAB, TEM1, T1, Y1, VAR
  real :: Z, tmp1, tmp2, step, deriv, tmpv

  real, parameter :: SG = 5.673d-5, A0 = 9./4.
  real, parameter :: delta = 1.d-6, delta1 = 1.d-7
  real, parameter :: due_rad_due = 2.82842712474619
  real,parameter :: rad2m1 = 0.707106781186547
  
  ! disabilito il calcolo delle quantita' non volute in state
  real,parameter :: dum = 1.d99
  character(len=15) :: caller="supera        "
  

  T1 = 1.01*T 
  PRE1 = (P-1.d9*arad_3*(T1**4))*1.d15 
  TEM1 = T1*1.d6 

  call STATE(caller,PRE1,TEM1,dum,dum,dum,PMU1) 

  Q = 1.-(PMU1-PMU)/(.01*PMU) 
  if(Q < 1.d-20) Q = 1.d-20 
  T1 = T*1.d6 
  P1 = P*1.d15 
  HP = P1/(RO*GI) 

  AML = ALFA*HP 
  STAB = RADIAT-ADIAB 
  tmpv = sqrt(Q*RO/P1)
  E = (CP*CAP*GI*RO**2*tmpv*rad2m1*AML**2)/(48.*SG*T1**3)
  B = (E**2*STAB/A0)**(1./3.) 

  tmpv = GI*tmpv*AML/due_rad_due

  if(B < 1.) then 
     tmp1 = A0*B**3
     tmp2 = A0*B**2
     Z = tmp2**3 * (1. - 3.*tmp1 + 9.*tmp1**2 - 3.*tmp2**3)
     Y1 = Z**(1./3.) 

     ! Metodo di Newton per trovare lo zero di GUESS1
     ! al variare di Y1
     stato = 0
     do L=1,100
        GUESS1 = Y1 + B*Y1**2 + tmp2*Y1**3 - tmp2 
        VAR = abs(GUESS1) 
        if(VAR < delta) then
           stato = 1
           exit
        endif
        deriv = 1. + 2.*B*Y1 + 3.*tmp2*Y1**2
        step = GUESS1/deriv
        if(abs(step) < delta1) then
           stato = 1
           exit
        endif
        Y1 = Y1-step
     end do
     if(stato == 1) then
        ! ho la soluzione, torno al chiamante
        Z = Y1**3 
        CV = B*Y1
        GRSAD = (1.-Z)*RADIAT+Z*ADIAB 
        ACCO = abs(CV/E)*tmpv
        return 
     else
        ! se sono qui non ho trovato lo zero, stampo errore e torno indietro
        write(*,100)
        return
     endif
  else 
     tmp2 = A0*B**2
     FAC1 = (1.+B)/tmp2
     FAC2 = (1.+2.*B)/3./(1.+B) 
     Z = 1.-FAC1*(1.-FAC1*FAC2+(1./9.)*(9.*FAC2**2-1)*FAC1**2) 
     Y1 = Z**(1./3.) 
     
     ! Metodo di Newton per trovare lo zero di GUESS1
     ! al variare di Y1
     stato = 0
     do l=1,100
        GUESS1 = 1. - Y1**3 - (Y1 + B*Y1**2)/tmp2
        VAR = abs(GUESS1) 
        if(VAR < delta) then
           stato = 1
           exit
        endif
        deriv = -3.*Y1**2 - (1. + 2.*B*Y1)/tmp2
        step = GUESS1/deriv
        if(abs(step) < delta1) then
           stato = 1
           exit
        endif
        Y1 = Y1 - step
     end do
     if(stato == 1) then
        ! ho la soluzione, torno al chiamante
        Z = Y1**3 
        CV = B*Y1
        GRSAD = (1.-Z)*RADIAT+Z*ADIAB 
        ACCO = abs(CV/E)*tmpv
        return 
     else
        ! se sono qui non ho trovato lo zero, stampo errore e torno indietro
        write(*,100)
        return
     endif
  endif

100 format(1X,'ATTENZIONE MIX-LEN NON CONVERGE') 
end subroutine SUPERA
