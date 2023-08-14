subroutine LOCAL(P,T,R,EL,EM,RO,GRAD,CP,PMU,ACCA,CAP,RADIAT,ADIAB,GI,Y,GRSAD)
  use interfaccia
  use costanti
  implicit none

  real :: P,T,R,EL,EM,RO,GRAD,CP,PMU,ACCA,CAP,RADIAT,ADIAB, &
       GI,Y,GRSAD

  real :: DGRAD, PRAD, PRE, TEM
  character(len=15) :: caller="local         "

  PRAD = 1.d9 * arad_3 * (T**4) 
  PRE = (P-PRAD)*1.d15 
  !!GI = 666.8*EM/(R*R)
  GI = 1.d10 * Ggrav * EM / (R * R)
  TEM = T*1.d6 
  call STATE(caller,PRE,TEM,RO,ADIAB,CP,PMU) 
  call KAPPA(RO,TEM,CAP,Y) 
  !!RADIAT = 3.94676*EL*P*CAP/(EM*(T**4))
  RADIAT = 1.d-1 * cte_grad * EL * P * CAP / (EM * (T**4)) 
  DGRAD = RADIAT-ADIAB 
  if(DGRAD < 0.) then 
     GRAD = RADIAT 
     ACCA = 0. 
     GRSAD = 0. 
  else 
     call SUPERA(RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCA,RADIAT) 
     GRAD = GRSAD 
  endif
  return 
end subroutine LOCAL
