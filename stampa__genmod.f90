        !COMPILER-GENERATED INTERFACE MODULE: Mon May 28 09:46:39 2012
        MODULE STAMPA__genmod
          INTERFACE 
            SUBROUTINE STAMPA(TEMPO,DMAS,NMOD,NABLA,KSD,MAXMEIN,MAXMV,  &
     &IFVA,IFMO,LUMIN,TMAX,NSMORZA,FASE,LOGCNO,KOVER)
              USE INTERFACCIA
              REAL(KIND=8) :: TEMPO
              REAL(KIND=8) :: DMAS
              INTEGER(KIND=4) :: NMOD
              INTEGER(KIND=4) :: NABLA
              INTEGER(KIND=4) :: KSD
              INTEGER(KIND=4) :: MAXMEIN
              INTEGER(KIND=4) :: MAXMV
              INTEGER(KIND=4) :: IFVA
              INTEGER(KIND=4) :: IFMO
              TYPE (LUMINOSITY) :: LUMIN
              REAL(KIND=8) :: TMAX
              INTEGER(KIND=4) :: NSMORZA
              INTEGER(KIND=4) :: FASE
              INTEGER(KIND=4) :: LOGCNO
              INTEGER(KIND=4) :: KOVER
            END SUBROUTINE STAMPA
          END INTERFACE 
        END MODULE STAMPA__genmod
