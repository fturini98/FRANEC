! La routine zone_convettive identifica le zone convettive etichettando i mesh
! nel modo seguente:
! 1 = inizio zona convettiva
! 2 = interno zona convettiva
! 3 = fine zona convettiva
! 4 = mesh isola convettivo

subroutine zone_convettive(MAX, do_shell)
  use parametri
  use fisica
  use overshoot
  use nummod
  use zone_conv

  implicit none

  integer :: MAX
  logical :: do_shell

  integer :: jf
  integer :: convettivo, elimina

  logical,parameter :: debug = .false., do_elimina = .false.
  integer,dimension(max+1):: mygrad

  convettivo = 0
  n_conv = 0

  ! max = maxme -1 in input
  call do_overshoot(1, max+1, mygrad)
  ! ora in mygrad ho 1 se il mesh è convettivo secondo G(6,) e overshooting
  ! ho 0 altrimenti

  do jf=1,MAX

     
     ! gestisco la fine della zona convettiva al mesh precedente
     if(convettivo == 3 .or. convettivo == 4) convettivo = 0
     
     ! identifica zona convettiva
     if(jf == 1) then
        convettivo = 0
        ! inizio zona
        if(mygrad(JF) > 0 .and. mygrad(JF+1) > 0) convettivo = 1
        ! mesh "isola" convettivo      
        if(mygrad(JF) > 0 .and. mygrad(JF+1) == 0) convettivo = 4
     else if(jf == MAX) then
        ! chiudo la zona convettiva di bordo
        if(vconv(MAX-1) == 2 .or. vconv(MAX-1) == 1) then
           convettivo = 3
        else
           convettivo = 0
        endif
     else 
        ! inizio zona
        if(do_shell .and. mygrad(JF) > 0 .and. mygrad(JF-1) == 0) convettivo = 1
        ! interno zona
        if(convettivo == 1 .and. mygrad(JF) > 0 .and. mygrad(JF-1) > 0 &
             .and. mygrad(JF+1) > 0.0) convettivo = 2
        ! fine zona
        if((convettivo == 1 .or. convettivo == 2) .and. mygrad(JF) > 0 &
             .and. mygrad(JF+1) == 0) convettivo = 3
           
        ! mesh "isola" convettivo 
        if(mygrad(JF) > 0 .and. mygrad(JF-1) == 0 .and. mygrad(JF+1) == 0) &
             convettivo = 4
     endif

     vconv(jf) = convettivo
     ! aggiorno numero zone e boundary delle zone
     if(convettivo == 1) then
        n_conv = n_conv + 1
        b_conv(1, n_conv) = jf
     else if(convettivo == 3) then
        b_conv(2, n_conv) = jf
     endif
  end do

  ! elimino possibili isole radiative <=2 mesh tra due convettive
  if(do_elimina) then
     elimina = 0
     do jf=n_conv-1,1,-1
        if(b_conv(1,jf+1) - b_conv(2,jf) <= 2) then
           elimina = elimina + 1
           b_conv(2,jf) = b_conv(2,jf+1)
        endif
     end do
     n_conv = n_conv - elimina 
  endif

  if(debug .eqv. .true.) then
     write(165,*) nmd, n_conv
     write(165, *) b_conv(1:2, 1:n_conv)
  endif
     
end subroutine zone_convettive


subroutine do_overshoot(INI, IFI, mygrad)
  use overshoot
  use fisica
  use chimic
  use mesh
  use nummod
  
  implicit none

  integer:: ini, ifi, ifin

  real :: ADMM_M, admm_MP1
  logical :: printout = .false.
  real :: rr, pp, em, pm, rm, emm, hpmm, hpm
  real :: GICOVER, GICOVERR, ERRO, DIFOV, GICO
  integer :: trovato, i, jf, m
  integer,dimension(maxme):: mygrad
  logical,parameter :: do_coreHe_over = .false.


  mygrad = 0
  ! overshooting core He
  if(xxx(1,1) <= 0.0 .and. xxx(3,1) > 0.0) then
     if(do_coreHe_over) then
        KOVER = 2
        par_OVER = 0.1
     else
        kover = 0
     endif
  else if(xxx(3,1) <= 0.0) then
     KOVER = 0
  endif

  ! Inserimento overshooting

  ! L3 = mesh fine core convettivo overshoot
  L3 = 0

  if(KOVER >= 1 .and. G(6,1) > 0.0) then
     call do_core_overshoot(ini, ifi, mygrad)

     gico = G(5,l3+1)
     if(kover == 2) call do_shell_overshoot(l3, ifi, mygrad)

  else
     do i=1,IFI
        if(G(6,i) < 0) then
           mygrad(i) = 0
        else
           mygrad(i) = 1
        endif
     end do
  endif

end subroutine do_overshoot


subroutine do_core_overshoot(ini, ifi, mygrad)
 use overshoot
  use fisica
  use chimic
  use mesh
  use nummod
  
  implicit none

  integer:: ini, ifi, ifin

  real :: ADMM_M, admm_MP1
  logical :: printout = .false.
  real :: rr, pp, em, pm, rm, emm, hpmm, hpm
  real :: GICOVER, GICOVERR, ERRO, DIFOV, GICO
  integer :: trovato, i, jf, m
  integer,dimension(maxme):: mygrad

  do jf=INI,IFI-1
     ! trovo il bordo del core convettivo
     if(G(6,JF) >= 0.0 .and. G(6,JF+1) < 0.0) exit
  end do
     
  if(jf < IFI-1) then
     ! ho trovato un core convettivo. Aggiorno il bordo con overshoot
     ! il nuovo bordo e' al mesh L3
     PM = G(3,JF)*1.d17
     RM = G(1,JF)*1.d10
     EMM = G(5,JF)*1.d33
     HPMM = 1.884d8*PM*(RM**4)/EMM
     HPM = HPMM*1.d-33
        
     do M = JF, ifi-1
        ADMM_m = G(5,M)-(G(5,JF)+par_OVER*HPM)
        ADMM_mp1 = G(5,M+1)-(G(5,JF)+par_OVER*HPM)
        if( ADMM_m <= 0.0 .and. ADMM_mp1 > 0.0 ) exit
     end do

     L3 = M
     mygrad(1:L3) = 1
     do i=l3+1,IFI
        if(G(6,i) < 0) then
           mygrad(i) = 0
        else
           mygrad(i) = 1
        endif
     end do
     
     if(printout) then
        GICOVER = G(5,L3+1)
        GICOVERR = G(5,JF)+par_OVER*HPM
        ERRO = G(5,L3+1)-(G(5,JF)+par_OVER*HPM)
        DIFOV = G(5,JF)+par_OVER*HPM-G(5,JF+1)
        GICO = G(5,JF+1)
        write(50,7362) NMD,GICOVER,GICOVERR,ERRO,GICO, &
             DIFOV,HPM,JF, L3
     endif
  endif
  
7362 format(I5,6e13.5,2I5)
end subroutine do_core_overshoot



subroutine do_shell_overshoot(ini, ifi, mygrad)
 use overshoot
  use fisica
  use chimic
  use mesh
  use nummod
  
  implicit none

  integer:: ini, ifi, ifin

  real :: ADMM_M, admm_MP1
  logical :: printout = .false.
  real :: rr, pp, em, pm, rm, emm, hpmm, hpm
  real :: GICOVER, GICOVERR, ERRO, DIFOV, GICO
  integer :: trovato, i, jf, m, b_int, b_ext, m1, m2
  integer,dimension(maxme):: mygrad

  b_int = 0
  b_ext = 0

  do jf=INI,IFI-1
     ! trovo il bordo interno/esterno della shell convettiva
     if(G(6,JF) <= 0.0 .and. G(6,JF+1) > 0.0 .and. xxx(1,jf) <= 0) then
        b_int = jf+1
     else if (b_int > 0 .and. G(6,JF) > 0.0 .and. G(6,JF+1) <= 0.0) then
        b_ext = jf
        exit
     end if
  end do

   write(65,*) "shell bordo", b_int, b_ext

  ! no shell di He
  if(b_int == 0 ) return
  
  ! ho trovato ua shell di He convettiva. Aggiorno il bordo con overshoot
  PM = G(3,b_int)*1.d17
  RM = G(1,b_int)*1.d10
  EMM = G(5,b_int)*1.d33
  HPMM = 1.884d8*PM*(RM**4)/EMM
  HPM = HPMM*1.d-33
  
  ! matt: controllare la logica per la ricerca del bordo...
  do M = b_int, ini, -1
     ADMM_m = G(5,M)-(G(5,b_int)-par_OVER*HPM)
     ADMM_mp1 = G(5,M-1)-(G(5,b_int)-par_OVER*HPM)
     write(65,*) "in", m, ADMM_m, ADMM_mp1
     if( ADMM_m > 0.0 .and. ADMM_mp1 <= 0.0 ) exit
  end do
  ! nuovo bordo interno
  m1 = m

  PM = G(3,b_ext)*1.d17
  RM = G(1,b_ext)*1.d10
  EMM = G(5,b_ext)*1.d33
  HPMM = 1.884d8*PM*(RM**4)/EMM
  HPM = HPMM*1.d-33
  
  do M = b_ext, ifi
     ADMM_m = G(5,M)-(G(5,b_ext)+par_OVER*HPM)
     ADMM_mp1 = G(5,M+1)-(G(5,b_ext)+par_OVER*HPM)
     if( ADMM_m <= 0.0 .and. ADMM_mp1 > 0.0 ) exit
  end do
  ! nuovo bordo esterno
  m2 = m

  mygrad(M1:M2) = 1
     
end subroutine do_shell_overshoot
