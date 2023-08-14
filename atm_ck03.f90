! CASTELLI & KURUCZ (2003)
! La routine non blocca l'esecuzione se Z e Log G non sono in tabella
! ma ritorna al chiamante, il quale procede con atmosfera grigia. 
                                         
subroutine ATM_CK03(GIINI,TE,TAU_RACC) 
  use interfaccia
  use parametri
  use atmosfere

  implicit none

  real :: GIINI,TE,TAU_RACC 

  real,parameter,dimension(NZ_CK) :: ZZ = (/ 0.5, 0.2, 0.0, &
       -0.5, -1.0, -1.5, -2.0, -2.5 /)

  real,parameter,dimension(NZ_CK) :: ZZZ =(/ 5.14d-2, 2.57d-2, 1.62d-2, &
       5.14d-3, 1.62d-3, 5.14d-4, 1.62d-4, 5.14d-5 /)

  real,parameter,dimension(NG_CK) :: GG = (/ 0.0, 0.5, 1.0, 1.5, &
       2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 /)

  real,dimension(NZ_CK,NG_CK,2) :: V_TE
  real,dimension(NZ_CK,2) :: V_AG
  real,dimension(NZ_CK,NG_CK) ::  TE_MIN, TE_MAX 
  integer,dimension(NZ_CK) :: NTE_GG

  integer :: iz, ig, it, ik, ll, nv_atm, ig00, ig01, it0, it1, ik_tau
  integer :: n, m, iop
  real :: ag, tau_r

  real,save, dimension(NT_CK,3,NG_CK,NZ_CK) :: y2tmps
  real,dimension(NT_CK) :: y2tmp
  real,dimension(NG_CK + 1) :: y2tmpa
  real,dimension(NZ_CK) :: y2tmpb
  real :: y, zv
  integer,save :: firsttime = 1
  real,dimension(1200) :: X, YY
  
  ! le tabelle in input
  real,save :: ZX_CK(NZ_CK), ALG_CK(NZ_CK,NG_CK), &
       TE_CK(NZ_CK,NG_CK,NT_CK), TAB_CK(NZ_CK,NG_CK,NT_CK,3,NE_CK)
  integer,save :: NPT_CK(NZ_CK,NG_CK,NT_CK), NTT_CK(NZ_CK,NG_CK)
  
  ! LETTURA TABELLA
  if(ICK03_LETTA /= 1) then
     write(*,*)' ' 
     write(*,'(" ATMOSFERA Castelli & Kurucz (2003)")') 

     open(7,file='ATMOS_CK03_TP.in') 
     ICK03_LETTA = 1 

     do iz=1,NZ_CK 
        read(7,'(8x,F4.1)')ZX_CK(iz) 

        do ig=1,NG_CK 
           read(7,'(9x,F5.2,9x,I3)') ALG_CK(iz,ig),NTT_CK(iz,ig) 

           do it=1,NTT_CK(iz,ig) 
              read(7,'(9x,F11.4,7x,I3)') TE_CK(iz,ig,it),NPT_CK(iz,ig,it) 
              read(7,*) 

              do ik=1,NPT_CK(iz,ig,it) 
                 read(7,*)(TAB_CK(iz,ig,it,ll,ik),ll=1,3) 
              end do
              ! ciclo te                                              
           end do
           ! ciclo g                                                
        end do
        ! ciclo metallicità                                      
     end do
     close(7) 

     write(*,*)' ' 
     write(*,'("  RACCORDO TAU = ",F8.4)')TAU_RACC 
     write(*,*)' ' 

     ! CONTROLLI SULLA LETTURA

     ! [Z/X]:                                                            
     do iz=1,NZ_CK 
        if(ZX_CK(iz) /= ZZ(iz)) then 
           write(*,'(" [Z/X] letto male!! ",2(F4.1,2x))') ZX_CK(iz),ZZ(iz) 
           write(*,'("  ATTENZIONE!! CASINO!!! ")')
           write(66,*)'20 - atm_ck03'
           write(66,'(" [Z/X] letto male!! ",2(F4.1,2x))') ZX_CK(iz),ZZ(iz) 
           write(66,'("  ATTENZIONE!! CASINO!!! ")')
           stop 
        end if
     end do

     write(*,'("-> Blocchi [Z/X] letti bene! ")') 
  
     ! LOG g:                                                            
     do iz=1,NZ_CK 
        do ig=1,NG_CK 
           if(ALG_CK(iz,ig) /= GG(ig)) then 
              write(*,'(" LOG g letto male!! ",2(F4.1,2x))')ALG_CK(iz,ig), &
                   GG(ig)
              write(*,'("  ATTENZIONE!! CASINO!!! ")') 
              write(66,*)'20 - atm_ck03'
              write(66,'(" LOG g letto male!! ",2(F4.1,2x))')ALG_CK(iz,ig), &
                   GG(ig)
              write(66,'("  ATTENZIONE!! CASINO!!! ")') 
              stop 
           end if
        end do
     end do
     write(*,'("-> Blocchi LOG g letti bene! ")') 
     write(*,*)' ' 

     do iz=1,NZ_CK 
        do ig=1,NG_CK 
           NV_ATM = NPT_CK(iz,ig,1) 
           do it=1,NTT_CK(iz,ig) 
              if(NPT_CK(iz,ig,it) == 0) then 
                 write(*,*)iz,ig,it
                 write(66,*)'21 - atm_ck03'
                 write(66,*)iz,ig,it
                 stop 
              end if
              if(NPT_CK(iz,ig,it) < NV_ATM) NV_ATM = NPT_CK(iz,ig,it) 
           end do
        end do
     end do

     write(*,*)'NV_ATM',NV_ATM 
  endif

  ! controllo se ho i valori in tabella 
  ! Z:                                                                   
  if(Z_SUP > ZZZ(1) .or. Z_SUP < ZZZ(NZ_CK)) then 
     return
  end if

  AG = log10(GIINI) 
  ! log G:                                                               
  if(AG > GG(NG_CK) .or. AG < GG(1)) then 
     return
  end if

  ! cerco i due valori più grande e più piccolo di Log G.              
  ig00 = 1 
  ig01 = 1 
  do ig=1,NG_CK 
     if(AG > ALG_CK(1,ig)) then 
        ig00 = ig 
     end if
     if(AG < ALG_CK(1,ig)) then 
        ig01 = ig 
        exit 
     end if
  end do

  ! Teff:                                                                
  do it=1,NTT_CK(1,ig00) 
     if(TE >= TE_CK(1,ig00,it)) exit 
  end do
  it0 = it 

  do it=1,NTT_CK(1,ig01) 
     if(TE >= TE_CK(1,ig01,it)) exit 
  end do
  it1 = it 

  if(it0 > NTT_CK(1,ig00) .or. it1 > NTT_CK(1,ig01)) then 
     write(*,'("  VALORI di TEFF non in tabella!!!!",F10.3)') TE 
     write(66,*)'24 - atm_ck03'
     write(66,'("  VALORI di TEFF non in tabella!!!!",F10.3)') TE 
     stop 
  end if

  ! definisco la TEFF minima e massima per ogni tablla [Z/X], log g:  
  do iz=1,NZ_CK 
     do ig=1,NG_CK 
        TE_MIN(iz,ig) = TE_CK(iz,ig,1) 
        TE_MAX(iz,ig) = TE_CK(iz,ig,NTT_CK(iz,ig)) 
     end do
  end do

  ! controllo quali LOG G hanno TE_MAX > TE                          
  do iz=1,NZ_CK 
     do ig=1,NG_CK 
        if(TE >= TE_MIN(iz,ig) .and. TE <= TE_MAX(iz,ig)) then 
           NTE_GG(iz) = ig 
           
           exit 
        end if
     end do
  end do

  ! cerco il punto in cui TAU = TAU_RACC                         
  do ik=1,NPT_CK(1,1,1) 
     if(TAB_CK(1,1,1,1,ik) >= TAU_RACC) exit 
  end do

  ik_tau = ik 
  ! tau di raccordo                
  TAU_R = TAB_CK(1,1,1,1,ik_tau) 

  ! interpolo prima le TE:                                           
  do iz=1,NZ_CK 
     !! parto dal LOG G più piccolo che ha TE 
     do ig=NTE_GG(iz),NG_CK 

        do it=1,NTT_CK(iz,ig) 
           ! vettore che contiene le TE          
           X(it) = TE_CK(iz,ig,it) 
        end do

        ! colonna temperatura e pressione                  
        do ll=2,3 
           do it=1,NTT_CK(iz,ig) 
              !! prendo T, P al tau vol
              YY(it) = TAB_CK(iz,ig,it,ll,ik_tau) 
           end do

           N = NTT_CK(iz,ig) 
           if(n == 1) return ! errore di interpolazione

           if(firsttime == 1) then
              call makespline(x,yy,n,y2tmp)
              do it=1,NTT_CK(iz,ig)
                 y2tmps(it,ll,ig,iz) = y2tmp(it)
              end do
           else
              do it=1,NTT_CK(iz,ig)
                 y2tmp(it) = y2tmps(it,ll,ig,iz)
              end do
           endif
           call usesplint(x,yy,y2tmp,n,TE,y)
           V_TE(iz,ig,ll-1) = y
        end do
     end do
  end do
  firsttime = 0

  ! interpolo log G:                                                 
  do iz=1,NZ_CK 
     !! solo i valori di LOG G che hanno le TE 
     do ig=NTE_GG(iz),NG_CK 
        !! vettore delle log G       
        X(ig+1-NTE_GG(iz)) = ALG_CK(iz,ig) 
     end do

     !! colonna temperatura e pressione                   
     do ll=1,2 
        do ig=NTE_GG(iz),NG_CK 
           YY(ig+1-NTE_GG(iz)) = V_TE(iz,ig,ll) 
        end do
       
        N = NG_CK - NTE_GG(iz) + 1 
        if(n == 1) return ! errore di interpolazione
        call makespline(x,yy,n,y2tmpa)
        call usesplint(x,yy,y2tmpa,n,AG,y)
        V_AG(iz,ll) = y
     end do
  end do

  ! interpolo in [Z/X]:                                              
  ! passo il logaritmo perchè la SPLINE ha dei grossi problemi con i numeri 
  ! piccoli e vicini
  Zv = log10(Z_SUP) 

  do iz=1,NZ_CK 
     X(iz) = log10(ZZZ(NZ_CK+1-iz)) 
  end do

  ! colonna temperatura e pressione                    
  do ll=1,2 
     do iz=1,NZ_CK 
        YY(iz) = V_AG(NZ_CK+1-iz,ll) 
     end do
     
     N = NZ_CK 
     call makespline(x,yy,n,y2tmpb)
     call usesplint(x,yy,y2tmpb,n,Zv,y)
     V_ATM(ll+1) = y
     if( .not. y > 0.0) then
        write(66,*)'22 - atm_ck03'
        write(*,*) 'atm_ck03: V_ATM negativa'
        stop 
     endif
  end do

  ! colonna di tau                               
  V_ATM(1) = TAU_R 
  N_CK03 = 1 
  
  return 
end subroutine ATM_CK03
