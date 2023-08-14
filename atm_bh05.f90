! *** BROTT & HAUSCHILDT (2005):                                       
subroutine ATM_BH05(GIINI,TE,TAU_RACC) 
  use interfaccia
  use parametri
  use atmosfere

  implicit none

  real :: GIINI, TE, TAU_RACC 

  real,parameter,dimension(NZ_BH) :: ZZ =(/ 0.5, 0.0, -0.5, -1.0, &
       -1.5,  -2.0, -2.5, -3.0, -3.5, -4.0 /)

  real,parameter,dimension(NZ_BH) :: ZZZ =(/ 5.29d-2, 1.73d-2, 5.60d-3, &
       1.76d-3, 5.60d-4, 1.76d-4, 5.60d-5, 1.76d-5, 5.60d-6, 1.76d-6 /)
  !!real,parameter,dimension(NZ_BH) :: ZZZ =(/ 5.28040022d-2, &
  !!     & 1.73236761d-2, 5.5438904d-3, 1.7597988d-3, &
  !!     & 5.571726d-4, 1.762590d-4, 5.57447d-5, &
  !!     & 1.76287d-5, 5.5747d-6, 1.7629d-6 /)
  real,parameter,dimension(NG_BH) :: GG = (/-0.5, 0.0, 0.5, 1.0, & 
       1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5 /)

  real,dimension(NZ_BH,NG_BH,2) ::  V_TE
  real,dimension(NZ_BH,2) :: V_AG 
  real,dimension(NZ_BH,NG_BH) :: TE_MIN, TE_MAX 
  integer,dimension(NZ_BH) ::  NTE_GG 

  integer :: iz, ig, it, ik, ll, nv_atm, ig00, ig01, it0, it1, ik_tau
  integer :: n, m, iop
  real :: ag, tau_r

  real,save, dimension(NT_BH,3,NG_BH,NZ_BH) :: y2tmps
  real,dimension(NT_BH) :: y2tmp
  real,dimension(NG_BH + 1) :: y2tmpa
  real,dimension(NZ_BH) :: y2tmpb
  real :: y, zv
  integer,save :: firsttime = 1
  real,dimension(1200) :: X, YY

  ! le tabelle in input
  real,save :: ZX_BH(NZ_BH), ALG_BH(NZ_BH,NG_BH), &
       TE_BH(NZ_BH,NG_BH,NT_BH), TAB_BH(NZ_BH,NG_BH,NT_BH,3,NE_BH)
  integer,save :: NPT_BH(NZ_BH,NG_BH,NT_BH), NTT_BH(NZ_BH,NG_BH)

  ! LETTURA TABELLA
  if(IBH05_LETTA /= 1) then
     write(*,*)' ' 
     write(*,'(" ATMOSFERA Brott & Hauschildt (2005)")') 

     open(7,file='ATMOS_BH05_TP.in') 
     IBH05_LETTA = 1 


     do iz=1,NZ_BH 
        read(7,'(8x,F4.1)') ZX_BH(iz) 
        do ig=1,NG_BH 
           read(7,'(9x,F5.2,9x,I3)') ALG_BH(iz,ig),NTT_BH(iz,ig) 
           do it=1,NTT_BH(iz,ig) 
              read(7,'(9x,F11.4,7x,I3)') TE_BH(iz,ig,it),NPT_BH(iz,ig,it) 
              read(7,*) 

              do ik=1,NPT_BH(iz,ig,it) 
                 read(7,*)(TAB_BH(iz,ig,it,ll,ik),ll=1,3) 
              end do
              ! ciclo te                                              
           end do
           ! ciclo g                                                
        end do
        ! ciclo metallicità                                      
     end do
     close(7) 

     write(*,*)' ' 
     write(*,'("  RACCORDO TAU = ",F8.4)') TAU_RACC 
     write(*,*)' ' 

     ! CONTROLLI SULLA LETTURA 

     ! [Z/X]:                                                            
     do iz=1,NZ_BH 
        if(ZX_BH(iz) /= ZZ(iz))then 
           write(*,'(" [Z/X] letto male!! ",2(F4.1,2x))') ZX_BH(iz),ZZ(iz) 
           write(*,'("  ATTENZIONE!! CASINO!!! ")') 
           write(66,*)'10 - atm_bh05 ' 
           write(66,'(" [Z/X] letto male!! ",2(F4.1,2x))') ZX_BH(iz),ZZ(iz) 
           write(66,'("  ATTENZIONE!! CASINO!!! ")') 
           stop
        end if
     end do

     write(*,'("-> Blocchi [Z/X] letti bene! ")') 

     ! LOG g:                                                            
     do iz=1,NZ_BH 
        do ig=1,NG_BH 
           if(ALG_BH(iz,ig) /= GG(ig))then 
              write(*,'(" LOG g letto male!! ",2(F4.1,2x))')ALG_BH(iz,ig), &
                   GG(ig)
              write(*,'("  ATTENZIONE!! CASINO!!! ")') 
              write(66,*)'10 - atm_bh05 ' 
              write(66,'(" LOG g letto male!! ",2(F4.1,2x))')ALG_BH(iz,ig), &
                   GG(ig)
              write(66,'("  ATTENZIONE!! CASINO!!! ")')
              stop
           end if
        end do
     end do

     write(*,'("-> Blocchi LOG g letti bene! ")') 
     write(*,*)' ' 

     do iz=1,NZ_BH 
        do ig=1,NG_BH 
           NV_ATM = NPT_BH(iz,ig,1) 
           do it=1,NTT_BH(iz,ig) 
              if(NPT_BH(iz,ig,it) == 0) then 
                 write(*,*) iz,ig,it 
                 write(66,*)'11 - atm_bh05 ' 
                 write(66,*) iz,ig,it 
                 stop 
              end if
              if(NPT_BH(iz,ig,it) < NV_ATM) then 
                 NV_ATM = NPT_BH(iz,ig,it) 
                 write(*,*)'   NPT non uguale nelle tabelle!!! ' 
                 write(66,*)'12 - atm_bh05 ' 
                 write(66,*)'   NPT non uguale nelle tabelle!!! ' 
                 stop 
              end if
           end do
        end do
     end do

     write(*,*)'NV_ATM',NV_ATM 
  endif

  ! Se Teff > 10000  esce e va nelle Castelli & Kurucz (2003)
  if(TE > 10000.0) return 

  AG = log10(GIINI) 

  ! controllo se ho i valori in tabella
  ! Z:                                                                   
  if(Z_SUP > ZZZ(1) .or. Z_SUP < ZZZ(NZ_BH))then 
     write(*,'(" Z non in tabella!!!! ",F4.1)')Z_SUP 
     write(66,*)'13 - atm_bh05 ' 
     write(66,'(" Z non in tabella!!!! ",F4.1)')Z_SUP 
     stop 
  end if

  ! log G:                                                               
  if(AG > GG(NG_BH) .or. AG < GG(1))then 
     write(*,'(" LOG G non in tabella!!!! ",F4.1)')AG 
     write(66,*)'14 - atm_bh05 ' 
     write(66,'(" LOG G non in tabella!!!! ",F4.1)')AG 
     stop 
  end if

  ! cerco i due valori più grande e più piccolo di Log G.              
  ig00 = 1 
  ig01 = 1 
  do ig=1,NG_BH  
     if(AG > ALG_BH(1,ig))then 
        ig00 = ig 
     end if
     if(AG < ALG_BH(1,ig))then 
        ig01 = ig 
        exit 
     end if
  end do

  do it=1,NTT_BH(1,ig00) 
     if(TE >= TE_BH(1,ig00,it)) exit 
  end do

  it0 = it 

  do it=1,NTT_BH(1,ig01) 
     if(TE >= TE_BH(1,ig01,it)) exit 
  end do

  it1 = it 

  if(it0 > NTT_BH(1,ig00) .or. it1 > NTT_BH(1,ig01)) then 
     write(*,'("  VALORI di TEFF non in tabella!!!!",F10.3)')TE 
     write(66,*)'15 - atm_bh05 ' 
     write(66,'("  VALORI di TEFF non in tabella!!!!",F10.3)')TE 
     stop 
  end if

  ! definisco la TEFF minima e massima per ogni tablla [Z/X], log G:  
  do iz=1,NZ_BH 
     do ig=1,NG_BH 
        TE_MIN(iz,ig) = TE_BH(iz,ig,1) 
        TE_MAX(iz,ig) = TE_BH(iz,ig,NTT_BH(iz,ig)) 
     end do
  end do

  ! controllo quali LOG G hanno TE_MAX > TE                          
  do iz=1,NZ_BH 
     do ig=1,NG_BH 
        if(TE >= TE_MIN(iz,ig) .and. TE <= TE_MAX(iz,ig)) then 
           NTE_GG(iz) = ig 
           exit 
        end if
     end do
  end do

  ! cerco il punto in cui TAU = TAU_RACC 
  do ik=1,NPT_BH(1,1,1) 
     if(TAB_BH(1,1,1,1,ik) >= TAU_RACC) exit 
  end do

  ik_tau = ik 
  ! tau di raccordo                
  TAU_R = TAB_BH(1,1,1,1,ik_tau) 

  ! comincio a interpolare                                           
  ! interpolo prima le TE:                                           
  do iz=1,NZ_BH 
     ! parto dal LOG G piu' piccolo che ha TE 
     do ig=NTE_GG(iz),NG_BH 

        do it=1,NTT_BH(iz,ig) 
           ! vettore che contiene le TE          
           X(it) = TE_BH(iz,ig,it) 
        end do

        ! colonna temperatura e pressione                  
        do ll=2,3 
           do it=1,NTT_BH(iz,ig) 
              ! prendo T, P al tau voluto
              YY(it) = TAB_BH(iz,ig,it,ll,ik_tau) 
           end do

           N = NTT_BH(iz,ig)
           if(n == 1) return ! errore di interpolazione
           ! la prima volta costruisco la spline
           ! le volte seguenti la valuto e basta
           if(firsttime == 1) then
              call makespline(x,yy,n,y2tmp)
              do it=1,NTT_BH(iz,ig)
                 y2tmps(it,ll,ig,iz) = y2tmp(it)
              end do
           else
              do it=1,NTT_BH(iz,ig)
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
  do iz=1,NZ_BH 

     ! solo i valori di LOG G che hanno le TE 
     do ig=NTE_GG(iz),NG_BH 
        ! vettore delle log G       
        X(ig+1-NTE_GG(iz)) = ALG_BH(iz,ig) 
     end do

     ! colonna temperatura e pressione                   
     do ll=1,2 
        do ig=NTE_GG(iz),NG_BH 
           YY(ig+1-NTE_GG(iz)) = V_TE(iz,ig,ll) 
        end do

        N = NG_BH - NTE_GG(iz) + 1 
        if(n == 1) return ! errore di interpolazione
        call makespline(x,yy,n,y2tmpa)
        call usesplint(x,yy,y2tmpa,n,AG,y)
        V_AG(iz,ll) = y
     end do
  end do

  ! interpolo in Z:                                                  

  ! passo il logaritmo perchè la SPLINE ha dei grossi problemi con i numeri
  ! piccoli e vicini                        
  do iz=1,NZ_BH 
     X(iz) = log10(ZZZ(NZ_BH+1-iz)) 
  end do

  Zv = log10(Z_SUP) 
  ! colonna temperatura e pressione                    
  do ll=1,2 
     do iz=1,NZ_BH 
        YY(iz) = V_AG(NZ_BH+1-iz,ll) 
     end do

     N = NZ_BH
     if(n == 1) return ! errore di interpolazione
     call makespline(x,yy,n,y2tmpb)
     call usesplint(x,yy,y2tmpb,n,Zv,y)
     V_ATM(ll+1) = y
     if( .not. y > 0.0) then
        write(66,*)'12 - atm_bh05'
        write(*,*) 'atm_bh05: V_ATM negativa'
     endif
  end do

  ! colonna di tau                               
  V_ATM(1) = TAU_R 

  N_BH05 = 1 
  return 

end subroutine ATM_BH05
