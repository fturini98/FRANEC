  !    QUESTA SUBROUTINE SERVE PER ALTERARE IL NUMERO DEI MESH            
subroutine OPTIM(MAXME,SCALA,PROVV,ECNO, URA,ULA,UPA,UTA,UMA,VMM, fase,T_DM,NMD) 
  use interfaccia
  use fisica
  use intero
  use strut
  use serchi
  use nuconv
  use chimic
  use chim
  use numer
  use Dark_Matter

  implicit none

  integer :: MAXME, fase, NMD
  real :: SCALA, PROVV, ECNO
  real,dimension(4) :: URA, ULA, UPA, UTA, UMA 
  real,dimension(5) :: VMM

  real :: PPROV, ECNOV
  real,dimension(MELE) :: PX
  real,dimension(6) :: PG
  integer :: n, i, j, k, l, kk, ll, nn, mam
  real :: ym1, ym2, dym2, gig, ul, up, ut, ur, um, p1, p2, p3, p4, p5, PHe3
  real :: gmx, xid, xel, uxh, uh4, wxh, wh4, dm , UHe3
  real :: XHe3_max, XHe3_min, PXHe3_tmp
  real, parameter :: UHe3_ori = 1.0d-2

  integer, parameter :: Nini = 8
  real, dimension(LIM) :: XHe3_tmp
  integer :: iskip_he3 = 1 !! 0/1 = uso/non-uso He3 nella OPTIM 

  !##################################
  !Variabili infittimento mesh per DM
  !##################################
  real :: Lumi_DM,&!Luminosità totale DM
         T_DM,&!Temperatura della dark matter
         Lumi_mesh,&!Luminosità mesh per DM
         Lumi_mesh_prec,&!Luminosità del mesh precedente
         Lumi_mesh_sommati!Luminosità dei mesh sommati(serve per sfittimento)
   
   integer :: infittisci_per_DM !Falg per decidere se infittire il mesh a causa della DM

  !! Ema & Matt: Massimo He3
  XHe3_max = maxval( XXX(2,1:MAXME-1) )
  XHe3_min = 5.d-2 * XHe3_max

  XHe3_tmp(1:MAXME) = XXX(2,1:MAXME)
  
  SCALA = 1. 
  PPROV = PROVV 
  ECNOV = ECNO 
  YM1 = 0.
  YM2 = 0.
  if(BB(1) > 1.d-2 .and. XXX(1,1) <= 0. .and. XXX(3,1) > 0.) then 
     YM1 = GICO/EMTOT-.01 
     if(YM1 < 0.) YM1 = 0. 
     DYM2 = (20.0-(EMTOT/1.989))*.005 
     if( DYM2 <= 0. )  DYM2 = 0. 
     if((EMTOT/1.989) > 2.) then 
        YM2 = YM1+DYM2 
     else 
        YM2 = YM1+0.05+DYM2 
     endif
  endif
  if(PROVV >= 2.) SCALA = 1.2*(int(PPROV/2.)) 
  if(ECNO >= 2.)   SCALA = 1.1*(int(ECNOV/5. )) 
  if(SCALA < 1.) SCALA = 1. 

  ! INIZIA PARTE INFITTIMENTO 
  N = Nini
  do while(MAXME <= LIM-5 .and. N < MAXME)
     N = N+1
      
      if (on_off_fine_mesh_DM==1 .and. N>Nini+N_mesh_aggiunti_DM) then!Conviene partire da 3 mesh più in la per l'infittimento se no infittisce troppo
         !al centro e quseto fa impazzire la ciacio
         Lumi_mesh_prec=(epsi_DM(N-1)*(G(5,N)-G(5,N-1))*1e33)!La luminosità della DM nel mesh precedente al mesh N(Si usa il mesh precedente 
         !perche l'intervallo in massa è tra N-1 e N)

         !Se nel mesh il contributo della DM alla luminsoità positiva/negativa è maggiore di una frazione definita
         !rispetta a quella totale positiva/negativa devo infittire il mesh
         if ( Lumi_mesh_prec>=0 ) then
            if ( Lumi_mesh_prec>Lumi_DM_positiva/N_min_mesh_DM ) then 
               infittisci_per_DM=1
            else
               infittisci_per_DM=0
            end if
         else
            if (abs(Lumi_mesh_prec)>abs(Lumi_DM_negativa)/N_min_mesh_DM) then
               infittisci_per_DM=1
            else
               infittisci_per_DM=0
            endif
         end if

      else
         infittisci_per_DM=0
      endif

     ! ricerca intervallo di VMM da utilizzare
     GIG = G(5,N)/EMTOT 
     do I=1,4 
        if(GIG >= VMM(I) .and. GIG < VMM(I+1)) exit
     end do
     if(I > 4) I = 4 
     UL = ULA(I)*SCALA 
     UP = UPA(I) 
     UT = UTA(I) 
     UR = URA(I) 
     UM = UMA(I)
     UHe3 = UHe3_ori
     
     P1 = abs(G(1,N)-G(1,N-1))/G(1,N-1) 
     P2 = abs(G(2,N)-G(2,N-1))/G(2,MAXME) 
     P3 = abs(G(3,N)-G(3,N-1))/G(3,N-1) 
     P4 = abs(G(4,N)-G(4,N-1))/G(4,N-1) 
     P5 = abs(G(5,N)-G(5,N-1))/EMTOT

     !! Ema & Matt: infittimento He3:
     PHe3 = 0.0d0
     if(XHe3_tmp(N-1) >= abs(XHe3_min) .and. iskip_He3 <= 0) then
        PHe3 = abs( (XHe3_tmp(N) - XHe3_tmp(N-1)) / XHe3_tmp(N-1) )
     endif
     
     GMX = G(5,N)/EMTOT 
     if(BB(1) > 0 .and. GMX >= YM1 .and. GMX <= YM2 .and. XXX(1,1) <= 0) then 
        UR = URA(I)/10. 
        UL = ULA(I)*SCALA/10.
        UP = UPA(I)/5.
        UT = UTA(I)/10. 
        UM = UMA(I)/10. 
     endif

     ! controllo se almeno un vincolo non e' soddisfatto
     if(P1 >= UR .or. P2 >= UL .or. P3 >= UP .or. P4 >= UT .or. P5 >= UM .or. &
          PHe3 >= UHe3 .or. infittisci_per_DM==1) then

        !!if(PHe3 > UHe3) then
        !!   write(*,*)'infittimento He3:', MAXME
        !!endif
        
        !  CALCOLO VALORI DI FISICA E CHIMICA DA ASSEGNARE AL NUOVO MESH   
        do J=1,6 
           PG(J) = (G(J,N-1)+G(J,N))/2. 
        end do
        do J=1,MELE 
           PX(J) = XXX(J,N-1) 
        end do
        PXHe3_tmp = (XHe3_tmp(N-1) + XHe3_tmp(N) ) / 2.
        
        XID = XSERV(1,N-1) 
        XEL = XSERV(2,N-1) 
        !  INSERIMENTO NUOVO MESH  
        do K=N,MAXME 
           L = MAXME-K+N 
           do I=1,6 
              G(I,L+1) = G(I,L) 
           end do
           do J=1,MELE 
              XXX(J,L+1) = XXX(J,L) 
           end do
           XHe3_tmp(L+1) = XHe3_tmp(L)
           
           IN(L+1) = IN(L) 
           XSERV(1,L+1) = XSERV(1,L) 
           XSERV(2,L+1) = XSERV(2,L) 
        end do
        do I=1,6 
           G(I,N) = PG(I) 
        end do
        do I=1,MELE 
           XXX(I,N) = PX(I) 
        end do
        XHe3_tmp(N) = PXHe3_tmp
        
        IN(N) = 0 
        XSERV(1,N) = XID 
        XSERV(2,N) = XEL 
        MAXME = MAXME+1
        
        if ( on_off_fine_mesh_DM==1 ) then !Se non sono interessato all'infittimento su criteri basati sulla DM non mi ricalcolo l'epsi, questo accellera il processo.
         call epsi_DM_routine(T_DM,Lumi_DM,NMD)!Calcolo nuovamente l'array dell'epsi_DM per la nuova struttura coi mesh infittiti
        end if
     endif
  end do
  ! FINE FASE INFITTIMENTO MESHPOINTS

  ! INIZIA FASE SFOLTIMENTO MESHPOINTS 
  N = Nini + 2
  do while(N < MAXME-1)
     N = N+1
     
      if (on_off_fine_mesh_DM==1 .and. N>Nini+N_mesh_aggiunti_DM) then!Conviene partire da 3 mesh più in la per l'infittimento se no infittisce troppo
         !al centro e quseto fa impazzire la ciacio

         Lumi_mesh=(epsi_DM(N)*(G(5,N+1)-G(5,N))*1e33)

         Lumi_mesh_prec=(epsi_DM(N-1)*(G(5,N)-G(5,N-1))*1e33)!La luminosità della DM nel mesh precedente al mesh N(Si usa il mesh precedente 
         !perche l'intervallo in massa è tra N-1 e N)

         !Se nel mesh il contributo della DM alla luminsoità positiva/negativa è maggiore di una frazione definita
         !rispetta a quella totale positiva/negativa devo infittire il mesh. Se al passo N esimo la somma della luminosità nel mesh tra N-1 N ed il mesh
         !N N+1 é maggiore del mio valore massimo di luminosità per singolo mesh non devo sfittire
     
         Lumi_mesh_sommati=Lumi_mesh+Lumi_mesh_prec

         if ( Lumi_mesh_sommati>=0 ) then
            if ( Lumi_mesh_sommati>Lumi_DM_positiva/N_min_mesh_DM ) then 
               infittisci_per_DM=1
            else
               infittisci_per_DM=0
            end if
         else
            if (abs(Lumi_mesh_sommati)>abs(Lumi_DM_negativa)/N_min_mesh_DM) then
               infittisci_per_DM=1
            else
               infittisci_per_DM=0
            endif
         end if
      else
         infittisci_per_DM=0
      endif
      
      ! ricerca intervallo di VMM da utilizzare
     GIG = G(5,N)/EMTOT 
     do I=1,4 
        if(GIG >= VMM(I) .and. GIG < VMM(I+1)) exit
     end do
     if(I > 4) I = 4 
     UL = ULA(I)*SCALA*0.8 
     UP = UPA(I)*0.8
     UT = UTA(I)*0.8 
     UR = URA(I)*0.8 
     UM = UMA(I)*0.8
     UHe3 = UHe3_ori * 0.8
     
     UXH = 0.001 
     UH4 = 0.001 
     if(XXX(1,1) <= 0. .and. XXX(3,1) > 1.d-8 .and. abs(G(6,N)) <= 0.02) &
          cycle
     if(XXX(1,N-1) <= 0. .and. XXX(1,N) > 0.) cycle
     if(XXX(1,N) <= 0. .and. XXX(1,N+1) > 0.) cycle
     GMX = G(5,N)/EMTOT 
     if(BB(1) > 0 .and. GMX >= YM1 .and. GMX <= YM2 .and. XXX(1,1) <= 0.) &
          cycle
     P1 = abs(G(1,N-1)-G(1,N+1))/G(1,N) 
     P2 = abs(G(2,N-1)-G(2,N+1))/G(2,MAXME) 
     P3 = abs(G(3,N-1)-G(3,N+1))/G(3,N) 
     P4 = abs(G(4,N-1)-G(4,N+1))/G(4,N) 
     P5 = abs(G(5,N-1)-G(5,N+1))/EMTOT

     !! Ema & Matt: sfoltimento He3 !!
     PHe3 = 0.0
     if(XHe3_tmp(N) >= XHe3_min .and. iskip_He3 <= 0) then
        PHe3 = abs( (XHe3_tmp(N-1) - XHe3_tmp(N+1)) / XHe3_tmp(N))
     endif
     
     WXH = abs(XXX(1,N-1)-XXX(1,N+1)) 
     WH4 = abs(XXX(3,N-1)-XXX(3,N+1)) 
     if(P1 > UR .or. P2 > UL .or. P3 > UP .or. P4 > UT .or. P5 > UM .or. &
          & PHe3 > UHe3 .or. infittisci_per_DM==1 ) cycle
     
     if( WXH > UXH .or. WH4 > UH4 ) cycle
     do J = 1,MELE 
        PX(J) = XXX(J,N-1)*(G(5,N)-G(5,N-1)) + XXX(J,N)*(G(5,N+1)-G(5,N)) 
     end do
     PXHe3_tmp = XHe3_tmp(N-1) * (G(5,N) - G(5,N-1)) + &
          & XHe3_tmp(N) * (G(5,N+1) - G(5,N))
     
     DM = G(5,N+1) - G(5,N-1) 
     do J = 1,MELE 
        XXX(J,N-1) = PX(J)/DM 
     end do
     XHe3_tmp(N-1) = PXHe3_tmp / dM
     
     KK = N+1 
     do K=KK,MAXME 
        do J=1,6 
           G(J,K-1) = G(J,K) 
        end do
        do J=1,MELE 
           XXX(J,K-1) = XXX(J,K) 
        end do
        XHe3_tmp(K-1) = XHe3_tmp(K)
        
        IN(K-1) = IN(K) 
        XSERV(1,K-1) = XSERV(1,K) 
        XSERV(2,K-1) = XSERV(2,K) 
     end do
     MAXME = MAXME-1
     if ( on_off_fine_mesh_DM==1) then !Se non sono interessato all'infittimento su criteri basati sulla DM non mi ricalcolo l'epsi, questo accellera il processo.
      call epsi_DM_routine(T_DM,Lumi_DM,NMD)!Calcolo nuovamente l'array dell'epsi_DM per la nuova struttura coi mesh sfittiti 
     end if    
  end do
  ! FINE FASE SFOLTIMENTO 

  !********** LEVA    MESH CHE DISTANO MENO DI 1.D-12 IN MASSA ***********
17 continue 
  MAM = MAXME-1 
  do K=1,MAM 
     if(((G(5,K+1)-G(5,K))/(G(5,K+1)+G(5,K))) > 1.d-12) cycle 
     do LL=K,MAM 
        do NN=1,6 
           G(NN,LL) = G(NN,LL+1) 
        end do
        do NN=1,MELE 
           XXX(NN,LL) = XXX(NN,LL+1) 
        end do
        XSERV(1,LL) = XSERV(1,LL+1) 
        XSERV(2,LL) = XSERV(2,LL+1) 
     end do
     MAXME = MAXME-1 
     goto 17 
  end do


  
  !****** IN CASO DI DEBUG STAMPA FISICA E CHIMICA MODIFICATA ************
  if(IPRALL == 0) return 
  write(2,100) 
  do K=1,MAXME 
     write(2,101) K,(G(J,K),J=1,6),(XXX(I,K),I=1,6) 
  end do
100 format(1X,'FISICA DOPO REZONING: K - R - L - P - T - M - DGRAD - H&
       & - HE3 - HE4 - C - N - O',/)                                       
101 format(1X,I4,1P,12E10.3) 
  return 
end subroutine OPTIM
