subroutine Cattura_DM_geom(C_geom,v_esc,xi_DM,M_sup,R_sup)
!Modulo per calcolarsi la cattura geometrica usando la definizione in https://arxiv.org/pdf/1702.02768.pdf

   use Dark_Matter!Modulo per le variabili relative alla DM, in particolare serve per salvarsi istante per istante N_DM_tot
   !Vale a dire il numero di particelle di dark matter peresenti nella stella in quel momento

   use costanti
   !Per usare Ggrav (costante di Grav. univ. in cgs) 
  
   
   use fisica
   !Il modulo fisica serve a definire la matrice fisica G dove:
   !     G(1,k) e'      r/10^10          in cgs
   !     G(2,k) e' la   l/10^32          in cgs
   !     G(3,k) e' la  (P_totale)/10^17  in cgs
   !     G(4,k) e' la  T/10^6            in cgs     
   !     G(5,k) e' la  m/10^33           in cgs
  
   
   use chimic
   !Il modulo chimic sereve a definire la matrice della chimica XXX(LIM,MAXME)

   use chim
   !Il modulo chim serve per definire l' array XX(MAXME)(vale a dire l'abbondanza totale dell'elemento i-esimo), il passo temporale HT1 etc..
   
   use mesh
   !Serve per trovare il MAXME

   implicit none
   
   integer :: k!indice per creare i cicli

   real :: xi_DM,v_esc, C_geom !fattore xi, velocità di fuga e Rate cattura geometrica

   real :: R_sup, M_sup !Raggio e Massa sulla superfice della stella
      
   real :: N_in_PASTEM !Numero di particelle catturate dalla struttura
   !che ha come rate di cattura quello catturato, e come passo temporale
   !quello nuovo calcolato dalla subroutine PASTEM


   !#################
   !Catura geometrica
   !#################
   !Calcolo il rate di cattura geometrico secondo formula 3.7 di https://arxiv.org/pdf/1702.02768.pdf
   !Prima mi serve il fattore di soppressione dovuto al movimento del sole all'interno della stella xi_DM
   !che è dipendente dalla velocità di fuga superficilae
   
   R_sup=G(1,MAXME)
   M_sup=G(5,MAXME)
   
   !Calcolo la velocità di fuga superficiale
   v_esc=sqrt(2*Ggrav*M_sup/R_sup*1e23)*1e-5 !la velocità di fuga é in Km/s per usarla in xi_DM che richiede le velcocità in km/s
   
   !Calcolo il fattore corettivo dovuto al moto della stella
   call xi_DM_sub(v_esc,xi_DM)

   !Calcolo il rate di cattura geometrico
   C_geom= pigre*(R_sup*1e10)**2*(rho_DM/mass_DM)*sqrt(8.0d0/(3.0d0*pigre))*&
   (v_disp*1e5)*&!si converte solo questa velocità da km/s a cm/s perché 
   !le altre unità di misura si semplificano tra di loro
   (1+(3.d0/2.d0)*((v_esc/v_disp)**2))*xi_DM

   !######################
   !Parte di test
   !######################
    
end subroutine Cattura_DM_geom

subroutine xi_DM_sub(v_esc,  xi_DM)
!Mi calcolo il fattore correttivo dovuto 
!al fatto che la stella è in moto sempre
!secondo articolo https://arxiv.org/pdf/1702.02768.pdf
   use costanti
   use Dark_Matter
   use FUNZIONI !Serve per l'erf
   real, intent(in) :: v_esc
   real, intent(out) ::  xi_DM

   xi_DM=(v_disp**2 * exp((-3.0d0/2.0d0) * (v_star/v_disp)**2) + &
   sqrt(pigre/6.0d0) * v_disp/v_star * &
   & (v_disp**2 + 3.0d0 * v_esc**2 + 3.0d0 * v_star**2) * &
   erf(sqrt(3.0d0/2.0d0) * (v_star/v_disp))) / &
   (2.0d0 * v_disp**2 + 3.0d0 * v_esc**2)

   
end subroutine xi_DM_sub

subroutine Chi_DM(a,  b, chi)
!Subroutine per calcolarsi il fattore chi dell'articolo di Gould1987
   use costanti
   use FUNZIONI !Serve per l'erf

   real, intent(in) :: a,b
   real, intent(out) ::  chi

   chi=(sqrt(pigre)/2)*(erf(b)-erf(a))
   
end subroutine Chi_DM

subroutine Cattura_DM(NMD,Tempo,C_tot)
!Calcolo il rate di cattura e trascurando l'evaporazione della DM assegno a N_tot_DM (definito in Dark_Matter) il numero totale di particelle
   
   use mesh
   !Serve per MAXME(numero massimo di mesh) 

   use Dark_Matter!Modulo per le variabili relative alla DM

   use Masse_ele!Modulo per le masse degli elementi (massa_ele, u_to_gramms)

   use costanti
   !Per usare Ggrav (costante di Grav. univ. in cgs) 
  
   !Il modulo fisica serve a definire la matrice fisica G dove:
   !     G(1,k) e'      r/10^10          in cgs
   !     G(2,k) e' la   l/10^32          in cgs
   !     G(3,k) e' la  (P_totale)/10^17  in cgs
   !     G(4,k) e' la  T/10^6            in cgs     
   !     G(5,k) e' la  m/10^33           in cgs
   use fisica
  
   !Il modulo chimic sereve a definire la matrice della chimica XXX
   use chimic
  
   use chim
   !Il modulo chim serve per definire l' array XX, il passo temporale HT1 etc..

   use Asplund_per_DM
   !Serve a caricare gli elementi e i valori relativi di log10(N_X/N_H)+12 di asplund 2009
   

   implicit none

   integer :: NMD !Il numero del modello corrente
   real :: Tempo !L'età del modello

   real:: C_tot,& !Rate di cattura totale
         C_geom,& !Rate di cattura geometrico
         v_esc,& !Velocità di fuga
         xi_DM,& !fattore correttivo superficiale dovuto al movimento della stella (dipende anche dalla v_fuga)
         xi_DM_mesh,&!fattore correttivo dovuto al movimento della stella mesh per mesh
         C_weak,& !Rate di cattura dovuto all'interazione debole
         N_in_PASTEM,& !Particelle catturate nell'intervallo di tempo HT1 considerando rate di cattura nell'intervallo C_tot
         sigma_ele,& !Sezione d'urto DM elemento
         M_sup,&!Massa totale della stella
         R_sup,&!Raggio stella
         fattore_mediato,&!Fattore mediato con dentro il potenziale grav. ridotto
         fattore_mediato_mesh,&!Fattore mediato con dentro il potenziale grav. ridotto (valore al mesh)
         v_esc_mesh,& !Velocità di fuga mesh per mesh
         A_mesh_ele_quad,&!Fattore A dlell'articolo di Gould
         A_p_mesh,A_m_mesh,&!Fattori A_+ e A_- dell'articolo di Gould
         eta,& !Fattore eta di Gould
         chi_A,chi_eta,& !Fattori Chi della formula di Gould
         mu, mu_men!Fattori mu dell'articolo di Gould

   real, dimension(MAXME):: phi!Array del potenziale gravitazionale ridotto
   integer :: ele,i!Indici per elementi e mesh

   real, dimension(MELE):: C_weak_ele !Array per i rate di cattura per i differenti elementi
   
   !###################################################################
   !Variabili relative a introduzione di elementi aggiuntivi da asplund
   !###################################################################
   !Considero che tutti gli elementi esclusi quelli già presenti e che non rientrano nelle reazioni nucleari diffondano come il ferro
   real :: C_weak_asplund,&!Cattura considerando rapporti tra Fe e X_ele come quelli in  asplund
            C_weak_ele_asplund,&!Rate di cattura dovuto all'interazione con l'elemento ele-esimo
            Abb_ele,&!Abbondanza elemento su tutta la stella calcolata da tabella aspludn partendo da abbondanza del ferro, considerando che il rapporto 
            !abbondanza tra ferro e elemento rimane costante
            Abb_ele_mesh,&!Stessa cosa della variabile precedente ma relativa al mesh
            C_tot_asplund!Cattura totale considerando tutti gli elementi in XXX e in asplund 2009


   call Cattura_DM_geom(C_geom,v_esc,xi_DM,M_sup,R_sup)

   eta=sqrt(3.0/2.0)*(v_star/v_disp)

   C_weak=0!Setto la cattura debole a 0 per sommarci sopra per ogni elemento


   !Calcolo il potenziale gravitazionale ridotto
   phi(MAXME)=1.00
   do  i= (MAXME-1),1 , -1
      ! ciclo dall'esterno all'interno per calcolarsi phi secondo 
      !https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1991ApJ...368..626D&defaultprint=YES&filetype=.pdf
      !Per denistà costante è come l'integrale di Gaus per una sfera uniformemente carica
      if ( i==1 ) then
         phi(i)=phi(i+1)+R_sup/M_sup*((G(5,2)+G(5,1))/2)/(((G(1,2)+G(1,1))/2)**2)*G(1,2)!Questo serve a evitare il diviso 0 al centro, uso come valori di fisica della shell quelli dell'intershell
      else
         phi(i)=phi(i+1)+R_sup/M_sup*(G(5,i)/(G(1,i)**2))*(G(1,i+1)-G(1,i))
      end if
      
   end do

   !Calcolo il rate di cattura debole
   do ele=1,MELE!Ciclo su gli elementi chimici

      !Sezione d'urto per l'elemento ele-esimo
      sigma_ele=sigma0_DM*(massa_ele(ele)/massa_ele(1))**4 *&
      (((mass_DM*GeV_grammi+massa_ele(1)*u_to_gramms)/(mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms))**2)!Sezione d'urto come  doi:10.1103/physrevd.66.053005

      !Mi calcolo mu per l'elemento ele-esimo
      mu=mass_DM*GeV_grammi/(massa_ele(ele)*u_to_gramms)
      mu_men=(mu-1)/2

      fattore_mediato=0!Setto a 0 il fattore mediato per sommarci sopra per ogni mesh dato che è un integrale in delta M
      do i = 1, MAXME

         if(xx(ele)/=0)then!Serve a evitare il divide by 0 quando xx(ele)==0
            v_esc_mesh=sqrt(phi(i))*v_esc !Mi calcolo la velocità di fuga mesh per mesh partendo dalla definizione di phi e invertendola. phi=v_esc^2/v_esc(R_sup)^2 é in km/s
            A_mesh_ele_quad=3/2*(v_esc_mesh/v_disp)**2*(mu/(mu_men**2))!Mi calcolo A^2 di Gould
            A_p_mesh=sqrt(A_mesh_ele_quad)+eta!Mi calcolo A_+
            A_m_mesh=sqrt(A_mesh_ele_quad)-eta!Mi calcolo A_-
            call Chi_DM(-eta,eta,chi_eta)
            call Chi_DM(A_m_mesh,A_p_mesh,chi_A)
         
            fattore_mediato_mesh=(phi(i)*&
                           (1./(2*eta*A_mesh_ele_quad))*&
                           ((A_p_mesh*A_m_mesh-1./2)*(chi_eta-chi_A)+&
                           1./2*A_p_mesh*exp(-(A_m_mesh)**2)-1./2*A_m_mesh*exp(-(A_p_mesh)**2)-eta*exp(-(eta)**2)))*&!devi sostituire l'1 con i vari elementi contenenti A.
                           (G(5,i+1)-G(5,i))*XXX(ele,i)&!Delat m per l'elemento ele-esimo
                           /(M_sup*xx(ele))!Normalizazione per tornare adimensionale della massa dell'elemento ele-esimo
            fattore_mediato=fattore_mediato+fattore_mediato_mesh
         else
            fattore_mediato=0
         endif
      end do


      C_weak_ele(ele)=sqrt(8/(3*pigre))*sigma_ele*&!Sto usando la definizione in https://ui.adsabs.harvard.edu/abs/1987ApJ...321..571G/abstract
                 rho_DM/mass_DM*&
                 v_disp*1e5*&!Lo metto in cm/s
                 xx(ele)*(M_sup*1e33)/(massa_ele(ele)*u_to_gramms)*&
                 3/2*(v_esc/v_disp)**2*&
                 fattore_mediato
      !write(*,*)"sigma_ele",sigma_ele,"fattore",fattore_mediato,"C_weak_ele(ele)",C_weak_ele(ele),"Abb",xx(ele)
      C_weak=C_weak+C_weak_ele(ele)!Somma dei rate per tutti gli elementi
   end do

   C_tot=C_weak*(1-exp(-(C_geom/C_weak)))!Tengo conto della saturazione del rate di cattura


   !####################################
   !Sezione aggiunta elementi da Asplund
   !####################################
   !Dato che franec non segue tutti gli elemnti, ma solo quelli più leggeri, considerando come buona approssimazione che gli elementi con
   !A>=31(Fosforo) diffondano come il ferro,e dato che non concorrono alle reazioni nulceari, mesh per mesh afremo X_fe/X_ele=cost. Il valore
   !del loro rapporto lo ricavo da Asplund
   if ( on_off_DM_Asplund==1) then!Posso attivare o disattivare questo meccanismo cambiano la flag che è definita in Asplund per DM
      C_weak_asplund=0!Setto il rate di cattura a 0 per ogni elemento così che da potermi calcolare quello totale come somma

      do ele=1,69!Ciclo su gli elementi chimici in asplund

         if ( numeri_massa_asplund_DM(ele)>=31 .and. numeri_massa_asplund_DM(ele)/=56) then!Cosi consideriamo che gli elementi in asplund dal P(A=31) in poi diffondono come Fe, in alcune
            !simulazioni si può mettere al posto di 31(A) 19, che corrisponde al F, ma gli elementi tra 19 e 31 sono già nella chimica e ogni tanto disattivati
            
            !Calcolo l'abbondanza totale dell'elemento considerando che il rapporto delle abbondanze Fe e elemento rimane costante
            Abb_ele=xx(21)*numeri_massa_asplund_DM(ele)/(56.)*&!Il valore 21 nell'array delle abbondaze corrisponde al Fe
                     10**(dati_asplund_DM(ele)-12.-log10(xx(21)/(xx(1)*56)))

            !Sezione d'urto per l'elemento ele-esimo
            sigma_ele=sigma0_DM*(numeri_massa_asplund_DM(ele)/(numeri_massa_asplund_DM(1)))**4 *&
            (((mass_DM*GeV_grammi+numeri_massa_asplund_DM(1)*p_mass)/(mass_DM*GeV_grammi+numeri_massa_asplund_DM(ele)*p_mass))**2)!Sezione d'urto come  doi:10.1103/physrevd.66.053005
   
            mu=mass_DM*GeV_grammi/(numeri_massa_asplund_DM(ele)*p_mass)
            mu_men=(mu-1)/2
            fattore_mediato=0
            do i = 1, MAXME         
               if(Abb_ele/=0)then!Serve a evitare il divide by 0 quando Abb_ele==0

                  Abb_ele_mesh=xxx(21,i)**numeri_massa_asplund_DM(ele)/(56.)*&!Il valore 21 nell'array delle abbondaze corrisponde al Fe
                     10**(dati_asplund_DM(ele)-12.-log10(xxx(21,i)/(xxx(1,i)*56)))

                  v_esc_mesh=sqrt(phi(i))*v_esc !Mi calcolo la velocità di fuga mesh per mesh partendo dalla definizione di phi e invertendola.
                  A_mesh_ele_quad=3/2*(v_esc_mesh/v_disp)**2*(mu/(mu_men**2))
                  A_p_mesh=sqrt(A_mesh_ele_quad)+eta!Mi calcolo A_+
                  A_m_mesh=sqrt(A_mesh_ele_quad)-eta!Mi calcolo A_-
                  call Chi_DM(-eta,eta,chi_eta)
                  call Chi_DM(A_m_mesh,A_p_mesh,chi_A)

                  fattore_mediato_mesh=(phi(i)*&
                           (1./(2*eta*A_mesh_ele_quad))*&
                           ((A_p_mesh*A_m_mesh-1./2)*(chi_eta-chi_A)+&
                           1./2*A_p_mesh*exp(-(A_m_mesh)**2)-1./2*A_m_mesh*exp(-(A_p_mesh)**2)-eta*exp(-(eta)**2)))*&!devi sostituire l'1 con i vari elementi contenenti A.
                           (G(5,i+1)-G(5,i))*Abb_ele_mesh&!Delat m per l'elemento ele-esimo
                           /(M_sup*Abb_ele)!Normalizazione per tornare adimensionale della massa dell'elemento ele-esimo
                  
                  fattore_mediato=fattore_mediato+fattore_mediato_mesh
               else
                  fattore_mediato=0
               endif
            end do
            
            C_weak_ele_asplund=sqrt(8/(3*pigre))*sigma_ele*&!Sto usando la definizione in https://ui.adsabs.harvard.edu/abs/1987ApJ...321..571G/abstract
                    rho_DM/mass_DM*&
                    v_disp*1e5*&
                    Abb_ele*(M_sup*1e33)/(numeri_massa_asplund_DM(ele)*p_mass)*&
                    3/2*(v_esc/v_disp)**2*&
                    fattore_mediato
            !write(*,*)"sigma_ele",sigma_ele,"fattore",fattore_mediato,"C_weak_ele_asplund",C_weak_ele_asplund,"elemento",char(9),elementi_asplund(ele),"Abb",Abb_ele,"A",numeri_massa_asplund_DM(ele),ele
            C_weak_asplund=C_weak_asplund+C_weak_ele_asplund!Somma dei rate per tutti gli elementi
         end if
      end do
   
      C_tot_asplund=(C_weak_asplund+C_weak)*(1-exp(-(C_geom/(C_weak+C_weak_asplund))))!Qui sommo i rate di cattura dovuti agli elementi in Asplund e quelli in XXX, quindi attenuo per la saturazione
      !write(*,*)"C_weak_asplund/C_weak",C_weak_asplund,"C_tot",C_tot,"C_tot con asplund",C_tot_asplund
      C_tot=C_tot_asplund
   end if
   !################################
   !Fine sezione elementi da Asplund
   !################################
   
   

   N_in_PASTEM=C_tot*(HT1*1e6)*sec_anno !Particelle catturate nel passo temporale in cui vale tale struttura.
   !HT1 è in milioni di anni.

   N_DM_tot=N_DM_tot+N_in_PASTEM !Particelle fino alla prossima struttura

   if ( ioDarkCattura_on_off==1 ) then !Abilito il file di scrittura relativo alla cattura della DM
   

      if (NMD==3)then!Solo al primo giro scrivo cosa sono le varie colonne nel file di output della cattura
         write(ioDarkCattura,'(A)')"NMD"//char(9)//char(9)//char(9)//char(9)//&
                                 "Tempo"//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_weak"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_geom"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_tot"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "HT1"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "N_in_pastem"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "N_DM_tot"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_H"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_HE3"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_HE4"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_C12"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_N14"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_O16"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_O18"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_NE20"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_NE22"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_MG24"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_MG25"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_MG26"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_SI28"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_NEU"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_NA23"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_C13"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_N15"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_O17"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_F19"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_NE21"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_Fe"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_Li6"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_Li7"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_Be"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_B"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                 "C_D"

      endif
      
         write(ioDarkCattura,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A)', advance='no')nmd,char(9),Tempo,char(9),C_weak,char(9),C_geom,char(9),C_tot,char(9),HT1,char(9),&
                                                                  N_in_PASTEM,char(9),N_DM_tot,char(9)!Scrivo le variabili della struttura totali
      do ele = 1, MELE!Scrivo i rate di cattura per ogni elemento
         if ( ele==MELE ) then!All'ultimo elemento avanzo di una riga così da registrare il passo temporale successivo
            write(ioDarkCattura,'(ES25.15)',advance='yes')C_weak_ele(ele)
         else
            write(ioDarkCattura,'(ES25.15, A)',advance='no')C_weak_ele(ele),char(9)
         end if
      
      end do
   
   else if (NMD==3) then!NMD 3 è la prima volta in cui abilito la cattura, questo serve per far scrivere una sola volta questo messaggio
      write(ioDarkCattura,*)"LA SCRITTURA DEL FILE  DarkMatterCattura.DAT E' DISABILITATA!"
      write(ioDarkCattura,*)"Per abilitarla impostare ioDarkCattura_on_off==1"
   end if 
end subroutine Cattura_DM
