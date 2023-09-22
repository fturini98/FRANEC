subroutine Cattura_DM()
   use Dark_Matter!Modulo per le variabili relative alla DM

   use costanti
   !Per usare Ggrav (costante di Grav. univ. in cgs) 
    

   !PER VEDERE LE FUNZIONI DEI MODULI GUARDA IL FILE moduli.f90
  
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
   

   implicit none
   
   integer :: k, counter
   real :: xi_DM,v_esc, C_geom !fattore xi, velocità di fuga e Rate cattura geometrica

   real :: R_sup, M_sup !Raggio e Massa sulla superfice della stella
      
   real :: N_in_PASTEM !Numero di particelle catturate da la struttura
   !che ha come rate di cattura quello catturato, e come passo temporale
   !quello nuovo calcolato dalla subroutine PASTEM


   !#################
   !Catura geometrica
   !#################
   !Calcolo il rate di cattura geometrico secondo formula 3.7 di https://arxiv.org/pdf/1702.02768.pdf
   !Prima mi serve il fattore di soppressione dovuto al movimento del sole all'interno della stella xi_DM
   !che è dipendente dalla velocità di fuga superficilae
   counter=0
   do k = 1, LIM !Loop per trovare i valori della struttura sulla superficie
      if (G(1,k)/=0) then 
         R_sup=G(1,k)
         M_sup=G(5,k)
         counter=counter+1
      endif
   end do
   
   !Calcolo la velocità di fuga superficiale
   v_esc=sqrt(2*Ggrav*M_sup/R_sup*1e23)*1e-5
   
   !Calcolo il fattore corettivo dovuto al moto della stella
   call xi_DM_sub(v_esc,xi_DM)

   !Calcolo il rate di cattura geometrico
   C_geom= pigre*(R_sup*1e10)**2*(rho_DM/mass_DM)*sqrt(8.0d0/(3.0d0*pigre))*&
   (v_disp*1e5)*&!si converte solo questa velocità da km/s a cm/s perché 
   !le altre unità di misura si semplificano tra di loro
   (1+(3.d0/2.d0)*((v_esc/v_disp)**2))*xi_DM

   N_in_PASTEM=C_geom*(HT1*1e6)*sec_anno !Particelle catturate nel passo temporale in cui vale tale struttura.
   !HT1 è in milioni di anni.

   N_DM_tot=N_DM_tot+N_in_PASTEM !Particelle fino alla prossima struttura

   !######################
   !Parte inutile di test
   !######################
    
   !write(*,*)"Counter:",counter,"Massa sup:",M_sup,"Raggio_sup:",R_sup, &
   !"Velocità di fuga:", v_esc,"HT1:",HT1,"HT1V",HT1v,"Xi:",xi_DM,"C_geom:",&
   !C_geom,"N chi totali fino a qui:",N_DM_tot

   !Controllo sole
   !R_sup=6.96
   !v_esc=sqrt(2*Ggrav*M_sup/R_sup*1e23)*1e-5
   !call xi_DM_sub(v_esc,xi_DM)
   !C_geom= pigre*(R_sup*1e10)**2*(rho_DM/mass_DM)*sqrt(8.0d0/(3.0d0*pigre))*&
   !(v_disp*1e5)*&!si converte solo questa velocità da km/s a cm/s perché 
   !le altre unità di misura si semplificano tra di loro
   !(1+(3.d0/2.d0)*((v_esc/v_disp)**2))*xi_DM
   !write(*,*)"V fuga sole:",v_esc,"Xi soloe",xi_DM ,"10^5",1e5,"C_geom Sole",C_geom 

end subroutine Cattura_DM

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