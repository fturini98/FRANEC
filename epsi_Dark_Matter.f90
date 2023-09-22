subroutine epsi_DM_routine(T_DM,epsi_DM_tot)

    use Dark_Matter!Modulo per le variabili relative alla DM

    use Masse_ele!Modulo per le masse degli elementi (massa_ele, u_to_gramms)

    use costanti
    !Per usare Ggrav (costante di Grav. univ. in cgs)
    !e la costante di boltzman in cgs 
    

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
    
    implicit none

    integer :: i, & !indice
               ele  !indice per elementi


    real,dimension(LIM) :: T_mesh,& !Variabile per salvare la temperatura della stella
                           M_mesh,& !Variabile di supporto per l'array dei raggi
                           exp_dist_DM !Fattore esponenziale della distirbuzione della DM

    real :: T_DM,& !Temperatura della DM per cui si calcola l'epsi_DM
            sigma_ele,& !Sezione d'urto per elemento
            norm_DM_distr !Costante di normalizzazione della distribuzione di densità

    real :: epsi_DM_tot !Integrale sull'epsi
    !##########################
    !Inizializzazione variabili
    !##########################


    !Ciclo per caricare l'array della temperatura dalla matrice G
    !a dipendenza dell'implementazione del metodo di Heneyey può essere tolta
    !Serve anche a calcolarsi i fattori della distribuzione della DM
    do i = 1,LIM 
        
        T_mesh(i)=G(4,i)
        M_mesh(i)=G(5,i)!Creo una variablie di suporto per la massa, a dipendenza dell'implementazione del metodo di Heneyey può essere tolta.

        !#################################################
        !Calcolo exp(-m_DM*V(r)/K_b T_DM) e la sua normal.
        !per calcolarmi l'array della densità della DM
        !#################################################

        if ( i==1 ) then !Questo if serve ad evitare problemi del divided by 0 del raggio nella shell più interna
            
            exp_dist_DM(i)=exp(-(mass_DM*GeV_grammi*&
            (Ggrav*(G(5,i)+G(5,i+1))/(G(1,i)+G(1,i+1))*1e23))/&!Ho fatto la media tra il mesh i e il mesh i+1 per le variabili fisiche
            (Kboltz*T_DM*1e6))

            norm_DM_distr=((G(1,i)+G(1,i+1))/2)**2*(G(1,i+1)-G(1,i))*exp_dist_DM(i)*1e30
        end if
        
        if (G(1,i)/=0) then

            exp_dist_DM(i)=exp(-(mass_DM*GeV_grammi*&
            (Ggrav*(G(5,i))/(G(1,i))*1e23))/&
            (Kboltz*T_DM*1e6))

            norm_DM_distr=((G(1,i))**2*(G(1,i+1)-G(1,i))*exp_dist_DM(i)*1e30)+norm_DM_distr

        else!Questo serve a prendere in considerazione quando il raggio è messo
            ! a 0 perché quell'elemento dell'array non ha senso fisico

            exp_dist_DM(i)=0!Setto a 0 la distribuzione di DM perché qui sono fuori dalla stella e non ho definito il potenziale grav.
        end if

    end do

    !#################################
    !Calcolo epsi con loop su elementi
    !#################################

    epsi_DM_tot=0!Inizializzo la funzione cumulartiva della epsi a 0

    do i = 1, LIM!Loop sui mesh
        epsi_DM(i)=0 !Lo setto a zero per poi sommarci su ad ogni elemento
        do ele = 1, MELE !Loop su tutti gli elementi
            sigma_ele=sigma0_DM*(massa_ele(ele)/massa_ele(1))**4 *&
            (((mass_DM*GeV_grammi+massa_ele(1)*u_to_gramms)/(mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms))**2)!Sezione d'urto come  doi:10.1103/physrevd.66.053005

            !Qua calcola l'epsi con la T_DM inserita. La calcola erg/grammi, quindi al posto della densità numerica dell'i-esimo elemento basta che usi la sua abbondanza
            epsi_DM(i)=epsi_DM(i)+(8*sqrt(2/pigre)*&
            N_DM_tot*exp_dist_DM(i)/norm_DM_distr*& !Densità numerica della DM mesh per mesh
            XXX(ele,i)/(massa_ele(ele)*u_to_gramms)*sigma_ele*& !Densità dell'elemento chimico * sezione d'urto
            (mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms)/((mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms)**2)*&!Fattore differenza massa
            sqrt(Kboltz*1e6*(massa_ele(ele)*u_to_gramms*T_DM+mass_DM*GeV_grammi*T_mesh(i))/(mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms))*& !Fattore somma temperature
            Kboltz*1e6*(T_mesh(i)-T_DM))!Differenza temperature

        end do

        !Mi calcolo anche l'integrale su epsi
        if ( epsi_DM(i)/=0 ) then !Questa condizione dovrebbe accellerare il calcolo e in più mi garantisce di evitare il problema del mesh i=LIM in cui R(i+1) non è
            !definito. Infatti epsi_DM(LIM)=0 per definizione in quanto le abbondanze in LIM sono 0 e epsi dipende da loro in maniera lineare.
            epsi_DM_tot=epsi_DM_tot+(epsi_DM(i)*(M_mesh(i+1)-M_mesh(i))*1e33)!Mi calcolo la lumi totale tramite integrazione in massa
        end if
    end do

    
    !#########
    !Zona Test
    !#########

end subroutine epsi_DM_routine

subroutine convergenza_epsi_DM()
    use Dark_Matter

    !Il modulo fisica serve a definire la matrice fisica G dove:
    !     G(1,k) e'      r/10^10          in cgs
    !     G(2,k) e' la   l/10^32          in cgs
    !     G(3,k) e' la  (P_totale)/10^17  in cgs
    !     G(4,k) e' la  T/10^6            in cgs     
    !     G(5,k) e' la  m/10^33           in cgs
    use fisica

    use costanti !Serve e per il pi greco

    implicit none

    integer :: i,&!indice
                estremo,& !Mi serve per memorizzare se l'estremo precedente dell'intervallo di temperatura era a destra o sinistra
                !estremo=1 Destra, quindi prima avevo calcolato per T_max
                !estremo=0 Sinistra, quindi prima avevo calcolato per T_min
                max_cicli=60!Numero massimo di cicli per la convergenza della T_DM.
    real :: epsi_DM_tot,& !Epsi cumulativa totale secondo la formula di Spergel and Press
            T_DM,& !Temperatura della Dark Matter
            epsi_DM_tot_vecchia,&!Variabile di supporto per salvarmi l'epsi tot durante la convergenza
            epsi_max_val !Il valore limite per cui si ritiene che l'errore che si commette è trascurabile

    real :: T_max=0,& !Variabile per registrare la temperatura massima della stella
            T_min !VAriabile per registrare la temperatura minima della stella
    !Queste due variabili permettono di avere un epsi cumulativa >0 o <0, il che garantisce
    !il funzionamento del metodo di bisezione per la convergenza.

    real,dimension(LIM) :: T_mesh

    !#########################
    !Inizializazione variabili
    !#########################

    !Ciclo per inizializzare:
    !Temperatura massima e minima della stella 
    do  i= 1, LIM 
        T_mesh(i)=G(4,i)
        !Al primo giro assegno a T_min la temperatura centrale, per dare un punto di partenza all'algoritmo.
        if ( i==1 ) then 
            T_min=T_mesh(i)
        end if

        !Mi trovo temperatura massima e minima all'interno della stella
        if ( T_mesh(i)>T_max ) then
            T_max=T_mesh(i)
        else if (T_mesh(i)<T_min .and. T_mesh(i)/= 0) then
            T_min=T_mesh(i)
        end if

    end do

    !Assegno come primo valore di test per la temperatura della DM la tmperatura massima della stella(che non per forza è quella centrale)
    T_DM=T_min

    !write(*,*)"La T_max è:",T_max,"La T_min è:",T_min

    !##############################
    !Calcolo dell'integrale di epsi
    !##############################
    
    !Mi calcolo la epsi_DM mesh per mesh a partire T_DM e mi calcolo il suo integrale sulla struttura
    call epsi_DM_routine(T_DM,epsi_DM_tot)
    estremo=0
    !write(*,*)"La T_DM è",T_DM,"e la epsi tot vale",epsi_DM_tot

    epsi_max_val=1e18!Setto il valore limite per l'epsi totale in erg/S

    !Uso metodo di bisezione per trovare la T_DM
    do i=1,max_cicli
        
        if(abs(epsi_DM_tot)<=epsi_max_val) then
            if ( i==1 ) then
                write(*,*)"La temperatura della DM è uguale a quella massima della stella e vale:",T_DM
            end if        
            exit
        end if

        epsi_DM_tot_vecchia=epsi_DM_tot !Salvo la epsi_tot precedente per compararla con quella nuova.
        
        T_DM=(T_max+T_min)/2!Mi calcolo la nuova temperatura della DM

        call epsi_DM_routine(T_DM,epsi_DM_tot)

        if ((epsi_DM_tot >= 0.0d0) == (epsi_DM_tot_vecchia >= 0.0d0))then
            if ( estremo==1 ) then
                T_max=T_DM
            else
                T_min=T_DM
            end if
        
        else
            if ( estremo==1 ) then
                T_min=T_DM
                estremo=0
            else
                T_max=T_DM
                estremo=1
            end if
        end if
        
        !write(*,*)"Ciclo",i,"integrale dell'epsi è:",epsi_DM_tot,"e la T_DM",T_DM

        if (abs(epsi_DM_tot)>epsi_max_val .and. i==max_cicli) then
            write(*,*)"La temperatura della DM non converge"
            exit
        end if
        
    end do
    !write(*,*)"La temperatura della DM vale:",T_DM,"L'epsi tot:",epsi_DM_tot
end subroutine convergenza_epsi_DM