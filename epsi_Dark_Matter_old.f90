subroutine epsi_DM_routine(T_DM,Lumi_DM,Delta_lumi_DM)
!Subroutine che si calcola la epsi DM in maniera vettoriale (slavandola nella variabile globale data dal modulo) 
!e restituisce Lumi_DM
    use Dark_Matter!Modulo per le variabili relative alla DM

    use Masse_ele!Modulo per le masse degli elementi (massa_ele, u_to_gramms)

    use costanti!Per usare Ggrav (costante di Grav. univ. in cgs)
                !e la costante di boltzman in cgs 
  
    
    use fisica!Il modulo fisica serve a definire la matrice fisica G dove:
        !     G(1,k) e'      r/10^10          in cgs
        !     G(2,k) e' la   l/10^32          in cgs
        !     G(3,k) e' la  (P_totale)/10^17  in cgs
        !     G(4,k) e' la  T/10^6            in cgs     
        !     G(5,k) e' la  m/10^33           in cgs
  
    use chimic!Il modulo chimic sereve a definire la matrice della chimica XXX

    use mesh !per usare MAXME
    
    use FUNZIONI !Serve per usare la funzione che calcola il logaritmo dell'integrale tramite il metodo motecarlo

    implicit none

    integer :: i,j, & !indice
               ele !indice per elementi
    
    real :: M_sup,&!Massa della stella
            R_sup!Raggio stella


    real,dimension(LIM) :: T_mesh,& !Variabile per salvare la temperatura della stella
                           M_mesh !Variabile di supporto per l'array dei raggi
                           
    real, dimension(MAXME) :: phi ,& !Array del potenziale gravitazionale ridotto
                            exp_dist_DM,& !Fattore esponenziale della distirbuzione della DM
                            esponenti

    real :: T_DM,& !Temperatura della DM per cui si calcola l'epsi_DM
            sigma_ele,& !Sezione d'urto per elemento
            fattore_norma_log,&
            norm_DM_distr,& !Costante di normalizzazione della distribuzione di densità
            epsi_mesh,&!L'epsi calcolata nel mesh i con la temoperatura T(i)
            epsi_mesh_p!L'epsi calcolata nel mesh i con la temperatura T(i+1)

    real :: Lumi_DM,& !Integrale sull'epsi(Luminosità DM)
            Lumi_mesh,& !La luminosità del mesh
            Delta_lumi_mesh,& !L'errore che si ottiene nella stima della luminosità della DM nel mesh
            Delta_lumi_DM !L'errore che si ottiene nella stima della luminosità totale della DM

    integer :: N_punti=10000 !Punti da generare per l'integrale col metodo montecarlo del fattore di normalizazione

    real, dimension(:), allocatable :: test_array
    allocate(test_array(N_punti))
    !##########################
    !Inizializzazione variabili
    !##########################
    !Mi trovo raggio e massa della stella
    R_sup=G(1,MAXME)
    M_sup=G(5,MAXME)


    !Calcolo il potenziale gravitazionale Ridotto        
    phi(MAXME)=1.00
    do  i= (MAXME-1),1 , -1! ciclo dall'esterno all'interno per calcolarsi phi secondo 
    !https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1991ApJ...368..626D&defaultprint=YES&filetype=.pdf
    !Per denistà costante è come l'integrale di Gaus per una sfera uniformemente carica
        if ( i==1 ) then
            phi(i)=phi(i+1)+R_sup/M_sup*((G(5,2)+G(5,1))/2)/(((G(1,2)+G(1,1))/2)**2)*G(1,2)!Questo serve a evitare il diviso 0 al centro, uso come valori di fisica della shell quelli dell'intershell
        else
            phi(i)=phi(i+1)+R_sup/M_sup*(G(5,i)/(G(1,i)**2))*(G(1,i+1)-G(1,i))
        end if
        
        esponenti(i)=((mass_DM*GeV_grammi*&!Nell'esponenziale ci va il + perché il potenziale gravitazionale è definito 
        !con il meno mentre io calcolo il suo modulo, quindi -*-=+
        ((Ggrav*M_sup/R_sup)*phi(i)*1e23))/&
        (Kboltz*T_DM*1e6))

    end do
    
    fattore_norma_log=log_int_exp(N_punti,esponenti)

    do i = 1,MAXME
    !Ciclo per caricare l'array della temperatura dalla matrice G
    !a dipendenza dell'implementazione del metodo di Heneyey può essere tolta
    !Serve anche a calcolarsi i fattori della distribuzione della DM
        
        T_mesh(i)=G(4,i)
        M_mesh(i)=G(5,i)!Creo una variablie di suporto per la massa, a dipendenza dell'implementazione del metodo di Heneyey può essere tolta.

        !########################
        !Questo serve a gestire il dleta M al MAXME che ha bisogno di Maxme+1
        !Se non si mette non converge perche L_DM(T_max) ha lo stesso segno di L_DM(T_min) per alcuni modelli
        !########################
        if (i==MAXME) then
            T_mesh(i+1)=G(4,i+1)
            M_mesh(i+1)=G(5,i+1)
        end if

        !#################################################
        !Calcolo exp(-m_DM*V(r)/K_b T_DM) e la sua normal.
        !per calcolarmi l'array della densità della DM
        !#################################################
        exp_dist_DM(i)=exp(esponenti(i)-fattore_norma_log)


        norm_DM_distr=((4*pigre*G(1,i))**2.*(G(1,i+1)-G(1,i))*exp_dist_DM(i)*1e30)+norm_DM_distr

        if (ABS(exp_dist_DM(i)) >= HUGE(exp_dist_DM(i))) then !Errore di Overflow nel calcolo della distribuzione radiale della DM
            write(ioDarkError,'(A)')"################################################################################################"
            write(ioDarkError,'(A)')"############################SIMULAZIONE BLOCCATA PER OVERFLOW###################################"
            write(ioDarkError,'(A)')"################################################################################################"
            write(ioDarkError,'(A)')"Impossibile calcolare la distribuzione Radiale della DM per overflow in xp(-m_DM*V(r)/K_b T_DM)"
            write(ioDarkError,*)"Il log del fattore di normalizazione è:",fattore_norma_log
            write(ioDarkError,'(ES25.15)')exp_dist_DM(1:Maxme)
            call system("killall -9 franec Lancia")
        endif
    end do

    !test_array= prob_r_quad(N_punti,R_sup)
    !write(ioDarkCattura,'(ES25.15)')test_array(1:N_punti)
    
    

    !write(ioDarkCattura,*)norm_DM_distr,fattore_norma_log
    !#################################
    !Calcolo epsi con loop su elementi
    !#################################

    Lumi_DM=0!Inizializzo la funzione cumulartiva della epsi a 0
    Delta_lumi_DM=0 !Inizializzo l'errore sulla stima della luminosita della DM a 0S
    Lumi_DM_positiva=0!Inizializo la funzione cumulativa della Luminosità positiva a 0
    Lumi_DM_negativa=0!Inizializo la funzione cumulativa della Luminosità negativa a 0
    Max_Lumi_DM=0!Inizializzo a o il massimo del modulo della luminosità a 0.
    
    do i = 1, MAXME!Loop sui mesh
        epsi_DM(i)=0 !Lo setto a zero per poi sommarci su ad ogni elemento
        Delta_epsi_DM(i)=0
        do ele = 1, MELE !Loop su tutti gli elementi
            sigma_ele=sigma0_DM*(massa_ele(ele)/massa_ele(1))**4 *&
            (((mass_DM*GeV_grammi+massa_ele(1)*u_to_gramms)/(mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms))**2)!Sezione d'urto come  doi:10.1103/physrevd.66.053005

            epsi_mesh=(8.*sqrt(2./pigre)*& !Mi calcolo l'epis nel mesh per le condizioni al bordo i
            N_DM_tot*exp_dist_DM(i)/norm_DM_distr*& !Densità numerica della DM mesh per mesh
            XXX(ele,i)/(massa_ele(ele)*u_to_gramms)*sigma_ele*& !Densità dell'elemento chimico * sezione d'urto
            (mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms)/((mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms)**2)*&!Fattore differenza massa
            sqrt(Kboltz*1e6*(massa_ele(ele)*u_to_gramms*T_DM+mass_DM*GeV_grammi*T_mesh(i))/(mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms))*& !Fattore somma temperature
            Kboltz*1e6*(T_DM-T_mesh(i)))!Differenza temperature


            epsi_mesh_p=(8.*sqrt(2./pigre)*& !Mi calcolo l'epsi nel mesh perl le condizioni al bordo (solo T)i+1 con la chimica dell'intermesh i
            N_DM_tot*exp_dist_DM(i)/norm_DM_distr*& !Densità numerica della DM mesh per mesh
            XXX(ele,i)/(massa_ele(ele)*u_to_gramms)*sigma_ele*& !Densità dell'elemento chimico * sezione d'urto
            (mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms)/((mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms)**2)*&!Fattore differenza massa
            sqrt(Kboltz*1e6*(massa_ele(ele)*u_to_gramms*T_DM+mass_DM*GeV_grammi*T_mesh(i+1))/(mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms))*& !Fattore somma temperature
            Kboltz*1e6*(T_DM-T_mesh(i+1)))!Differenza temperature

            !Qua calcola l'epsi con la T_DM inserita. La calcola erg/grammi, quindi al posto della densità numerica dell'i-esimo elemento basta che usi la sua abbondanza
            !Si considera l'epsi come la media aritmetica tra quella calcolata al bordo i e quella al bordo i+1
            epsi_DM(i)=epsi_DM(i)+(1./2.)*(epsi_mesh+epsi_mesh_p)

            !Si considera come errore nella stima dell'epsi la semi differenza tra l'epsi calcolata con la T al bordo i e quella al bordo i+1
            Delta_epsi_DM(i)=Delta_epsi_DM(i)+(1./2.)*(epsi_mesh-epsi_mesh_p)
           
        end do

        !Mi calcolo anche l'integrale su epsi e l'errore sull'integrale
        if ( epsi_DM(i)/=0 .or. i==MAXME) then !Questa condizione dovrebbe accellerare il calcolo e in più mi garantisce di evitare il problema del mesh i=LIM in cui R(i+1) non è
            !definito. Infatti epsi_DM(LIM)=0 per definizione in quanto le abbondanze in LIM sono 0 e epsi dipende da loro in maniera lineare.
            Lumi_mesh=(epsi_DM(i)*(M_mesh(i+1)-M_mesh(i))*1e33)!La luminosità della DM nel mesh
            Delta_lumi_mesh=(Delta_epsi_DM(i)*(M_mesh(i+1)-M_mesh(i))*1e33)

            Delta_lumi_DM=Delta_lumi_DM+(Delta_lumi_mesh)**2 !Errore in quadratura sulla luminosità della DM
            Lumi_DM=Lumi_DM+Lumi_mesh!Mi calcolo la lumi totale tramite integrazione in massa
            
            if ( Lumi_mesh>=0 ) then!Ogni volta che la luminosità è positiva/negativa la addizziono alla funzione cumulativa corrispondente (Questo serve per l'optim e l'infittimento dei mesh)
                Lumi_DM_positiva=Lumi_DM_positiva+Lumi_mesh
            else
                Lumi_DM_negativa=Lumi_DM_negativa+Lumi_mesh
            end if

            if ( abs(Lumi_DM)>Max_Lumi_DM ) then!Cerco il valore massimo del modulo della luminosità dovuta alla DM
                Max_Lumi_DM=abs(Lumi_DM)       
            end if
            
        end if
    end do

    Delta_lumi_DM=sqrt(Delta_lumi_DM)!Faccio la radice quadrta della luminosità della DM per avere il Delta e non il delta quad, visto che prima ho sommato in quadratura
end subroutine epsi_DM_routine


subroutine convergenza_epsi_DM(Lumi_DM,Delta_lumi_DM,T_DM,NMD,tempo)
    
    use Dark_Matter!Carica le variabili della DM

    use fisica!Il modulo fisica serve a definire la matrice fisica G dove:
    !     G(1,k) e'      r/10^10          in cgs
    !     G(2,k) e' la   l/10^32          in cgs
    !     G(3,k) e' la  (P_totale)/10^17  in cgs
    !     G(4,k) e' la  T/10^6            in cgs     
    !     G(5,k) e' la  m/10^33           in cgs

    use costanti !Serve e per il pi greco

    use mesh !per mettere il massimo del mesh

    use tempe !Serve per ellog(luminosità totale log(stella/L_sun))

    implicit none

    integer :: i,j,&!indice
                estremo,& !Mi serve per memorizzare se l'estremo precedente dell'intervallo di temperatura era a destra o sinistra
                !estremo=1 Destra, quindi prima avevo calcolato per T_max
                !estremo=0 Sinistra, quindi prima avevo calcolato per T_min
                max_cicli=21,&!Numero massimo di cicli per la convergenza della T_DM.(Con venti ho una precisone di circa 1e-6 in T_DM)
                NMD !Servwe per quando ho l'errore di convergenza per stamaprimi il modello in DarkMatterERROR.DAT
                
    real :: Lumi_DM,& !Epsi cumulativa totale secondo la formula di Spergel and Press(Luminosità DM)
            Delta_lumi_DM,& !L'errore che si commette nello stimare la luminosità della DM discretizzando l'integrale
            T_DM,& !Temperatura della Dark Matter
            Lumi_DM_vecchia,&!Variabile di supporto per salvarmi la luminosità della DM durante la convergenza
            errore_max=1e-2,& !Il valore limite per cui si ritiene che l'errore che si commette è trascurabile(Rapporto tra epsi e epsi tot)
            Rapporto_Lumi_DM,&
            Tempo,&
            epsi_tot_min

    real :: T_max=0,& !Variabile per registrare la temperatura massima della stella che poi diventerà la temperatura massima del range in cui si trova quella della DM
            T_min !Variabile per registrare la temperatura minima della stella che poi diventerà la temperatura minima del range in cui si trova quella della DM
    !Queste due variabili permettono di avere un epsi cumulativa >0 o <0, il che garantisce
    !il funzionamento del metodo di bisezione per la convergenza.

    real,dimension(LIM) :: T_mesh,epsi_DM_min

    !##############################################
    !Variabili per killare franec in caso di errore
    !##############################################
    ! Comando e argomenti
    character(100) :: command
    command="killall -9 franec Lancia"
    !#########################
    !Inizializazione variabili
    !#########################

    !Ciclo per inizializzare:
    !Temperatura massima e minima della stella 
    do  i= 1, MAXME 
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
    T_DM=T_max



    !##############################
    !Calcolo dell'integrale di epsi
    !##############################
    
    if (controllo_estremi_DM==1) then
        !Mi calcolo la epsi_DM mesh per mesh a partire T_DM e mi calcolo il suo integrale sulla struttura
        call epsi_DM_routine(T_min,Lumi_DM,Delta_lumi_DM)
        epsi_tot_min=Lumi_DM
        do i=1,MAXME
            epsi_DM_min(i)=epsi_DM(i)!Mi salvo l'array per l'epsi DM per la temperatura minima
            !Questo mi serve in caso di errore
        end do
    endif
    call epsi_DM_routine(T_DM,Lumi_DM,Delta_lumi_DM)
    estremo=1
    
    !Mi calcolo il rapporto tra la luminosità della DM al bordo e quella massima
    Rapporto_Lumi_DM=abs(Lumi_DM)/Max_Lumi_DM
    
    !##############################
    !#      ERRORE GRAVE          #
    !##############################
    if ((epsi_tot_min*Lumi_DM)>=0 .and. controllo_estremi_DM==1) then !Se la luminosità totale della DM per il massimo e il minimo delle Temperature hanno lo stesso
        !segno non è garantito lo zero, ma ciò non è possibile, vuol dire che c'è un errore, qunidi lo scrivo nel file DarkMatterError.DAT
        write(ioDarkError,*)"#####################################################"
        write(ioDarkError,*)"#           ATTENZIONE al seguente modello          #"
        write(ioDarkError,*)"#####################################################"
        !Mi salvo numero modello e età del modello che mi da errore
        write(ioDarkError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),T_DM,char(9),Lumi_DM,char(9),&
        10.00**(ELLOG)*38.27*1e32,char(9),Rapporto_Lumi_DM
        write(ioDarkError,*)"#####################################################"
        write(ioDarkError,*)"# Estremi con lo stesso segno, lo 0 non è garantito #"
        write(ioDarkError,*)"#####################################################"
       !Mi salvo le info per il modello che mi da errore
        write(ioDarkError,*)"La luminosità max è", Lumi_DM,"con una T_max di:",T_max,"la luminosità min è:", epsi_tot_min,"con una T_min di",T_min
        
        if ( ioDarkError_on_off==1 ) then !Attivo la scrittura della scrittura, disattivarla salva spazio su disco
            
        
            write(ioDarkError,*)"###############################################"
            write(ioDarkError,*)"# Per ogni shell valgono le seguenti quantità #"
            write(ioDarkError,*)"###############################################"
            write(ioDarkError,'(A,A,A,A,A,A,A,A,A,A,A)')"#Epsi max",char(9)//char(9)//char(9)//char(9)//char(9)//char(9),&
                            "Epsi cumulativa max",char(9)//char(9)//char(9)//char(9)//char(9),&
                            "Epsi min",char(9)//char(9)//char(9),&
                            "Epsi cumulativa min",char(9)//char(9)//char(9)//char(9),&
                            "Temperatura mesh",char(9)//char(9)//char(9)//char(9),&
                            "Delta m"
        
            !Risetto le variabili della luminosità totale a 0 per poter riscrivere punto per punto nel file di errore il loro andamento
            Lumi_DM=0
            epsi_tot_min=0

            !Mi salvo i vari valori per ogni shell nel file di errore
            do i=1,MAXME
                Lumi_DM=Lumi_DM+(epsi_DM(i)*(G(5,i+1)-G(5,i))*1e33)
                epsi_tot_min=epsi_tot_min+(epsi_DM_min(i)*(G(5,i+1)-G(5,i))*1e33)
                write(ioDarkError,'(ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')epsi_DM(i),char(9),&
                                                                              Lumi_DM,char(9),&
                                                                              epsi_DM_min(i),char(9),&
                                                                              epsi_tot_min,char(9),&
                                                                              T_mesh(i),char(9),&
                                                                              G(5,i+1)-G(5,i)

            
                if ( i==maxme ) then
                    write(ioDarkError,*)"####################################"
                    write(ioDarkError,*)"# Siamo arriviati infonod al MAXME #"
                    write(ioDarkError,*)"####################################"
                end if
            end do
        end if

        call system(command)
    endif

    !Uso metodo di bisezione per trovare la T_DM
    do i=1,max_cicli+N_cicli_extra_convergenza_DM
        
        !Mi calcolo il rapporto tra la luminosità della DM al bordo e quella massima
        Rapporto_Lumi_DM=abs(Lumi_DM)/Max_Lumi_DM

        if(Rapporto_Lumi_DM<=errore_max .and. abs(Lumi_DM)<Delta_lumi_DM) then !Se il rapporto tra la luminosità DM al bordo e quella massima è sotto una certa soglia esco dal ciclo
            if ( i==1 ) then
                write(*,*)"La temperatura della DM è uguale a quella massima della stella e vale:",T_DM,"err_max",errore_max,"Rapporto",Rapporto_Lumi_DM,"Lumi_DM",Lumi_DM,"Max_lumi",Max_Lumi_DM
            end if
            
            if ( i>max_cicli .and. on_off_Error_DM_aggiuntivi==1 .and. ioDarkError_on_off==1 ) then
                write(ioDarkError,'(A)')"CONVERGE"
                write(ioDarkError, '(A, A, A, A, A, A, A, A, A, A, A)')"#_Modello",char(9)//char(9),&
                "Tempo_yr",char(9)//char(9)//char(9)//char(9),&
                "T_DM-1e-6kel",char(9)//char(9)//char(9)//char(9),&
                "Luminosità_DM",char(9)//char(9)//char(9)//char(9),&
                "Lum_tot",char(9)//char(9)//char(9)//char(9),"L_DM(R)/max(L_DM)"
                !Mi salva le varie info per quando la T_DM non converge
                write(ioDarkError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),T_DM,char(9),Lumi_DM,char(9),&
                10.00**(ELLOG)*38.27*1e32,char(9),Rapporto_Lumi_DM
                write(ioDarkError,'(A)')"##################"
                write(ioDarkError,'(A)')"# Andamento  epsi#"
                write(ioDarkError,'(A)')"##################"
                write(ioDarkError,'(A)')"mesh"//char(9)//char(9)//&
                                        "epsi[erg/(gs)]"//char(9)//char(9)//char(9)//char(9)//&
                                        "R"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                        "deltaR"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                        "m"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                        "deltaM"
                do j = 1, MAXME
                    write(ioDarkError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')j,char(9),epsi_DM(j),char(9),G(1,j)*1e10,&
                                                                                                            char(9),(G(1,j+1)-G(1,j))*1e10,char(9),G(5,j)*1e33,&
                                                                                                            char(9),(G(5,j+1)-G(5,j))*1e33
                end do
                
                ! Esegui il comando per bloccare franec utilizzando execute_command_line
                call system(command)
            end if
            
            exit

        !Se entro i max_cicli non converge la T_DM mi salvo le varie variabili su un file di errore
        else if ((Rapporto_Lumi_DM>=errore_max .or. abs(Lumi_DM)>Delta_lumi_DM) .and. i>=max_cicli ) then
            write(*,*)"La temperatura della DM non converge" !Stampa a schermo

            write(ioDarkError,*)"###########################################"
            write(ioDarkError,*)"#  Errore NON di segno degli estremi      #"
            write(ioDarkError,*)"###########################################"
            write(ioDarkError, '(A, A, A, A, A, A, A, A, A, A, A)')"#_Modello",char(9)//char(9),&
            "Tempo_yr",char(9)//char(9)//char(9)//char(9),&
            "T_DM-1e-6kel",char(9)//char(9)//char(9)//char(9),&
            "Luminosità_DM",char(9)//char(9)//char(9)//char(9),&
            "Lum_tot",char(9)//char(9)//char(9)//char(9),"L_DM(R)/max(L_DM)"
            !Mi salva le varie info per quando la T_DM non converge
            write(ioDarkError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),T_DM,char(9),Lumi_DM,char(9),&
            10.00**(ELLOG)*38.27*1e32,char(9),Rapporto_Lumi_DM
            if ( ioDarkError_on_off==1 ) then !Abilito la scrittura delle variabili mesh per mesh, disattivare per salvare spazio nel disco
                write(ioDarkError,'(A)')"##################"
                write(ioDarkError,'(A)')"# Andamento  epsi#"
                write(ioDarkError,'(A)')"##################"
                write(ioDarkError,'(A)')"mesh"//char(9)//char(9)//&
                                    "epsi[erg/(gs)]"//char(9)//char(9)//char(9)//char(9)//&
                                    "R"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                    "deltaR"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                    "m"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                                    "deltaM"
                do j = 1, MAXME
                    write(ioDarkError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')j,char(9),epsi_DM(j),char(9),G(1,j)*1e10,&
                                                                                                        char(9),(G(1,j+1)-G(1,j))*1e10,char(9),G(5,j)*1e33,&
                                                                                                        char(9),(G(5,j+1)-G(5,j))*1e33
                end do
            end if
            
            if ( on_off_Error_DM_aggiuntivi/=1 ) then !Se non voglio informazioni aggiuntive, tipo se converge ai passi successivi, al primo errore esco, altrimenti no
                
                exit 
            end if

        end if

        Lumi_DM_vecchia=Lumi_DM !Salvo la epsi_tot precedente per compararla con quella nuova.
        
        T_DM=(T_max+T_min)/2!Mi calcolo la nuova temperatura della DM
        
        call epsi_DM_routine(T_DM,Lumi_DM,Delta_lumi_DM)

        if ((Lumi_DM >= 0.0d0) == (Lumi_DM_vecchia >= 0.0d0))then
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
        
        
        
        
    end do

end subroutine convergenza_epsi_DM



contains

  function logsumexp(nums) result(result)
    !Funzione per calcolarsi il logaritmo della somma degli esponenziali
    !dell'array nums
    real, dimension(:), intent(in) :: nums
    integer :: ct
    real :: max_exp, sum
    real :: result  ! Dichiarazione del tipo di ritorno
    integer :: i
        
    ct=size(nums)
    max_exp = nums(1)
    sum = 0.0d0
  
    do i = 2, ct
      if (nums(i) > max_exp) max_exp = nums(i)
    end do
  
    do i = 1, ct
      sum = sum + exp(nums(i) - max_exp)
    end do
  
    result = log(sum) + max_exp
  end function logsumexp

  function prob_r_quad(N_punti,R_max) result(punti_generati)
    !Genera un array distribuito  con una distribuzione di probabilità  che scala come r^2
    integer :: N_punti,index
    real,dimension(N_punti) :: punti_generati
    real numero_casuale, r_casuale,&
        R_max

    !Inizializza il generatore di numeri casuali
    call random_seed()

    do index = 1, N_punti
      call random_number(numero_casuale)

      r_casuale=R_max*(numero_casuale)** (1.0/3.0)

      punti_generati(index)=r_casuale
    end do

  end function prob_r_quad

  function log_int_exp(N_punti,esponenti) result(log_int)
    !Mi calcola il logaritmo dell'integrale di 4*pi*r^2*exp(esponente) tramite metodo montecarlo
    !Sfruttando il fatto che l'integrale di f(x)*g(x) è uguale a 1/N* somma di f(x_i) distribuite come g(x)
    use mesh
    use fisica
    use costanti !Serve per il pi greco

    integer :: N_punti,mesh_val,indice
    real :: R_max, log_int,R_mesh
    real,dimension(:) :: esponenti
    real,dimension(:),allocatable :: esponenti_con_R_distribuito,r_quad

    allocate(esponenti_con_R_distribuito(N_punti))
    allocate(r_quad(N_punti))

    R_max=G(1,MAXME)!Prendo il massimo del raggio della struttura

    r_quad=prob_r_quad(N_punti,R_max)!Genero un set di valori per il raggio distribuiti secondo r^2
    do indice=1,N_punti
      do mesh_val = 1, MAXME

        R_mesh=G(1,mesh_val)
        R_mesh_p=G(1,mesh_val+1)

        if ( r_quad(indice)>=R_mesh .and. r_quad(indice)<R_mesh_p ) then !Controllo se l'R genearto casualmente sta in un certo range, e se sta in tale
          !range gli assegno l'esponente in quel range, cosi genero l'array f(x_i) con x_i distribuito come g(x).
          esponenti_con_R_distribuito(indice)=esponenti(mesh_val)
          exit ! Esci dal ciclo quando hai trovato il valore dell'esponente per
          !per il raggio indice     
        end if
      
      end do
    end do
    log_int=logsumexp(esponenti_con_R_distribuito)-log(real(N_punti))+log((4.0/3.)*pigre*((R_max*1e10)**(3)))!La divisione per N porta il fattore  -log(N)
    !Mentre il fattore di volume deriva dal fatto che la distribuzione g(x) era stata normalizzata a 1, quindi per correttezza dimensionale bisogna rimoltiplicare
    !Per il volume

  end function log_int_exp