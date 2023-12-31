subroutine epsi_DM_routine(T_DM,Lumi_DM,NMD)
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
    
    use FUNZIONI !Serve per usare la erf
    
    implicit none

    integer :: i,j, & !indice
               ele,& !indice per elementi
               NMD !Numero del modello
    
    real :: M_sup,&!Massa della stella
            R_sup!Raggio stella


    real,dimension(LIM) :: T_mesh,& !Variabile per salvare la temperatura della stella
                           M_mesh !Variabile di supporto per l'array delle masse
            
    real, dimension(MAXME):: phi!Array del potenziale gravitazionale ridotto

    real :: T_DM,& !Temperatura della DM per cui si calcola l'epsi_DM
            sigma_ele,& !Sezione d'urto per elemento
            V_centrale,& !potenziale gravitazionale al centro della stella
            densita_cost,& !Denistà per modello stella a densità costante (calcolata mettoendo che i potenziali reale e parametrizzato metchino al centro)
            r_scala_DM_quad,& !Quadrato raggio scala DM
            exp_dist_DM !distiribuzione maxwell Boltzman della DM (la parametrizo con -r^2)
    
    real :: Lumi_DM,& !Integrale sull'epsi(Luminosità DM)
            Lumi_mesh !La luminosità del mesh
         
    !##########################
    !Inizializzazione variabili
    !##########################

    !Mi trovo raggio e massa della stella
    R_sup=G(1,MAXME)
    M_sup=G(5,MAXME)


    !Calcolo il potenziale gravitazionale Ridotto,(Per l'ultima implementazione serve solo quello al centro)        
    phi(MAXME)=1.00
    do  i= (MAXME-1),1 , -1! ciclo dall'esterno all'interno per calcolarsi phi secondo 
    !https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1991ApJ...368..626D&defaultprint=YES&filetype=.pdf
    !Per denistà costante è come l'integrale di Gaus per una sfera uniformemente carica
        if ( i==1 ) then
            phi(i)=phi(i+1)+R_sup/M_sup*((G(5,2)+G(5,1))/2)/(((G(1,2)+G(1,1))/2)**2)*G(1,2)!Questo serve a evitare il diviso 0 al centro, uso come valori di fisica della shell quelli dell'intershell
        else
            phi(i)=phi(i+1)+R_sup/M_sup*(G(5,i)/(G(1,i)**2))*(G(1,i+1)-G(1,i))
        end if
               
    end do

    V_centrale=Ggrav*(M_sup/R_sup)*phi(1)*1e23 !Mi calcolo il potenziale gravitazionale centrale

    densita_cost=V_centrale/(2*Ggrav*pigre*(R_sup*1e10)**2) !Mi calcolo la densità per la parametrizazione a densità costante

    r_scala_DM_quad=(Kboltz*T_DM)/(2*pigre*Ggrav*densita_cost*(mass_DM*GeV_grammi))!Mi calcolo il quadrato del raggio scala della DM

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

        
    end do

    !#################################
    !Calcolo epsi con loop su elementi
    !#################################

    Lumi_DM=0!Inizializzo la funzione cumulartiva della epsi a 0
    Lumi_DM_positiva=0!Inizializo la funzione cumulativa della Luminosità positiva a 0
    Lumi_DM_negativa=0!Inizializo la funzione cumulativa della Luminosità negativa a 0
    Max_Lumi_DM=0!Inizializzo a o il massimo del modulo della luminosità a 0.
    
    do i = 1, MAXME!Loop sui mesh
        epsi_DM(i)=0 !Lo setto a zero per poi sommarci su ad ogni elemento
        
        exp_dist_DM=exp(-(G(1,i)**2/r_scala_DM_quad)-log(pigre*r_scala_DM_quad*(sqrt(pigre*r_scala_DM_quad)*erf((R_sup*1e10)/(sqrt(r_scala_DM_quad)))-2*(R_sup*1e10)*exp(-((R_sup*1e10)**2)/r_scala_DM_quad))))

        if (exp_dist_DM>1 .or. isnan(exp_dist_DM) ) then
            write(ioDarkError,*)"Errore calcolo densità DM per il modello:",NMD,"con la densità della DM normalizzata=",exp_dist_DM,"al mesh:",i
            !call system("killall -9 franec Lancia")
        end if
        do ele = 1, MELE !Loop su tutti gli elementi
            sigma_ele=sigma0_DM*(massa_ele(ele)/massa_ele(1))**4 *&
            (((mass_DM*GeV_grammi+massa_ele(1)*u_to_gramms)/(mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms))**2)!Sezione d'urto come  doi:10.1103/physrevd.66.053005

            !Mi calcolo l'epis nel mesh per le condizioni al bordo i
            epsi_DM(i)=epsi_DM(i)+(8.*sqrt(2./pigre)*& 
            N_DM_tot*exp_dist_DM*& !Densità numerica della DM mesh per mesh(la normalizazione è gia inserita)
            XXX(ele,i)/(massa_ele(ele)*u_to_gramms)*sigma_ele*& !Densità dell'elemento chimico * sezione d'urto
            (mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms)/((mass_DM*GeV_grammi+massa_ele(ele)*u_to_gramms)**2)*&!Fattore differenza massa
            sqrt(Kboltz*1e6*(massa_ele(ele)*u_to_gramms*T_DM+mass_DM*GeV_grammi*T_mesh(i))/(mass_DM*GeV_grammi*massa_ele(ele)*u_to_gramms))*& !Fattore somma temperature
            Kboltz*1e6*(T_DM-T_mesh(i)))!Differenza temperature
           
        end do

        !Mi calcolo anche l'integrale su epsi
        if ( epsi_DM(i)/=0 .or. i==MAXME) then !Questa condizione dovrebbe accellerare il calcolo e in più mi garantisce di evitare il problema del mesh i=LIM in cui R(i+1) non è
            !definito. Infatti epsi_DM(LIM)=0 per definizione in quanto le abbondanze in LIM sono 0 e epsi dipende da loro in maniera lineare.
            Lumi_mesh=(epsi_DM(i)*(M_mesh(i+1)-M_mesh(i))*1e33)!La luminosità della DM nel mesh
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

end subroutine epsi_DM_routine


subroutine convergenza_epsi_DM(Lumi_DM,T_DM,NMD,tempo)
        
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
                i_max,i_min,&!Indice per il mesh dove c'è temperatura massima e la temperatura minima
                i_DM,& !Indice del mesh in cui la T_mesh=T_DM
                estremo,& !Mi serve per memorizzare se l'estremo precedente dell'intervallo di temperatura era a destra o sinistra
                !estremo=1 Destra, quindi prima avevo calcolato per T_max
                !estremo=0 Sinistra, quindi prima avevo calcolato per T_min
                max_cicli=21,&!Numero massimo di cicli per la convergenza della T_DM.(Con venti ho una precisone di circa 1e-6 in T_DM)
                NMD !Servwe per quando ho l'errore di convergenza per stamaprimi il modello in DarkMatterERROR.DAT
                
    real :: Lumi_DM,& !Epsi cumulativa totale secondo la formula di Spergel and Press(Luminosità DM)
            T_DM,& !Temperatura della Dark Matter
            Lumi_DM_vecchia,&!Variabile di supporto per salvarmi la luminosità della DM durante la convergenza
            Rapporto_Lumi_DM,&
            Tempo,&
            epsi_tot_min

    real :: T_max=0,& !Variabile per registrare la temperatura massima della stella che poi diventerà la temperatura massima del range in cui si trova quella della DM
            T_min !Variabile per registrare la temperatura minima della stella che poi diventerà la temperatura minima del range in cui si trova quella della DM
    !Queste due variabili permettono di avere un epsi cumulativa >0 o <0, il che garantisce
    !il funzionamento del metodo di bisezione per la convergenza.

    real,dimension(LIM) :: T_mesh,epsi_DM_min

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
            i_max=i !Prendo l'indice del mesh in cui ho la T_max
        else if (T_mesh(i)<T_min .and. T_mesh(i)/= 0) then
            T_min=T_mesh(i)
            i_min=i !Prendo l'indice del mesh in cui ho la T_min
        end if

    end do

    !Assegno come primo valore di test per la temperatura della DM la tmperatura massima della stella(che non per forza è quella centrale)
    T_DM=T_max



    !##############################
    !Calcolo dell'integrale di epsi
    !##############################
    
    if (controllo_estremi_DM==1) then
        !Mi calcolo la epsi_DM mesh per mesh a partire T_DM e mi calcolo il suo integrale sulla struttura
        call epsi_DM_routine(T_min,Lumi_DM,NMD)
        epsi_tot_min=Lumi_DM
        do i=1,MAXME
            epsi_DM_min(i)=epsi_DM(i)!Mi salvo l'array per l'epsi DM per la temperatura minima
            !Questo mi serve in caso di errore
        end do
    endif
    call epsi_DM_routine(T_DM,Lumi_DM,NMD)
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
    endif

    !Uso metodo di bisezione per trovare la T_DM
    do i=1,max_cicli
        
        !Mi calcolo il rapporto tra la luminosità della DM al bordo e quella massima
        Rapporto_Lumi_DM=abs(Lumi_DM)/Max_Lumi_DM

        if(i_max + 1  ==i_min) then !Se la temperatura della DM conerge a lla temperatura di un mesh esco dal ciclo
            if ( i==1 ) then
                write(*,*)"La temperatura della DM è uguale a quella massima della stella e vale:",T_DM,"Rapporto",Rapporto_Lumi_DM,"Lumi_DM",Lumi_DM,"Max_lumi",Max_Lumi_DM
            end if

            exit

        !Se entro i max_cicli non converge la T_DM mi salvo le varie variabili su un file di errore
        else if (i_min/=i_max + 1  .and. i>=max_cicli) then
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
            
            exit 

        end if
    
        Lumi_DM_vecchia=Lumi_DM !Salvo la epsi_tot precedente per compararla con quella nuova.

        i_DM=(i_max+i_min)/2.
        
        T_DM=T_mesh(i_DM)!Mi calcolo la nuova temperatura della DM

        call epsi_DM_routine(T_DM,Lumi_DM,NMD)

        if ((Lumi_DM >= 0.0d0) == (Lumi_DM_vecchia >= 0.0d0))then
            if ( estremo==1 ) then
                i_max=i_DM
            else
                i_min=i_DM
            end if
        
        else
            if ( estremo==1 ) then
                i_min=i_DM
                estremo=0
            else
                i_max=i_DM
                estremo=1
            end if
        end if
        
        


    end do
        
end subroutine convergenza_epsi_DM