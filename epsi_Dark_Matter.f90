subroutine epsi_DM_routine(T_DM,epsi_DM_tot)
!Subroutine che si calcola la epsi DM in maniera vettoriale (slavandola nella variabile globale data dal modulo) 
!e restituisce epsi_DM_tot
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
    
    implicit none

    integer :: i,j, & !indice
               ele !indice per elementi


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

        if ( i==1 .or. G(1,i)==0 ) then !Questo if serve ad evitare problemi del divided by 0 del raggio nella shell più interna
            
            exp_dist_DM(i)=exp(-(mass_DM*GeV_grammi*&
            (Ggrav*(G(5,i)+G(5,i+1))/(G(1,i)+G(1,i+1))*1e23))/&!Ho fatto la media tra il mesh i e il mesh i+1 per le variabili fisiche
            (Kboltz*T_DM*1e6))

            norm_DM_distr=((G(1,i)+G(1,i+1))/2)**2*(G(1,i+1)-G(1,i))*exp_dist_DM(i)*1e30

            if (G(1,i+1)==0 ) then!Questo serve a prendere in considerazione quando il raggio è messo
            ! a 0 sia nel mesh i che nel mesh i+1

            exp_dist_DM(i)=0!Setto a 0 la distribuzione di DM perché qui non ho definito il potenziale grav.

            endif
        
        else if (G(1,i)/=0) then

            exp_dist_DM(i)=exp(-(mass_DM*GeV_grammi*&
            (Ggrav*(G(5,i))/(G(1,i))*1e23))/&
            (Kboltz*T_DM*1e6))

            norm_DM_distr=((G(1,i))**2*(G(1,i+1)-G(1,i))*exp_dist_DM(i)*1e30)+norm_DM_distr
        end if

    end do

    !#################################
    !Calcolo epsi con loop su elementi
    !#################################

    epsi_DM_tot=0!Inizializzo la funzione cumulartiva della epsi a 0

    do i = 1, MAXME!Loop sui mesh
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
            Kboltz*1e6*(T_DM-T_mesh(i)))!Differenza temperature

        end do

        !Mi calcolo anche l'integrale su epsi
        if ( epsi_DM(i)/=0 .or. i==MAXME) then !Questa condizione dovrebbe accellerare il calcolo e in più mi garantisce di evitare il problema del mesh i=LIM in cui R(i+1) non è
            !definito. Infatti epsi_DM(LIM)=0 per definizione in quanto le abbondanze in LIM sono 0 e epsi dipende da loro in maniera lineare.
            epsi_DM_tot=epsi_DM_tot+(epsi_DM(i)*(M_mesh(i+1)-M_mesh(i))*1e33)!Mi calcolo la lumi totale tramite integrazione in massa
            !epsi_DM_tot=epsi_DM_tot+(epsi_DM(i)*(G(5,i+1)-G(5,i))*1e33)

            !Errore a schermo per quando il delta m al maxme è minore di 0
            if((M_mesh(i+1)-M_mesh(i))<0 .and. i/=MAXME )then 
                write(*,*)"Delta m minore di 0, mesh "
                do j=1,(MAXME+10)
                    write(*,*)M_mesh(j)
                    if ( j==MAXME ) then
                        write(*,*)"SIAMO AL MAXME"
                    end if
                end do
                stop
            endif
        end if
    end do

end subroutine epsi_DM_routine


subroutine convergenza_epsi_DM(epsi_DM_tot,T_DM,NMD,tempo)
    
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

    integer :: i,&!indice
                estremo,& !Mi serve per memorizzare se l'estremo precedente dell'intervallo di temperatura era a destra o sinistra
                !estremo=1 Destra, quindi prima avevo calcolato per T_max
                !estremo=0 Sinistra, quindi prima avevo calcolato per T_min
                max_cicli=20,&!Numero massimo di cicli per la convergenza della T_DM.(Con venti ho una precisone di circa 1e-6 in T_DM)
                NMD !Servwe per quando ho l'errore di convergenza per stamaprimi il modello in DarkMatterERROR.DAT

    real :: epsi_DM_tot,& !Epsi cumulativa totale secondo la formula di Spergel and Press(Luminosità DM)
            T_DM,& !Temperatura della Dark Matter
            epsi_DM_tot_vecchia,&!Variabile di supporto per salvarmi l'epsi tot durante la convergenza
            errore_max=1e-4,& !Il valore limite per cui si ritiene che l'errore che si commette è trascurabile(Rapporto tra epsi e epsi tot)
            Rapporto_Lumi_DM,&
            Tempo,&
            epsi_tot_min

    real :: T_max=0,& !Variabile per registrare la temperatura massima della stella
            T_min !VAriabile per registrare la temperatura minima della stella
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
        else if (T_mesh(i)<T_min .and. T_mesh(i)/= 0) then
            T_min=T_mesh(i)
        end if

    end do

    !Assegno come primo valore di test per la temperatura della DM la tmperatura massima della stella(che non per forza è quella centrale)
    T_DM=T_max



    !##############################
    !Calcolo dell'integrale di epsi
    !##############################
    
    !Mi calcolo la epsi_DM mesh per mesh a partire T_DM e mi calcolo il suo integrale sulla struttura
    call epsi_DM_routine(T_min,epsi_DM_tot)
    epsi_tot_min=epsi_DM_tot
    do i=1,MAXME
        epsi_DM_min(i)=epsi_DM(i)!Mi salvo l'array per l'epsi DM per la temperatura minima
        !Questo mi serve in caso di errore
    end do
    call epsi_DM_routine(T_DM,epsi_DM_tot)
    estremo=1
    
    
    !##############################
    !#      ERRORE GRAVE          #
    !##############################
    if ((epsi_tot_min*epsi_DM_tot)>=0) then !Se la luminosità totale della DM per il massimo e il minimo delle Temperature hanno lo stesso
        !segno non è garantito lo zero, ma ciò non è possibile, vuol dire che c'è un errore, qunidi lo scrivo nel file DarkMatterError.DAT
        write(ioDrakError,*)"#####################################################"
        write(ioDrakError,*)"#           ATTENZIONE al seguente modello          #"
        write(ioDrakError,*)"#####################################################"
        !Mi salvo numero modello e età del modello che mi da errore
        write(ioDrakError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),T_DM,char(9),epsi_DM_tot,char(9),&
        10.00**(ELLOG)*38.27*1e32,char(9),Rapporto_Lumi_DM
        write(ioDrakError,*)"#####################################################"
        write(ioDrakError,*)"# Estremi con lo stesso segno, lo 0 non è garantito #"
        write(ioDrakError,*)"#####################################################"
       !Mi salvo le info per il modello che mi da errore
        write(ioDrakError,*)"L'epsi max è", epsi_DM_tot,"con una T_max di:",T_max,"l'epsi min è:", epsi_tot_min,"con una T_min di",T_min
        
        write(ioDrakError,*)"###############################################"
        write(ioDrakError,*)"# Per ogni shell valgono le seguenti quantità #"
        write(ioDrakError,*)"###############################################"
        write(ioDrakError,'(A,A,A,A,A,A,A,A,A,A,A)')"#Epsi max",char(9)//char(9)//char(9)//char(9)//char(9)//char(9),&
                            "Epsi cumulativa max",char(9)//char(9)//char(9)//char(9)//char(9),&
                            "Epsi min",char(9)//char(9)//char(9),&
                            "Epsi cumulativa min",char(9)//char(9)//char(9)//char(9),&
                            "Temperatura mesh",char(9)//char(9)//char(9)//char(9),&
                            "Delta m"
        
        !Risetto le variabili della luminosità totale a 0 per poter riscrivere punto per punto nel file di errore il loro andamento
        epsi_DM_tot=0
        epsi_tot_min=0

        !Mi salvo i vari valori per ogni shell nel file di errore
        do i=1,MAXME
            epsi_DM_tot=epsi_DM_tot+(epsi_DM(i)*(G(5,i+1)-G(5,i))*1e33)
            epsi_tot_min=epsi_tot_min+(epsi_DM_min(i)*(G(5,i+1)-G(5,i))*1e33)
            write(ioDrakError,'(ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')epsi_DM(i),char(9),&
                                                                              epsi_DM_tot,char(9),&
                                                                              epsi_DM_min(i),char(9),&
                                                                              epsi_tot_min,char(9),&
                                                                              T_mesh(i),char(9),&
                                                                              G(5,i+1)-G(5,i)

            
            if ( i==maxme ) then
                write(ioDrakError,*)"####################################"
                write(ioDrakError,*)"# Siamo arriviati infonod al MAXME #"
                write(ioDrakError,*)"####################################"
            end if
        end do
        
        stop
    endif

    !Uso metodo di bisezione per trovare la T_DM
    do i=1,max_cicli
        
        !Mi calcolo il rapporto tra la luminosità della DM e quella totale
        Rapporto_Lumi_DM=abs(epsi_DM_tot)*1e-32/&
                        (10.0**(ELLOG)*&
                        38.27)!Luminosità solare *10^-32 erg/s

        if(Rapporto_Lumi_DM<=errore_max) then !Se il rapporto tra la luminosità DM e quella totale è sotto una certa soglia esco dal ciclo
            if ( i==1 ) then
                write(*,*)"La temperatura della DM è uguale a quella massima della stella e vale:",T_DM,"err_max",errore_max,"Rapporto",Rapporto_Lumi_DM
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
        
        
        !Se entro i max_cicli non converge la T_DM mi salvo le varie variabili su un file di errore.
        if (Rapporto_Lumi_DM>=errore_max .and. i==max_cicli) then
            write(*,*)"La temperatura della DM non converge" !Stampa a schermo

            write(ioDrakError,*)"###########################################"
            write(ioDrakError,*)"#  Errore NON di segno degli estremi      #"
            write(ioDrakError,*)"###########################################"
            write(ioDrakError, '(A, A, A, A, A, A, A, A, A, A, A)')"#_Modello",char(9)//char(9),&
            "Tempo_yr",char(9)//char(9)//char(9)//char(9),&
            "T_DM-1e-6kel",char(9)//char(9)//char(9)//char(9),&
            "Luminosità_DM",char(9)//char(9)//char(9)//char(9),&
            "Lum_tot",char(9)//char(9)//char(9)//char(9),"L_DM/L_tot"
            !Mi salva le varie info per quando la T_DM non converge
            write(ioDrakError,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),T_DM,char(9),epsi_DM_tot,char(9),&
            10.00**(ELLOG)*38.27*1e32,char(9),Rapporto_Lumi_DM
            exit
        end if
        
    end do

end subroutine convergenza_epsi_DM

subroutine stampa_epsi_DM(NMD,epsi_DM_tot,T_DM,Tempo,C_tot)
!Questa subroutine serve a stamaprsi i vari valori dell'epsi totale relativa alla DM.
!Il file relativo si chiama epsi_DM e la sua unit vale unit=ioDark==411(è definita all'interno del moduleo Dark_Matter)
    
    use Dark_Matter

    use tempe !Serve per ellog(luminosità totale stella)

    implicit none

    integer :: NMD
    real :: epsi_DM_tot,T_DM,Tempo,Rapporto_Lumi_DM !La luminosità totale della DM e la temperatura della DM
    real :: C_tot !Rate di cattura nell'intervallo di tempo

    if (first_DM_write==1) then!Serve a segnare cosa sono le varie colonne
        write(ioDark, '(A, A, A, A, A, A, A, A, A, A, A, A, A)')"#_Modello",char(9)//char(9),&
                    "Tempo_yr",char(9)//char(9)//char(9)//char(9),&
                    "T_DM-1e-6kel",char(9)//char(9)//char(9)//char(9),&
                    "Luminosità_DM",char(9)//char(9)//char(9)//char(9),&
                    "Lum_tot",char(9)//char(9)//char(9)//char(9)//char(9),&
                    "L_DM/L_tot",char(9)//char(9)//char(9)//char(9)//char(9),&
                    "Teff",char(9)//char(9)//char(9)//char(9)//char(9),&
                    "C_tot [HZ]"
        first_DM_write=0
    endif
    Rapporto_Lumi_DM=abs(epsi_DM_tot)*1e-32/&
                        (10.0**(ELLOG)*&
                        38.27)!Luminosità solare *10^-32 erg/s
    write(ioDark,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),T_DM,char(9),epsi_DM_tot,char(9),10.0**(ELLOG)*38.27*1e32,char(9),Rapporto_Lumi_DM,char(9),10.0**TEFF,char(9),C_tot
end subroutine stampa_epsi_DM