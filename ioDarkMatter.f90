subroutine stampa_epsi_DM(NMD,epsi_DM_tot,T_DM,Tempo,C_tot)
!Questa subroutine serve a stamaprsi i vari valori dell'epsi totale relativa alla DM.
!Il file relativo si chiama epsi_DM e la sua unit vale unit=ioDark==411(è definita all'interno del moduleo Dark_Matter)
        
        use Dark_Matter
    
        use tempe !Serve per ellog(luminosità totale stella)
    
        implicit none
    
        integer :: NMD !Il numero del modello corrente

        real :: epsi_DM_tot,&!La luminosità totale della DM
                T_DM,&!temperatura della DM
                Tempo,&!Età della stella
                Rapporto_Lumi_DM,&! Rapporto tra luminosità totale e della DM
                C_tot !Rate di cattura nell'intervallo di tempo
    
        if (first_DM_write==1) then!Serve a segnare cosa sono le varie colonne
            write(ioDark, '(A)')"# Modello"//char(9)//char(9)//&
                        "Tempo_yr"//char(9)//char(9)//char(9)//char(9)//&
                        "Teff"//char(9)//char(9)//char(9)//char(9)//char(9)//&
                        "Lum_tot"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                        "T_DM-1e-6kel"//char(9)//char(9)//char(9)//char(9)//&
                        "C_tot_[HZ]"//char(9)//char(9)//char(9)//char(9)//char(9)//char(9)//&
                        "Luminosità_DM"//char(9)//char(9)//char(9)//char(9)//&
                        "L_DM/L_tot"//char(9)//char(9)//char(9)//char(9)//char(9)//&
                        "N_DM"
                        
            first_DM_write=0
        endif
        Rapporto_Lumi_DM=abs(epsi_DM_tot)*1e-32/&
                            (10.0**(ELLOG)*&
                            38.27)!Luminosità solare *10^-32 erg/s
        write(ioDark,'(I0, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15, A, ES25.15)')NMD,char(9),Tempo,char(9),10.0**TEFF,char(9),10.0**(ELLOG)*38.27*1e32,char(9),T_DM,char(9),C_tot,char(9),epsi_DM_tot,char(9),Rapporto_Lumi_DM,char(9),N_DM_tot
    end subroutine stampa_epsi_DM