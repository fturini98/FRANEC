# Dark Matter
clicca per tornare a:
- [README](../../README.md)
- [Documentazione FRANEC](../Documentazione.md)

#
E'  stato modificato il codice standard per includere l'effetto di un'ipotetica particella WIMP catturata all'interno della stella.

Le basi da cui si è partiti per sviluppare il codice sono:

- La WIMP è abbastanza massiva da non evaporare, quindi si accumula nella stella in funzione del flusso

- La WIMP non annichila

- Si trascura la self-capture

- La WIMP ha una sezione d'urto molto bassa tale che siamo nel regime di **Kundsen** (Vale a dire si considera isoterma su tutta la struttura.)

## Indice
- [Modifiche rispetto a FRANEC standard](#modifiche-rispetto-franec-standard)

- [Subroutine per la DM](#subroutine-per-la-dm)

## Modifiche rispetto FRANEC  standard
- Creazione dei file:

    - **Cattura_DM.f90**: Questo file contiene le routine per calcolarsi la cattura della DM da parte della stella.

    - **epsi_Dark_Matter.f90**: Questo file contiene routine per calcolarsi l'array mesh per mesh con i valori dell'epsi dovuto alla DM

- Modifica dei file per includere le routine:

    - **Makefile**: 
        - sono stati aggiunti nella sezione OBJ il file:
            - **Cattura_DM.o**, che legge e interpreta il file catturaDM.f90 e genera il file catturaDM.o da usare nel programma.

            - **epsi_Dark_Matter.o**, che fa a stessa cosa per il file relativo alla epsi della DM

    - **franec.f90**:

- Sono stati modificati i programmi
    - [**lancia.c**](./../../Gestione_lancia/lancia.c)

    - [**driver.c**](./../../Gestione_driver/driver.c)

    per poter genereare dei file in output (*DarkMatter.DAT* e *DarMatterERROR.DAT*) dove registrare varie variabili tramite la subroutine [stampa_epsi_DM](Subroutine/epsi_Dark_Matter.md) (definita in [epsi_Dark_Matter.f90](./../../epsi_Dark_Matter.f90)). Per poter generare questi file di putput è stato necessario creare nelle cartelle tools/driver e tools/lancia i rispettivi link simbolici tramite:
        
        ln -s ./../../DMsrc.90/Gestione_driver/driver.c driver.c

        ln -s ./../../DMsrc.90/Gestione_lancia/lancia.c lancia.c

    *Attenzione:* ogni volta che vengono fatte delle modifiche al file sorgente bisogna andare nelle cartelle tools/driver e tools/lancia e lanciare il make per generare i rispettivi programmi.


## Subroutine per la DM:

- [Cattura_DM.f90](../../Cattura_DM.f90)

    - **Cattura_DM( out: C_tot)** : La routine cattura_DM gestisce la cattura della DM da parte della stella. Restituisce il rate di cattura totale (*C_tot*) e aggiorna il numero totale di particelle di DM nella stella nella variabile globale *N_tot_DM* (definita nel modulo [Dark_Matter](./../../moduli.f90)).

    - **Cattura_DM_geometrica** : La routine si calcola il rate di cattura geometrico della stella, necessario per calcolarsi l'effetto di saturazione della cattura

    - **xi_DM_sub**: La routine si calcola il fattore $\xi$ di correzzione dovuto al movimento della stella nel rest frame della DM. Tale fattore dipende dalla velocità della stella, dalla velocità di dispersione delle stelle e dalla velocità di fuga, quindi può essere definito shell per shell.



