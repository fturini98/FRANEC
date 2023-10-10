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

    - **ioDarkMatter.f90**: Questo file contiene la routine per riempire il file di output DarkMatter.DAT chiamata *stampa_epsi_DM()*

- Modifica dei file per includere le routine:

    - **Makefile**: 
        - sono stati aggiunti nella sezione OBJ il file:
            - **Cattura_DM.o**, che legge e interpreta il file catturaDM.f90 e genera il file catturaDM.o da usare nel programma.

            - **epsi_Dark_Matter.o**, che fa a stessa cosa per il file relativo alla epsi della DM

            - **ioDarkMatter.o**, he legge e interpreta il file ioDarkMatter.f90 e genera il file ioDarkMatter.o da usare nel programma.

    - **franec.f90**:

- Sono stati modificati i programmi
    - [**lancia.c**](./../../Gestione_lancia/lancia.c)

    - [**driver.c**](./../../Gestione_driver/driver.c)

    per poter genereare dei file in output (*DarkMatter.DAT*, *DarMatterERROR.DAT*,*DarkMatterCattura.DAT*) dove registrare varie variabili tramite la subroutine [stampa_epsi_DM](Subroutine/epsi_Dark_Matter.md) (definita in [ioDarkMatter.f90](./../../ioDarkMatter.f90)) e altri comanndi nelle subroutini **Cattura_DM** e **epsi_DM_routine**(file di errore). Per poter generare questi file di putput è stato necessario creare nelle cartelle tools/driver e tools/lancia i rispettivi link simbolici tramite:
        
        ln -s ./../../DMsrc.90/Gestione_driver/driver.c driver.c

        ln -s ./../../DMsrc.90/Gestione_lancia/lancia.c lancia.c

    *Attenzione:* ogni volta che vengono fatte delle modifiche al file sorgente bisogna andare nelle cartelle tools/driver e tools/lancia e lanciare il make per generare i rispettivi programmi.


## Subroutine per la DM:

- [Cattura_DM.f90](../../Cattura_DM.f90)

    - **Cattura_DM( out: C_tot)** : La routine cattura_DM gestisce la cattura della DM da parte della stella. Restituisce il rate di cattura totale (*C_tot*) e aggiorna il numero totale di particelle di DM nella stella nella variabile globale *N_tot_DM* (definita nel modulo [Dark_Matter](./../../moduli.f90)). Questa subrutine se la flag *ioDarkCattura_on_off* è impostata ad 1 scrive nel file DarkMatterCattura.DAT i seguenti valori per ogni passo temporale:

        - *Modello*(NMD)
        - *Età della stella [yr]*
        - *Rate di cattura dovuto all'interazione [Hz]*
        - *Rate di catura gemoetrico [Hz]*
        - *Rate di cattura totale [Hz]* (Cioé quello di cattura dovuto all'interazione attenuato dalla saturazione)
        - *Passo temporale*
        - *Numero di particelle di Dark Matter catturate nel passo temporale*
        - *Numero totale di particelle catturate fino a quell'istante di tempo dalla struttura*
        - *Rate di cattura deboli per ogni singolo elemento*

    - **Cattura_DM_geometrica** : La routine si calcola il rate di cattura geometrico della stella, necessario per calcolarsi l'effetto di saturazione della cattura

    - **xi_DM_sub**: La routine si calcola il fattore $\xi$ di correzzione dovuto al movimento della stella nel rest frame della DM. Tale fattore dipende dalla velocità della stella, dalla velocità di dispersione delle stelle e dalla velocità di fuga, quindi può essere definito shell per shell.

- [ioDarkMatter.f90](../../ioDarkMatter.f90)

    - **stampa_epsi_DM(NMD,epsi_DM_tot,T_DM,Tempo,C_tot)**: Questa routine scrive nel file DarkMatter.DAT(*La scrittura di questo file può essere disabilitata per salvare spazio sul disco impostando* **ioDark_on_off/=1** ) i seguenti valori per ogni step temporale:

        - *Numero modello*(NMD)
        - *Età della stella[yr]*
        - *Temperatura effettiva [$k$]*
        - *Luminosità della stella [$\frac{erg}{s}$]*
        - *Temepratura Dark Matter [$\frac{k}{10^6}$]*
        - *Rate di cattura Dark Matter [$Hz$]*
        - *Luminosità Dark Matter [$\frac{erg}{s}$]*
        - *Rapporto tra luminosità della DM e quella totale della stella*

            - Questo rapporto deve stare sotto una soglia di $10^{-4}$ se no viene considerato che la temperatura della Dark Matter non è giustoa, in quanto la luminosità totale della Dark Matter per la temperatura corretta deve essere nulla.

        - *Numero totale di particelle di Dark Matter all'interno della stella*

- [epsi_Dark_Matter.f90](../../epsi_Dark_Matter.f90)

    - **epsi_DM_routine(T_DM,epsi_DM_tot)**: Questa routine si occupa di calcolare l'array su tutta la struttura del'$\epsilon_{DM}$ data una temperatura della Dark Matter

    - **convergenza_epsi_DM(epsi_DM_tot,T_DM,NMD,tempo)**: Questa routine, sfruttando il metodo di bisezione e la routine precedente trova la temepratura della Dark Matter ottimale, infatti la luminosità della Dark Matter totale deve essere nulla. Dato che non è possibile trovare uno 0 preciso si considera nulla la luminosità quando è minore di una parte su $10^4$ di quella totale; in caso tale condizione non venga raggiunta dopo un numero di iterazioni (solitamente impostato a 20, che impone una precisone di circa $10^{-6}$ sulla temperatura della dark matter) si scrive nel file **DarkMatterERROR.DAT** quale modello ha dato problemi, e se la flag *ioDarkError_on_off* è impostata a 1 viene scritta l'epsi_DM mesh per mesh e altre informazioni mesh per mesh di tale modello. Altro errore che viene scritto in tale file è se la luminosità totale per la temperatura minima e massima all'interno della struttura non sono di segno discordi, infatti tale condizione è necessaria per il funzionamento del metodo di bisezione. 


