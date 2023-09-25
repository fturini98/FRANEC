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
- Come prima cosa è stato creato il file:

    - **Cattura_DM.f90**: Questo file contiene la routine per calcolarsi la cattura della DM da parte della stella.

- Modifica dei file per includere le routine:

    - **Makefile**: 
        - è stato aggiunto nella sezione OBJ il file:
            - **Cattura_DM.o**, che legge e interpreta il file catturaDM.f90 e genera il file catturaDM.o da usare nel programma.

    - **franec.f90**:
        - è stata aggiunta dopo la routine __ciaccioLi__ la chiamata all routine **Cattura_DM** 

- Sono stati modificati i programmi
    - [**lancia.c**](./../../Gestione_lancia/lancia.c)

    - [**driver.c**](./../../Gestione_driver/driver.c)

    per poter genereare dei file in output (*DarkMatter.DAT* e *DarMatterERROR.DAT*) dove registrare varie variabili tramite la subroutine [stampa_epsi_DM](Subroutine/epsi_Dark_Matter.md) (definita in [epsi_Dark_Matter.f90](./../../epsi_Dark_Matter.f90)). Per poter generare questi file di putput è stato necessario creare nelle cartelle tools/driver e tools/lancia i rispettivi link simbolici tramite:
        
        ln -s ./../../DMsrc.90/Gestione_driver/driver.c driver.c

        ln -s ./../../DMsrc.90/Gestione_lancia/lancia.c lancia.c

    *Attenzione:* ogni volta che vengono fatte delle modifiche al file sorgente bisogna andare nelle cartelle tools/driver e tools/lancia e lanciare il make per generare i rispettivi programmi.


## Subroutine per la DM:

- [Cattura_DM](Subroutine/Cattura_DM.md) : La routine cattura_DM gestisce la cattura della DM da parte della stella. (**Si trova nel file [Cattura_DM.f90](../../Cattura_DM.f90)**)   



