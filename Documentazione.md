# Franec DM
## Indice

- [Gestione della compilazione del codice](#gestione-della-compilazione-del-codice)

- [Modifiche rispetto FRANCEC standard](#modifiche-rispetto-franec-standard)

- [File catturaDM.f90](#catturadm)

## Gestione della compilazione del codice
Per evitare di danneggiare il codice standard è stato fatto un beckup con la cartella *src90-std*, dopop di che è stata creata la cartella *DMsrc90* dove verranno effettuate le modifiche del codice.

### Compilazione
Per la compilazione è stato effettuato un link simbolico tramite il comando:

    ln -s DMsrc.90 src.90

Il che permette di compilare il codice mediante il comando:

    make SM=AS09a0
Se si vuole tornare a usare il codice standard basta modificare il link simbolico alla cartella standard mediante:

    rm -r src.90
    ln -s src.90-std src.90
## Modifiche rispetto FRANEC  standard
Come prima cosa è stato creato il file:

- **CatturaDM.f90**: Questo file contiene la routine per calcolarsi la cattura della DM da parte della stella.

Modifica dei file per includere le routine:

- **Makefile**: 
    - è stato aggiunto nella sezione OBJ il file:
        - **catturaDM.o**, che legge e interpreta il file catturaDM.f90 e genera il file catturaDM.o da usare nel programma.

- **franec.f90**:
    - è stata aggiunta dopo la routine __ciaccioLi__ la chiamata all routine **catturaDM** (per ora è disabilitata mediante commento)

## CatturaDM
La routine catturaDM gestisce la cattura della DM da parte della stella.    
    
    