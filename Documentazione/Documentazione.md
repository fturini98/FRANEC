# Documentazione FRANEC
Clicca per tornare al [README](../README.md)


## Indice

- [Gestione della compilazione del codice](#gestione-della-compilazione-del-codice)

- [Suggerimenti per il test del codice](#suggerimenti-per-il-test-del-codice)

- [Sezione relativa alla Dark Matter](Dark_Matter/Dark_Matter.md)

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
### In caso di errori di compilazione:
In generale i problemi di compilazione sono dati da un errata scrittura del codice sorgente, ma ogni tanto si possono corrompere i vari file, si può risolvere tale problema con:
    
    make clean
e poi ricompilando.

## Suggerimenti per il test del codice
Per testare il codice, dopo la compilazione, basta posizionarsi nel main folder **FRANEC-net** e lanciare il programma tramite

    ./lancia.sh
Tale procedura permetterà di visualizzare i vari *write* sulla shell. La normale procedura se no scriverebbe gli output su dei file di log.

