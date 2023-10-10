# Documentazione FRANEC
Clicca per tornare al [README](../README.md)


## Indice

- [Gestione della compilazione del codice](#gestione-della-compilazione-del-codice)

- [Suggerimenti per il test del codice](#suggerimenti-per-il-test-del-codice)

- [Cose da fare](#cose-da-fare)

- [Struttura FRANEC (franec.f90)](Franec_standard/Franec_standard.md)

    - [subroutines](Franec_standard/franec_subroutines/franec_indice_subroutines.md)

    - [moduli](Franec_standard/moduli.md)

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
Per testare il codice, dopo la compilazione, si possono utilizzare 3 prinicpali metodi partendo dal main folder **FRANEC-net**:
    
- **./franec.sh**, legge il Modstart.in (*è stato creato un link simbolico tra il **Modstrat.in in DMsrc.90/Modstart_Backup** e quello nel main folder*) nel main folder e lancia i processi collegati alla shell, ciò permette di:
    - leggere i *write* a schermo
    - stoppare il programma con **stop**

- **./lancia.sh**, questo programma mi permette di lanciare i calcoli sempre collegando i processi alla shell, ma a differenza di *./franec.sh* non si ferma con *stop*. Questo succede perché lancia gestisce i vari problemi delle ripartenze.

- **./Driver-run "n-processi-in parallelo" lancia.sh tools/driver/out**, Il metodo classico di lanciare il programma.
    - Legge il file [**prove-pre-driver.dat**](./../../prove-per-driver.dat) dove ha la lista dei vari modstart da lanciare.
    - Legge i modstart in [tools/driver/start](./../../tools/driver/start) elencati in pove-pre-driver e ne genera uno shedule per lanciarli in parallelo.
        - in fase di test è sato modificato **M1.00_Z0.01760_He0.2734_ML1.90** ([vai a](./../../tools/driver/start/M1.00_Z0.01760_He0.2734_ML1.90)), questo file è stato linkato tramite link simbolico tra Modstart_Backup/M1.00_Z0.01760_He0.2734_ML1.90 e tools/driver/start/M1.00_Z0.01760_He0.2734_ML1.90
    - Lancia i vari *./lancia.sh* e registra i vari file di output nella cartella [*tools/driver/out*](./../../tools/driver/out) e quella lincata simbolicamente ad essa.

    Il vantaggio di lanciare il prgramma così è che i processi non sono legati alla shell, quindi è possibile sologgare dalla macchina senza interrompere i processi.

**Attenzione**:

- Il programma legge i parametri iniziali da **parametri.in** nel main folder, questo però è linkato simbolicamente con il file in *DMsrc.90/File_aggiuntivi/parametri.in*, quindi va cambiato questo.
## Cose da fare
- **[X]** Sistemare M_mesh e G che da problemi al bordo MAXME+1 è diverso
- **[X]** cambiare variabile *modello* con *NMD* in epsi_DM.f90
- **[]** Sistemare convergenza
- **[X]** LIM com MAXME<---------Controllali
- **[]** Commentare
- **[X]** Passare N_DM a ogni ripartenza del lancia tramite modifica di Modstart ripartenza
- **[X]** Fare un print nel file->Subroutine
    - **[]** Cattura
    - **[]** epsi
- **[]** aggiungere epsi come con matteo
- **[]** Guardare il fine mesh
- **[]** Mass_DM.in rho_DM.in
- **[]** Scrivere nel file run.log l'on_off_DM si fa tramite la subroutine in io.f90
- **[]** Mettere apposto l'elemento NEU (massa?)
- **[X]** Fare il file out di cattura
- **[X]** Sistemare il file Dark_Matter.DAT
- **[X]** Sistemare il rapporto epsi_DM/epsi_tot, deve venire uguale su i due file


