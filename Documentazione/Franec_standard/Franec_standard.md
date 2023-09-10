# Struttura FRANEC

Clicca per tornare a:
- [README](../../README.md)
- [Documentazione FRANEC](../Documentazione.md)

#

# Main (STELEV)

Il programma base di FRANEC, in cui non si considera la presenza di WIMPs, è strutturato a partire dal main che si trova nell file [franec.f90](../../franec.f90), utilizza delle variabili globali definite in [moduli.f90](../../moduli.f90), e chiama diverse subrutine definite all'interno dei vari file **f90** aggiuntivi.

Per approfondire  le variabili definite in [moduli](moduli.md) e le varie [subrutines](franec_subroutines/franec_indice_subroutines.md) guarda le sezioni della documentazione relative ad esse.

## Funzionamento

Dopo aver importato i vari moduli e definito le varie variabili il programma chiama:

**inout**: 
tale subrutine apre i file e gli assegna dei "puntatori", per una lista copleta vedi [qui](franec_subroutines/inout.md).

Successivamente viene letto il **Modstart.in**
- in particolare viene letta la variabile **IREAD**, che corrisponde al punto di partenza della simulazione.
    
    - **IREAD=1** Modello precedente (ripartenza)
    - **IREAD=2** Presequenza
    - **IREAD=3** Sequenza principale
    - **IREAD=4** Horizontal Branch

A dipendenza della scelta vengono generate le condizioni di start del modello utilizzando o il vecchio status per il *Modello precedente* o chiamando la subroutine [INNES](franec_subroutines/INNES.md).

Vinene aperta il file **scenario.log** e registrate le varie variabili da cui parte l'evoluzione, tra cui:

- La massa della stella in masse solari
- L'abondanza di idrogeno
- La metallicità
- Alfa
- la mistura

A questo punto si inizia con l'evoluzione della stella nel **mainloop**.

### Mainloop
#
Inanzi tutto vi sono due possibilità a dipendenza se siamo nel primo giro di iterazione, quindi se siamo al primo passo temporale, oppure no. Se non siamo al primo passo temporale ma siamo al primo tentativo di ottenere la convergenza vengono chiamate:

- **MIXING**, che si occupa di mescolare gli elementi all'interno delle varie zone dove è presente la convezione.

- **EQLB**

Dopo di che, indipendemente dal numero di tentativi di convergenza, si cerca di trovare le soluzioni della struttura per tramite il metodo di *Henyey*. Vengono chiamate quindi:

- **QUATM**, che gestisce l'atmosfera (per ottenere le soluzioni al contorno?)

- **HENYEY**, che si occupa di trovare una soluzone mesh per mesh e dare una stima dell'errore che si effettua nel far evolvere la struttura.

A questo punto viene controllato se l'errore trovato da HENYEY è piccolo abbastanza, ciò dipende a quale iterazione per la convergenza siamo (variabile: *MAIS*), se è sufficente si continua, altrimenti si prova nuovamente con nuovi parametri.

Tutta questa parte per il passo temporale 1 non viene fatto in quanto le condizioni di partenza sono già state trovate mediante INNES.

Successivamente si passa a un rilassamento della subatmosfera se necessario (call **rilassamento_subatm**)

Si salvano i dati mediante la subroutine **STAMPA**.

A questo punto viene calcolato il passo temporale successivo (HT1) tramite la funzione **PASTEM**.

Dopo di che vengono chiamate le subroutine:

- **OPTIM**
- **zone_convettive**
- **ciacioLi** che calcola la diffuzione degli elementi chimici
- **EVOLUT**  al suo interno chiama la
    
    - **STATE**
    - **epsi**

- **MASLOS** 

Quindi a questo punto si riparte con il loop passando al passo temporale successivo.





