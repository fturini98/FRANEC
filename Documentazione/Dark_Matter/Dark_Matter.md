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

## Modifiche rispetto FRANEC  standard per includere DM
Come prima cosa è stato creato il file:

- **Cattura_DM.f90**: Questo file contiene la routine per calcolarsi la cattura della DM da parte della stella.

Modifica dei file per includere le routine:

- **Makefile**: 
    - è stato aggiunto nella sezione OBJ il file:
        - **Cattura_DM.o**, che legge e interpreta il file catturaDM.f90 e genera il file catturaDM.o da usare nel programma.

- **franec.f90**:
    - è stata aggiunta dopo la routine __ciaccioLi__ la chiamata all routine **Cattura_DM** 

## Subroutine per la DM:

- [Cattura_DM](Subroutine/Cattura_DM.md) : La routine cattura_DM gestisce la cattura della DM da parte della stella. (**Si trova nel file [Cattura_DM.f90](../../Cattura_DM.f90)**)   
