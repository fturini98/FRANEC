# inout
Clicca per tornare a:
- [README](../../../README.md)
- [Documentazione FRANEC](../../Documentazione.md)
- [Struttura_FRANEC](../Franec_standard.md)
- [Subroutines](franec_indice_subroutines.md)
#
Tale subrutine si occupa di aprire i file necessari, sia per lettura che scrittura e assegna loro un "puntatore".

I file caricati sono:

- Modstart.in(unit=4,**lettura**)
- TABEL.in (unit=25,**lettura**)
- EOS.in (unit=8,**lettura**)
- EOSTOTALE.in (unit=5,**lettura**)
- NUCL.DAT (unit=39,**scrittura**)
- OUT.DAT (unit=ioOUT,**scrittura**)
- BIGTAB.DAT (unit=ioBIGTAB,**scrittura**)
- CHIMICA.DAT (unit=ioCHIMICA,**scrittura**)
- FISICA.DAT (unit=ioFISICA,**scrittura**)
- SUBAATM.DAT (unit=23,**scrittura**)
- CHUMICA_SUP.DAT (ioCHIMICASUP)
- abbondanze.in (unit=81, **lettura**)

    - Questa viene chiusa dopo aver letto da tale file i valori per :
        
        - tipo di mistura

        - Abbondanze: He,C,N,O,Fe,Li6,Li7,Be,B,D

- misture.in (unit=81,**lettura**) 
    - ! attenzione unit=81 prima era le abbondanze.
- PRINT.DAT (unit=2,**scrittura**)
- GRAFI.DAT (unit=10,**scrittura**)
- SEGNALI.DAT (unit=15,**scrittura**)
- fasievolutive.log (unit=68,**scrittura**)
- mixcno.dat (unit=18,**scrittura**)
- error.log (unit=66,**scrittura**)
- run.log (unit=iolog,**scrittura**)
- PROFILI_H.DAT (unit=70,**scrittura**)
- PROFILI_HE.DAT (unit=71,**scrittura**)
- PROFILI_C.DAT (unit=72,**scrittura**)
- PROFILI_N.DAT (unit=73,**scrittura**)
- PROFILI_O.DAT (unit=74,**scrittura**)
- CHECHK.DAT (unit=50,**scrittura**)
- parametri.in (unit=80,**lettura**)
- OVERSHOOT-CENTRAL.DAT (unit=85,**scrittura**)
    - solo se è presente l'overshoot

