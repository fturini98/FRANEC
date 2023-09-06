# FRANEC
Codice sorgente FRANEC modificato per considerare la cattura di WIMP da parte della stella.

- esite un file [Documentazione.md](Documentazione/Documentazione.md) che contiene una sottospecie di documentazione per le modifiche al codice fatte riguardanti la DM.

## Git
I comandi git funzionano solo sul server, è possibile clonare la cartella in locale, ma una volta fatto non è possibile più caricarare le modifiche nuovamente sul server a causa dei firewall (motivazione plausibile data da Matteo).

### GitHub 
Esiste una repository privata (utente: fturini98) su GitHub per tenere meglio traccia dei cambiamenti effettuati.
>Per uplodare i cambiamenti su GitHub:
    
    git push github
>Per vedere quali altri link remoti sono presenti:
    
    git remote -v
>Per aggiungere altri link remoti:
    
    git remote add "nome" "link"

### Salvare e fare i commit
Anche usando l'estensione ssh di VSCode **è necessario** chiamare:

    git add .
ma vengono aggiunti tutti i file salvati con *Ctrl+S*.
Successivamnete è sufficente chiamare:

    git commit -m "(descrizione commit)"

# Python3
 Per far funzionare le estensioni di Vscode sul server serve python3 con le librerie aggiuntive:

    
    pip install fortls
    pip install fortran-language-server