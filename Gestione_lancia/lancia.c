#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct 
{
  int diff;
  int atm_bh05, atm_ck03, overshoot;
  double div_pastem;
  double max_pastem, par_OVER, maxAGE;
  double URA[4],  ULA[4], UPA[4], UTA[4], UMA[4];
  double VMM[5];
  int stdchim;
} PARAM;

void leggi_parametri(PARAM *st);

void scrivi_parametri(PARAM st, FILE *fz, int cnt);

int  mod_parametri(PARAM *st, int cod_err, int cnt, int *rep, FILE *fz);

void unisci(int cnt, FILE *fz);

void salvamodello(int cnt);

void make_modstart_pepper();

#define MAXITER 511

int nummod[MAXITER];
double pastem[MAXITER];
static int corto = 1; // disattiva la modifica di UMA se corto = 1

main(int argc, char **argv)
{
  FILE *fp, *fz, *ft;
  char com1[256], com2[256], com3[256];
  char eseguibile[256], dummy[256];
  int diff, cod_err, cnt, mod, parte1, parte2, retry, i, numarg;
  PARAM P;
  
  int MAXCNT = 15;  // massimo numero di lanci del franec

  if(argc != 2) 
    {
      printf("Usage: %s <eseguibile_franec>\n", argv[0]);
      exit(-1);
    }

  retry = 1;
  // leggo da riga di comando il nome dell'eseguibile del franec
  sprintf(eseguibile, "%s", argv[1]);   

  // comandi di lancio FRANEC
  sprintf(com3, "./cancella; ./%s", eseguibile);  // serve per il primo giro
  sprintf(com1, "cp MOD1.DAT Modstart.in; ./cancella; ./%s", eseguibile);
  sprintf(com2, "cp MOD2.DAT Modstart.in; ./cancella; ./%s", eseguibile);

  // copia di backup file di parametri e pulizia
  system("cp parametri.in parametri.in-bak; cp Modstart.in Modstart.in-bak");
  system("rm BIGTAB.DAT_* OUT.DAT_* MOD1.DAT_* MOD2.DAT_*  DarkMAtter.DAT_* DarkMatterERROR.DAT_* >& /dev/null");
  system("rm PRINT.DAT_* GRAFI.DAT_* >& /dev/null");

  // apertura file di log per controllare i passi fatti in caso di errore
  // del franec
  fz = fopen("lancia.log", "w");  

  // lettura parametri dinamici
  leggi_parametri(&P);

  cnt = 1;
  retry = 1;
  scrivi_parametri(P, fz, cnt);
  do
    {
      // lancio FRANEC
      if(cnt == 1) 
	{ // prima esecuzione
	  system(com3);
	  salvamodello(cnt);
	}
      else
	{ // dalla seconda esecuzione in poi
	  if(mod > 0)
	    {  // il modello e' buono, riparto dal MOD<1,2> piu' vecchio
	      retry = 1; // resetto il contatore di tentativi per questo start
	      parte1 = parte2 = -1;
	      if((ft = fopen("MOD1.DAT", "r"))  != NULL)
		{
		  for(i=0;i<20;i++)
		    fgets(dummy,255,ft);
		  fscanf(ft,"%d", &parte1);
		  fclose(ft);
		}

	      if((ft = fopen("MOD2.DAT", "r")) != NULL)
		{
		  for(i=0;i<20;i++)
		    fgets(dummy,255,ft);
		  fscanf(ft,"%d", &parte2);
		  fclose(ft);
		  fprintf(fz,"partenze: MOD1 %d  MOD2 %d\n", parte1,parte2);
		  fflush(fz);
		}

	      if((parte2 < 0 && parte1 > 0) || (parte1 < parte2 && parte1 > 0))
		{ // MOD1 e' piu' vecchio, uso questo
		  system(com1);
		  salvamodello(cnt);
		}
	      else if(parte2 > 0)
		{ // MOD2 e' piu' vecchio, uso questo
			  system(com2);
		  salvamodello(cnt);
		}
	      else
		{
		  // qui gestisco i casi in cui non ho i MOD1/2
		  // perche' la simulazione muore precocemente
		  system(com3);
		  salvamodello(cnt);
		}
	    }
	  else
	    { // il modello non e' buono, riparto con il vecchio modstart
	      system(com3);
	      salvamodello(cnt);
	    }
	}
      
      // controllo il codice di uscita del franec
      fp = fopen("error.log", "r");
      numarg = fscanf(fp,"%d", &cod_err);
      fclose(fp);

      if(numarg <= 0) // probabilmente interrotto franec con CTRL+c
	cod_err = 1;

      if(cod_err > 0) 
	{
	  // il codice esce con errore...
	  // ================================================================
	  // ATTENZIONE: modificare i valori sotto se cambio il set di regole
	  // ================================================================
	  // 24 = numero di regole che conosco in mod_parametri + 1 
	  if((corto == 0 && retry == 24) || (corto == 1 && retry == 17)) 
	    {
	      // sono in un loop che non so trattare. Esco con errore
	      printf("ERRORE: sono in un loop da cui non esco.\n");
	      fprintf(fz,"ERRORE: sono in un loop da cui non esco.\n");
	      unisci(cnt-1, fz);
	      system("cp parametri.in-bak parametri.in");
	      system("cp Modstart.in-bak Modstart.in");
	      exit(-1);
	    }
	  // modifiche ai parametri in input
	  mod = mod_parametri(&P, cod_err, cnt, &retry, fz);
	 
	  if(mod < 0)  // scarto il run
	    // e aumento il contatore di tentativi con lo stesso start
	    {
	      cnt--;
	      retry ++;
	    }

	  cnt++; // aumento il contatore di lancio franec

	  // troppe iterazioni... esecuzione sospesa
	  if(cnt > MAXCNT) 
	    {
	      printf("Numero di iterazioni massime raggiunto.\n\n");
	      fprintf(fz, "Troppe iterazioni. Modello sospeso.\n");
	      fflush(fz);
	      unisci(cnt-1, fz);
	      // ripristino i parametri di partenza
	      system("cp parametri.in-bak parametri.in");
	      system("cp Modstart.in-bak Modstart.in");
	      exit(-1);
	    }
	  
	  // scrittura parametri modificati
	  scrivi_parametri(P, fz, cnt);
	}
      
    } while(cod_err > 0);

  unisci(cnt, fz);
  // ripristino i parametri di partenza
  system("cp parametri.in-bak parametri.in");
  system("cp Modstart.in-bak Modstart.in");

  // Il flash di He ha interrotto la simulazione.
  // Creazione Modstart.in-pepper per pepper...
  if(cod_err < 0)
    make_modstart_pepper();
 } 


void make_modstart_pepper()
{
  // Attenzione! Funziona se lanciato dalla dir radice.
  // Se lanciato mediante driver viene aggiustato automaticamente il path.

  // la copia non dovrebbe servire (viene fatta da Creamod-pepper) 
  // ma e' un sistema di sicurezza dato che ci sono comportamenti non 
  // ancora debuggati
  system("mkdir -p safe; cp BIGTAB.DAT OUT.DAT GRAFI.DAT PRINT.DAT* FISICA.DAT CHIMICA.DAT CHIMICA_SUP.DAT mixcno.dat scenario.log run.log error.log fasievolutive.log safe DarkMatter.DAT DarkMatterERROR.DAT >& /dev/null; cp Modstart.in-bak Modstart.in-preflash");

  if(fopen("run-by-driver","r") != NULL)
    // lanciato mediante driver
      system("../../tools/per_pepper/Creamod-pepper -m");
  else
    system("tools/per_pepper/Creamod-pepper -m");
  
  // ripristino gli out
  system("mv safe/* .; rmdir safe");
}

/* ========================================================= */

/* Questa routine mette da parte gli out di ogni run, etichettandoli
   con un indice progressivo  */
void salvamodello(int cnt)
{
  char com[256];
  
  sprintf(com,"cp BIGTAB.DAT BIGTAB.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp OUT.DAT OUT.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp MOD1.DAT MOD1.DAT_%d >& /dev/null", cnt);
  system(com);
  sprintf(com,"cp MOD2.DAT MOD2.DAT_%d  >& /dev/null", cnt);
  system(com);
  sprintf(com,"cp GRAFI.DAT GRAFI.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp PRINT.DAT PRINT.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp CHIMICA.DAT CHIMICA.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp FISICA.DAT FISICA.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp SUBATM.DAT SUBATM.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp CHIMICA_SUP.DAT CHIMICA_SUP.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp run.log run.log_%d", cnt);
  system(com);
  sprintf(com,"cp error.log error.log_%d", cnt);
  system(com);
  sprintf(com,"cp mixcno.dat mixcno.dat_%d", cnt);
  system(com);
  sprintf(com,"cp fasievolutive.log fasievolutive.log_%d", cnt);
  system(com);
  sprintf(com,"cp NUCL.DAT NUCL.DAT_%d", cnt);
  system(com);
  sprintf(com,"cp DarkMatter.DAT DarkMatter.DAT_%d",cnt);
  system(com);
  sprintf(com,"cp DarkMatterERROR.DAT DarkMatterERROR.DAT_%d",cnt);
  system(com);
}

/* Routine finale, che "cuce" gli output prodotti dai vari step, fino a
   generare i file finali */
void unisci(int cnt, FILE *fz)
{
  int partenze_run[MAXITER];
  int i, f, m;
  FILE *fp, *fpp, *fpp1;
  char filename[256], comando[256], dummy[256];

  fprintf(fz,"\n ====== MODELLI DI AVVIO ======\n");
  for(i=0;i<cnt;i++)
    {  // legge i modelli di avvio dei vari run
      sprintf(filename,"BIGTAB.DAT_%d", i+1);
      fp = fopen(filename, "r");
      fgets(dummy, 255, fp);
      fscanf(fp,"%d", &partenze_run[i]);

      fprintf(fz,"run %d: %d\n", i+1, partenze_run[i]);
      fclose(fp);
    }
  fflush(fz);
  // aggiusto l'offset del primo valore 
  partenze_run[0]--;
  
  // unione dei BIGTAB
  unlink("tmpBIGTAB.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d BIGTAB.DAT_%d >> tmpBIGTAB.DAT",
		    partenze_run[i+1]-partenze_run[i], i+1); 
	  else
	    {
	      sprintf(comando,
		      "grep -v \"#\" BIGTAB.DAT_%d | head -%d>> tmpBIGTAB.DAT",
		      i+1, partenze_run[i+1]-partenze_run[i]); 
	    }
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" BIGTAB.DAT_%d >> tmpBIGTAB.DAT", i+1);
      else 
	sprintf(comando,"cp BIGTAB.DAT_1 tmpBIGTAB.DAT");
      system(comando);
    }

  // unione degli OUT
  unlink("tmpOUT.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d OUT.DAT_%d >> tmpOUT.DAT",
		    partenze_run[i+1]-partenze_run[i], i+1); 
	  else
	    sprintf(comando,
		    "grep -v \"#\" OUT.DAT_%d | head -%d>> tmpOUT.DAT",
		    i+1, partenze_run[i+1]-partenze_run[i]); 
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" OUT.DAT_%d >> tmpOUT.DAT", i+1);
      else
	sprintf(comando,"cp OUT.DAT_1 tmpOUT.DAT");
      system(comando);
    }

  // unione di CHIMICA_SUP
  unlink("tmpCS.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d CHIMICA_SUP.DAT_%d >> tmpCS.DAT",
		    partenze_run[i+1]-partenze_run[i], i+1); 
	  else
	    sprintf(comando,
		    "grep -v \"#\" CHIMICA_SUP.DAT_%d | head -%d>> tmpCS.DAT",
		    i+1, partenze_run[i+1]-partenze_run[i]); 
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" CHIMICA_SUP.DAT_%d >> tmpCS.DAT", i+1);
      else
	sprintf(comando,"cp CHIMICA_SUP.DAT_1 tmpCS.DAT");
      system(comando);
    }

  // unione di NUCL
  unlink("tmpNC.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d NUCL.DAT_%d >> tmpNC.DAT",
		    partenze_run[i+1]-partenze_run[i]-1, i+1); 
	  else
	    sprintf(comando,
		    "grep -v \"#\" NUCL.DAT_%d | head -%d>> tmpNC.DAT",
		    i+1, partenze_run[i+1]-partenze_run[i]); 
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" NUCL.DAT_%d >> tmpNC.DAT", i+1);
      else
	sprintf(comando,"cp NUCL.DAT_1 tmpNC.DAT");
      system(comando);
    }

  // unione di CHIMICA
  unlink("tmpC.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	sprintf(comando,"cp CHIMICA.DAT_%d CHIMICA.DAT; ./Estrai -c %d %d >> tmpC.DAT",
		i+1, partenze_run[i], partenze_run[i+1]-1); 
      else if(cnt>1)
	sprintf(comando,"cp CHIMICA.DAT_%d CHIMICA.DAT; ./Estrai -c %d %d >> tmpC.DAT",
		i+1, partenze_run[i], 100000);
      else
	sprintf(comando,"cp CHIMICA.DAT_1 tmpC.DAT");
      system(comando);
    }

 // unione di FISICA
  unlink("tmpF.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	sprintf(comando,"cp FISICA.DAT_%d FISICA.DAT; ./Estrai -f %d %d >> tmpF.DAT",
		i+1, partenze_run[i], partenze_run[i+1]-1); 
      else if(cnt>1)
	sprintf(comando,"cp FISICA.DAT_%d FISICA.DAT; ./Estrai -f %d %d >> tmpF.DAT",
		i+1, partenze_run[i], 100000);
      else
	sprintf(comando,"cp FISICA.DAT_1 tmpF.DAT");
      system(comando);
    }

  // unione di SUBATM
  unlink("tmpS.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	sprintf(comando,"cp SUBATM.DAT_%d SUBATM.DAT; ./Estrai -a %d %d >> tmpS.DAT",
		i+1, partenze_run[i], partenze_run[i+1]-1); 
      else if(cnt>1)
	sprintf(comando,"cp SUBATM.DAT_%d SUBATM.DAT; ./Estrai -a %d %d >> tmpS.DAT",
		i+1, partenze_run[i], 100000);
      else
	sprintf(comando,"cp SUBATM.DAT_1 tmpS.DAT");
      system(comando);
    }

 // unione dei GRAFI
  unlink("tmpGRAFI.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d GRAFI.DAT_%d >> tmpGRAFI.DAT",
		    12*(partenze_run[i+1]-1-partenze_run[i]), i+1); 
	  else
	    sprintf(comando,
		    "head -%d GRAFI.DAT_%d >> tmpGRAFI.DAT",
		    12*(partenze_run[i+1]-partenze_run[i]), i+1); 
	}
      else if(cnt>1)
	sprintf(comando,"cat GRAFI.DAT_%d >> tmpGRAFI.DAT", i+1);
      else
	sprintf(comando,"cp GRAFI.DAT_1 tmpGRAFI.DAT");
      system(comando);
    }

  // unione dei mixcno.dat
  unlink("tmpmixcno.dat");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d mixcno.dat_%d >> tmpmixcno.dat",
		    partenze_run[i+1]+1-partenze_run[i], i+1);
	  else
	    sprintf(comando,
		    "grep -v \"#\" mixcno.dat_%d | head -%d>> tmpmixcno.dat",
		    i+1, partenze_run[i+1]-partenze_run[i]); 
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" mixcno.dat_%d >> tmpmixcno.dat", i+1);
      else
	sprintf(comando,"cp mixcno.dat_1 tmpmixcno.dat");
      system(comando);
    }

// unione dei fasievolutive.dat
  unlink("fasievolutive.log");
  fpp = fopen("tmpfasievolutive.log", "w");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  sprintf(filename,"fasievolutive.log_%d", i+1);
	  fpp1 = fopen(filename, "r");
	  while(fscanf(fpp1, "%d %d", &f, &m) != EOF)
	    {
	      if(m < partenze_run[i+1])
		fprintf(fpp,"%d %d\n", f,m);
	    }
	  fclose(fpp1);
	}
      else if(i==cnt-1)
	{
	  sprintf(filename,"fasievolutive.log_%d", i+1);
	  fpp1 = fopen(filename, "r");
	  while(fscanf(fpp1, "%d %d", &f, &m) != EOF)
	    {
	      fprintf(fpp,"%d %d\n", f,m);
	    }
	  fclose(fpp1);
	}
    }
  fclose(fpp);

   // unione degli DarkMatter
  unlink("tmpDarkMatter.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d DarkMatter.DAT_%d >> tmpDarkMatter.DAT",
		    partenze_run[i+1]-partenze_run[i], i+1); 
	  else
	    sprintf(comando,
		    "grep -v \"#\" DarkMatter.DAT_%d | head -%d>> tmpDarkMatter.DAT",
		    i+1, partenze_run[i+1]-partenze_run[i]); 
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" DarkMatter.DAT_%d >> tmpDarkMatter.DAT", i+1);
      else
	sprintf(comando,"cp DarkMatter.DAT_1 tmpDarkMAtter.DAT");
      system(comando);
    }

     // unione degli DarkMatterERROR
  unlink("tmpDarkMatterERROR.DAT");
  for(i=0;i<cnt;i++)
    {
      if(i<cnt-1)
	{
	  if(i==0) 
	    sprintf(comando,"head -%d DarkMatterERROR.DAT_%d >> tmpDarkMatterERROR.DAT",
		    partenze_run[i+1]-partenze_run[i], i+1); 
	  else
	    sprintf(comando,
		    "grep -v \"#\" DarkMatterERROR.DAT_%d | head -%d>> tmpDarkMatterERROR.DAT",
		    i+1, partenze_run[i+1]-partenze_run[i]); 
	}
      else if(cnt>1)
	sprintf(comando,"grep -v \"#\" DarkMatterERROR.DAT_%d >> tmpDarkMatterERROR.DAT", i+1);
      else
	sprintf(comando,"cp DarkMatterERROR.DAT_1 tmpDarkMAtterERROR.DAT");
      system(comando);
    }
  
  // unione di run.log
  unlink("run.log");
  unlink("error.log");
  for(i=0;i<cnt;i++)
    {
      sprintf(comando, 
	      "echo \"@@@@@@@@@@@@@@ %2d @@@@@@@@@@@@@@@\" >> run.log", i+1);
      system(comando);                     
      sprintf(comando,"cat run.log_%d >> run.log", i+1);
      system(comando);

      sprintf(comando, 
	      "echo \"@@@@@@@@@@@@@@ %2d @@@@@@@@@@@@@@@\" >> error.log", i+1);
      system(comando);
      sprintf(comando,"cat error.log_%d >> error.log", i+1);
      system(comando);
    }

  system("mv tmpOUT.DAT OUT.DAT; mv tmpBIGTAB.DAT BIGTAB.DAT; mv tmpGRAFI.DAT GRAFI.DAT; mv tmpmixcno.dat mixcno.dat; mv tmpDarkMAtter.DAT DarkMatter.DAT;mv tmpDarkMAtterERROR.DAT DarkMatterERROR.DAT;");
  system("mv tmpCS.DAT CHIMICA_SUP.DAT; mv tmpC.DAT CHIMICA.DAT; mv tmpF.DAT FISICA.DAT; mv tmpS.DAT SUBATM.DAT; mv tmpfasievolutive.log fasievolutive.log");
  system("mv tmpNC.DAT NUCL.DAT");
}

/* Questa routine legge da file parametri.in il valore dei parametri
   dinamici   */
void leggi_parametri(PARAM *p)
{
  FILE *fp;
  char linea[256];

  fp = fopen("parametri.in", "r");
  fgets(linea, 255, fp);
  sscanf(linea, "%d", &(p->diff));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf", &(p->div_pastem));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf", &(p->max_pastem));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf %lf %lf %lf", 
	 &(p->URA[0]), &(p->URA[1]), &(p->URA[2]), &(p->URA[3]));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf %lf %lf %lf", 
	 &(p->ULA[0]), &(p->ULA[1]), &(p->ULA[2]), &(p->ULA[3]));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf %lf %lf %lf", 
	 &(p->UPA[0]), &(p->UPA[1]), &(p->UPA[2]), &(p->UPA[3]));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf %lf %lf %lf", 
	 &(p->UTA[0]), &(p->UTA[1]), &(p->UTA[2]), &(p->UTA[3]));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf %lf %lf %lf", 
	 &(p->UMA[0]), &(p->UMA[1]), &(p->UMA[2]), &(p->UMA[3]));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf %lf %lf %lf %lf", 
	 &(p->VMM[0]), &(p->VMM[1]), &(p->VMM[2]), &(p->VMM[3]), &(p->VMM[4]));
  fgets(linea, 255, fp);
  sscanf(linea, "%d", &(p->atm_bh05));
  fgets(linea, 255, fp);
  sscanf(linea, "%d", &(p->atm_ck03));
  fgets(linea, 255, fp);
  sscanf(linea, "%d", &(p->overshoot));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf", &(p->par_OVER));
  fgets(linea, 255, fp);
  sscanf(linea, "%lf", &(p->maxAGE));
  fgets(linea, 255, fp);
  sscanf(linea, "%d", &(p->stdchim));

  fclose(fp);
 }

/* Questa routine scrive i valori dei parametri correnti,
   sia su file parametri.in, sia su file lancia.log       */
void scrivi_parametri(PARAM p, FILE *fz, int cnt)
{
  FILE *fp;
  
  fp = fopen("parametri.in", "w");

  fprintf(fp,"%d\t\t\t\t! diffusione: (1/0) = (on/off)\n", p.diff);
  fprintf(fp,"%.4f\t\t\t\t! divisore passo temporale\n", p.div_pastem);
  fprintf(fp,"%.4e\t\t\t! massimo valore passo temporale\n", p.max_pastem);
  fprintf(fp,"%.3e %.3e %.3e %.3e\t! URA (optim)\n", 
	  p.URA[0],p.URA[1],p.URA[2],p.URA[3]);
  fprintf(fp,"%.3e %.3e %.3e %.3e\t! ULA (optim)\n", 
	  p.ULA[0],p.ULA[1],p.ULA[2],p.ULA[3]);
  fprintf(fp,"%.3e %.3e %.3e %.3e\t! UPA (optim)\n", 
	  p.UPA[0],p.UPA[1],p.UPA[2],p.UPA[3]);
  fprintf(fp,"%.3e %.3e %.3e %.3e\t! UTA (optim)\n", 
	  p.UTA[0],p.UTA[1],p.UTA[2],p.UTA[3]);
  fprintf(fp,"%.3e %.3e %.3e %.3e\t! UMA (optim)\n", 
	  p.UMA[0],p.UMA[1],p.UMA[2],p.UMA[3]);
  fprintf(fp,"%.3f %.3f %.3f %.3f %.3f\t! VMM (optim)\n", 
	  p.VMM[0],p.VMM[1],p.VMM[2],p.VMM[3],p.VMM[4]);
  fprintf(fp,"%d\t\t\t\t! (1/0) = ON/OFF atmos. BH05 \n", p.atm_bh05);
  fprintf(fp,"%d\t\t\t\t! (1/0) = ON/OFF atmos. CK03 \n", p.atm_ck03);
  fprintf(fp,"%d\t\t\t\t! (1/0) = ON/OFF overshooting \n", p.overshoot);
  fprintf(fp,"%.4f\t\t\t\t! parametro di overshooting \n", p.par_OVER);
  fprintf(fp,"%.4f\t\t\t\t! max age in Gyr (< 0 per disabilitare)\n", p.maxAGE);
  fprintf(fp,"%d\t\t\t\t! STDCHIM: (1/0) = (on/off)\n", p.stdchim);  
  fclose(fp);

  fprintf(fz, "=========================================\n");
  fprintf(fz, "========== NUOVA ITERAZIONE === %.2d ======\n", cnt);
  fprintf(fz, "=========================================\n");
  fprintf(fz,"%d\t\t\t\t! diffusione: (1/0) = (on/off)\n", p.diff);
  fprintf(fz,"%.4f\t\t\t\t! divisore passo temporale\n", p.div_pastem);
  fprintf(fz,"%.4e\t\t\t! massimo valore passo temporale\n", p.max_pastem);
  fprintf(fz,"%.3e %.3e %.3e %.3e\t! URA (optim)\n", 
	  p.URA[0],p.URA[1],p.URA[2],p.URA[3]);
  fprintf(fz,"%.3e %.3e %.3e %.3e\t! ULA (optim)\n", 
	  p.ULA[0],p.ULA[1],p.ULA[2],p.ULA[3]);
  fprintf(fz,"%.3e %.3e %.3e %.3e\t! UPA (optim)\n", 
	  p.UPA[0],p.UPA[1],p.UPA[2],p.UPA[3]);
  fprintf(fz,"%.3e %.3e %.3e %.3e\t! UTA (optim)\n", 
	  p.UTA[0],p.UTA[1],p.UTA[2],p.UTA[3]);
  fprintf(fz,"%.3e %.3e %.3e %.3e\t! UMA (optim)\n", 
	  p.UMA[0],p.UMA[1],p.UMA[2],p.UMA[3]);
  fprintf(fz,"%.3f %.3f %.3f %.3f %.3f\t! VMM (optim)\n", 
	  p.VMM[0],p.VMM[1],p.VMM[2],p.VMM[3],p.VMM[4]);
  fprintf(fz,"%d\t\t\t\t! (1/0) = ON/OFF atmosfera BH05\n", p.atm_bh05);
  fprintf(fz,"%d\t\t\t\t! (1/0) = ON/OFF atmosfera CK03\n", p.atm_ck03);
  fprintf(fz,"%d\t\t\t\t! (1/0) = ON/OFF overshooting\n", p.overshoot);
  fprintf(fz,"%.4f\t\t\t\t! parametro di overshooting \n", p.par_OVER);
  fprintf(fz,"%.4f\t\t\t\t! max age in Gyr (< 0 per disabilitare)\n", p.maxAGE);
  fprintf(fz,"%d\t\t\t\t! STDCHIM: (1/0) = (on/off)\n", p.stdchim);
  fflush(fz);
}



/* =================================================================== 
   Questa routine e' responsabile delle modifiche ai parametri
   nei casi in cui il codice si blocca con errore.
   Come prima soluzione si prova a variare il passo temporale.
   Se la soluzione non funziona, si alterano i vettori della optim.
   =================================================================== */
int mod_parametri(PARAM *p, int cod_err, int cnt, int *retry, FILE *fz)
{
  double eps = 1e-3;
  int nmd, modifica = 1; 
  FILE *fp;

  system("tail -1 OUT.DAT | grep -v '#' > nmod");

  fp = fopen("nmod", "r");
  if( fscanf(fp,"%d", &nmd) == EOF )
    nmd = nummod[cnt-2];
  fclose(fp);

  fprintf(fz,"MODELLO n. %d\n", nmd);

  nummod[cnt-1] = nmd;
  pastem[cnt-1] = p->div_pastem;

  if(cnt > 1 && nmd - nummod[cnt-2] < 20) 
    modifica = -1;  // ho fatto troppo pochi passi, non considero il run

  fprintf(fz,"PASSI FATTI: %d\n", nmd - nummod[cnt-2]);

  if(modifica == 1) 
    *retry = 1;
  
  // il run precedente aveva passtem std, quindi lo devo cambiare
  if(*retry == 1 && fabs(p->div_pastem-1.0) < eps) 
    {
      //      p->UMA[0] = 0.01;
      *retry = 2;
    }

  if(*retry == 1)
    {
      //      p->UMA[0] = 0.01;
      p->div_pastem = 1.0;
    }
  else if(*retry == 2)
    p->div_pastem = 0.9090909;
  else if(*retry == 3)
    p->div_pastem = 0.8333333;
  else if(*retry == 4)
    p->div_pastem = 1.1;
  else if(*retry == 5)
    p->div_pastem = 1.2;
  else if(*retry == 6)
    p->div_pastem = 1.5;
  else if(*retry == 7)
    p->div_pastem = 0.6666667;
    // ci sono problemi, provo a diminuire parecchio il passo 
  // per uscire da una zona difficile da trattare
  else if(*retry == 8)
    p->div_pastem = 2.0;
  else if(*retry == 9)
    p->div_pastem = 3.0;
  else if(*retry == 10)
    p->div_pastem = 4.0;
  else if(*retry == 11)
    p->div_pastem = 8.0;
  else if(*retry == 12)
    p->div_pastem = 16.0;
  else if(*retry == 13)
    p->div_pastem = 30.0;
  else if(*retry == 14)
    p->div_pastem = 40.0;
  else if(*retry == 15)
    p->div_pastem = 75.0;
  else if(*retry == 16)
    p->div_pastem = 100.0;

  if(corto == 1)
    return modifica;

  // a questo punto prova a cambiare il grigliato della optim
  // girando fra due valori possibili
  if(*retry == 17)
    {
      p->div_pastem = 1.0;
      p->UMA[0] = 0.005;
    }
  else if(*retry == 18)
    p->div_pastem = 0.9090909;
  else if(*retry == 19)
    p->div_pastem = 0.8333333;
  else if(*retry == 20)
    p->div_pastem = 1.1;
  else if(*retry == 21)
    p->div_pastem = 1.2;
  else if(*retry == 22)
    p->div_pastem = 1.5;
  else if(*retry == 23)
    p->div_pastem = 0.6666667;

    
  return modifica;
}
