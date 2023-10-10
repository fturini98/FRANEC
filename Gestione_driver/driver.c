#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <unistd.h>
#include <time.h>

#define MAX 10000

main(int argc, char **argv)
{
  FILE *fp, *fpz, *fout, *fpm;
  char prova[MAX][256], PR[MAX][256], com[1256], nmods[64];
  char eseguibile[256], zs[256], nomeout[256], filest[256], linea[256];
  char nomeout1[256], nomeout2[256], nomelog[256], pathout[1024], nOVER[1024];
  char *mistura, *numdrv;
  int i, j, n, k, ip=1, PAR, SAVEHB, km, iread, nmix, old, ln;
  int io, OVER_attivo, im, shf;
  double par_OVER;
  double dum, met;
  char mix[128][16];
  time_t t;

  int NUMP=8; // numero perturbazioni
  char listP[8][32] = {"! pp", "! C12", "! N14", "! 3alpha", "! neutr", 
		       "! OP_rad", "! OP_cond", "! diffusion"};
  
  // errore di argomenti... uscita
  if(argc != 7) 
    {
      printf("Usage: %s <eseguibile_franec> <prove in input> <PAR=0/1> <COD> <SAVEHB=0/1> <path out>\n", 
	     argv[0]);
      exit(-1);
    }

  // processo gli argomenti da linea di comando
  sprintf(eseguibile, "%s", argv[1]);  // eseguibile
  PAR = atoi(argv[3]);                 // per fare prove con vari parametri.in 
  fp = fopen(argv[2], "r");
  // cosa fare con eventuali cartelle gia' presenti
  old = atoi(argv[4]);
  SAVEHB = atoi(argv[5]);
  sprintf(pathout, "%s", argv[6]);  // path out

  // apro il file di log in cui scrivo i modelli terminati
  numdrv = strrchr(argv[2], '_');
  sprintf(nomelog, "driver.log%s", numdrv);
  fout = fopen(nomelog, "w");
  
  // lettura modelli da simulare
  i = 0;
  while(fscanf(fp, "%s", prova[i]) != EOF)
    i++;
  fclose(fp);

  // lettura mistura
  fp = fopen("out/PAR/misture.in", "r");
  j = 0;
  while(fscanf(fp,"%s", mix[j]) != EOF)
    j++;
  fclose(fp);
  fp = fopen("out/PAR/abbondanze.in", "r");
  fscanf(fp, "%d", &nmix);
  mistura = &mix[nmix-1][0];
  fclose(fp);

  if(PAR == 1) 
    {
      ip = 0;
      // file di gestione variazioni parametri per lo stesso modello
      fp = fopen("prove_parametri.dat", "r");
      while(fscanf(fp, "%s", PR[ip]) != EOF)
	ip++;
      fclose(fp);
    }

  for(j=0;j<i;j++) // ciclo sui modelli
    {
      for(k=0;k<ip;k++) // ciclo sui parametri.in
	{
	  // creo e popolo la cartella di lavoro temporanea...
	  sprintf(com, "mkdir -p out/%s", prova[j]);
	  system(com);
	  sprintf(com, "cd out/PAR; cp * ../%s; cd ../..; touch out/%s/run-by-driver", prova[j], prova[j]);
	  system(com);
	  
	  if(PAR == 1) // se devo gestire vari file di parametri ...
	    {
	      sprintf(com, "cp tools/driver/param/%s out/%s/parametri.in",
		      PR[k], prova[j]);
	      system(com);
	    }

	  // metto al suo posto il Modstart corretto...
	  sprintf(nmods, "%s", prova[j]);
	  if(nmods[30]=='_')
	    {
	      nmods[30] = '\0';
	      shf = 0;
	    }
	  else
	    {
	      nmods[31] = '\0';
	      shf = 1;
	    }
	  
	  sprintf(com, "out/%s/perturbazioni.in", prova[j]);
	  fpm = fopen(com, "w");
	  for(im=0;im<NUMP;im++)
	    fprintf(fpm, "%c\t\t%s\n", prova[j][31+im+shf], listP[im]);
	  fclose(fpm);

	  sprintf(com, "cp tools/driver/start/%s out/%s/Modstart.in", 
		  nmods, prova[j]);
	  system(com);
	  
	  // questo blocco serve per leggere eventuale overshoot
	  sprintf(com, "out/%s/parametri.in", prova[j]);
	  fpz = fopen(com, "r");
	  //	  fpz = fopen("out/PAR/parametri.in", "r");
	  for(io=0;io<12;io++)
	    fgets(linea,255,fpz);
	  sscanf(linea,"%d", &OVER_attivo);
	  fgets(linea,255,fpz);
	  sscanf(linea,"%lf", &par_OVER);
	  fclose(fpz);
	  sprintf(nOVER,"");
	  if(OVER_attivo == 1) 
	    sprintf(nOVER,"_OV%.3f", par_OVER);

	  // questo blocco serve per recuperare la metallicita'
	  // da utilizzare.
	  sprintf(filest, "out/%s/Modstart.in", prova[j]);
	  fpz = fopen(filest, "r");
	  fgets(linea,255,fpz);
	  sscanf(linea,"%d", &iread);
	  if(iread == 1) // devo leggere il MOD.DAT
	    {
	      for(km=1;km<26;km++)
		fgets(linea,255,fpz);
	      for(km=0;km<112;km++)
		fscanf(fpz,"%lf", &dum);
	      fscanf(fpz,"%lf", &met);
	      sprintf(zs,"%.5f",met);
	    }
	  else
	    {
	      for(km=1;km<19;km++)
		fgets(linea,255,fpz);
	      fscanf(fpz,"%lf", &met);
	      sprintf(zs,"%.5f",met);
	    }
	  fclose(fpz);
	  // fine lettura metallicita' -> zs

	  // metto al loro posto le tabelle EOS
	  sprintf(com, 
		  "cp archivioEOS/EOS-Z%s.DAT out/%s/EOS.in", zs, prova[j]);
	  system(com);
	  sprintf(com, 
		  "cp archivioEOS/EOSTOTALE-Z%s.DAT out/%s/EOSTOTALE.in", 
		  zs, prova[j]);
	  system(com);
	  
	  if(PAR == 1) // gestisco il nome di out nel caso di multiparametri
	    sprintf(nomeout,"%s%s_%s_%s",prova[j], nOVER, mistura, PR[k]);
	  else
	    sprintf(nomeout,"%s%s_%s", prova[j], nOVER, mistura);
	  
	  // gestisco il salvataggio di casi di HB morti anzitempo
	  // va recuperata la cartella in cui inserire i file di out
	  if(SAVEHB == 1)
	    {
	      strcpy(nomeout1,prova[j]);
	      ln = strlen(nomeout1)-9;
	      nomeout1[ln] = '\0';
	      sprintf(nomeout2,"%s/HB/%s", nomeout1, &prova[j][ln+3]);
	    }

	  if(old != 0)
	    {
	      if(old == 1)
		{
		  time(&t);
		  // backup vecchi out
		  sprintf(com, "mkdir -p  tools/driver/old-out/%s-%ld", 
			  nomeout, t);
		  system(com);
		  sprintf(com, "cp %s/%s/* tools/driver/old-out/%s-%ld >& /dev/null", 
			  pathout, nomeout, nomeout, t);
		  system(com);
		}
	       // rimuovo i vecchi out
	      sprintf(com, "rm -f %s/%s/*", pathout, nomeout);
	      system(com);
	    }

	  // lancio simulazione
	  sprintf(com, "cd out/%s; ./%s > out.dat", prova[j], eseguibile);
	  system(com);
	  
	  // crazione cartella di output
	  sprintf(com, "mkdir -p %s/%s", pathout, nomeout);
	  system(com);
	  
	  // compressione degli out della simulazione nella cartella creata
	  sprintf(com, "cd out/%s; bzip2 GRAFI.DAT  CHIMICA.DAT FISICA.DAT", prova[j]);
	  system(com);

	  // copia degli out della simulazione nella cartella creata
	  if(pathout[0] != '/')
	    {
	      sprintf(com, "cd out/%s; cp BIGTAB.DAT OUT.DAT GRAFI.DAT.bz2  CHIMICA.DAT.bz2 FISICA.DAT.bz2 error.log lancia.log run.log Modstart.in-* scenario.log mixcno.dat abbondanze.in abbondanze-HB.in* scansioneHB.in flash.log fasievolutive.log parametri.in perturbazioni.in lastmass DarkMatter.DAT DarkMatterERROR.DAT DarkMatterCattura.DAT ../../tools/driver/out/%s >& /dev/null", 
		      prova[j], nomeout); 
	    }
	  else
	    {
	      sprintf(com, "cd out/%s; cp BIGTAB.DAT OUT.DAT GRAFI.DAT.bz2  CHIMICA.DAT.bz2 FISICA.DAT.bz2 error.log lancia.log run.log Modstart.in-* scenario.log mixcno.dat abbondanze.in abbondanze-HB.in* scansioneHB.in flash.log fasievolutive.log parametri.in perturbazioni.in lastmass DarkMatter.DAT DarkMatterERROR.DAT DarkMatterCattura.DAT %s/%s >& /dev/null", 
		      prova[j], pathout, nomeout); 
	  
	    }
	  system(com);

	  // nel caso di salvataggio di HB sposto i file di out nella
	  // cartella corretta
	  if(SAVEHB == 1)
	    {

	      sprintf(com,"rm -f %s/%s/*", pathout, nomeout2);
	      system(com);
	      
	      sprintf(com,"cd %s/%s; mv BIGTAB.DAT OUT.DAT GRAFI.DAT.bz2 CHIMICA.DAT.bz2 FISICA.DAT.bz2 run.log error.log lancia.log Modstart.in-bak lastmass DarkMatter.DAT DarkMatterERROR.DAT DarkMatterCattura.DAT ../%s", 
		      pathout, nomeout, nomeout2);
	      system(com);
	      sprintf(com,"rm -rf %s/%s", pathout, nomeout);
	      system(com);
	    }
	  
	  // log del modello terminato
	  fprintf(fout, "%s\n", prova[j]);
	  fflush(fout);

	  // rimozione cartella di lavoro temporanea
	  sprintf(com, "rm -rf out/%s", prova[j]);
	  system(com);
	}
     }
  // rimozione file lista parziale
  unlink(argv[2]);

  fclose(fout);
}
