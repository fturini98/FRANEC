subroutine inout(IDIFF, div_pastem, max_pastem, URA, ULA, UPA, UTA, UMA, VMM, &
     mix, Kover, par_OVER, maxAGE)
  use mistura
  use atmosfere
  use sceltachim
  use def_io

  implicit none

  integer :: idiff, kover
  real :: div_pastem, max_pastem, par_OVER, maxAGE
  real,dimension(4) :: URA, ULA, UPA, UTA, UMA 
  real,dimension(5) :: VMM
  character(len=16) :: mix(32)

  integer :: i, stat
  integer,parameter :: fileaggiuntivi = 0

  ! file in lettura
  ! i file di opacita' sono nella kappa 

  open(unit=4, FILE='Modstart.in',STATUS='OLD')
  open(unit=25, FILE='TABEL.in',STATUS='OLD')
  open(unit=8, FILE='EOS.in',STATUS='OLD')
  open(unit=5, file='EOSTOTALE.in',STATUS='OLD')
  open(unit=39,FILE='NUCL.DAT',STATUS='UNKNOWN')
  open(unit=ioOUT,FILE='OUT.DAT',STATUS='UNKNOWN')
  open(unit=ioBIGTAB,FILE='BIGTAB.DAT',STATUS='UNKNOWN')
  open(unit=ioCHIMICA,FILE='CHIMICA.DAT',STATUS='UNKNOWN')
  open(unit=ioFISICA,FILE='FISICA.DAT',STATUS='UNKNOWN')
  open(unit=23,FILE='SUBATM.DAT',STATUS='UNKNOWN')
  write(ioBIGTAB,103)
  write(ioOUT,108)

  ! CHIMICA SUPERFICIALE
  open(ioCHIMICASUP,file='CHIMICA_SUP.DAT')
  write(ioCHIMICASUP,1980)

  ! file con le abbondanze
  open(unit=81, file='abbondanze.in',STATUS='OLD')
  read(81,*) misturaop
  read(81,*) DEFAUHE3
  read(81,*) DEFAUC
  read(81,*) DEFAUN
  read(81,*) DEFAUO
  read(81,*) DEFAUFe
  read(81,*) DEFAULi6
  read(81,*) DEFAULi7
  read(81,*) DEFAUBe
  read(81,*) DEFAUB
  read(81,*) DEFAUD
  close(81) 

  ! file con le etichette delle misture note al sistema
  open(unit=81, file='misture.in',STATUS='OLD')
  i = 1
  do
     read(81, '(a16)', iostat=stat) mix(i)
     if (stat /= 0) exit ! se sono alla fine del file esco dal ciclo
     i = i + 1
  end do
  close(81) 

  ! file prodotti:
  open(UNIT=2 ,FILE='PRINT.DAT',STATUS='new')
  open(UNIT=10,FILE='GRAFI.DAT',STATUS='new')
  open(UNIT=15,FILE='SEGNALI.DAT',STATUS='new')
  open(UNIT=68,FILE='fasievolutive.log',status='unknown')
  open(unit=18,file='mixcno.dat',status='unknown')
  write(18,'("#        |--------------                core He", &
       "                  ----------!------------   zona esterna core He")') 
  write(18,'("# NMD    C             N             O             Fe        ", &
       "    z             C             N             O             Fe",&
       "            He3           z")')

  ! file per messaggi di errore
  open(66, FILE='error.log', status='unknown')

  ! file per messaggi standard di funzionamento
  open(iolog, FILE='run.log', status='unknown')

  ! log delle abbondanze
  write(iolog,*) "@ ================================"
  write(iolog,*) "@             Abbondanze          "
  write(iolog,*) "@ ================================"
  write(iolog,'(" @ ",a4,1x,e12.5)') "He3", DEFAUHE3
  write(iolog,'(" @ ",a4,1x,e12.5)') "C", DEFAUC
  write(iolog,'(" @ ",a4,1x,e12.5)') "N", DEFAUN
  write(iolog,'(" @ ",a4,1x,e12.5)') "O", DEFAUO
  write(iolog,'(" @ ",a4,1x,e12.5)') "Fe", DEFAUFe
  write(iolog,'(" @ ",a4,1x,e12.5)') "Li6", DEFAULi6
  write(iolog,'(" @ ",a4,1x,e12.5)') "Li7", DEFAULi7
  write(iolog,'(" @ ",a4,1x,e12.5)') "Be", DEFAUBe
  write(iolog,'(" @ ",a4,1x,e12.5)') "B", DEFAUB
  write(iolog,'(" @ ",a4,1x,e12.5)') "D", DEFAUD
  write(iolog,*) "@ ================================"
  write(iolog,*)

  open(UNIT=70,FILE='PROFILI_H.DAT',STATUS='new')     
  open(UNIT=71,FILE='PROFILI_HE.DAT',STATUS='new')     
  open(UNIT=72,FILE='PROFILI_C.DAT',STATUS='new')     
  open(UNIT=73,FILE='PROFILI_N.DAT',STATUS='new')     
  open(UNIT=74,FILE='PROFILI_O.DAT',STATUS='new')     
  open(UNIT=50,FILE='CHECK.DAT',STATUS='new')

  !     leggo il file di parametri dinamici
  open(80, FILE='parametri.in', status='old')
  read(80, *) IDIFF
  read(80, *) div_pastem
  read(80, *) max_pastem
  read(80,*) (URA(i), i=1,4)
  read(80,*) (ULA(i), i=1,4)
  read(80,*) (UPA(i), i=1,4)
  read(80,*) (UTA(i), i=1,4)
  read(80,*) (UMA(i), i=1,4)
  read(80,*) (VMM(i), i=1,5)
  read(80,*) do_bh_05
  read(80,*) do_ck_03
  read(80,*) kover
  read(80,*) par_OVER
  read(80,*) maxAGE
  read(80,*) STDCHIM    ! scelta chimica
  close(80)
  if(IDIFF  ==  0) then
     write(*,*) 'DIFFUSIONE OFF'
  else
     write(*,*) 'DIFFUSIONE ON'
  endif

  if(kover == 1) open(85,FILE='OVERSHOOT-CENTRAL.DAT',STATUS='unknown')

  ! log dei parametri dinamici su run.log
  write(iolog,*) "@ ================================"
  write(iolog,*) "@       Parametri dinamici        "
  write(iolog,*) "@ ================================"
  write(iolog,'(" @ diffusione",1x,i4)') IDIFF
  write(iolog,'(" @ div_pastem",1x,f10.6)') div_pastem
  write(iolog,'(" @ max_pastem",1x,e10.4)') max_pastem
  write(iolog,'(" @ URA",1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)') (URA(i), i=1,4)
  write(iolog,'(" @ ULA",1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)') (ULA(i), i=1,4)
  write(iolog,'(" @ UPA",1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)') (UPA(i), i=1,4)
  write(iolog,'(" @ UTA",1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)') (UTA(i), i=1,4)
  write(iolog,'(" @ UMA",1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)') (UMA(i), i=1,4)
  write(iolog,'(" @ VMM",1x,f7.5,1x,f7.5,1x,f7.5,1x,f7.5,1x,f7.5)') &
       (VMM(i), i=1,5)
  write(iolog,'(" @ atmos. BH05",1x,i4)') do_bh_05
  write(iolog,'(" @ atmos. CK03",1x,i4)') do_ck_03
  write(iolog,'(" @ overshooting",1x,i4)') kover
  write(iolog,'(" @ chimica standard (1=on)",1x,i4)') stdchim
  write(iolog,'(" @ parametro overshooting",1x,f8.5)') par_OVER
  write(iolog,'(" @ max age",1x,f9.5)') maxAGE
  write(iolog,*) "@ ================================"
  write(iolog,*) 

1980 format('#NMD',3x,'log(AGE)',9x,'H',12x,'H2',12x,'He3',10x,'He4', &
       12x,'C',13x,'N',13x,'O',12x,'Fe',12x,'Li6',11x,'Li7',12x,'Be9', &
       11x,'B11')

103 format('# ',1X,'NMOD   LOG(T)        H/HE    LOG L   LOG TE   LOG TC  ', &
         'LOG RHOc',&
         3x,' TMAX    MTMAX     M-CC    M-HEC       M-CE       L-PP', &
         '      L-CNO      L-3A      L-GRA   HE_SUP  LOGTce   GICO       ', &
         'R/Rsun      LOG(g)      Eneu')
108 format('# NMOD   LOG(T)          H/HE    LOG L         LOG TE', &
         '        MASS', '      L-GRA      L-3A      log(Fe/H)     [Fe/H]' &
         '      R         Logg       Dni        nimax')

end subroutine inout


! questa routine viene usata per leggere i file di abbondanze modificati per 
! la partenza da HB (pepper)
! Vengono rimpiazzati i valori delle variabili gia' lette e vengono assegnati
! NCo, OCO, XMEin 
subroutine inout_hb(NCO,OCO,XMEin)
  use mistura
  use def_io

  implicit none
  
  real :: NCO,OCO,XMEin

  open(unit=81, file='abbondanze-HB.in',STATUS='OLD')
  read(81,*) misturaop
  read(81,*) DEFAUHE3
  read(81,*) DEFAUC
  read(81,*) DEFAUN
  read(81,*) DEFAUO
  read(81,*) DEFAUFe
  read(81,*) DEFAULi6
  read(81,*) DEFAULi7
  read(81,*) DEFAUBe
  read(81,*) DEFAUB
  read(81,*) DEFAUD

  read(81,*) NCO
  read(81,*) OCO
  read(81,*) XMEin
  close(81) 

  write(iolog,*) "@ ================================"
  write(iolog,*) "@    Abbondanze HB (sostituzione) "
  write(iolog,*) "@ ================================"
  write(iolog,'(" @ ",a6,1x,e12.5)') "He3", DEFAUHE3
  write(iolog,'(" @ ",a6,1x,e12.5)') "C", DEFAUC
  write(iolog,'(" @ ",a6,1x,e12.5)') "N", DEFAUN
  write(iolog,'(" @ ",a6,1x,e12.5)') "O", DEFAUO
  write(iolog,'(" @ ",a6,1x,e12.5)') "Fe", DEFAUFe
  write(iolog,'(" @ ",a6,1x,e12.5)') "Li6", DEFAULi6
  write(iolog,'(" @ ",a6,1x,e12.5)') "Li7", DEFAULi7
  write(iolog,'(" @ ",a6,1x,e12.5)') "Be", DEFAUBe
  write(iolog,'(" @ ",a6,1x,e12.5)') "B", DEFAUB
  write(iolog,'(" @ ",a6,1x,e12.5)') "D", DEFAUD
  write(iolog,'(" @ ",a6,1x,e12.5)') "N core", NCO
  write(iolog,'(" @ ",a6,1x,e12.5)') "O core", OCO
  write(iolog,'(" @ ",a6,1x,e12.5)') "z core", XMEin
  write(iolog,*) "================================"
  write(iolog,*)

  NCO = NCO*XMEin
  OCO = OCO*XMEin
  
end subroutine inout_hb
