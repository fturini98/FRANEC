subroutine FATO(INDU,AS,AT,AU,AV,IPP, EL0,EL1,EL2,P0,P1,P2,R0,R1,R2,T0,T1,T2, &
          DPC,DTC,DELL,DTEF,EL0P,EL1P,EL2P,P0P,P1P,P2P,R0P,R1P,R2P, &
          T0P,T1P,T2P,CORZ) 
  use strut
  use fitt
  use tempe
  use numer

  implicit none
  
  real :: AS, AT, AU, AV
  integer :: INDU, IPP
  real :: EL0, EL1, EL2, P0, P1, P2, R0, R1, R2, T0, T1, T2, DPC, DTC, &
       DELL, DTEF, EL0P, EL1P, EL2P, P0P, P1P, P2P, R0P, R1P, R2P, &
       T0P, T1P, T2P, CORZ

  real :: aa, ab, ac, ad, ae, af, ag, ah, ai, al, am, an, ao, ap, aq, ar 
  real :: alf1, alf2, alf3, bet1, bet2, bet3, gam1, gam2, gam3, eta1
  real :: eta2, eta3, zot1, zot2, zot3, zot4, zot5, zot6
  real :: disc1, disc2, errel, errte, errpc, errtc, zum1, zum2, zum3, zum4
  real :: emmu, fat1, fat2, fat3, fat4, fat5, fat6
  
  real,parameter :: prec = 0.1
  

  INDU = 0 
  write(2,555) 
  AA = (EL2-EL0)/DPC 
  AB = (R2-R0)/DPC 
  AC = (P2-P0)/DPC 
  AD = (T2-T0)/DPC 
  if(.not. (IREAD == 2 .and. IPP == 0)) then 
     AE = (EL1-EL0)/DTC 
     AF = (R1-R0)/DTC 
     AG = (P1-P0)/DTC 
     AH = (T1-T0)/DTC 
  endif
  AI = (EL2P-EL0P)/DELL 
  AL = (R2P-R0P)/DELL 
  AM = (P2P-P0P)/DELL 
  AN = (T2P-T0P)/DELL 
  AO = (EL1P-EL0P)/DTEF 
  AP = (R1P-R0P)/DTEF 
  AQ = (P1P-P0P)/DTEF 
  AR = (T1P-T0P)/DTEF 
  if(.not. (IREAD == 2 .and. IPP == 0)) then
     ALF1 = AL-AF*AI/AE 
     ALF2 = AM-AG*AL/AF 
     ALF3 = AN-AH*AM/AG 
     BET1 = AP-AF*AO/AE 
     BET2 = AQ-AG*AP/AF 
     BET3 = AR-AH*AQ/AG 
     GAM1 = AF*AA/AE-AB 
     GAM2 = AG*AB/AF-AC 
     GAM3 = AH*AC/AG-AD 
     ETA1 = AT-AF*AS/AE 
     ETA2 = AU-AG*AT/AF 
     ETA3 = AV-AH*AU/AG 
     ZOT1 = ALF2-GAM2*ALF1/GAM1 
     ZOT2 = BET2-GAM2*BET1/GAM1 
     ZOT3 = ETA2-GAM2*ETA1/GAM1 
     ZOT4 = ALF3-GAM3*ALF2/GAM2 
     ZOT5 = BET3-GAM3*BET2/GAM2 
     ZOT6 = ETA3-GAM3*ETA2/GAM2 
     DISC1 = ZOT1*ZOT5-ZOT2*ZOT4 
     ERREL = (ZOT5*ZOT3-ZOT2*ZOT6)/DISC1 
     ERRTE = (ZOT1*ZOT6-ZOT3*ZOT4)/DISC1 
     ERRPC = ETA1/GAM1-ALF1*ERREL/GAM1-BET1*ERRTE/GAM1 
     ERRTC = AI*ERREL/AE+AO*ERRTE/AE-AA*ERRPC/AE-AS/AE 

     write(2,333) ERREL,ERRTE,ERRPC,ERRTC 

     ZUM1 = ERREL/ELL 
     ZUM2 = ERRTE/TEF 
     ZUM3 = ERRPC/PCEN 
     ZUM4 = ERRTC/TCEN 

     write(2,333) ZUM1,ZUM2,ZUM3,ZUM4 
     write(2,789) 

     if(abs(ZUM1) > 2. .or. abs(ZUM2) > 2. .or. abs(ZUM3) > 2. &
          .or. abs(ZUM4) > 2.) then
        INDU = 1
        return
     endif

     if(abs(ERREL)-prec*ELL > 0) then
        ERREL = prec*ELL*(ERREL/abs(ERREL))
     endif
     ELL = ELL+ERREL*CORZ 
     if(abs(ERRTE)-prec*TEF > 0) then 
        ERRTE = prec*TEF*(ERRTE/abs(ERRTE)) 
     endif
     TEF = TEF+ERRTE*CORZ 
     if(abs(ERRPC)-prec*PCEN > 0) then 
        ERRPC = prec*PCEN*(ERRPC/abs(ERRPC))
     endif
     PCEN = PCEN+ERRPC*CORZ 
     if(abs(ERRTC)-prec*TCEN > 0) then 
        ERRTC = prec*TCEN*(ERRTC/abs(ERRTC))
     endif
     TCEN = TCEN+ERRTC*CORZ 
     ELLOG = log10(ELL/3.826E+33) 
     TEFF = log10(TEF) 
     EMMU = EMTOT/1.989E+33 
     write(2,547) EMMU,ELLOG,TEFF 
     return 
  endif

  FAT1 = AC*AL/AB-AM 
  FAT2 = AC*AP/AB-AQ 
  FAT3 = AC*AT/AB-AU 
  FAT4 = AD*AL/AB-AN 
  FAT5 = AD*AP/AB-AR 
  FAT6 = AD*AT/AB-AV 
  DISC2 = FAT1*FAT5-FAT2*FAT4 
  ERREL = (FAT3*FAT5-FAT2*FAT6)/DISC2 
  ERRTE = (FAT1*FAT6-FAT3*FAT4)/DISC2 
  ERRPC = AL*ERREL/AB+AP*ERRTE/AB-AT/AB 
  ZUM1 = ERREL/ELL 
  ZUM2 = ERRTE/TEF 
  ZUM3 = ERRPC/PCEN 

  write(2,333) ZUM1,ZUM2,ZUM3 
  write(2,789) 

  if(abs(ZUM1) > 2 .or. abs(ZUM2) > 2 .or. abs(ZUM3) > 2 ) then
     INDU = 1 
     return
  endif

  if(abs(ERREL)-prec*ELL > 0) then 
     ERREL = prec*ELL*(ERREL/abs(ERREL)) 
  endif
  ELL = ELL+ERREL*CORZ 
  if(abs(ERRTE)-prec*TEF > 0) then 
     ERRTE = prec*TEF*(ERRTE/abs(ERRTE)) 
  endif
  TEF = TEF+ERRTE*CORZ 
  if(abs(ERRPC)-prec*PCEN > 0) then 
     ERRPC = prec*PCEN*(ERRPC/abs(ERRPC))
  endif
  PCEN = PCEN+ERRPC*CORZ 
  ELLOG = log10(ELL/3.826E33) 
  TEFF = log10(TEF) 
  EMMU = EMTOT/1.989E33 

  write(2,547) EMMU,ELLOG,TEFF 
  return 

547 format(3F8.4) 
789 format(///) 
555 format(15X,'DL',7X,'DTE',7X,'DPC',7X,'DTC') 
333 format(11X,10E10.3) 

end subroutine FATO
