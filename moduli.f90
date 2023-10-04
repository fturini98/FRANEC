module parametri
  integer,parameter :: LIM = 3000
  integer,parameter :: MELE = 26
  !     parametri atmosfera:
  integer,parameter :: NZ_BH = 10, NG_BH = 13, NT_BH = 1000, NE_BH = 80
  integer,parameter :: NZ_CK = 8, NG_CK = 11, NT_CK = 600, NE_CK = 80
  integer,parameter :: NN_ATM = 72
  ! numero reazioni seguite da cross
  integer,parameter :: nsezi = 40
end module parametri


module def_io
  integer,parameter :: ioOUT = 37, ioBIGTAB = 36, iolog = 67, &
       ioCHIMICASUP = 38, ioFISICA = 22, ioCHIMICA = 21
end module def_io


module costanti 
  real,parameter :: ssb = 5.6704d-5, &
       pigre = 3.1415926535d0, &
       Ggrav = 6.674d-8, &
       MeV_erg = 1.602177d-6, arad_3 = 2.52192d-15, &
       Kboltz = 1.38065d-16, p_mass = 1.672622d-24, &
       sec_anno = 3.1536d+7, cluce = 2.99792458d+10, &
       AVN = 6.022141d23, sec_giorno = 86400.0d0
  !! sigma (erg cm-2 s-1 K-4), G (cm3 gr-1 s-2), conversione MeV->erg,
  !! sigma*4/cluce, costante boltzmann (erg K^-1), massa protone (gr),
  !! secondi in un anno,
  !! velocita' della luce (cm s^-1), numero di Avogadro (mol^-1), 
  !! secondi in un giorno
  real, parameter :: cte_grad = 3.d-8 / (64.d0 * pigre * ssb * Ggrav)
  real, parameter :: cte_rad = 3.d10 / (64.d0 * pigre * ssb)
  real, parameter :: Rsun = 6.9598d0, Lsun = 3.8418d0, Msun = 1.989d0
  !! Rsun(10^10 cm), Msun(10^33 gr), Lsun(10^33 erg)
  real :: Teff_sun = ( 1.0d13 * Lsun / (4.0d0 * pigre * ssb * (Rsun**2)) )**0.25
  real, parameter :: temposole = 4.57d9 !! eta' sole anni
  
end module costanti

module fasievolut
  ! Punti di stop definiti:
  ! fase = -10 : MS
  ! fase = -11 : Sole
  ! fase = -30 : esaurimento H centrale
  ! fase = -50 : RGB
  ! fase = -60 : bump RGB
  ! fase = -90 : flash He
  ! fase = -91 : innesco He
  ! fase = -100 : ripartenza dopo pepper
  ! fase = -110 : esaurimento He centrale
  ! fase = -150: pulsi
  integer,parameter :: fase_MS = -10, fase_Sole = -11
  integer,parameter :: fase_exH = -30, fase_RGB = -50
  integer,parameter :: fase_bump1 = -60, fase_bump2 = -65
  integer,parameter :: fase_FlashHe = -90,  fase_innHe = -91
  integer,parameter :: fase_postPepper = -100
  integer,parameter :: fase_exHe = -110, fase_pulsiT = -150
  integer :: idfase
end module fasievolut

module sceltachim
  integer :: STDCHIM 
end module sceltachim

module overshoot
  ! eqlb, idro, epsi, stampa, mixing
  integer :: L3, KOVER
  real :: par_OVER
end module overshoot

module zone_conv
  ! evolut, epsi, zone_convettive
  use parametri
  integer,dimension(lim) :: vconv
  integer :: n_conv
  integer,dimension(2, 100) :: b_conv
end module zone_conv

module intero 
  ! main, atmos, optim
  use parametri
  integer,dimension(LIM) :: IN
end module intero

module mistura 
! main, io, innes, kappa
 real :: DEFAUC, DEFAUN, DEFAUO
 real :: DEFAUFe, DEFAUHe3, DEFAULi6, DEFAULi7, DEFAUBe, DEFAUB, DEFAUD
 integer :: misturaop
end module mistura

module fisica 
  ! main, veiove, atmos, carbon, elio, epsi, epsig, eqlb, evolut, fitta,
  ! hbrezo, henyey, idreli, idro, innes, maslos, mixing, optim,
  ! pastem, plotta, quatm, stampa, risulta, risultati, ciacioLi
  use parametri
  real,dimension(6,LIM) :: G, GG
end module fisica

module preco 
  ! main, atmos(KSF), quatm(MAYA)
  integer :: MAYA, KSF
end module preco

module nummod  
  ! main, epsi, epsig, henyey, stampa, eqlb
  integer :: NMD
end module nummod

module chequ 
  ! stampa, idro, epsi
  use parametri
  integer,dimension(lim) :: KECUIL
end module chequ

module gratm 
  ! main, quatm
  real,dimension(3,4,4) :: U, V 
  real,dimension(4) :: EL, TE 
  real :: PMA1, PMA2
end module gratm

module strut 
  ! main, ciacioLi, risulta, fato, innes, maslos, mixing, optim
  ! pastem, plotta, quatm, stampa, veiove
  real :: EMTOT, ROCEN, ETAFPR
end module strut

module mesh 
  ! atmos, quatm, stampa 
  integer :: MAXME, MAXMIV
end module mesh

module serchi 
  ! evolut, optim, stampa
  use parametri
  real,dimension(2,LIM) :: XSERV
end module serchi

module nuconv 
  ! stampa, mixing, optim, pastem
  real,dimension(12) :: BB 
  real :: GSEMI, GICO
end module nuconv
 
module atmosfere 
  ! main, atm_ck03, atm_bh05, atmos, stampa
  real,dimension(3) :: V_ATM 
  real :: X_SUP, Y_SUP, Z_SUP
  integer :: IBH05_LETTA, ICK03_LETTA, N_BH05, N_CK03, NMOD_PRE
  integer :: do_bh_05, do_ck_03
end module atmosfere

module equil 
  ! epsi, eqlb, idro
  real,dimension(5) :: COEFF
  integer :: ISHEL
end module equil

module fitt  
  ! innes, fato, veiove
  real :: ELL,TEF,PCEN,TCEN
end module fitt

module neutri 
  ! epsi, stampa
 real,dimension(8) :: PPNEU
 real :: TCRIT 
end module neutri

module second  
  ! fitta, henyey, resnuc, stampa 
  use parametri
  real,dimension(4,LIM) :: ALF, BET, GAM
  real,dimension(5,LIM) :: DG
end module second

module sistema 
  ! henyey, resnuc
  real,dimension(4) :: ALF1, BET1, GAM1, E
  real,dimension(4,4) :: B, C
end module sistema

module atm 
  ! fitta, innes, quatm
  real :: RB, PB, TB, DPL, DPT, DRL, DRT, DTL, DTT
end module atm

module chimic 
  ! main, eqlb, optim, atmos, ciacioLi, decidi, risulta, risultati, epsig,
  ! evolut, hbrezo, henyey, innes, mixing, pastem, plotta, stampa, veiove
  ! stop_evolut, epsi
  use parametri
  real,dimension(MELE,LIM) :: XXX, XXV
  real,dimension(MELE,2) :: XCE
end module chimic

module indi 
  ! innes, veiove, hbrezo
  integer :: LAXME, LAST, NABLA, IOTA
end module indi

module chim 
  ! main, hbrezo, optim, atmos, cross, epsi, screen, maslos, pastem, stampa,
  ! fitta, evolut, henyey, veiove, kappa, mixing, neutr, epsig, state,
  ! state_free, innes, ciacioLi, risultati
  use parametri
  real,dimension(MELE) :: XX
  real :: HT1, HT1V, EFFEP
  !  1=XH,     2=XHE3,   3=XHE4,    4=XC12,   5=XN14,
  !  6=XO16,   7=XO18,   8=XNE20,   9=XNE22, 10=XMG24,
  ! 11=XMG25, 12=XMG26, 13=XSI28,  14=XNEU,  15=XNA23,
  ! 16=XC13,  17=XN15,  18=XO17,   19=XF19,  20=XNE21,
  ! 21=XFE,   22=XLI6,  23=XLI7,   24=XBE,   25=XB
  ! 26=XD
end module chim

module tempe  
  ! main, maslos, pastem, stampa, ciacioLi, risultati,
  ! fato, fitta, innes, quatm
  real :: ELLOG, TEFF, TEFFV, TEFFVV
end module tempe

module numer 
  ! main, fitta, fato, veiove, innes, quatm, optim, hbrezo, atmos, resnuc,
  ! stampa, henyey, kappa, epsig
  integer :: MAIS, IPRALL, ISUB, IBAT, IREAD, IQUAT, IHP, IFAST
end module numer

module print1
  ! main, stampa
  real :: AHMAX,AHEMAX,AGRMAX,ANUMAX, FG1,FG2,FG3,FG4,FG5,FG6
end module print1

module varie 
  ! main, atmos, epsi, innes, kappa, maslos, quatm, stampa, supera, stop_evolut
  real :: XH, XME, HECE, ALFA, BMAG, EMMA, FRAZ, EMAXH, EMAHE
end module varie

! =================================================
!              moduli per xcotrin
! =================================================

module parametri_xcotrin
  integer,parameter :: mx=1, mc=8, mo=8, nrm=19, nrb=1, nre=19
  integer,parameter :: nr=nre+1-nrb, ntabs=60, ntm=70, ntb=1, nt=ntm+1-ntb
  integer,parameter :: IPR=20
end module parametri_xcotrin

module aaa 
  ! opac, cointerp, t6rinterp
  use parametri_xcotrin
  real,dimension(mx,mc) :: oxf, cxf, xcdf, xodf, cxdf, oxdf
  real,dimension(mx,nt,nr) :: opl
  integer,dimension(mx) :: itime = (/0/)                
end module aaa

module aa 
  ! opac, cointerp, t6rinterp, readco
  use parametri_xcotrin
  real,dimension(mc) :: xcd, xod, xc, cxd, cx 
  real,dimension(mo) :: xo, oxd, ox
  real,parameter,dimension(mc) :: xcs = (/0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/)
  real,parameter,dimension(mc) :: xos = (/0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/)
  real ::  zzz, xxh
  real,dimension(mx) :: xx
  real,dimension(4) :: h, q
  integer :: nc,no
end module aa

module a1 
   ! opac, cointerp, t6rinterp, readco
  use parametri_xcotrin
  real,dimension(mx,mc,mo,nt,nr,10) :: covett 
  real,dimension(mx,mc,nt,nr,10) :: diagvett
  real,dimension(nt) :: t6list, alt, dfs 
  real, dimension(nr) :: alr, dfsr

  real,dimension(mx,mo,nt,nr,10) :: diagovett
  real,dimension(nt,nr) :: opk
  real,dimension(3,mx) :: a
  real,dimension(3) :: b
  
  real,parameter,dimension(8) :: xa = (/0.0,0.03,0.1,0.35,0.7,0.,0.,0./)
  
  real,dimension(nrm) :: alrf 
  real,dimension(ntm,nrm) :: cof 
  real,dimension(ntm) :: t6listf
  real,dimension(nt,nr) :: opk2 
  real,dimension(mx) :: dfsx

  integer,parameter,dimension(101) :: index1 = (/1,2,2,3,3,3,3,3,3,3,4,4, &
       4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6, &
       6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, &
       7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/)

  integer,dimension(mx,mc) :: n
  integer :: m, mf  
end module a1

module alink  
  ! readco, opaltab
  use parametri_xcotrin
  integer :: NTEMP, NSM, nrlow, nrhigh
  real :: RLE
  real,dimension(100) :: t6arr 
  real,dimension(100,nr) ::  coff
end module alink

module cst  
  ! fity, fitx, interp, smoothing, opaltab, readco
  real :: tmax, RLS
  integer :: NRL, nset
end module cst

module d 
  ! opac, t6rinterp, quad
  real :: dkap
end module d

module bb  
  ! opac, cointerp, t6rinterp
  real :: xodp, xcdp, xxco, cxx, oxx
  integer :: l1, l2, l3, l4, k1, k2, k3, k4, ip, iq
end module bb

module b1  
  ! opac, readco
  use parametri_xcotrin
  real,dimension(mx,ntabs) :: x, y, zz, xca, xoa
  integer,dimension(mx,ntabs) :: itab
  integer,parameter,dimension(nrm) :: nta = (/70,70,70,70,70,70,70,70,70,70, &
       70,70,70,70,69,64,60,58,57/)
end module b1

module cf  
  ! smoothing, opaltab, fity, fitx, interp
  use parametri_xcotrin
  real,dimension(85,IPR) :: F, FX, FY, FXY
end module cf

module e  
  ! opac, t6rinterp, readco, kappa
  real :: opact,dopact,dopacr,dopactd
end module e

module indicefile 
  ! kappa, opac, cointerp, readco
  integer :: indice, icntrl, ireadco
end module indicefile

module recoin  
  ! opac, readco
  integer :: itimeco = 0
  integer :: mxzero
end module recoin


module strutture
  type pow_of_T
     real :: T9, T92, T93, T912, T913, T923, T932, T6, T612
  end type pow_of_T

  type luminosity
     real :: Lpp, Lcno, L3a, Lgra, L
     real :: teff
     real :: M_bordo_conv, nuclear
     real,dimension(3) :: pn
  end type luminosity
end module strutture

module interfaccia
  use strutture
  interface
     function SK(RO,T,X,Z1,Z2, cache) 
       use parametri
       real :: sk,RO,T,X(MELE)
       integer :: Z1, Z2, cache
     end function SK

     real function f(y)
       real :: y
     end function f
  end interface
end module interfaccia

module interfaccia2
  use strutture
  ! external procedures
  interface
     subroutine ATM_BH05(GIINI,TE,TAU_RACC) 
       real :: GIINI, TE, TAU_RACC 
     end subroutine ATM_BH05
     
     subroutine ATM_CK03(GIINI,TE,TAU_RACC) 
       real :: GIINI,TE,TAU_RACC
     end subroutine ATM_CK03
     
     subroutine ATMOS(FLUMIN,TEFF,P,T,R,EMASS,ISCR,IRA) 
       real :: FLUMIN, TEFF, P, T, R, EMASS
       integer :: ISCR, IRA
     end subroutine ATMOS
     
     subroutine CARBON (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S) 
       use parametri
       real :: RO, T, PAS
       integer :: JF, INI, IFI, N
       real,dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real,dimension(32) :: S
     end subroutine CARBON
     
     subroutine CROSS (RO,T,XXn,IR,ISHORT,S,powt,nabla)
       use strutture
       use parametri
       real :: RO, T
       real,dimension(mele) :: xxn
       integer :: IR, ISHORT, nabla
       real,dimension(32) :: S
       type(pow_of_T) :: powt
     end subroutine CROSS

     subroutine CUB(A,B,ALFA,BETA,GAMMA,DELTA) 
       real,dimension(4) :: A, B
       real :: ALFA, BETA, GAMMA, DELTA
     end subroutine CUB

     subroutine CUB2dr(A, B, idx, xABGD) 
       real,dimension(4) :: B
       real,dimension(4,4) :: A
       real,dimension(320) :: xABGD
       integer :: idx
     end subroutine CUB2dr

     subroutine CUB2dt(A, B, xtABGD, is, ie) 
       real,dimension(4) :: B
       real,dimension(80) :: xtABGD, A
       integer :: is, ie
     end subroutine CUB2dt

     subroutine ELIO (RO,T,PAS,JF,INI,IFI,N,ISHORT, YV, YF, NEQUI, S) 
       use parametri
       real :: RO, T, PAS
       integer :: JF, INI, IFI, N, ISHORT
       real,dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real,dimension(32) :: S
     end subroutine ELIO

     subroutine EPSI (RO,T,JF,NABLA,HT,EPS,EPSA, ECAR, DELTA) 
       real :: RO, T, HT, EPS, EPSA, ECAR, DELTA
       integer :: JF, NABLA
     end subroutine EPSI

     subroutine EPSIG(P,T,RO,CP,HT,JF,GRAVI,MAXMV,DAD,IPP) 
       real :: P, T, RO, CP, HT, GRAVI, DAD
       integer :: JF, MAXMV, IPP
     end subroutine EPSIG
     
     subroutine EQLB(MAXME) 
       integer :: maxme
     end subroutine EQLB

     subroutine EVOLUT(INI, IFI) 
       integer :: INI, IFI
     end subroutine EVOLUT

     subroutine akappetta(Y,xC,xO,tlog,rlog,cap) 
       real :: Y,xC,xO,tlog,rlog,cap
     end subroutine akappetta

     subroutine FATO(INDU,AS,AT,AU,AV,IPP, EL0,EL1,EL2,P0,P1,P2,R0,R1,R2, &
          T0,T1,T2,DPC,DTC,DELL,DTEF,EL0P,EL1P,EL2P,P0P,P1P,P2P,R0P,R1P,R2P, &
          T0P,T1P,T2P,CORZ) 
       real :: AS, AT, AU, AV
       integer :: INDU, IPP
       real :: EL0, EL1, EL2, P0, P1, P2, R0, R1, R2, T0, T1, T2, DPC, DTC, &
       DELL, DTEF, EL0P, EL1P, EL2P, P0P, P1P, P2P, R0P, R1P, R2P, &
       T0P, T1P, T2P, CORZ
     end subroutine FATO

     subroutine FITTA(M,DAM,NUM1,NUM2,FRAT, ierr)
       integer :: M, NUM1, NUM2, ierr
       real :: DAM, FRAT
     end subroutine FITTA

     subroutine HBREZO (AMCO, YCO, CCO, NCO, OCO)
       real :: AMCO, YCO, CCO, NCO, OCO
     end subroutine HBREZO
     
     subroutine HENYEY(MAXME, MAXMV, ERROR, INUM, IERR) 
       real :: ERROR
       integer :: MAXME, MAXMV, INUM, IERR
     end subroutine HENYEY

     subroutine KAPPA(RHO,TEMP,CAP,Y) 
       real :: RHO, TEMP, CAP, Y
     end subroutine KAPPA

     subroutine IDRELI (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S) 
       use parametri
       real :: RO, T, PAS
       integer :: JF, INI, IFI, N
       real,dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real,dimension(32) :: S
     end subroutine IDRELI

     subroutine IDRO (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S, oish) 
       use parametri
       real :: RO, T, PAS
       integer :: JF, INI, IFI, N, oish
       real,dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real,dimension(32) :: S
     end subroutine IDRO

     subroutine INNES(ITCC,IPP,MAXMEin,ILEG,LBA, YB,ZB, He, Zeta, Alpha, do_relax) 
       real :: YB, ZB, He, Zeta, Alpha
       integer :: ITCC, IPP, MAXMEin, ILEG, LBA 
       logical :: do_relax
     end subroutine INNES
     
     subroutine LOCAL(P,T,R,EL,EM,RO,GRAD,CP,PMU,ACCA,CAP,RADIAT,ADIAB, &
          GI,Y,GRSAD)
       real :: P,T,R,EL,EM,RO,GRAD,CP,PMU,ACCA,CAP,RADIAT,ADIAB, &
            GI,Y,GRSAD
     end subroutine LOCAL

     subroutine MASLOS(MAXME,DM,nmd,iread, fase, fase_start, fase_stop) 
       integer :: MAXME,nmd,iread, fase, fase_start, fase_stop
       real :: DM
     end subroutine MASLOS

     subroutine MIXING(MAXME,iread)
       integer :: MAXME,iread
     end subroutine MIXING

     subroutine NEUTR(RO,T,QT) 
       real :: RO, T, QT
     end subroutine NEUTR

     subroutine OPTIM(MAXME,SCALA,PROVV,ECNO, URA,ULA,UPA,UTA,UMA,VMM,fase) 
       integer :: MAXME,fase
       real :: SCALA, PROVV, ECNO
       real,dimension(4) :: URA, ULA, UPA, UTA, UMA
       real,dimension(5) :: VMM
     end subroutine OPTIM

     subroutine PASTEM(MAXME,MAXMV,PROV,PROVV,ECNO,ECNOV,FG6,FG6V,TMAX, &
          div_pastem,max_pastem)
       integer :: MAXME, MAXMV
       real :: PROV, PROVV, ECNO, ECNOV, FG6, FG6V, div_pastem, max_pastem, &
            TMAX
     end subroutine PASTEM

     subroutine PLOTTA(MAXME,RTOT) 
       integer :: MAXME
       real :: RTOT
     end subroutine PLOTTA

     subroutine QUATM(NABLA,INI,MAXMEin,MAXMV,EMTOV) 
       integer :: NABLA, INI, MAXMEin, MAXMV
       real :: EMTOV
     end subroutine QUATM

     subroutine RESNUC(JFF)                                             
       integer :: jff
     end subroutine RESNUC

     subroutine SIMQ(A,B,N,KS) 
       real,dimension(361) :: A
       real,dimension(19) :: B
       integer :: N, KS
     end subroutine SIMQ

     subroutine STAMPA(TEMPO,DMAS,NMOD,NABLA,KSD,MAXMEin,MAXMV,IFVA,IFMO, &
          lumin,TMAX,nsmorza,fase,logCNO,kover, nsole)
       use strutture
       real :: TEMPO, DMAS, TMAX
       integer :: NMOD,NABLA,KSD,MAXMEin,MAXMV,IFVA,IFMO,nsmorza,fase
       integer :: logCNO, kover, nsole
       type(luminosity) :: lumin
     end subroutine STAMPA

     subroutine STATE(caller,PR,TE,RHO,GRAD,CSPE,PMOL) 
       character(len=15) :: caller
       real :: PR, TE, RHO, GRAD, CSPE, PMOL
     end subroutine STATE

     subroutine stop_evolut(fase, lumin, nmd, age, logCNO, bump, tei, tef, &
          He, shellH, nsole)
       use strutture
       integer :: fase, nmd, logCNO, bump, shellH, nsole 
       type(luminosity) :: lumin
       real :: age, tei, tef, He
     end subroutine stop_evolut

     subroutine SUPERA(RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCO,RADIAT)
       real :: RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCO,RADIAT
     end subroutine SUPERA

     function SK(RO,T,X,Z1,Z2, cache) 
       use parametri
       real :: sk,RO,T,X(MELE)
       integer :: Z1, Z2, cache
     end function SK

     subroutine VEIOVE(EM,R,EL,P,T,RO,NUMERO,IPP) 
       real :: EM, R, EL, P, T, RO
       integer :: NUMERO, IPP
     end subroutine VEIOVE

     ! ciacioLi
     subroutine ciacioLi(tempo,nmd,maxme,nsmorza, nsole, nfine_diff)
       real :: tempo
       integer :: nmd, maxme, nsmorzansole, nfine_diff
     end subroutine ciacioLi

     real function f(y)
       real :: y
     end function f

     subroutine deriva(max,x,y,dy,nsmo)
       integer :: max, nsmo
       real,dimension(max) :: x, y, dy
     end subroutine deriva

     subroutine smooth(y,npti)
       integer :: npti
       real,dimension(npti) ::  y
     end subroutine smooth

     subroutine smoothb(y,npti)
       integer :: npti
       real,dimension(npti) ::  y
     end subroutine smoothb

     subroutine smoothc(y,npti)
       integer :: npti
       real,dimension(npti) ::  y
     end subroutine smoothc

     subroutine smdia(y,npti)
       integer :: npti
       real,dimension(npti) ::  y
     end subroutine smdia

     subroutine SPLININT(x,y,n,yp1,ypn,num,xa,ya,niscri,erin, upsilon, delta)
       use parametri
       integer :: n, num, niscri
       real,dimension(n) :: x, y
       real,dimension(num) :: xa, ya 
       real :: yp1,ypn
       real,dimension(lim) :: erin, upsilon
       real :: delta
     end subroutine SPLININT
     
     subroutine rnmdna(y,max)
       integer :: max
       real :: y(max)
     end subroutine rnmdna

     subroutine cambia(X,Y)
       use parametri
       real,dimension(lim) :: X,Y
     end subroutine cambia
     
     subroutine integra(dm_dr_v,dm_dr_n,deltadiff,maxme, r, vl, dv_dr)
       use parametri
       real,dimension(lim) :: dm_dr_v, dm_dr_n, r, vl, dv_dr
       real :: deltadiff
       integer :: maxme
     end subroutine integra
     
     subroutine decidi(nmd,maxme,xsave, nsole, nscr, tempomod, dm_dr_v, &
          emutot, r, vl, dv_dr, m, ele)
       use parametri
       integer :: nmd, maxme, nsole, nscr, m
       real :: tempomod
       real,dimension(lim) :: xsave, dm_dr_v, emutot, r, vl, dv_dr
       integer,dimension(m) :: ele
     end subroutine decidi

     subroutine check(nmd,maxme,nscr,tempomod, dm_dr_v, emutot, r, vl, &
          dv_dr, m)
       use parametri
       integer :: nmd, maxme, nscr, m
       real :: tempomod
       real,dimension(lim) :: dm_dr_v, emutot, r, vl, dv_dr
     end subroutine check

     subroutine risulta(nmd,tempo,maxme)
       integer :: nmd, maxme
       real :: tempo
     end subroutine risulta

     subroutine risultati(nmd,tempo,maxme)
       integer :: nmd, maxme
       real :: tempo
     end subroutine risultati
     ! end ciacioLi
 
     ! condint    
     subroutine potekhin(Zion,RLG,TLG,CK, index) 
       real :: Zion, RLG, TLG, CK
       integer :: index
     end subroutine potekhin

     subroutine CONINTER(Zion,TLG,RLG,CK, index) 
       real :: Zion, TLG, RLG, CK
       integer :: index
     end subroutine CONINTER

     subroutine CINTERP3(ZM,Z0,Z1,ZP,Z,N0,MXNV,VM,V0,V1,VP,VF)
       real :: ZM,Z0,Z1,ZP,Z,VM,V0,V1,VP,VF
       integer :: N0, MXNV
     end subroutine CINTERP3

     subroutine HUNT(XX,N,X,JLO) 
       real,dimension(*) :: XX
       integer :: N,JLO
       real :: X
     end subroutine HUNT
     ! end condint
     
     ! routinedoc     
     subroutine DIFFUSION(M,A,Z,X,CL,AP,AT,AX) 
       integer :: M
       real,dimension(M) :: A,Z,X,AP,AT
       real,dimension(M,M) :: AX, CL 
     end subroutine DIFFUSION

     subroutine LUBKSB(A,N,NP,INDX,B) 
       integer :: N, NP
       real,dimension(NP,NP) :: A
       real,dimension(N) :: B 
       integer,dimension(N) :: INDX 
     end subroutine LUBKSB

     subroutine LUDCMP(A,N,NP,INDX,D) 
       real :: D                                                           
       integer :: N, NP
       real,dimension(NP,NP) :: A
       integer,dimension(N) :: INDX
     end subroutine LUDCMP
     ! end routinedoc      

     ! xcotrin
     subroutine opac (z,xh,xxc,xxo,t6,r,err) 
       real :: z, xh, xxc, xxo, t6, r
       integer :: err
     end subroutine opac

     subroutine cointerp(xxc,xxo) 
       real :: xxc, xxo
     end subroutine cointerp

     subroutine t6rinterp(slr,slt) 
       real :: slr, slt
     end subroutine t6rinterp

     subroutine readco
     end subroutine readco

     subroutine SPLINE(X,Y,N,Y2) 
       integer :: N
       real,dimension(N) :: X,Y,Y2
     end subroutine SPLINE

     subroutine SPLINT(XA,YA,N,Y2A,X,Y,YP) 
       integer :: N
       real,dimension(N) :: XA, YA, Y2A 
       real :: X, Y, YP
     end subroutine SPLINT

     subroutine FITY 
     end subroutine FITY

     subroutine FITX
     end subroutine FITX

     subroutine GETD(F,N,D,FP1,FPN) 
       integer :: N
       real,dimension(N) ::  F, D
       real :: FP1, FPN
     end subroutine GETD

     subroutine INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR) 
       real :: FLT, FLRHO, G, DGDT, DGDRHO
       logical :: IERR
     end subroutine INTERP

     subroutine SMOOTHING
     end subroutine SMOOTHING

     subroutine opaltab
     end subroutine opaltab
     ! end xcotrin     

  end interface
end module interfaccia2

!########################################
!Moduli per le variabili realtive alla DM
!########################################

module Dark_Matter
!Modulo per la gestione delle varie variabili relative alla DM
  use parametri !Serve per avere la variabile LIM


  !###############################
  !Per la subroutine di scrittura
  !###############################
  integer,parameter :: ioDark=411 !La uint per il file di log delle epsi della Dark matter
  integer :: first_DM_write=1!Indice per la scrittura in ioDark 

  integer,parameter :: ioDrakError=412 !La uniti per il file che mi salva i vari modelli in cui la T_DM non converge

  real :: rho_DM !GeV/(c^2*cm^3)
  real :: sigma0_DM !Sezione d'urto per DM-idrogeno secondo upper limit Xenon-1ton 
  real :: N_DM_tot !Numero di particelle all'intermo della struttura per passo temporale
  real :: mass_DM !GeV/c^2 Massa della particella di DM
  real,dimension(LIM) :: epsi_DM !Array per contenere l'epsi dovuta alla DM
  real,parameter :: GeV_grammi=1.783e-24 !Costante di conversione da GeV/c^2 a grammi

  
  !Scelgo le velocità media e la velocità di dispersione della stella rispetto al rest frame dell DM
  !Questo serve per la funzione di cattura
  real :: v_star=220d0 !km/s In questo caso scelgo la velocità del Sole
  real :: v_disp=270d0 !Km/s Uso sempre il Sole
end module Dark_Matter

module FUNZIONI
!Modulo per alcune funzioni matematiche necessarie per il calcolo della DM

  !######
  !ERF(x)
  !######
  interface
  real function erf(x) bind(c, name="erf")
    import
    real, value :: x
  end function erf
  end interface
end module FUNZIONI

module Masse_ele
!Modulo per avere a disposizione le varie masse degli elementi all'interno delle varie subroutine
  implicit none

  !Definizione masse elementi, la metto in masse atomiche, la sequenza degli elementi la prendo da chim
    !  1=XH,     2=XHE3,   3=XHE4,    4=XC12,   5=XN14,
    !  6=XO16,   7=XO18,   8=XNE20,   9=XNE22, 10=XMG24,
    ! 11=XMG25, 12=XMG26, 13=XSI28,  14=XNEU,  15=XNA23,
    ! 16=XC13,  17=XN15,  18=XO17,   19=XF19,  20=XNE21,
    ! 21=XFE,   22=XLI6,  23=XLI7,   24=XBE,   25=XB
    ! 26=XD
  real, parameter :: massa_ele(26) = &
  [1.00784, 3.0160293, 4.002602, 12.0, 14.003074, 15.994914, 17.9991610, 19.992440, &
   21.9913855, 23.98504170, 24.9858370, 25.9825930, 27.9769265, 1.0, 22.989770, 13.00335, &
   15.000108899, 16.999131, 18.998403, 20.9938467, 55.845, 6.01512, 7.01600, 9.012182, &
   10.81, 2.014]

  ! Elenco dei vari elementi chimici
  ! H, He3, He4, C12, N14, O16, O18, Ne20, Ne22, Mg24, Mg25, Mg26, 
  ! Si28, Booooo 14, Na23, C 13, N15, O17, F19, Ne21, Fe,
  ! Li6, Li7, Be, B, D

  real,parameter :: u_to_gramms=1.66054e-24!Conversione Unità di massa a grammi
  
end module Masse_ele

module Asplund_per_DM
!Modulo per caricare i valori log10(X_ele/X_H*m_H/m_ele)+12 da asplund 2009

!Viene usato in: Cattura_DM.f90 (sub: Cattura_DM)

  implicit none

  integer :: on_off_DM_Asplund=0!Flag per attivare o disattivare gli elementi più pesanti guardando le abbondanze da asplung

  character(2), dimension(69) :: elementi_asplund=&
    ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', &
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', &
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
    'Ga', 'Ge', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', &
    'Rh', 'Pd', 'Ag', 'In', 'Sn', 'Xe', 'Ba', 'La', 'Ce', 'Pr', &
    'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
    'Lu', 'Hf', 'W', 'Os', 'Ir', 'Au', 'Tl', 'Pb', 'Th']

    real,dimension(69):: dati_asplund_DM=& !log10(N_ele/N_H)+12
    [12.00, 10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93, &
    6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34, &
    3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56, &
    3.04, 3.65, 3.25, 2.52, 2.87, 2.21, 2.58, 1.46, 1.88, 1.75, &
    0.91, 1.57, 0.94, 0.80, 2.04, 2.24, 2.18, 1.10, 1.58, 0.72, &
    1.42, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84, &
    0.10, 0.85, 0.85, 1.40, 1.38, 0.92, 0.90, 1.75, 0.02]

    real,dimension(69) :: numeri_massa_asplund_DM=&! A degli elementi
    [1, 4, 7, 9, 11, 12, 14, 16, 19, 20, &
    23, 24, 27, 28, 31, 32, 35, 40, 39, 40, &
    45, 48, 51, 52, 55, 56, 59, 58, 63, 64, &
    69, 74, 84, 85, 88, 89, 94, 93, 98, 102, &
    103, 106, 107, 115, 120, 131, 136, 139, 140, 141, &
    142, 152, 157, 156, 159, 164, 165, 166, 169, 174, &
    175, 180, 184, 190, 193, 197, 204, 207, 232]
end module Asplund_per_DM