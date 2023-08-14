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

module precision
  integer, parameter :: WP = kind(1.0D0)
end module precision

module fasievolut
  ! Punti di stop definiti:
  ! fase = -10 : MS
  ! fase = -30 : esaurimento H centrale
  ! fase = -50 : RGB
  ! fase = -60 : bump RGB
  ! fase = -90 : flash He
  ! fase = -91 : innesco He
  ! fase = -100 : ripartenza dopo pepper
  ! fase = -110 : esaurimento He centrale
  ! fase = -150: pulsi
  integer,parameter :: fase_MS = -10, fase_exH = -30, fase_RGB = -50
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
  use precision
  ! eqlb, idro, epsi, stampa, mixing
  integer :: L3, KOVER
  real(WP) :: par_OVER
end module overshoot

module intero 
  ! main, atmos, optim
  use parametri
  integer,dimension(LIM) :: IN
end module intero

module mistura
  use precision
! main, io, innes, kappa
 real(WP) :: DEFAUC, DEFAUN, DEFAUO
 real(WP) :: DEFAUFe, DEFAUHe3, DEFAULi6, DEFAULi7, DEFAUBe, DEFAUB, DEFAUD
 integer :: misturaop
end module mistura

module fisica 
  ! main, veiove, atmos, carbon, elio, epsi, epsig, eqlb, evolut, fitta,
  ! hbrezo, henyey, idreli, idro, innes, maslos, mixing, optim,
  ! pastem, plotta, quatm, stampa, risulta, risultati, ciacioLi
  use parametri
  use precision
  real(WP),dimension(6,LIM) :: G, GG
  real(WP),dimension(5,LIM) :: GGV
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
  use precision
  ! main, quatm
  real(WP),dimension(3,4,4) :: U, V 
  real(WP),dimension(4) :: EL, TE 
  real(WP) :: PMA1, PMA2
end module gratm

module strut 
  ! main, ciacioLi, risulta, fato, innes, maslos, mixing, optim
  ! pastem, plotta, quatm, stampa, veiove
  use precision
  real(WP) :: EMTOT, ROCEN, ETAFPR
end module strut

module mesh 
  ! atmos, quatm, stampa 
  integer :: MAXME, MAXMIV
end module mesh

module serchi 
  ! evolut, optim, stampa
  use parametri
  use precison
  real(WP),dimension(2,LIM) :: XSERV
end module serchi

module nuconv 
  ! stampa, mixing, optim, pastem
  use precison
  real(WP),dimension(12) :: BB 
  real(WP) :: GSEMI, GICO
end module nuconv
 
module atmosfere 
  ! main, atm_ck03, atm_bh05, atmos, stampa
  use precision
  real(WP),dimension(3) :: V_ATM 
  real(WP) :: X_SUP, Y_SUP, Z_SUP
  integer :: IBH05_LETTA, ICK03_LETTA, N_BH05, N_CK03, NMOD_PRE
  integer :: do_bh_05, do_ck_03
end module atmosfere

module equil 
  ! epsi, eqlb, idro
  use precision
  real(WP),dimension(5) :: COEFF
  integer :: ISHEL
end module equil

module fitt  
  ! innes, fato, veiove
  use precision
  real(WP) :: ELL,TEF,PCEN,TCEN
end module fitt

module neutri 
  ! epsi, stampa
  use precision
 real(WP),dimension(8) :: PPNEU
 real(WP) :: TCRIT 
end module neutri

module second  
  ! fitta, henyey, resnuc, stampa 
  use parametri
  use precision
  real(WP),dimension(4,LIM) :: ALF, BET, GAM
  real(WP),dimension(5,LIM) :: DG
end module second

module sistema 
  ! henyey, resnuc
  use precision
  real(WP),dimension(4) :: ALF1, BET1, GAM1, E
  real(WP),dimension(4,4) :: B, C
end module sistema

module atm 
  ! fitta, innes, quatm
  use precision
  real(WP) :: RB, PB, TB, DPL, DPT, DRL, DRT, DTL, DTT
end module atm

module chimic 
  ! main, eqlb, optim, atmos, ciacioLi, decidi, risulta, risultati, epsig,
  ! evolut, hbrezo, henyey, innes, mixing, pastem, plotta, stampa, veiove
  ! stop_evolut, epsi
  use parametri
  use precision
  real(WP),dimension(MELE,LIM) :: XXX, XXV
  real(WP),dimension(MELE,2) :: XCE
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
  use precision
  real(WP),dimension(MELE) :: XX
  real(WP) :: HT1, HT1V, EFFEP
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
  use precision
  real(WP) :: ELLOG, TEFF, TEFFV, TEFFVV
end module tempe

module numer 
  ! main, fitta, fato, veiove, innes, quatm, optim, hbrezo, atmos, resnuc,
  ! stampa, henyey, kappa, epsig
  integer :: MAIS, IPRALL, ISUB, IBAT, IREAD, IQUAT, IHP, IFAST
end module numer

module print1
  ! main, stampa
  use precision
  real(WP) :: AHMAX,AHEMAX,AGRMAX,ANUMAX, FG1,FG2,FG3,FG4,FG5,FG6
end module print1

module varie 
  ! main, atmos, epsi, innes, kappa, maslos, quatm, stampa, supera, stop_evolut
  use precision
  real(WP) :: XH, XME, HELI, HECE, ALFA, BMAG, EMMA, FRAZ, EMAXH, EMAHE
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
  use precision
  real(WP),dimension(mx,mc) :: oxf, cxf, xcdf, xodf, cxdf, oxdf
  real(WP),dimension(mx,nt,nr) :: opl
  integer,dimension(mx) :: itime = (/0/)                
end module aaa

module aa 
  ! opac, cointerp, t6rinterp, readco
  use parametri_xcotrin
  use precision
  real(WP),dimension(mc) :: xcd, xod, xc, cxd, cx 
  real(WP),dimension(mo) :: xo, oxd, ox
  real(WP),parameter,dimension(mc) :: xcs = (/0.0,0.01,0.03,0.1,0.2,0.4, &
       0.6,1.0/)
  real(WP),parameter,dimension(mc) :: xos = (/0.0,0.01,0.03,0.1,0.2,0.4, &
       0.6,1.0/)
  real(WP) ::  zzz, xxh
  real(WP),dimension(mx) :: xx
  real(WP),dimension(4) :: h, q
  integer :: nc,no
end module aa

module a1 
   ! opac, cointerp, t6rinterp, readco
  use parametri_xcotrin
  use precision
  real(WP),dimension(mx,mc,mo,nt,nr,10) :: covett 
  real(WP),dimension(mx,mc,nt,nr,10) :: diagvett
  real(WP),dimension(nt) :: t6list, alt, dfs 
  real(WP),dimension(nr) :: alr, dfsr

  real(WP),dimension(mx,mo,nt,nr,10) :: diagovett
  real(WP),dimension(nt,nr) :: opk
  real(WP),dimension(3,mx) :: a
  real(WP),dimension(3) :: b
  
  real(WP),parameter,dimension(8) :: xa = (/0.0,0.03,0.1,0.35,0.7,0.,0.,0./)
  
  real(WP),dimension(nrm) :: alrf 
  real(WP),dimension(ntm,nrm) :: cof 
  real(WP),dimension(ntm) :: t6listf
  real(WP),dimension(nt,nr) :: opk2 
  real(WP),dimension(mx) :: dfsx

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
  use precision
  integer :: NTEMP, NSM, nrlow, nrhigh
  real(WP) :: RLE
  real(WP),dimension(100) :: t6arr 
  real(WP),dimension(100,nr) ::  coff
end module alink

module cst  
  ! fity, fitx, interp, smoothing, opaltab, readco
  use precision
  real(WP) :: tmax, RLS
  integer :: NRL, nset
end module cst

module d 
  ! opac, t6rinterp, quad
  use precision
  real(WP) :: dkap
end module d

module bb  
  ! opac, cointerp, t6rinterp
  use precision
  real(WP) :: xodp, xcdp, xxco, cxx, oxx
  integer :: l1, l2, l3, l4, k1, k2, k3, k4, ip, iq
end module bb

module b1  
  ! opac, readco
  use parametri_xcotrin
  use precision
  real(WP),dimension(mx,ntabs) :: x, y, zz, xca, xoa
  integer,dimension(mx,ntabs) :: itab
  integer,parameter,dimension(nrm) :: nta = (/70,70,70,70,70,70,70,70,70,70, &
       70,70,70,70,69,64,60,58,57/)
end module b1

module cf  
  ! smoothing, opaltab, fity, fitx, interp
  use parametri_xcotrin
  use precision
  real(WP),dimension(85,IPR) :: F, FX, FY, FXY
end module cf

module e  
  ! opac, t6rinterp, readco, kappa
  use precision
  real(WP) :: opact,dopact,dopacr,dopactd
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
  use precision
  type pow_of_T
     real(WP) :: T9, T92, T93, T912, T913, T923, T932, T6, T612
  end type pow_of_T

  type luminosity
     real(WP) :: Lpp, Lcno, L3a, Lgra, L
     real(WP) :: teff
     real(WP) :: M_bordo_conv, nuclear
     real(WP),dimension(3) :: pn
  end type luminosity
end module strutture

module interfaccia
  use strutture
  ! external procedures
  interface
     subroutine ATM_BH05(GIINI,TE,TAU_RACC) 
       use precision
       real(WP) :: GIINI, TE, TAU_RACC 
     end subroutine ATM_BH05
     
     subroutine ATM_CK03(GIINI,TE,TAU_RACC) 
       use precision
       real(WP) :: GIINI,TE,TAU_RACC
     end subroutine ATM_CK03
     
     subroutine ATMOS(FLUMIN,TEFF,P,T,R,EMASS,ISCR,IRA) 
       use precision
       real(WP) :: FLUMIN, TEFF, P, T, R, EMASS
       integer :: ISCR, IRA
     end subroutine ATMOS
     
     subroutine CARBON (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S) 
       use parametri
       use precision
       real(WP) :: RO, T, PAS
       integer :: JF, INI, IFI, N
       real(WP),dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real(WP),dimension(32) :: S
     end subroutine CARBON
     
     subroutine CROSS (RO,T,XXn,IR,ISHORT,S,powt,nabla)
       use strutture
       use parametri
       use precision
       real(WP) :: RO, T
       real(WP),dimension(mele) :: xxn
       integer :: IR, ISHORT, nabla
       real(WP),dimension(32) :: S
       type(pow_of_T) :: powt
     end subroutine CROSS

     subroutine CUB(A,B,ALFA,BETA,GAMMA,DELTA) 
       use precision
       real(WP),dimension(4) :: A, B
       real(WP) :: ALFA, BETA, GAMMA, DELTA
     end subroutine CUB

     subroutine CUB2dr(A, B, idx, xABGD) 
       use precision
       real(WP),dimension(4) :: B
       real(WP),dimension(4,4) :: A
       real(WP),dimension(320) :: xABGD
       integer :: idx
     end subroutine CUB2dr

     subroutine CUB2dt(A, B, xtABGD, is, ie) 
       use precision
       real(WP),dimension(4) :: B
       real(WP),dimension(80) :: xtABGD, A
       integer :: is, ie
     end subroutine CUB2dt

     subroutine ELIO (RO,T,PAS,JF,INI,IFI,N,ISHORT, YV, YF, NEQUI, S) 
       use parametri
       use precision
       real(WP) :: RO, T, PAS
       integer :: JF, INI, IFI, N, ISHORT
       real(WP),dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real(WP),dimension(32) :: S
     end subroutine ELIO

     subroutine EPSI (RO,T,JF,NABLA,HT,EPS,EPSA, ECAR, DELTA) 
       use precision
       real(WP) :: RO, T, HT, EPS, EPSA, ECAR, DELTA
       integer :: JF, NABLA
     end subroutine EPSI

     subroutine EPSIG(P,T,RO,CP,HT,JF,GRAVI,MAXMV,DAD,IPP) 
       use precision
       real(WP) :: P, T, RO, CP, HT, GRAVI, DAD
       integer :: JF, MAXMV, IPP
     end subroutine EPSIG
     
     subroutine EQLB(MAXME) 
       integer :: maxme
     end subroutine EQLB

     subroutine EVOLUT(INI, IFI) 
       integer :: INI, IFI
     end subroutine EVOLUT

     subroutine akappetta(Y,xC,xO,tlog,rlog,cap) 
       use precision
       real(WP) :: Y,xC,xO,tlog,rlog,cap
     end subroutine akappetta

     subroutine FATO(INDU,AS,AT,AU,AV,IPP, EL0,EL1,EL2,P0,P1,P2,R0,R1,R2, &
          T0,T1,T2,DPC,DTC,DELL,DTEF,EL0P,EL1P,EL2P,P0P,P1P,P2P,R0P,R1P,R2P, &
          T0P,T1P,T2P,CORZ) 
       use precision
       real(WP) :: AS, AT, AU, AV
       integer :: INDU, IPP
       real(WP) :: EL0, EL1, EL2, P0, P1, P2, R0, R1, R2, T0, T1, T2, DPC,DTC, &
       DELL, DTEF, EL0P, EL1P, EL2P, P0P, P1P, P2P, R0P, R1P, R2P, &
       T0P, T1P, T2P, CORZ
     end subroutine FATO

     subroutine FITTA(M,DAM,NUM1,NUM2,FRAT)
       use precision
       integer :: M, NUM1, NUM2
       real(WP) :: DAM, FRAT
     end subroutine FITTA

     subroutine HBREZO (AMCO, YCO, CCO, NCO, OCO)
       use precision
       real(WP) :: AMCO, YCO, CCO, NCO, OCO
     end subroutine HBREZO
     
     subroutine HENYEY(MAXME, MAXMV, ERROR, INUM, IERR) 
       use precision
       real(WP) :: ERROR
       integer :: MAXME, MAXMV, INUM, IERR
     end subroutine HENYEY

     subroutine KAPPA(RHO,TEMP,CAP,Y) 
       use precision
       real(WP) :: RHO, TEMP, CAP, Y
     end subroutine KAPPA

     subroutine IDRELI (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S) 
       use parametri
       use precision
       real(WP) :: RO, T, PAS
       integer :: JF, INI, IFI, N
       real(WP),dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real(WP),dimension(32) :: S
     end subroutine IDRELI

     subroutine IDRO (RO,T,PAS,JF,INI,IFI,N, YV, YF, NEQUI, S, oish) 
       use parametri
       use precision
       real(WP) :: RO, T, PAS
       integer :: JF, INI, IFI, N, oish
       real(WP),dimension(MELE) :: YV, YF 
       integer,dimension(MELE) :: NEQUI
       real(WP),dimension(32) :: S
     end subroutine IDRO

     subroutine INNES(ITCC,IPP,MAXMEin,ILEG,LBA, YB,ZB, He, Zeta, Alpha, &
          do_relax) 
       use precision
       real(WP) :: YB, ZB, He, Zeta, Alpha
       integer :: ITCC, IPP, MAXMEin, ILEG, LBA
       logical :: do_relax
     end subroutine INNES
     
     subroutine LOCAL(P,T,R,EL,EM,RO,GRAD,CP,PMU,ACCA,CAP,RADIAT,ADIAB, &
          GI,Y,GRSAD)
       use precision
       real(WP) :: P,T,R,EL,EM,RO,GRAD,CP,PMU,ACCA,CAP,RADIAT,ADIAB, &
            GI,Y,GRSAD
     end subroutine LOCAL

     subroutine MASLOS(MAXME,DM,nmd,iread, fase, fase_start, fase_stop) 
       use precision
       integer :: MAXME,nmd,iread, fase, fase_start, fase_stop
       real(WP) :: DM
     end subroutine MASLOS

     subroutine MIXING(MAXME,iread)
       integer :: MAXME,iread
     end subroutine MIXING

     subroutine NEUTR(RO,T,QT) 
       use precision
       real(WP) :: RO, T, QT
     end subroutine NEUTR

     subroutine OPTIM(MAXME,SCALA,PROVV,ECNO, URA,ULA,UPA,UTA,UMA,VMM,fase) 
       use precision
       integer :: MAXME,fase
       real(WP) :: SCALA, PROVV, ECNO
       real(WP),dimension(4) :: URA, ULA, UPA, UTA, UMA
       real(WP),dimension(5) :: VMM
     end subroutine OPTIM

     subroutine PASTEM(MAXME,MAXMV,PROV,PROVV,ECNO,ECNOV,FG6,FG6V,TMAX, &
          div_pastem,max_pastem)
       use precision
       integer :: MAXME, MAXMV
       real(WP) :: PROV, PROVV, ECNO, ECNOV, FG6, FG6V, div_pastem, max_pastem,&
            TMAX
     end subroutine PASTEM

     subroutine PLOTTA(MAXME,RTOT) 
       use precision
       integer :: MAXME
       real(WP) :: RTOT
     end subroutine PLOTTA

     subroutine QUATM(NABLA,INI,MAXMEin,MAXMV,EMTOV) 
       use precision
       integer :: NABLA, INI, MAXMEin, MAXMV
       real(WP) :: EMTOV
     end subroutine QUATM

     subroutine RESNUC(JFF)                                             
       integer :: jff
     end subroutine RESNUC

     subroutine SIMQ(A,B,N,KS) 
       use precision
       real(WP),dimension(361) :: A
       real(WP),dimension(19) :: B
       integer :: N, KS
     end subroutine SIMQ

     subroutine STAMPA(TEMPO,DMAS,NMOD,NABLA,KSD,MAXMEin,MAXMV,IFVA,IFMO, &
          lumin,TMAX,nsmorza,fase,logCNO,kover)
       use strutture
       use precision
       real(WP) :: TEMPO, DMAS, TMAX
       integer :: NMOD,NABLA,KSD,MAXMEin,MAXMV,IFVA,IFMO,nsmorza,fase
       integer :: logCNO, kover
       type(luminosity) :: lumin
     end subroutine STAMPA

     subroutine STATE(PR,TE,RHO,GRAD,CSPE,PMOL) 
       use precision
       real(WP) :: PR, TE, RHO, GRAD, CSPE, PMOL
     end subroutine STATE

     subroutine stop_evolut(fase, lumin, nmd, age, logCNO, bump, tei, tef, &
          He, shellH)
       use strutture
       use precision
       integer :: fase, nmd, logCNO, bump, shellH 
       type(luminosity) :: lumin
       real(WP) :: age, tei, tef, He
     end subroutine stop_evolut

     subroutine SUPERA(RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCO,RADIAT)
       use precision
       real(WP) :: RO,P,T,CAP,PMU,CP,ADIAB,GI,GRSAD,ACCO,RADIAT
     end subroutine SUPERA

     function SK(RO,T,X,Z1,Z2, cache) 
       use parametri
       use precision
       real(WP) :: sk,RO,T,X(MELE)
       integer :: Z1, Z2, cache
     end function SK

     subroutine VEIOVE(EM,R,EL,P,T,RO,NUMERO,IPP) 
       use precision
       real(WP) :: EM, R, EL, P, T, RO
       integer :: NUMERO, IPP
     end subroutine VEIOVE

     ! ciacioLi
     subroutine ciacioLi(tempo,nmd,maxme,nsmorza, nsole, nfine_diff)
       use precision
       real(WP) :: tempo
       integer :: nmd, maxme, nsmorza, nsole, nfine_diff
     end subroutine ciacioLi

     real function f(y)
       use precision
       real(WP) :: y
     end function f

     subroutine deriva(max,x,y,dy,nsmo)
       use precision
       integer :: max, nsmo
       real(WP),dimension(max) :: x, y, dy
     end subroutine deriva

     subroutine smooth(y,npti)
       use precision
       integer :: npti
       real(WP),dimension(npti) ::  y
     end subroutine smooth

     subroutine smoothb(y,npti)
       use precision
       integer :: npti
       real(WP),dimension(npti) ::  y
     end subroutine smoothb

     subroutine smoothc(y,npti)
       use precision
       integer :: npti
       real(WP),dimension(npti) ::  y
     end subroutine smoothc

     subroutine smdia(y,npti)
       use precision
       integer :: npti
       real(WP),dimension(npti) ::  y
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
end module interfaccia
