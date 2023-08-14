!      subroutine instruct
!-----VERSION of November 20, 1995-----------------------------------------
!-----------------------------------------------------------------------
!
!          This subroutine contains instructions for using the subroutine
!    OPAC( z, xh, xxc, xxo, t6, r ).
!          The purpose of the subroutine OPAC is to perform 4 or 5 variable
!    interpolation on  log10(kappa).  The opacity tables to be interpolated
!    are known to have somewhat random numerical errors of a few percent.
!    Consequently, adjusting the data prior to performing the interpolation
!    is justified at this level.  The code first reads the original(unsmoothed)
!    tabular data, this data is then passed through a smoothing filter; using
!     a set of routines developed by Mike Seaton (see M.J. Seaton,MNRAS 265,
!  L25,1993). It is the adjusted data that is actually used in carrying out the
!     interpolations in OPAC.  Steps are taken, as described below to insure
!     that the interpolated data is also smooth. The initial adjustment step
!     helps improve the smoothness of the OPAC output,particularly at the
!     smallest values of R. The medium to large R output is only slightly
!   effected by this step. It takes only a few seconds to carryout the initial
!     data smoothing step.  This step can be skipped by setting the parameter
!     ismdata=1 in subroutine readco.

!          The interpolation variables are :

!          xh     The hydrogen mass fraction, X
!          xxc    The enhanced carbon mass fraction, delta Xc.  
!                 The total C mass fraction, Xc, is the sum of the initial 
!                 amount included in the metal mass fraction, Z, and 
!                 delta Xc
!          xxo    The enhanced oxygen mass fraction, delta Xo
!          t6     The temperature in millions of degrees Kelvin, T6
!          r      =rho(gm/cc)/T6**3, R

!     Additional input to OPAC is:

!           z      The metallicity, Z

!          An interpolation between overlapping quadratics is used to obtain
!     smoothed results.  A 4x4 grid in logT6 and logR is used to interpolate
!     in four different 3x3 sub-grids. Linear interpolation between quadratic
!     fits in these different  sub-grids gives smoothed results in both log T6
!     and Log R directions. This procedure produces results that are similar
!     to bicubic spline interpolation, but require storage of only local 
!     information.  
!     Each of the individual tables in a file Gx**x**z covers 70 temperatures
!     in the range logT=3.75[T6=0.0056341325]( referred to as temperature 1) to
!     logT=8.7[T6=501.187] and 19 values of log R in the range -8 (referred to 
!     as 1) to +1; at half-integer steps.  
!     (NOTE: earlier tables were explicitly in
!     terms of T6. For convenience the present tables tabulate log Kappa vs 
!     logT. The
!     interpolation however still uses T6 for the temperature variable)
!     For specialized problems, if storage space is a problem, a reduced set of
!     data can be input.  This requires a recompilation with altered parameter
!     values. In order to limit the range of R, set the parameter nrb= index of
!     lowest value of log R to use(count from log R=-8). Then set the parameter
!     nre to the index of the largest value of log R to use. 
!     (NOTE: nre-nrb must
!     be at least 5). To ignore the first ntb-1 values of T6 (starting from
!     temperature 1) set the parameter ntb to desired value. 
!     No other parameters should be modified.
!     A five variable interpolation is done when X is greater than 0.  In this
!     case it is assumed that variable values of X are also needed.  
!     We have provided sets of tables for X=0, 0.03, 0.1, 0.35 and 0.7, 
!     for each of the metallicities 0.0,0.001, 0.004,0.01,0.02,0.03,05 and .1.
!     The five sets of tables associated with
!     a particular value of Z (eg. Gx0z01, Gx03x01,Gx10z01, Gx35z01,
!     Gx70z01) should be placed in files named codataa, codatab, codatac, 
!     codatad, and codatae, respectively.  To create storage for this data 
!     the routines must be recompiled with the parameter mx=5.  
!     Again if storage is a problem the T6
!     and log R ranges can be restricted.  This version of the interpolation 
!     routine does not interpolate in Z.  Values of Z not in the table 
!     can be obtained by interpolating the existing tables to produce 
!     similar tables for the Z of interest.
!     Interpolation in xh,xxo,xxc,t6, and r can be obtained as just described 
!     A 4 variable interpolation in xxc, xxo, t6 and r is performed in 
!     the special case  X=0.  
!     The set of 60 data tables in (xxc,xxo) for a given Z that have
!     been provided for X=0 should be placed in a file called 'codataa'.  
!     This file
!     will be read from unit 2 in the subroutine readco.  
!     In this special case the
!     set of routines provided should be compiled with the parameter mx=1 
!     (there are 4 occurances).  
!     (NOTE: The version of the code described above, intended
!     for X> 0, also  handles this case, but takes more storage space 
!     since mx=5).
!     If you want to work with a single value of X which is not zero, i.e.,
!     X=0.03, 0.1, 0.35, or 0.70, then compile the code with mx=1 but change 
!     the statement

!          data (xa(i), i=1,5 )/0.0, 0.03, 0.1, 0.35, 0.7/

!     If for example you want to use only the table for X=.70, then

!          data (xa(i), i=1,5 )/0.70, 0.03, 0.1, 0.35, 0.0/

!     You also need to place the tables for X=0.7 into a file named
!     codataa.
!     ***CAUTION***
!       As a result of the mixing procedure used to calculate the data a few
!     X=0.0, low T-small R, table values fell outside the range of T and R 
!     accessible
!     from the X=0.35 data directly calculated for this purpose.  These T-R 
!     locations are filled in with 9.99 (or for diagnostic purposes in some 
!     cases larger values.  At the same locations the derivatives are set 
!     to 99.9.  When T-R falls in this region a message is issued by 
!     the interpolation code to 
!     inform the user of this situation.  Presumable very few users will have
!     applications that take them into this region.

!          Your routine that performs the call to OPAC should include the
!      statement:

!         use e
!         which contains:      opact,dopact,dopacr,dopactd
!         These variables have the following meanings:

!         OPACT       Is the Log of the Rosseland mean opacity: Log(kappa)  
!         DOPACT      Is Log(kappa)/Log(T6)   ! at constant R
!         DOPACTD     Is Log(kappa)/Log(T6)   ! at constant density
!         DOPACR      Is Log(kappa)/Log(R),   ! at constant T6

!***********************************************************************

subroutine opac (z,xh,xxc,xxo,t6,r, err) 
  use interfaccia
  use aaa
  use aa
  use a1
  use d
  use bb
  use b1
  use e
  use indicefile
  use recoin

  !..... The purpose of this subroutine is to interpolate log kappa       
  !          and obtain smooth derivatives.                               
  !      in C/O abundance and in T6, R,i.e. (xc,xo,T6,R)                  
  !      xxc=Xc=carbon mass fraction                                      
  !      xxo=Xo=oxygen mass abundance                                     
  !      t6=T6=temperature in millions of degrees kelvin                  
  !      r=R=density(g/cm**3)/T6**3                                       
  !..... to use opac insert   use e    in the calling routine.             
  !      This module contains interpolated values for kappa and its       
  !      first derivities.                                                
  !                                                                       

  implicit none

  real :: z, xh, xxc, xxo, t6, r
  integer :: err
  save 

  !..... OPACT- opacity obtained from a quadraric interpolation at        
  !      fixed log T6 at three values of log R; followed by quadratic     
  !      interpolation along log T6. Results smoothed by mixing           
  !      overlapping quadratics.                                          
  !..... DOPACT- is Log(k)/Log(T6) smoothed by mixing quadratics.       
  !..... DOPACR- is  Log(k)/Log(R) smoothed by mixing quadratics.       

  real :: xxx, slt, slr, t6o, xho, ro

  integer :: iop, ilo, ihi, mg, mh, mi, mf2,  ntd, k3s, l3s, imd
  integer :: kmin, k1in, iadvance, mfin, is, iw, ir, it, i, isset
  integer :: istep1, ntlimit
  real :: xxci, xxoi, xxi, t6i, ri, cmod, xhe, dixr, quad, xxcs, xxos
  integer :: firsttime = 1

  if(indice == 1)   z = 0.0 
  if(indice == 2)   z = 1.d-4 
  if(indice == 3)   z = 3.d-4 
  if(indice == 4)   z = 1.d-3 
  if(indice == 5)   z = 2.d-3 
  if(indice == 6)   z = 4.d-3 
  if(indice == 7)   z = 1.d-2 
  if(indice == 8)   z = 2.d-2 
  if(indice == 9)   z = 3.d-2 
  if(indice == 10)  z = 4.d-2 

  !..... INDEX refers to the index,i, of xc(i) or xo(i); xc(i) and xo(i)  
  !      are the abundance grid points.                                   
  ! provides smoothed interpolations; iop=0 gives no smoothing
  iop = 1
  if(nr < 6) then
     write(*,'("Too few R values; NRE+1-NRB < 6")') 
     write(66,*)'500 - opac'
     write(66,'("Too few R values; NRE+1-NRB < 6")') 
     stop 
  endif

  if(xh > 1.d-6 .and. mx < 4) then
     write(*,'(" X not equal to zero: To run this case it &
          &is necessary"/ "to recompile with parameter (mx=5)")')  
     write(66,*)'501 - opac'
     write(66,'(" X not equal to zero: To run this case it &
          &is necessary"/ "to recompile with parameter (mx=5)")')
     stop 
  endif

  !..... set-up C/O axis points                                           
  xxco = xxc+xxo 

  if(z+xh+xxco-1.d-6 > 1  .and. ireadco == 12345678) then
     ! restores input value; required if stop replaced with a return     
     xxc = xxci 
     xxo = xxoi 
     !controllo per la kappa                                
     icntrl = 0 
     return 
  endif

  ! giada
  if (firsttime == 1) then 
     alr(1) = -8.0+(nrb-1)*0.5 

     do i=2,nr 
        alr(i) = alr(i-1)+0.5 
     enddo
     alt(1) = -2.25+(ntb-1)*0.05 

     do i=ntb+1,46 
        alt(i) = alt(i-1)+.05 
     enddo
     ntd = 47 
     if (ntb+1 > 47) ntd = ntb+1 
     do i=ntd,68 
        alt(i) = alt(i-1)+.1 
     enddo
     do i=68,70 
        alt(i) = alt(i-1)+.2 
     enddo

     do i=1,nt 
        t6list(i) = 10.**alt(i) 
     enddo
     firsttime = 0
  endif

  zzz = z+0.001 
  xxh = xh 
  xxci = xxc 
  xxoi = xxo 
  xxi = xh 
  t6i = t6 
  ri = r 
  !                                                                       
  !..... convert xxc and xxo to logarithmic shifted by Z                  
  cxx = log10(zzz+xxc) 
  oxx = log10(zzz+xxo) 
  ! giada                                                                 
  if(xh /= xho)  then 
     xxx = log10(0.005+xh) 
     xho = xh 
  endif
  if(t6 /= t6o) then 
     slt = log10(t6) 
     t6o = t6 
  endif
  if(r /= ro) then 
     slr = log10(r) 
     ro = r 
  endif

  !..... set X indices                                                    
  ilo = 2 
  ihi = mx 
  do while(ihi-ilo > 1) 
     imd = (ihi+ilo)/2 
     if(xh <= xa(imd)+1.d-7) then 
        ihi = imd 
     else 
        ilo = imd 
     endif
  end do
  i = ihi 
  mf = i-2 
  mg = i-1 
  mh = i 
  mi = i+1 
  mf2 = mi 
  istep1=1 
  if (mx > 1) then 
     istep1 = mx-1 
     if((xh <= xa(2)+1.d-7) .or. (xh >= xa(istep1)-1.d-7)) mf2 = mh 
  endif

  if (mx == 1 .or. xh < 1.d-6) then 
     mf = 1 
     mg = 1 
     mh = 1 
     mi = 2 
     mf2 = 1 
  endif

  ilo = 2 
  ihi = nr 
  do while(ihi-ilo > 1) 
     imd = (ihi+ilo)/2 
     if(slr <= alr(imd)+1.d-7) then 
        ihi = imd 
     else 
        ilo = imd 
     endif
  end do
  !ha trovato le 4 densita`                      
  i = ihi 
  l1 = i-2 
  l2 = i-1 
  l3 = i 
  l4 = l3+1 
  !                                                                       
  ilo = 2 
  ihi = nt 
  do while(ihi-ilo > 1)
     imd = (ihi+ilo)/2 
     if(t6 <= t6list(imd)+1.d-7) then 
        ihi = imd 
     else 
        ilo = imd 
     endif
  end do
  !ha trovato le 4 temperature                    
  i = ihi 
  k1 = i-2 
  k2 = i-1 
  k3 = i 
  k4 = k3+1 
  k3s = k3+ntb-1 
  l3s = l3+nrb-1 

  !-----set-up indices to skip when low T&R data missing for X=0.         
  kmin = 0 
  k1in = k1 
  iadvance = 0 
  mfin = mf 
  if ((mfin  ==  1) .and. (covett(1,1,1,k1,l1,indice) > 9.)) then 
     do i=1,6 
        if (covett(1,1,1,i,l1,indice) > 9.)  then 
           if (xh < .1) then 
              kmin = i+1 
           else 
              ! sfift X index to avoid X=0.    
              if (iadvance  ==  0) then 
                 iadvance = iadvance+1 
                 mf = mf+1 
                 mg = mg+1 
                 mh = mh+1 
                 mi = mi+1 
                 mf2 = mf2+1 
              endif
           endif
        endif
     enddo
     if (iadvance  ==  0 .and. k1 <= kmin .and. slt <= alt(kmin)) then
        k1 = kmin 
        if (covett(1,1,1,kmin,l1+1,indice) < 9. .and.                &
             (slr+.01) > alr(l1+1)) then                         
           l1 = l1+1 
           kmin = 0 
           k1 = k1in 
           do i=1,6 
              if (covett(1,1,1,i,l1,indice) > 9.) kmin = i+1 
           enddo
           if (kmin /= 0 .and. k1in < kmin) k1 = kmin 
        endif
     endif
     if ((slt+.001) < alt(k1)) then 
        write (*,'("OPAL data not available for X=", f7.5," logT6=", f7.3, &
             &" logR=",f7.3)') xh,slt,slr                            
        opact = 30. 
        dopact = 99. 
        dopacr = 99. 
        dopactd = 99. 
        return 
     endif
     l2 = l1+1 
     l3 = l2+1 
     l4 = l3+1 
     l3s = l3+nrb-1 
     k2 = k1+1 
     k3 = k2+1 
     k4 = k3+1 
     k3s = k3+ntb-1 
  endif
  !-----end of check for missing data                                     
  do m=mf,mf2 
     if(mx >= 4) then 
        ! C and O  fractions determined by the ray through the origin that
        ! also passes through the point (Xc,Xo). Specific interpolation   
        ! values determined by tabulated X values;i.e. xa(m).  Inter-     
        ! polation along the ray gives log (kappa(Xc,Xo)).  (Advantage    
        ! of method: keeps indices within table boundaries)               
        ! Subtract Z to prevent out-of-range C+O values for small X        
        if(1.-xh-z > 1.d-6) then 
           cmod = (1.-xa(m)-z)/(1.-xh-z) 
        else 
           cmod = 0. 
        endif
        xxc = cmod*xxci 
        xxo = cmod*xxoi 
        cxx = log10(zzz+xxc) 
        oxx = log10(zzz+xxo) 
     endif

     do i=1,mc 
        xhe = 1.-xa(m)-z 
        nc = i 
        no = i 
        xc(i) = xcs(i) 
        xo(i) = xos(i) 
        if(xcs(i) > xhe) then 
           xc(i) = xhe 
           xo(i) = xhe 
           exit
        endif
     enddo

     if(itime(m) /= 12345678) then 
        itime(m) = 12345678 
        mxzero = 0 
        if(isset /= 1) then
           do i=1,mx 
              xx(i) = log10(0.005+xa(i)) 
              if(xa(i)  ==  0.0) mxzero = i 
           end do
           isset = 1
        endif
        ! this is the first time throught this m. Calculate the decadic    
        ! log of the perimeter points shifted by Z+0.001(to avoid divergenc
        ! at origin); m refers to xa(m); the hydrogen table value.         
        ! note that the nc-th elements are sometimes needed!               
        do i=1,nc 
           oxf(m,i) = log10(zzz+xo(i)) 
           ! giada 
           if(xc(i) /= xo(i)) then 
              cxf(m,i) = log10(zzz+xc(i)) 
           else 
              cxf(m,i) = oxf(m,i) 
           endif
           xcdf(m,i) = -xo(i)+xo(no) 
           xodf(m,i) = -xc(i)+xc(nc) 
           cxdf(m,i) = log10(zzz+xcdf(m,i)) 
           ! giada     
           if(xcdf(m,i) /= xodf(m,i)) then 
              oxdf(m,i) = log10(zzz+xodf(m,i)) 
           else 
              oxdf(m,i) = cxdf(m,i) 
           endif
        enddo

        !      note that the nc-th elements are sometimes needed! 
        do i=1,nc 
           ox(i) = oxf(m,i) 
           cx(i) = cxf(m,i) 
           xcd(i) = xcdf(m,i) 
           xod(i) = xodf(m,i) 
           cxd(i) = cxdf(m,i) 
           oxd(i) = oxdf(m,i) 
        enddo

        !.....read the data files                                             
        call readco 
        !aggiunto per poter chiamare piu' volte la readco    
        itime(m) = 0 
     endif

     do i=1,nc 
        ox(i) = oxf(m,i) 
        cx(i) = cxf(m,i) 
        xcd(i) = xcdf(m,i) 
        xod(i) = xodf(m,i) 
        cxd(i) = cxdf(m,i) 
        oxd(i) = oxdf(m,i) 
     enddo

     !..... Determine log R and log T6 grid points to use in the             
     !      interpolation.                                                   
     if(slt < alt(1) .or. slt > alt(nt)) then
        ! restores input value; required if stop replaced  with a return    
        xxc = xxci 
        xxo = xxoi 
        !nelle zone in cui la OPAC va fuori dalle tabelle domina
        ! comunque l'opacita` conduttiva                        
        OPACT = 0.0 
        err = 1
        return 
     endif

     if(slr < alr (1) .or. slr > alr(nr))  then
        ! restores input value; required if stop replaced  with a return    
        xxc = xxci 
        xxo = xxoi 
        !nelle zone in cui la OPAC va fuori dalle tabelle domina
        ! comunque l'opacita` conduttiva                        
        OPACT = 0.0 
        err = 1
        return 
     endif

     !  calculate table indices                   
     if (m  ==  mf) then 
        if(mf2 /= mxzero .and. k3s > ntm)  then
           ! restores input value; required if stop replaced  with a return    
           xxc = xxci 
           xxo = xxoi 
           !nelle zone in cui la OPAC va fuori dalle tabelle domina
           ! comunque l'opacita` conduttiva                        
           OPACT = 0.0 
           err = 1
           return 
        endif
        do i=14,18 
           if(l3s > i .and. k3s > nta(i+1))  then
              ! restores input value; required if stop replaced with a return 
              xxc = xxci 
              xxo = xxoi 
              !nelle zone in cui la OPAC va fuori dalle tabelle domina
              ! comunque l'opacita` conduttiva                        
              OPACT = 0.0 
              err = 1
              return 
           endif
        enddo
        ip = 3 
        iq = 3 
        ntlimit = nta(l3s) 
        if(k3s  ==  ntlimit .or. iop  ==  0) then 
           ip = 2 
           iq = 2 
        endif
        if(t6 <= t6list(2)+1.d-7 .or. iop  ==  0) ip = 2 

        ! right edge of full table
        if(l3  ==  nr .or. iop  ==  0) then 
           iq = 2 
           ip = 2 
        endif
        if(slr <= alr(2)+1.d-7 .or. iop  ==  0) iq = 2 
     endif

     xodp = max(-xxc+xc(nc),0.) 
     xcdp = max(-xxo+xo(no),0.) 
     is = 0 

     call cointerp(xxc,xxo) 
  end do

  if((zz(mg,1) /= zz(mf,1)) .or. (zz(mh,1) /= zz(mf,1))) then 
     write(*,'("Z does not match Z in codata* files you are using")')   
     write(*,*) "sono nell'iffone" 
     write(66,*)'502 - opac'
     write(66,'("Z does not match Z in codata* files you are using")')
     stop 
  endif
  if(z /= zz(mf,1)) then
     write(*,'(" Z does not match Z in codata* files you are using")') 
     write(66,*)'503 - opac'
     write(66,'("Z does not match Z in codata* files you are using")')
     stop
  endif

  ! restores input value; necessary if stop replaced with return       
  xxc = xxci 
  xxo = xxoi 
  is = 0 
  iw = 1 
  do ir=l1,l1+iq 
     do it=k1,k1+ip 
        if(mx == 1 .or. mf2 == 1) then 
           opk(it,ir) = opl(mf,it,ir) 
        else
           opk(it,ir) = quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir)          &
                ,opl(mh,it,ir),xx(mf),xx(mg),xx(mh))      
           is = 1 
        endif
     end do
  end do

  ! interpolate between quadratics           
  if (mi  ==  mf2) then 
     is = 0 
     iw = 1 
     dixr = (xx(mh)-xxx)*dfsx(mh) 
     do ir=l1,l1+iq 
        do it=k1,k1+ip 
           opk2(it,ir) = quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir)       &
                ,opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
           opk(it,ir) = opk(it,ir)*dixr+opk2(it,ir)*(1.-dixr) 
           is = 1 
        enddo
     end do
     !     interpolate X between two overlapping quadratics          
  endif
  is = 0 

  !..... completed H,C,O interpolation. Now interpolate T6 and log R on a 
  !      4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure 
  !      mixes overlapping quadratics to obtain smoothed derivatives.     
  !                                                                       
  call t6rinterp(slr,slt) 
  return 

end subroutine opac

!***********************************************************************
subroutine cointerp(xxc,xxo) 
  use interfaccia
  use aaa
  use aa
  use a1
  use bb
  use indicefile
  
  implicit none
  !     The purpose of this subroutine is to interpolate in C and O
  !     abundances.     

  real :: xxc, xxo

  save 
  integer :: w 

  integer :: is, ir, it, iw, i, m1, m2, ie, iei, i1, i2, i3, iej
  integer :: j1, j2, j3, jx, ix
  real :: oxdp, fac, quad, cxdp

  is = 0 
  if(xxco < 1.d-6) then 
     do ir=l1,l1+iq 
        do it=k1,k1+ip 
           opl(m,it,ir) = covett(m,1,1,it,ir,indice) 
        enddo
     enddo
     is = 1 
     return 
  endif
  !     include boundaries that could later cause division by 0!          
  if(xxc > xcd(3)-1.d-6) then 
     oxdp = log10(zzz+xodp) 
     !     handle possibility that xodp=0                                    
     fac = max(min((oxx-ox(1))/max(oxdp-ox(1),1.d-6),1.),0.) 
     do it=k1,k1+ip    ! giada: loop interchange
        do ir=l1,l1+iq 
           ! interpolation in region c1                         
           ! include boundaries that could later cause division by 0! 
           if(xxc > xcd(2)-1.d-6) then 
              iw = 1  
              a(1,m) = quad(is,iw,cxx,covett(m,nc-2,1,it,ir,indice),  &
                   covett(m,nc-1,1,it,ir,indice),diagvett(m,1,it,ir, &
                   indice), cx(nc-2),cx(nc-1),cx(nc))                   
              iw = iw+1 
              a(2,m) = quad(is,iw,cxx,diagvett(m,1,it,ir,indice), &
                   diagvett(m,2,it, ir,indice),diagvett(m,3,it,ir, &
                   indice),cxd(1),cxd(2),cxd(3))     
              do w=1,2 
                 b(w) = a(w,m) 
              enddo
              !     handle possibility that xodp=0  
              opl(m,it,ir) = b(1)+(b(2)-b(1))*fac 
              is = 1 
              cycle
           endif
           ! interpolation in region c2  
           iw = 1 
           a(1,m) = quad(is,iw,cxx,covett(m,nc-2,1,it,ir,indice), &
                covett(m,nc-1,1,it,ir,indice),diagvett(m,1,it,ir, &
                indice),cx(nc-2),cx(nc-1),cx(nc)) 
           iw = iw+1 
           a(2,m) = quad(is,iw,cxx,covett(m,n(m,2)-2,2,it,ir,indice), &
                covett(m, n(m,2)-1,2,it,ir,indice),diagvett(m,2,it,ir, &
                indice),cx(n(m,2)-2),cx(n(m,2)-1),cxd(2))
           iw = iw+1 
           a(3,m) = quad(is,iw,cxx,diagvett(m,1,it,ir,indice), &
                diagvett(m,2,it, ir,indice),diagvett(m,3,it,ir,indice), &
                cxd(1),cxd(2),cxd(3))   
           do w=1,3 
              b(w) = a(w,m) 
           enddo
           iw = iw+1 
           opl(m,it,ir) = quad(is,iw,oxx,b(1),b(2),b(3),ox(1),ox(2),oxdp) 
           is = 1 
        end do
     end do
     if(is  ==  1) return 
  endif
  ! interpolation in region c3 to c6                   
  is = 0 
  if(nc >= 5) then 
     do i=4,nc-1 
        ! do not go beyond middle (where c3-c6 overlaps o3-o6), and         
        if((xxc > xcd(i)-1.d-6) .and. (xxo > xo(i-1)-1.d-6) .and. &
             (xcd(i-1) > xc(i-1))) then     
           oxdp = log10(zzz+xodp)  ! giada: estratto dal loop
           do it=k1,k1+ip          
              do ir=l1,l1+iq 
                 iw = 1 
                 m1 = i-1 
                 m2 = i-2 
                 a(1,m) = quad(is,iw,cxx,covett(m,n(m,m2)-2,m2,it,ir,indice), &
                      covett(m,n(m,m2)-1,m2,it,ir,indice),diagvett(m,m2,it, &
                      ir,indice), cx(n(m,m2)-2),cx(n(m,m2)-1),cxd(m2)) 
                 iw = iw+1 
                 a(2,m) = quad(is,iw,cxx,covett(m,n(m,m1)-2,m1,it,ir,indice), &
                      covett(m,n(m,m1)-1,m1,it,ir,indice),diagvett(m,m1,it, &
                      ir,indice), cx(n(m,m1)-2),cx(n(m,m1)-1),cxd(m1))
                 iw = iw+1 
                 a(3,m) = quad(is,iw,cxx,diagvett(m,m2,it,ir,indice), &
                      diagvett(m,m1,it,ir,indice),diagvett(m,i,it,ir, &
                      indice),cxd(m2),cxd(m1),cxd(i)) 
                 do w=1,3 
                    b(w) = a(w,m) 
                 enddo
                 iw = iw+1 
                 opl(m,it,ir) = quad(is,iw,oxx,b(1),b(2),b(3),ox(i-2), &
                      ox(i-1),oxdp) 
                 is = 1 
              end do
           end do
           if (is  ==  1) return 
        endif
     end do
  endif

  if (is  ==  1) return 
  !     include boundaries that could later cause division by 0!          
  if(xxo > xod(3)-1.d-6) then 
     cxdp = log10(zzz+xcdp) 
     !     handle possibility that xcdp=0                                    
     fac = max(min((cxx-cx(1))/max(cxdp-cx(1),1.d-6),1.),0.) 
     do ir=l1,l1+iq 
        do it=k1,k1+ip 
           ! interpolation in region  o1                        
           ! include boundaries that could later cause division by 0!
           if(xxo > xod(2)-1.d-6) then 
              iw = 1 
              a(1,m) = quad(is,iw,oxx,covett(m,1,no-2,it,ir,indice), &
                   covett(m,1,no-1,it,ir,indice),diagovett(m,no-1,it,ir, &
                   indice),ox(no-2),ox(no-1),ox(no)) 
              iw = iw+1 
              a(2,m) = quad(is,iw,oxx,diagovett(m,no-1,it,ir,indice), &
                   diagovett(m,no-2,it,ir,indice),diagovett(m,no-3,it,ir, &
                   indice),oxd(1),oxd(2),oxd(3))
              do w=1,2 
                 b(w) = a(w,m) 
              enddo
              !     handle possibility that xcdp=0     
              opl(m,it,ir) = b(1)+(b(2)-b(1))*fac 
              is = 1 
              cycle
           endif
           ! interpolation in region  o2  
           iw = 1 
           a(1,m) = quad(is,iw,oxx,covett(m,1,no-2,it,ir,indice),covett(m,1, &
                no-1,it,ir,indice),diagovett(m,no-1,it,ir,indice),ox(no-2),  &
                ox(no-1),ox(no))
           iw = iw+1 
           a(2,m) = quad(is,iw,oxx,covett(m,2,n(m,2)-2,it,ir,indice), &
                covett(m,2,n(m,2)-1,it,ir,indice),diagovett(m,no-2,it,ir, &
                indice),ox(n(m,2)-2),ox(n(m,2)-1),oxd(2))
           iw = iw+1 
           a(3,m) = quad(is,iw,oxx,diagovett(m,no-1,it,ir,indice), &
                diagovett(m,no-2,it,ir,indice),diagovett(m,nc-3,it,ir, &
                indice),oxd(1),oxd(2),oxd(3)) 
           do w=1,3 
              b(w) = a(w,m) 
           enddo
           iw = iw+1 
           opl(m,it,ir) = quad(is,iw,cxx,b(1),b(2),b(3),cx(1),cx(2),cxdp) 
           is = 1 
        end do
     end do
     if(is  ==  1) return 
  endif
  ! interpolation in region  o3 to o6                  
  is = 0 
  if(no >= 5) then 
     do i=4,no-1 
        ! do not go beyond middle (where o3-o6 overlaps c3-c6), and         
        if((xxo > xod(i)-1.d-6) .and. (xxc > xc(i-1)-1.d-6) .and.   &
             (xod(i-1) > xo(i-1)-1.d-6)) then 
           cxdp = log10(zzz+xcdp) 
           do ir=l1,l1+iq 
              do it=k1,k1+ip 
                 iw = 1 
                 m2 = i-2 
                 m1 = i-1 
                 a(1,m) = quad(is,iw,oxx,covett(m,m2,n(m,m2)-2,it,ir, &
                      indice),covett(m,m2,n(m,m2)-1,it,ir,indice), &
                      diagovett(m,no-m2,it,ir,indice),ox(n(m,m2)-2), &
                      ox(n(m,m2)-1),oxd(m2))
                 iw = iw+1 
                 a(2,m) = quad(is,iw,oxx,covett(m,m1,n(m,m1)-2,it,ir, &
                      indice),covett(m,m1,n(m,m1)-1,it,ir,indice), &
                      diagovett(m,no-m1,it,ir,indice),ox(n(m,m1)-2), &
                      ox(n(m,m1)-1),oxd(m1))
                 iw = iw+1 
                 a(3,m) = quad(is,iw,oxx,diagovett(m,no-m2,it,ir,indice), &
                      diagovett(m,no-m1,it,ir,indice),diagovett(m,no-i,it, &
                      ir,indice),oxd(m2),oxd(m1),oxd(i))
                 do w=1,3 
                    b(w) = a(w,m) 
                 enddo
                 iw = iw+1 
                 opl(m,it,ir) = quad(is,iw,cxx,b(1),b(2),b(3),cx(m2), &
                      cx(m1),cxdp) 
                 is = 1 
              end do
           end do
           if (is  ==  1) return 
        endif
     end do
  endif

  if (is  ==  1) return 

  !.....find index1 of C grid.                                             
  ie = 100*xxc+1 
  iei = index1(ie)+1 
  !     must also allow index1 = nc, to avoid extrapolation                
  if (iei > nc) iei = nc 

  if(iei > 3) then 
     i1 = iei-2 
     i2 = iei-1 
     i3 = iei 
  else 
     i1 = 1 
     i2 = 2 
     i3 = 3 
  endif
  !.....find index1 of O grid                                              
  ie = 100*xxo+1 
  iej = index1(ie)+1 
  !     must also allow index1 = no, to avoid extrapolation                
  if (iej > no) iej = no 

  if(iej > 3) then 
     j1 = iej-2 
     j2 = iej-1 
     j3 = iej 
  else 
     j1 = 1 
     j2 = 2 
     j3 = 3 
  endif
  !     lower-O part of grid: interpolate C before O                      
  if(j3 < no .and. i3 <= n(m,j3) .and.                             &
       (xxc < xcd(j3)+1.d-6 .or. xxc >= xxo)) then  
     do ir=l1,l1+iq 
        do it=k1,k1+ip 
           iw = 0 
           do jx=j1,j1+2 
              iw = iw+1 
              !     if i3=n(m,jx), then must replace cx(i3) with cxd(jx)
              a(iw,m) = quad(is,iw,cxx,covett(m,i1,jx,it,ir,indice),covett(m, &
                   i2,jx,it,ir,indice),covett(m,i3,jx,it,ir,indice),cx(i1), &
                   cx(i2),min(cx(i3),cxd(jx))) 
           enddo
           do w=1,3 
              b(w) = a(w,m) 
           enddo
           iw = iw+1 
           opl(m,it,ir) = quad(is,iw,oxx,b(1),b(2),b(3),ox(j1),ox(j2),ox(j3)) 
           is = 1
        end do
     end do
     !     else: high-O part of grid: must interpolate O before C            
  else 
     do ir=l1,l1+iq 
        do it=k1,k1+ip 
           iw = 0 
           do ix=i1,i1+2 
              iw = iw+1 
              if(j3 < n(m,ix)) then 
                 a(iw,m) = quad(is,iw,oxx,covett(m,ix,j1,it,ir,indice), &
                      covett(m,ix,j2,it,ir,indice),covett(m,ix,j3,it,ir, &
                      indice),ox(j1),ox(j2), ox(j3))
              else 
                 a(iw,m) = quad(is,iw,oxx,covett(m,ix,j1,it,ir,indice), &
                      covett(m,ix,j2,it,ir,indice),diagovett(m,no-ix,it,ir, &
                      indice),ox(j1),ox(j2),oxd(ix)) 
              endif
           enddo
           do w=1,3 
              b(w) = a(w,m) 
           enddo
           iw = iw+1 
           opl(m,it,ir) = quad(is,iw,cxx,b(1),b(2),b(3),cx(i1),cx(i2),cx(i3)) 
           is = 1 
        enddo
     enddo
  endif
  return
end subroutine cointerp

!***********************************************************************
subroutine t6rinterp(slr,slt) 
  use interfaccia
  use aaa
  use aa
  use a1
  use d
  use bb
  use e

  implicit none
  !     The purpose of this subroutine is to interpolate in logT6 and logR
  
  real :: slr, slt

  save 
  integer :: is, iu, kx, iw, lx, idx
  real :: quad, dkap1, dkap2, opactq, opactq2, dkapq1, dkapq2, opact2, dix
  real :: dopactq, opacr, opacrq, dopacrq, dix2,  dopactr

  real,dimension(4) :: x
  real,dimension(4*(ip+1)) :: y

  is = 0 
  iu = 0 

  do kx=k1,k1+ip 
     iw = 1 
     iu = iu+1 
     h(iu) = quad(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3),          &
          alr(l1),alr(l2),alr(l3))         
     if(iq ==  3) then 
        iw = 2 
        q(iu) = quad(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4),      &
             alr(l2),alr(l3),alr(l4))   
     endif
     is = 1 
  enddo

  is = 0 
  iw = 1 
  !..... k and Log(k)/log(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)   
  opact = quad(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3)) 
  dopact = dkap 
  dkap1 = dkap 
  if (iq == 3) then 
     !.....    k and Log(k)/Log(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2) 
     opactq = quad(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3)) 
     dkapq1 = dkap 
  endif
  if(ip == 3) then 
     !.....    k and Log(k)/Log(T6) in lower-left 3x3.                     
     opact2 = quad(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4)) 
     dkap2 = dkap 
     !.....    k and Log(k)/Log(T6) smoothed in left 3x4                   
     dix = (alt(k3)-slt)*dfs(k3) 
     dopact = dkap1*dix+dkap2*(1.-dix) 
     opact = opact*dix+opact2*(1.-dix) 
     if(iq == 3) then 
        !.....    k and Log(k)/Log(T6) in upper-right 3x3.
        opactq2 = quad(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4)) 
        dkapq2 = dkap 
        dopactq = dkapq1*dix+dkapq2*(1.-dix) 
        opactq = opactq*dix+opactq2*(1.-dix) 
     endif
  endif

  iu = 0 
  do lx=l1,l1+iq 
     iw = 1 
     iu = iu+1 
     h(iu) = quad(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx),          &
          alt(k1),alt(k2),alt(k3))             
     if(ip == 3) then 
        iw = 2 
        q(iu) = quad(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx),      &
             alt(k2),alt(k3),alt(k4))                
     endif
     is = 1 
  end do
  
  is = 0 
  iw = 1 
  !..... k and Log(k)/Log(R) in lower-left 3x3                          
  opacr = quad(is,iw,slr,h(1),h(2),h(3),alr(l1),alr(l2),alr(l3)) 
  dopacr = dkap 
  if(ip == 3) then 
     opacrq = quad(is,iw,slr,q(1),q(2),q(3),alr(l1),alr(l2),alr(l3)) 
     !.....    k and Log(k)/Log(R) in upper-left 3x3.                      
     dopacrq = dkap 
  endif
  if(iq == 3) then 
     !.....    k and Log(k)/Log(R) in lower-right 3x3.                     
     opact2 = quad(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4)) 
     dix2 = (alr(l3)-slr)*dfsr(l3) 
     dopacr = dopacr*dix2+dkap*(1.-dix2) 
     !.....        k and Log(k)/Log(T6) smoothed in both log(T6) and log(R)
     dopact = dopact*dix2+dopactq*(1.-dix2) 
     opact = opact*dix2+opactq*(1.-dix2) 
  endif
  if(ip == 3) then 
     if(iq == 3) then 
        !.....    k and Log(k)/Log(R) in upper-right 3x3.   
        opacrq = quad(is,iw,slr,q(2),q(3),q(4),alr(l2),alr(l3),alr(l4)) 
        !.....        Log(k)/Log(R) smoothed in both log(T6) and Log(R). 
        dopacrq = dopacrq*dix2+dkap*(1.-dix2) 
     endif
     dopacr = dopacr*dix+dopacrq*(1.-dix) 
  endif
  dopactd = dopact-3.*dopacr 
  if (opact > 1.d15) then 
     write(*,'("Interpolation indices out of range; please report &
          &conditions.")')   
     write(66,*)'510 - t6rinterp'
     write(66,'("Interpolation indices out of range",&
          &";please report conditions.")') 
     stop 
  endif

  if (opact > 9) then 
     opact = 30. 
     dopact = 99. 
     dopactr = 99. 
     dopactd = 99. 
  endif
  return 
end subroutine t6rinterp

!***********************************************************************
subroutine readco 
  use interfaccia
  use aa
  use a1
  use alink
  use cst
  use b1
  use e
  use indicefile
  use recoin
  
  implicit none
  !..... The purpose of this subroutine is to read the data tables modified 
  save 

  integer,parameter :: ismdata=0
  real,dimension(mx,mc,mo,nt,nr):: co
  real,dimension(mx,mc,nt,nr) :: diag
  real,dimension(mx,mo,nt,nr) :: diago
  real,dimension(mx,ntabs,10) :: zzvett

  integer :: i, j, k, l, mq, int, kk, ll, isett6, istep
  real :: dum, altin

  if (ireadco /= 12345678) then 
     write(*,*)'Entro in readco' 

     if (itimeco /= 12345678) then 
        do i=1,mx 
           do j=1,mc 
              do k=1,mo 
                 do l=1,nt 
                    do mq=1,nr 
                       co(i,j,k,l,mq) = 1.d35 
                    end do
                 end do
              end do
           end do
        end do
        do i=1,mx 
           do j=1,mc 
              do l=1,nt 
                 do mq=1,nr 
                    diag(i,j,l,mq) = 1.d35 
                    diago(i,j,l,mq) = 1.d35 
                 end do
              end do
           end do
        end do
        itimeco = 12345678 
     endif
     do j=1,nc-1 
        do i=1,nc 
           if(xcd(j) >= xc(i)) then 
              n(m,j) = i+1 
              if(xcd(j) < xc(i)+1.d-6) n(m,j) = i 
           endif
        end do
     end do
     n(m,nc) = 0 

     !..... read X=0.0 tables                                                
     if(m  ==  1) then 
        ! Nel mio caso devo leggere uno alla volta i 10 files   
        ! con x=0, ciascuno con 60 tabelle   
        if(indice == 1)  open(7,FILE='X00.Z0.0000',status='old') 
        if(indice == 2)  open(7,FILE='X00.Z0.0001',status='old') 
        if(indice == 3)  open(7,FILE='X00.Z0.0003',status='old') 
        if(indice == 4)  open(7,FILE='X00.Z0.0010',status='old') 
        if(indice == 5)  open(7,FILE='X00.Z0.0020',status='old') 
        if(indice == 6)  open(7,FILE='X00.Z0.0040',status='old') 
        if(indice == 7)  open(7,FILE='X00.Z0.0100',status='old') 
        if(indice == 8)  open(7,FILE='X00.Z0.0200',status='old') 
        if(indice == 9)  open(7,FILE='X00.Z0.0300',status='old') 
        if(indice == 10) open(7,FILE='X00.Z0.0400',status='old') 
     endif
     
     !..... read X=0.03 tables                                               
     if(m  ==  2) open(7, FILE='codatab') 
     !..... read X=0.10 tables                                               
     if(m  ==  3) open(7, FILE='codatac') 
     !..... read X=0.35 tables                                               
     if(m  ==  4) open(7, FILE='codatad') 
     !.....read X=0.70 tables                                                
     if(m  ==  5) open(7, FILE='codatae') 

     rewind(7) 

     int = 0 
     do j=1,no-1 
        do i=1,n(m,j) 
           int = int+1 
           read(7,'(a)') dum 
           read (7,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)')     &
                itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
           xca(m,int) = min(xca(m,int),1.-x(m,int)-zz(m,int)-xoa(m,int)) 

           read(7,'(a)') dum,dum,dum 
           read(7,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm) 
           read(7,'(a)') dum 
           do k=1,ntm 
              read(7,'(f4.2,19f7.3)') altin,(cof(k,l), l=1,nrm) 
              ! modified                                    
              do ll=1,nrm 
                 coff(k,ll) = cof(k,ll) 
              end do
           end do
           if (isett6 /= 1234567) then 
              do k=1,ntm 
                 t6arr(k) = t6list(k) 
              end do
           endif
           isett6 = 1234567 
           zzvett(m,int,indice) = zz(m,int) 

           if (ismdata  ==  0) then 
              if (nrm /= nr .or. ntm /= nt) then 
                 write (*,'("Not set up to smooth data with reduced ",&
                      &"T-Rho grid")') 
                 write(66,*)'520 - readco'
                 write (66,'("Not set up to smooth data with reduced ",&
                      &"T-Rho grid")') 
                 stop 
              endif
              ! modified                                           
              tmax = 10. 
              ! 65 in earlier version                              
              nset = 67 
              RLS = -8. 
              nsm = 1 
              RLE = 1.0 
              nrlow = 1 
              nrhigh = 2*(RLE-RLS)+1 
              !modified                                       
              call opaltab 
           endif

           ll=1 
           do kk=nrb,nre 
              do k=1,nt 
                 if (ismdata  ==  0) then 
                    if (m  == 1 .and. k <= 9) then 
                       co(m,i,j,k,ll) = cof(k+ntb-1,kk) 
                    else 
                       co(m,i,j,k,ll) = coff(k+ntb-1,kk) 
                    endif
                 else 
                    co(m,i,j,k,ll) = coff(k+ntb-1,kk) 
                 endif
              end do
              ll = ll+1 
           end do
        end do
     end do

     if(x(m,1) /= xa(m)) then 
        write(*,'(" X in the codata? file does not match xa(m)")') 
        write(66,*)'521 - readco'
        write(66,'(" X in the codata? file does not match xa(m)")')
        stop 
     endif
     do i=1, nc-1 
        do k=1,nt 
           do l=1,nr 
              diag(m,i,k,l) = co(m,n(m,i),i,k,l) 
           end do
        end do
     end do

     do j=1,no-1 
        int = int+1 
        read(7,'(a)') dum 
        read (7,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)')     &
             itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int), &
             xoa(m,int)  
        read(7,'(a)') dum,dum,dum 
        read(7,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm) 
        read(7,'(a)') dum 
        do k=1,ntm 
           read (7,'(f4.2,19f7.3)') dum, (cof(k,l), l=1,nrm) 
           !     set up to smooth final "diago" opacity tables
           do l=1,nrm 
              coff(k,l)=cof(k,l) 
           end do
        enddo
        !     smooth final "diago" opacity tables too!  
        if (ismdata  ==  0) then 
           ! modified                                           
           tmax = 10. 
           !65 in earlier version                                  
           nset = 67 
           RLS = -8. 
           nsm = 1 
           RLE = 1.0 
           nrlow = 1 
           nrhigh = 2*(RLE-RLS)+1 
           !modified                                       
           call opaltab 
           do k=3,NTEMP-2 
              do ll=nrlow,nrhigh 
                 ! Following skip required because, due to missing data,
                 ! the X=0  low T data cannot be smoothed
                 if (m == 1 .and. k <= 9) then 
                    cof(k,ll) = cof(k,ll) 
                 else 
                    cof(k,ll) = coff(k,ll) 
                 endif
              end do
           end do
           ll = 1 
           do kk=nrb,nre 
              do k=1,nt 
                 diago(m,j,k,ll) = cof(k+ntb-1,kk) 
              end do
              ll = ll+1 
           end do
        endif
     end do

     do i=2,nt 
        dfs(i) = 1./(alt(i)-alt(i-1)) 
     enddo
     do i=2,nr 
        dfsr(i) = 1./(alr(i)-alr(i-1)) 
     enddo
     istep = -1 
     if (mx > 1 ) then 
        istep = 1 
        do i=2,mx,istep 
           dfsx(i) = 1./(xx(i)-xx(i-1)) 
        end do
     endif
     close(7) 

     do i=1,mx 
        do j=1,mc 
           do k=1,mo 
              do l=1,nt 
                 do mq=1,nr 
                    covett(i,j,k,l,mq,indice) = co(i,j,k,l,mq) 
                 enddo
              enddo
           enddo
        enddo
     enddo
     do i=1,mx 
        do j=1,mc 
           do l=1,nt 
              do mq=1,nr 
                 diagvett(i,j,l,mq,indice) = diag(i,j,l,mq) 
                 diagovett(i,j,l,mq,indice) = diago(i,j,l,mq) 
              end do
           end do
        end do
     end do
  else 
     do int=1,ntabs 
        zz(m,int) = zzvett(m,int,indice) 
     end do
  end if

  return 
end subroutine readco

!***********************************************************************
real function quad(ic,i,x,y1,y2,y3,x1,x2,x3) 
  use d
  implicit none
  !..... this function performs a quadratic interpolation.                
  
  integer :: ic, i
  real :: x, y1, y2, y3, x1, x2, x3
  save 
  
  real :: xx(3), yy(3), xx12(30), xx13(30), xx23(30), xx1sq(30), xx1pxx2(30)
  real :: c1, c2, c3, c3a
  
  xx(1) = x1 
  xx(2) = x2 
  xx(3) = x3 
  yy(1) = y1 
  yy(2) = y2 
  yy(3) = y3 
  if(ic  ==  0) then 
     xx12(i) = 1./(xx(1)-xx(2)) 
     xx13(i) = 1./(xx(1)-xx(3)) 
     xx23(i) = 1./(xx(2)-xx(3)) 
     xx1sq(i) = xx(1)*xx(1) 
     xx1pxx2(i) = xx(1)+xx(2) 
  endif
  
  c3a = (yy(1)-yy(2))*xx12(i) 
  c3 = (c3a-(yy(2)-yy(3))*xx23(i))*xx13(i) 
  c2 = c3a-xx1pxx2(i)*c3 
  c1 = yy(1)-xx(1)*c2-xx1sq(i)*c3 
  dkap = c2+(x+x)*c3 
  quad = c1+x*(c2+x*c3) 
  
end function quad


subroutine SPLINE(X,Y,N,Y2) 
  implicit none
  integer,parameter :: NMAX=100

  integer :: N
  real,dimension(N) :: X,Y,Y2

  real,dimension(NMAX) :: U 
  integer :: i, k
  real :: yp1, ypn, sig, p, qn, un
  

  !     FIRST DERIVATIVES AT END POINTS USING CUBIC FIT                   
  YP1 = ((Y(3)-Y(1))*(X(2)-X(1))**2                                &
       -(Y(2)-Y(1))*(X(3)-X(1))**2)/                               &
       ((X(3)-X(1))*(X(2)-X(1))*(X(2)-X(3)))                          
  YPN = ((Y(N-2)-Y(N))*(X(N-1)-X(N))**2                            &
       -(Y(N-1)-Y(N))*(X(N-2)-X(N))**2)/                           &
       ((X(N-2)-X(N))*(X(N-1)-X(N))*(X(N-1)-X(N-2)))              

  Y2(1) = -0.5 
  U(1) = (3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1) 
  do I=2,N-1 
     SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1)) 
     P = SIG*Y2(I-1)+2. 
     Y2(I) = (SIG-1.)/P 
     U(I) = (6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))             &
          /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P               
  end do
  QN = 0.5 
  UN = (3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1))) 
  Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.) 
  do K=N-1,1,-1 
     Y2(K) = Y2(K)*Y2(K+1)+U(K) 
  end do
  return 
end subroutine SPLINE

!*********************************************************************  
subroutine SPLINT(XA,YA,N,Y2A,X,Y,YP) 
 implicit none
 
 integer :: N
 real,dimension(N) :: XA, YA, Y2A 
 real :: X, Y, YP

 integer :: klo, khi, k
 real :: h, a, b

 KLO = 1 
 KHI = N 
 do while(KHI-KLO > 1)
    K = (KHI+KLO)/2 
    if(XA(K) > X)then 
       KHI = K 
    else 
       KLO = K 
    endif
 end do

 H = XA(KHI)-XA(KLO) 
 if (H == 0.) then 
    write(*,*) 'Bad XA input.' 
    write(66,*)'530 - SPLINT'
    write(66,*) 'Bad XA input.' 
    stop 
 endif

 A = (XA(KHI)-X)/H 
 B = (X-XA(KLO))/H 
 Y = A*YA(KLO)+B*YA(KHI)+                                            &
      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.             
 YP = 0.05*((-YA(KLO)+YA(KHI))/H + (-(3*A**2-1)*Y2A(KLO)             &
      +(3*B**2-1)*Y2A(KHI) )*H/6. )                         
 return 
end subroutine SPLINT

!*********************************************************************  
subroutine FITY 
  use cst
  use cf

  implicit none
  !  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS             
  !  FY AND FXY                                                           

  real,dimension(IPR) :: A1, B1, AD, BD 

  integer :: i, j
  real :: ap1, apn, bp1, bpn

  do I=1,nset 
     do J=1,NRL 
        A1(J) = F(I,J) 
        B1(J) = FX(I,J) 
     end do

     call GETD(A1,NRL,AD,AP1,APN) 
     call GETD(B1,NRL,BD,BP1,BPN) 

     FY(I,1) = AP1 
     FY(I,NRL) = APN 
     FXY(I,1) = BP1 
     FXY(I,NRL) = BPN 
     do J=2,NRL-1 
        FY(I,J) = -A1(J)+A1(J+1)-2.*AD(J)-AD(J+1) 
        FXY(I,J) = -B1(J)+B1(J+1)-2.*BD(J)-BD(J+1) 
     end do
  end do
  return 
end subroutine FITY

!*********************************************************************  
subroutine FITX 
  use cst
  use cf

  implicit none
  !  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.                           
  !  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.                    
  
  real,dimension(85) ::  A, D

  integer :: i, j
  real ::  ap1, apn

  do J=1,NRL 
     ! modified                                      
     do I=1,nset 
        A(I) = F(I,J) 
     end do
     ! modified                        
     call GETD(A,nset,D,AP1,APN) 
     FX(1,J) = AP1 
     ! modified                                    
     FX(nset,J) = APN 
     ! modified                                   
     do I=2,nset-1 
        FX(I,J) = -A(I)+A(I+1)-2.*D(I)-D(I+1) 
     end do
  end do
  return 
end subroutine FITX

!***********************************************************************
subroutine GETD(F,N,D,FP1,FPN) 
  implicit none
  !  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS       
  !  OF UNITY.                                                            

  integer :: N
  real,dimension(N) ::  F, D
  real :: FP1, FPN

  real,dimension(85) :: T 
  integer :: j

  FP1 = (-11.*F(1)+18.*F(2)-9.*F(3)+2.*F(4))/6. 
  FPN = (11.*F(N)-18.*F(N-1)+9.*F(N-2)-2.*F(N-3))/6. 

  D(1) = -.5 
  T(1) = .5*(-F(1)+F(2)-FP1) 

  do J=2,N-1 
     D(J) = -1./(4.+D(J-1)) 
     T(J) = -D(J)*(F(J-1)-2.*F(J)+F(J+1)-T(J-1)) 
  end do

  D(N) = (FPN+F(N-1)-F(N)-T(N-1))/(2.+D(N-1)) 

  do J=N-1,1,-1 
     D(J) = D(J)*D(J+1)+T(J) 
  end do
  return 
end subroutine GETD

!*********************************************************************  
subroutine INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR) 
  use cst
  use cf

  implicit none

  real :: FLT, FLRHO, G, DGDT, DGDRHO
  logical :: IERR

  !  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE               
  !  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN              
  !  Numerical Recipes, PP. 118 TO 120                                    

  real,dimension(16) :: B

  real :: ff, ffx, ffy
  integer :: i, j
  real :: x, flr, y, u, v, s, t

  !  EXTREME LIMITS ALLOWED ARE:-                                         
  !     (3.800-0.0125) TO (8.000+0.0125) FOR LOG10(T)                     
  !     (RLS-0.125) TO (RLE+0.1254) FOR LOG10(R)                          
  !     (ALLOWING FOR SMALL EXTRAPOLATIONS BEYOND TABULAR VALUES)         
  !                                                                       
  !  FUNCTION DEFINITIONS FOR CUBIC EXPANSION                             
  !                                                                       
  FF(S,T) = B( 1)+T*(B( 2)+T*(B( 3)+T*B( 4)))                    &
       +S*(B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))                     &
       +S*(B( 9)+T*(B(10)+T*(B(11)+T*B(12)))                     &
       +S*(B(13)+T*(B(14)+T*(B(15)+T*B(16))) )))                

  FFX(S,T) =  B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))                   &
       +S*(2*(B( 9)+T*(B(10)+T*(B(11)+T*B(12))))                  &
       +S*(3*(B(13)+T*(B(14)+T*(B(15)+T*B(16)))) ))                 

  FFY(S,T) =  B( 2)+S*(B( 6)+S*(B(10)+S*B(14)))                   &
       +T*(2*(B( 3)+S*(B( 7)+S*(B(11)+S*B(15))))                  &
       +T*(3*(B( 4)+S*(B( 8)+S*(B(12)+S*B(16)))) ))                 

!!$  FFXY(S,T) = B( 6)+T*(2*B( 7)+3*T*B( 8))                          &
!!$       +S*(2*B(10)+T*(4*B(11)+6*T*B(12))                           &
!!$       +S*(3*B(14)+T*(6*B(15)+9*T*B(16)) ))                        

  IERR = .false. 

  X = 20.*(FLT-3.800)+1 
  FLR = FLRHO+18.-3.*FLT 
  Y = 2*( FLR - RLS )+1 

  if(X < 2.) then 
     if(X < 0.75) then 
        IERR = .true. 
     else 
        I = 1 
     endif
  else if(X > 84) then 
     if(X > 85.25) then 
        IERR = .true. 
     else 
        I = 84 
     endif
  else 
     I = X 
  endif
  U = X-I 

  if(Y < 2.) then 
     if(Y < 0.75) then 
        IERR = .true. 
     else 
        J = 1 
     endif
  elseif(Y > NRL-1) then 
     if(Y > NRL+.25) then 
        IERR = .true. 
     else 
        J = NRL-1 
     endif
  else 
     J = Y 
  endif
  V = Y-J 

  if(IERR) then 
     G = 9.999 
     DGDT = 9.999 
     DGDRHO = 9.999 
     return 
  endif

  !  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
  B(1) = F(I,J) 
  B(2) = FY(I,J) 
  B(3) = 3*(-F(I,J)+F(I,J+1))-2*FY(I,J)-FY(I,J+1) 
  B(4) = 2*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1) 

  B(5) = FX(I,J) 
  B(6) = FXY(I,J) 
  B(7) = 3*(-FX(I,J)+FX(I,J+1))-2*FXY(I,J)-FXY(I,J+1) 
  B(8) = 2*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1) 

  B(9) = 3*(-F(I,J)+F(I+1,J))-2*FX(I,J)-FX(I+1,J) 
  B(10) = 3*(-FY(I,J)+FY(I+1,J))-2*FXY(I,J)-FXY(I+1,J) 
  B(11) = 9*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))                     &
       +6*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))+4*FXY(I,J)            &
       +3*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))               &
       +2*(FXY(I,J+1)+FXY(I+1,J))+FXY(I+1,J+1)                          
  B(12) = 6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))                    &
       +4*(-FX(I,J)+FX(I,J+1))                                        &
       +3*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))                  &
       +2*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))                &
       -FXY(I+1,J)-FXY(I+1,J+1)                                         

  B(13) = 2*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J) 
  B(14) = 2*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J) 
  B(15) = 6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))                    &
       +4*(-FY(I,J)+FY(I+1,J))                                        &
       +3*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))                  &
       +2*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))                 &
       -FXY(I+1,J+1)-FXY(I,J+1)                                         
  B(16) = 4*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))                     &
       +2*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)                    &
       +FY(I,J)-FY(I+1,J)-FY(I+1,J+1)+FY(I,J+1))                      &
       +FXY(I,J)+FXY(I+1,J)+FXY(I+1,J+1)+FXY(I,J+1)                    
  
  !  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),                    
  !      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)                                
  G = FF(U,V) 
  DGDT = 20.*FFX(U,V)-6.*FFY(U,V) 
  DGDRHO = 2.*FFY(U,V) 
  return 
end subroutine INTERP

!*********************************************************************  
subroutine SMOOTHING 
  use cst
  use cf
  
  implicit none
  !  THIS SUBROUTINE USES A 2-DIMENSIONAL GENERALISATION OF THE SMOOTHING 
  !  TECHNIQUES DESCRIBED ON PP. 644 TO 649 OF Numerical Recipes.         
  !                                                                       
  !  CONSIDER THE 25 POINTS DEFINED BY                                    
  !       I+n, n=-2,-1,0,1,2 AND J+m, m=-2,-1,0,1,2.                      
  !  THE FUNCTION TO BE SMOOTHED IS FITTED TO A BI-CUBIC, INVOLVING       
  !  16 COEFFICIENTS, USING TECHNIQUES OF LEAST-SQUARES. THE SMOOTHED     
  !  FUNCTION (TEMPORARILY STORED IN FXY) IS GIVEN BY THE FITTED VALUE    
  !  AT THE POINT I AND J.                                                

  !  THE FITTING IS SHIFTED FOR POINTS CLOSE TO BOUNDARIES.               

  real,parameter,dimension(6) :: GAM = (/0.0073469388,-0.0293877551, &
       -0.0416326531,0.1175510204,0.1665306122,0.2359183673/)
  real,parameter,dimension(11) :: BET =(/-0.0048979592,-0.0661224490, &
       -0.0293877551,0.0195918367,0.2644897959,0.1175510204,-0.0783673469, &
       0.0277551020,0.3746938776,0.1665306122,-0.1110204082/)
  real,parameter,dimension(11) :: ALP =(/-0.0844897959,-0.0048979592, &
       0.0073469388,0.0012244898,0.3379591837,0.0195918367,-0.0293877551, &
       0.4787755102,0.0277551020,-0.0416326531,-0.0069387755/)

  integer :: i, j
 
  do I=3,nset-2
     J = 1 
     FXY(I,J) =                                                          &
          ALP(1)*( F(I-2,J  )+F(I+2,J  ) )                               &
          +ALP(2)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J+3)+F(I+2,J+3)          &
          +F(I-1,J+4)+F(I+1,J+4) )                                       &
          +ALP(3)*( F(I-2,J+2)+F(I+2,J+2) )                              &
          +ALP(4)*( F(I-2,J+4)+F(I+2,J+4) )                              &
          +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )                              &
          +ALP(6)*( F(I-1,J+1)+F(I+1,J+1)+F(I-1,J+3)+F(I+1,J+3) )        &
          +ALP(7)*( F(I-1,J+2)+F(I+1,J+2) )                              &
          +ALP(8)*  F(I  ,J  )                                           &
          +ALP(9)*( F(I  ,J+1)+F(I  ,J+3) )                              &
          +ALP(10)* F(I  ,J+2) +ALP(11)*F(I  ,J+4)                       

     J = 2 
     FXY(I,J)=                                                           &
          BET(1)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J+3)+F(I+2,J+3) )         &
          +BET(2)*( F(I-2,J  )+F(I+2,J  ) )                              &
          +BET(3)*( F(I-2,J+1)+F(I+2,J+1) )                              &
          +BET(4)*( F(I-2,J+2)+F(I+2,J+2)+F(I-1,J-1)+F(I+1,J-1)          &
          +F(I-1,J+3)+F(I+1,J+3) )                                       &
          +BET(5)*( F(I-1,J  )+F(I+1,J  ) )                              &
          +BET(6)*( F(I-1,J+1)+F(I+1,J+1) )                              &
          +BET(7)*( F(I-1,J+2)+F(I+1,J+2) )                              &
          +BET(8)*( F(I  ,J-1)+F(I  ,J+3) )                              &
          +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J+1) +BET(11)*F(I  ,J+2)     

     do J=3,NRL-2 
        FXY(I,J)=                                                      &
             GAM(1)*( F(I-2,J-2)+F(I-2,J+2)+F(I+2,J-2)+F(I+2,J+2) )    &
             +GAM(2)*( F(I-2,J+1)+F(I-2,J-1)+F(I-1,J-2)+F(I-1,J+2)     &
             +F(I+1,J-2)+F(I+1,J+2)+F(I+2,J-1)+F(I+2,J+1) )            &
             +GAM(3)*( F(I-2,J  )+F(I+2,J  )+F(I  ,J-2)+F(I  ,J+2) )   &
             +GAM(4)*( F(I-1,J-1)+F(I-1,J+1)+F(I+1,J-1)+F(I+1,J+1) )   &
             +GAM(5)*( F(I-1,J  )+F(I  ,J-1)+F(I  ,J+1)+F(I+1,J  ) )   &
             +GAM(6)*  F(I  ,J  )                                      
     end do

     J = NRL-1 
     FXY(I,J)=                                                          &
          BET(1)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J-3)+F(I+2,J-3) )        &
          +BET(2)*( F(I-2,J  )+F(I+2,J  ) )                             &
          +BET(3)*( F(I-2,J-1)+F(I+2,J-1) )                             &
          +BET(4)*( F(I-2,J-2)+F(I+2,J-2)+F(I-1,J+1)+F(I+1,J+1)         &
          +F(I-1,J-3)+F(I+1,J-3) )                                      &
          +BET(5)*( F(I-1,J  )+F(I+1,J  ) )                             &
          +BET(6)*( F(I-1,J-1)+F(I+1,J-1) )                             &
          +BET(7)*( F(I-1,J-2)+F(I+1,J-2) )                             &
          +BET(8)*( F(I  ,J+1)+F(I  ,J-3) )                             &
          +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J-1) +BET(11)*F(I  ,J-2)    

     J = NRL 
     FXY(I,J)=                                                          &
          ALP(1)*( F(I-2,J  )+F(I+2,J  ) )                              &
          +ALP(2)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J-3)+F(I+2,J-3)         &
          +F(I-1,J-4)+F(I+1,J-4) )                                      &
          +ALP(3)*( F(I-2,J-2)+F(I+2,J-2) )                             &
          +ALP(4)*( F(I-2,J-4)+F(I+2,J-4) )                             &
          +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )                             &
          +ALP(6)*( F(I-1,J-1)+F(I+1,J-1)+F(I-1,J-3)+F(I+1,J-3) )       &
          +ALP(7)*( F(I-1,J-2)+F(I+1,J-2) )                             &
          +ALP(8)*  F(I  ,J  )                                          &
          +ALP(9)*( F(I  ,J-1)+F(I  ,J-3) )                             &
          +ALP(10)* F(I  ,J-2) +ALP(11)*F(I  ,J-4)                      
  end do
 
  ! modified                                     
  do I=3,nset-2 
     do J=1,NRL 
        F(I,J) = FXY(I,J) 
     end do
  end do
  return 
end subroutine SMOOTHING

!*********************************************************************  
subroutine opaltab 
  use alink
  use cst
  use cf

  implicit none

  !  CODE FOR FITTING AND SMOOTHING OPAL DATA. ADAPTED FROM A CODE        
  !     WRITTEN BY MIKE SEATON(obtained june 1993)                        
  !                                                                       
  !     OPAL DATA.                                                        
  !     ASSUMES FIRST T6=0.006, LAST T6=10.OR 0.04). Depending on position
  !     in the table.                                                     
  !     USES RECTANGULAR ARRAY FOR VARIABLES T6 AND LOG10(R)              
  !                                                                       
  !     (1) NSM=NUMBER OF PASSES THROUGH SMOOTHING FILTER.                
  !     USE OF NSM=1 OR 2 IS RECOMMENDED.                                 
  !     NO SMOOTHING WITH NSM=0                                           
  !     (2) RANGE FOR LOG10(R),                                           
  !     RLS=FIRST VALUE, RLE=LAST VALE                                    
  !     (RLS MUST BE FIRST VALUYE IN TABLE)                               
  !                                                                       
  !  SUBROUTINE INTERP                                                    
  !     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF                   
  !               SUBROUTINE INTERP                                       
  !     WHICH GIVES LOG(ROSS) AND TWO FIRST DERIVATIVES FOR ANY           
  !     VALUES OF LOG(T) AND LOG(RHO). SEE BELOW FOR FURTHER              
  !     EXPLANATION.                                                      
  !                                                                       
  !  OUTPUT FOR THE CASE OF NSM>0.                                     
  !     INTERP IS USED TO OBTAIN SMOOTHED DATA INTERPOLATED               
  !     BACK TO THE ORIGINAL OPAL MESH. TWO FILES ARE WRITTEN.            
  !                                                                       
  !                                                                       
  !  THE SUBROUTINES SPLINE AND SPLINT ARE ADAPTED FROM THOSE GIVE BY     
  !  W.H. Press, S.A. Teulolsky, W.T. Vettering and B.P. Flannery,        
  !  "Numerical Recipes in FORTRAN", 2nd edn., 1992, C.U.P.               
  !  OTHER REFERENCES ARE MADE TO METHODS DESCRIBED IN THAT BOOK.         

  integer,parameter :: IP = 100
  real,dimension(IP) :: U, V, V2 
  real,dimension(IP,IPR) :: ROSSL
  logical :: IERR 

  integer :: i, j, ns, l, k
  real :: t6, flt, flr, flrho, g, dgdt, dgdrho

  NRL = 2*(RLE-RLS)+1 
  !     STORE LOG10(T) IN U AND LOG10(ROSS) IN ROSSL                      
  !     CHECK FIRST VALUE OF T6                                           
  T6 = t6arr(1) 
  do j=1,NRL 
     ROSSL(1,j) = coff(1,j) 
  enddo

  if( abs(T6 -.0056341325) < 1.d-8) then 
     U(1) = 6.+log10(T6) 
  endif
  !     SET ROSSL UP TO T6=t6arr(nset)                                    
  I = 1       
  do
     I = I+1                                                             
     T6 = t6arr(I)                                                       
     do j=1,NRL                                                        
        ROSSL(I,j) = coff(I,j)                                              
     enddo
     U(I) = 6.+log10(T6)                                              
     if(T6 >= tmax) exit
  end do
  NTEMP = I                                                               
  if(NTEMP > IP)then                                                   
     print*,' REQUIRE PARAMETER IP OF AT LEAST ',NTEMP    
     write(66,*)'540 - opaltab'
     write(66,*)' REQUIRE PARAMETER IP OF AT LEAST ',NTEMP
     stop                                                           
  endif

  !     DEFINE VARIABLES                                                  
  !         X=20.0*(LOG10(T)-3.80)+1                                      
  !         Y=2.0*(LOG10(R)-RLS)+1                                        
  !     USE INDICES I=1 TO nset AND J=1 TO NRL                            
  !     X AND Y ARE SUCH THAT, ON MESH-POINT (I,J), X=I AND Y=J           
  !     OBTAIN:-                                                          
  !         F(I,J)=LOG10(ROSS)                                            
  !         FX(I,J)=dF/dX                                                 
  !         FY(I,J)=dF/dY                                                 
  !         FXY(I,J)=ddF/dXdY                                             
  !                                                                       
  !     FIRST GET F AND FX, INTERPOLATING FROM OPAL T6 TO                 
  !     INTERVAL OF 0.05 IN LOG10(T).                                     
  do J=1,NRL                                                     
     !        FOR EACH LOG10(R), STORE LOG10(ROSS) IN V(I)                   
     do I=1,NTEMP                                                    
        V(I) = ROSSL(I,J)                                             
     end do
     !        GET FIRST DERIVATIVES AT END POINTS                            
     !        GET SECOND DERIVATIVES FOR SPLINE FIT                          
     call SPLINE(U,V,NTEMP,V2)                                          
     !        INTERPOLATE TO LOG10(T)=FLT, FLT=3.8(0.05)8.0                  
     do I=1,nset                                                 
        FLT = 3.75+0.05*I                                             
        call SPLINT(U,V,NTEMP,V2,FLT,F(I,J),FX(I,J))   
     end do
  end do

  !  OPTION FOR SMOOTHING                                                 
  if(NSM > 0)then                                                  
     do NS=1,NSM                                                 
        call SMOOTHING 
     end do
     call FITX                                                      
  endif

  !  GET FY AND FXY                                                       
  call FITY                                                         
  !  THE ARRAYS F, FX, FY AND FXY ARE NOW STORED                          
  !                                                                       
  !  CAN NOW DO INTERPOLATIONS USING                                      
  !       CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)                       
  !       INPUT IS FLT=LOG10(T), FLRHO=LOG10(RHO)                         
  !       OUTPUT IS G=LOG10(ROSS)                                         
  !              DGDT=dG/d(LOG10(T))                                      
  !            DGDRHO=dG/d(LOG10(RHO))                                    
  !              IERR=.TRUE. IF INPUT FLT, FLRHO ARE OUT-OF-RANGE,        
  !                          ELSE IERR=.FALSE.                            
  !                                                                       
  ! INTERPOLATE BACK TO OPAL POINTS                                       
  if(NSM > 0)then                                                  
     do l=1,NRL                                                     
        coff(1,l) = ROSSL(1,l)                                           
     enddo

     do K=2,NTEMP                                                    
        FLT = U(K)                                                    
        do L=nrlow,nrhigh                                        
           FLR = RLS+.5*(L-1)                                         
           FLRHO = FLR-18.+3.*FLT                                     
           call INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)                
           if(IERR)then                                             
           endif
           V(L) = G   
        end do
        T6 = t6arr(K)                                                 
        do l=nrlow,nrhigh                                           
           coff(K,l) = V(l)                                              
        enddo
     end do
  endif
  return                                                            

end subroutine opaltab
