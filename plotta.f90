subroutine PLOTTA(MAXME,RTOT) 
  use interfaccia
  use fisica
  use strut
  use chimic

  implicit none
  integer :: MAXME
  real :: RTOT

  real :: VA 
  character(len=1),dimension(51,120) :: A
  real,dimension(10) :: VAR, VARM, IVAV 
  character(len=1) :: BLANK = ' ', STAR = '*', SII = 'I', POINT = '.'
  character(len=1),dimension(10) :: SIMB =(/'R','L','P','T','H','3','4',&
       'C','N','O'/)
  character(len=1),dimension(11) :: PNUM =(/'0','9','8','7','6','5','4',&
       '3','2','1','0'/)

  integer :: k, l, n, j, jj, jjj, ll, m, mm, iva, mac, mu
  real :: pap, em, rat

  write(2,555) 
  do K=1,51 
     do L=1,120 
        A(K,L) = BLANK 
     end do
     A(K,12) = SII 
  end do
  N = 1 
  do K=1,11 
     A(N,7) = PNUM(1) 
     A(N,8) = POINT 
     A(N,9) = PNUM(K) 
     N = N+5 
  end do
  A(1,7) = PNUM(10) 
  VARM(1) = RTOT 
  do J=2,4 
     VARM(J) = G(J,1) 
     do JJ=2,MAXME 
        PAP = G(J,JJ) 
        VARM(J) = max(PAP,VARM(J)) 
     end do
  end do
  do J=1,6 
     JJJ = J+4 
     VARM(JJJ) = XXX(J,1) 
     do JJ=2,MAXME 
        VARM(JJJ) = max(VARM(JJJ),XXX(J,JJ)) 
     end do
  end do
  N = 0 
  do K=1,101 
     EM = .01* real(K-1) 
     EM = EM*EMTOT 
     if(EM > G(5,MAXME)) cycle
     N = N+1 
     do L=2,MAXME 
        if(EM-G(5,L) <= 0) then 
           exit
        else 
           cycle
        endif
     end do
     LL = L-1 
     if(LL > (MAXME-1)) LL = MAXME-1 
     RAT = (EM-G(5,LL))/(G(5,LL+1)-G(5,LL)) 
     do M=1,4 
        VAR(M) = G(M,LL)+RAT*(G(M,LL+1)-G(M,LL)) 
     end do
     do M=1,6 
        MM = M+4 
        VAR(MM) = XXX(M,LL)+RAT*(XXX(M,LL+1)-XXX(M,LL)) 
     end do
     MAC = K+11 
     do M=1,10 
        if(VARM(M) == 0) then 
           VA = 50.26 
        else 
           VA = (VAR(M)/VARM(M))*50.+.26 
        end if
        IVA = 51-int(VA) 
        if(IVA < 1) IVA = 1 
        if(M == 1) then
           A(IVA,MAC) = SIMB(M)
           IVAV(M) = IVA
           cycle
        endif
        MU = M-1 
        do MM=1,MU 
           if(IVA == IVAV(MM)) goto 14 
        end do
        A(IVA,MAC) = SIMB(M) 
        goto 13 
14      A(IVA,MAC) = STAR 
13      IVAV(M) = IVA 
     end do
  end do
  N = N+11 
  write(2,333) 
  write(2,222) 
  do M=1,51 
     write(2,111) (A(M,K),K=1,N) 
  end do
  write(2,222) 
  write(2,333) 
  VARM(1) = VARM(1)*1.d10 
  VARM(2) = VARM(2)*1.d32 
  VARM(3) = VARM(3)*1.d17 
  VARM(4) = VARM(4)*1.d+6 
  write(2,666)(VARM(K),K=1,10) 
  write(2,444) 
  return 

111 format(120A1) 
222 format(11X,'I---------I---------I---------I---------I---------I---&
       &------I---------I---------I---------I---------I')
333 format(10X,'0.0',7X,'0.1',7X,'0.2',7X,'0.3',7X,'0.4',7X,'0.5',7X, &
         '0.6',7X,'0.7',7X,'0.8',7X,'0.9',7X,'1.0',3X,'M/MTOT') 
444 format(//////) 
555 format(1H1) 
666 format(//,1X,'R=',1PE9.3,2X,'L=',E9.3,2X,'P=',E9.3,2X,'T=',E9.3,2X,&
       'H=',E9.3,2X,'3=',E9.3,2X,'4=',E9.3,2X,'C=',E9.3,2X,'N=',E9.3,2X,&
       'O=',E9.3)
end subroutine PLOTTA
