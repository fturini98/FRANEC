! Questa e` la subroutine che interviene quando la OPAC va fuori dalle t
! La subroutine ha bisogno dell' abbondanza in massa di He, C ed O.     
! Inoltre serve Log_{10}(T) e Log_{R} dove R=(rho(g/cm3))/((T_6)**3)    

subroutine akappetta(Y,xC,xO,tlog,rlog,akapparad) 

  implicit none
  real :: Y,xC,xO,tlog,rlog,akapparad

  integer,save ::  ileggi = 0 
  real,parameter,dimension(29) :: rval = (/-8.0,-7.5,-7.0,-6.5,-6.0,     &
       -5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,  &
       1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0/)

  real,save,dimension(70,29) :: tabelio, tabcarb, tabossi 
  real,save,dimension(70) :: tval
  real,dimension(3) :: optot
  real,dimension(4,4,3) :: TR
  real,dimension(4,3) :: trop
  real :: rhigh, rlow, thigh, tlow 
  integer :: i, k, j
  integer :: indthigh, indtlow, indrlow,  indrhigh
  real :: ALFA, BETA, GAMMA, DELTA
  real,dimension(4) :: A,B

  !****************  Lettura e memorizzazione iniziali   **************   

  if (ileggi /= 123456) then 

     write(*,*)'entro nella kappetta' 

     open(40,file='extelio.in',status='old') 
     open(41,file='extcarbonio.in',status='old') 
     open(42,file='extossigeno.in',status='old') 


     do i=1,70 
        read(40,fmt='(F6.4,29(1X,F6.3))')tval(i),(tabelio(i,j),j=1,29) 
     end do
     do i=1,70 
        read(41,fmt='(F6.4,29(1X,F6.3))')tval(i),(tabcarb(i,j),j=1,29) 
     end do
     do i=1,70 
        read(42,fmt='(F6.4,29(1X,F6.3))')tval(i),(tabossi(i,j),j=1,29) 
     end do

     close(40) 
     close(41) 
     close(42) 
     ileggi = 123456 

     write(*,*)'Ho letto i files di puro He, C ed O' 

  end if

  !*********** Ricerca delle T e R prima e dopo Tlog e Rlog **************
  !************  e degli indici di tabella corrispondenti  ***************

  ! matt
  ! modificato schema di ricerca indici per renderlo piu' rapido
  if(tlog < 6.) then
     indtlow = int((tlog-3.75)/0.05)+1
  else if(tlog < 8.1) then
     indtlow = int((tlog-6.0)/0.1)+46
  else 
     indtlow = int((tlog-8.1)/0.2)+67
  endif
  indthigh = indtlow+1
  thigh = tval(indthigh) 
  tlow = tval(indtlow) 

  indrlow = int((rlog+8.)/0.5)+1 
  indrhigh = indrlow+1  
  rhigh = rval(indrhigh) 
  rlow = rval(indrlow) 

  !*******************  Interpolazione cubica    *************************
  ! per evitare uscite dalle      
  if (indtlow == 69) indtlow=68 
  ! tabelle                       
  if (indrlow == 28) indrlow=27 

  do i=1,4 
     do j=1,4 
        TR(j,i,1) = tabelio(indtlow+j-2,indrlow+i-2) 
     end do
  end do

  do i=1,4 
     do j=1,4 
        TR(j,i,2) = tabcarb(indtlow+j-2,indrlow+i-2) 
     end do
  end do

  do i=1,4 
     do j=1,4 
        TR(j,i,3) = tabossi(indtlow+j-2,indrlow+i-2) 
     end do
  end do

  do k=1,3 
     do j=1,4 
        B(1) = rval(indrlow-1) 
        B(2) = rval(indrlow) 
        B(3) = rval(indrlow+1) 
        B(4) = rval(indrlow+2) 

        A(1) = TR(j,1,k) 
        A(2) = TR(j,2,k) 
        A(3) = TR(j,3,k) 
        A(4) = TR(j,4,k) 

        call CUB(A,B,ALFA,BETA,GAMMA,DELTA) 

        trop(j,k) = ALFA*(rlog**3)+BETA*(rlog**2)+GAMMA*rlog+DELTA 
     end do
  end do

  do k=1,3 
     B(1) = tval(indtlow-1) 
     B(2) = tval(indtlow) 
     B(3) = tval(indtlow+1) 
     B(4) = tval(indtlow+2) 

     A(1) = trop(1,k) 
     A(2) = trop(2,k) 
     A(3) = trop(3,k) 
     A(4) = trop(4,k) 

     call CUB(A,B,ALFA,BETA,GAMMA,DELTA) 

     optot(k) = ALFA*(tlog**3)+BETA*(tlog**2)+GAMMA*tlog+DELTA 
  end do

  akapparad = (optot(1)*Y+optot(2)*xC+optot(3)*xO)/(Y+xC+xO) 

  return 

end subroutine akappetta
