subroutine edoradmodif(iz1,nel,nmix,iz,fa,TEMP,flt,flrho,ri,grl)
  use opserver_dp
  
  implicit none
  
  integer,parameter :: nch=5
  real,parameter :: chil0=-0.1, dchil=0.05
  integer :: nel, nmix
  integer :: iz1
  integer,dimension(nel,nmix) :: iz
  real,dimension(nel,nmix) :: fa
  real :: TEMP
  real,dimension(nmix) :: flt, flrho, ri
  real,dimension(nmix,nch) :: grl
  
  integer :: j, k
  integer,dimension(nmix) :: imix
  
  real,dimension(nmix,nch) :: zet, g, gp, fross, fp
  integer :: ierr
  
  do j=1,nmix
     imix(j)=j
  enddo
  
  call op_ax(iz1, nel, nmix, nch, nmix, iz, fa, chil0, dchil, TEMP, &
       flt, flrho, ri, imix, zet, g, gp, fross, fp, grl, ierr)
  write(*,*)"controllo routine g_rad: ierr =",ierr
  return
end subroutine edoradmodif
