MODULE VCA_IO
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_INPUT_VARS
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  !USE SF_MISC,only : assert_shape
  implicit none
  private



  !Retrieve self-energy through routines:
  interface vca_get_sigma_matsubara
     module procedure vca_get_sigma_matsubara_full
     module procedure vca_get_sigma_matsubara_ij
  end interface vca_get_sigma_matsubara

  interface vca_get_sigma_realaxis
     module procedure vca_get_sigma_real_full
     module procedure vca_get_sigma_real_ij
  end interface vca_get_sigma_realaxis


  !Retrieve imp GF through routines.
  interface vca_get_gimp_matsubara
     module procedure vca_get_gimp_matsubara_full
     module procedure vca_get_gimp_matsubara_ij
  end interface vca_get_gimp_matsubara

  interface vca_get_gimp_realaxis
     module procedure vca_get_gimp_real_full
     module procedure vca_get_gimp_real_ij
  end interface vca_get_gimp_realaxis


  interface vca_gf_cluster
     module procedure :: vca_gf_cluster_scalar
     module procedure :: vca_gf_cluster_array
  end interface vca_gf_cluster

  
  public :: vca_gf_cluster

  !
  public :: vca_get_sigma_matsubara
  public :: vca_get_sigma_realaxis
  !
  public :: vca_get_gimp_matsubara
  public :: vca_get_gimp_realaxis
  !
  public :: vca_print_impSigma
  public :: vca_print_impG
  public :: vca_print_impG0
  !
  public :: vca_read_impSigma
  !
  public :: vca_get_dens
  public :: vca_get_mag
  public :: vca_get_docc
  public :: vca_get_sft_potential


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,wr



contains

  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the cluster green's functions 
  !+-----------------------------------------------------------------------------+!
  include "VCA_IO/gf_cluster.f90"

  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!
  include "VCA_IO/get_sigma_matsubara.f90"
  include "VCA_IO/get_sigma_realaxis.f90"


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  include "VCA_IO/get_gimp_matsubara.f90"
  include "VCA_IO/get_gimp_realaxis.f90"




  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  include "VCA_IO/get_dens.f90"
  include "VCA_IO/get_mag.f90"
  include "VCA_IO/get_docc.f90"
  include "VCA_IO/get_sft_potential.f90"






  !+------------------------------------------------------------------+
  !                         PRINT SIGMA:
  !+------------------------------------------------------------------+  


  subroutine vca_print_impSigma
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,jorb,ispin
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !Print the impurity Sigma:
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
           do jorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                call splot("impSigma"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impSmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                call splot("impSigma"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impSreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
             enddo
            enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine vca_print_impSigma


  !+------------------------------------------------------------------+
  !                         PRINT G
  !+------------------------------------------------------------------+  


  subroutine vca_print_impG
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,ispin,jorb
    !
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb  
            do jorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                call splot("impG"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impGmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                call splot("impG"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impGreal(ilat,jlat,ispin,ispin,iorb,jorb,:))                
             enddo
            enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine vca_print_impG

  !+------------------------------------------------------------------+
  !                         PRINT G0
  !+------------------------------------------------------------------+  

  subroutine vca_print_impG0
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,ispin,jorb
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                call splot("impG0"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impG0mats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                call splot("impG0"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impG0real(ilat,jlat,ispin,ispin,iorb,jorb,:))
             enddo
            enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine vca_print_impG0

   ! PURPOSE: Read self-energy function(s) - also for inequivalent sites.
   !+-----------------------------------------------------------------------------+!
   subroutine vca_read_impSigma
     integer                                           :: i,ispin,isign,unit(2),iorb,jorb,ilat,jlat
     character(len=30)                                 :: suffix
     !
     if(.not.allocated(wm))allocate(wm(Lmats))
     if(.not.allocated(wr))allocate(wr(Lreal))
     wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
     wr     = linspace(wini,wfin,Lreal)
     !
     !!
     !Print the impurity functions:
     do ispin=1,Nspin
      do ilat=1,Nlat
        do jlat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                call sread("impSigma"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impSmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                call sread("impSigma"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impSreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
            enddo
          enddo
        enddo
      enddo
     enddo
     !
     if(allocated(wm))deallocate(wm)
     if(allocated(wr))deallocate(wr)
     !
   end subroutine vca_read_impSigma

END MODULE VCA_IO




