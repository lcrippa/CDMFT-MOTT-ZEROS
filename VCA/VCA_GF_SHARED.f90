MODULE VCA_GF_SHARED
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_IO                     !< this contains the routine to print GF,Sigma and G0
  USE VCA_EIGENSPACE
  USE VCA_SETUP
  USE VCA_BATH_FUNCTIONS
  USE VCA_HAMILTONIAN
  USE VCA_AUX_FUNX
  !
  implicit none



  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer                :: state_cvec
  real(8)                                        :: state_e,max_exc

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm


  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:)       :: auxGmats,auxGreal 



contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    !print*,beta
    wm     = pi/beta*(2*arange(1,Lmats)-1)
    do i=0,Lmats
       vm(i) = pi/beta*2*i
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids


END MODULE VCA_GF_SHARED