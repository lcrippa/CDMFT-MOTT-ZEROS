! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT
MODULE VCA_HAMILTONIAN_DIRECT_HxV
  USE VCA_HAMILTONIAN_COMMON
  implicit none
  private

  integer                              :: iiup,iidw
  integer                              :: iud,jj
  integer                              :: i,iup,idw
  integer                              :: j,jup,jdw
  integer                              :: m,mup,mdw
  integer                              :: ishift
  integer                              :: isector,jsector
  integer                              :: ms
  integer                              :: impi
  integer                              :: ilat,jlat,iorb,jorb,ispin,jspin,is,js,ibath
  integer                              :: kp,k1,k2,k3,k4
  integer                              :: ialfa,ibeta
  real(8)                              :: sg1,sg2,sg3,sg4
  complex(8)                           :: htmp,htmpup,htmpdw
  logical                              :: Jcondition
  integer                              :: Nfoo,Nfoo2
  real(8),dimension(:,:,:,:),allocatable :: diag_hybr ![Nlat,Nspin,Norb,Nbath]
  real(8),dimension(:,:,:,:),allocatable :: bath_diag ![Nlat,Nspin,Norb/1,Nbath]





  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  !public  :: directMatVec_orbs

#ifdef _MPI
  public  :: directMatVec_MPI_main
  !public  :: directMatVec_MPI_orbs
#endif


contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    integer,dimension(Ns)               :: ibup,ibdw
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath] 
    integer,dimension(Nlat,Norb)        :: Nup,Ndw 
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    if(Nbath>0)then
      include "VCA_HAMILTONIAN/diag_hybr_bath.f90"
    endif
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "VCA_HAMILTONIAN/direct/HxV_local.f90" 
    !
    !UP HAMILTONIAN TERMS
    include "VCA_HAMILTONIAN/direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "VCA_HAMILTONIAN/direct/HxV_dw.f90"
    !-----------------------------------------------!
    !
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    return
  end subroutine directMatVec_main

#ifdef _MPI
  subroutine directMatVec_MPI_main(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: Vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    integer,dimension(Ns)               :: ibup,ibdw
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Nlat,Norb)        :: Nup,Ndw 
    !
    integer,allocatable,dimension(:)    :: Counts,Displs
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    !Get diagonal hybridization, bath energy
    if(Nbath>0)then
      include "VCA_HAMILTONIAN/diag_hybr_bath.f90"
    endif
    !
    !    
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "VCA_HAMILTONIAN/direct_mpi/HxV_local.f90"
    !
    !NON-LOCAL TERMS:
    ! 
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "VCA_HAMILTONIAN/direct_mpi/HxV_up.f90"    
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw)) ;vt=zero
    allocate(Hvt(mpiQup*DimDw));Hvt=zero
    call vector_transpose_MPI(DimUp,MpiQdw,Vin,DimDw,MpiQup,vt) !Vin^T --> Vt
    include "VCA_HAMILTONIAN/direct_mpi/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=zero         !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
    Hv = Hv + Vt
    !-----------------------------------------------!
    !
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    return
  end subroutine directMatVec_MPI_main

#endif

END MODULE VCA_HAMILTONIAN_DIRECT_HxV
