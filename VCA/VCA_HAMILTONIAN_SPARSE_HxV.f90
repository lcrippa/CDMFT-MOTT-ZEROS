! > BUILD SPARSE HAMILTONIAN of the SECTOR
MODULE VCA_HAMILTONIAN_SPARSE_HxV
  USE VCA_HAMILTONIAN_COMMON
  implicit none
  private

  !>build sparse hamiltonian of the sector
  !public  :: buildH_c
  public :: vca_buildh_main
  !public :: vca_buildh_orbs
  !
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_main
  !public :: spMatVec_orbs
#ifdef _MPI
  public  :: spMatVec_MPI_main
  !public  :: spMatVec_MPI_orbs
#endif
  !
  !
  !> Related auxiliary routines:
  !public  :: setup_Hv_sector
  !public  :: delete_Hv_sector

  integer                                :: i,iup,idw
  integer                                :: j,jup,jdw
  integer                                :: m,mup,mdw
  integer                                :: ms,iud
  integer                                :: impi
  integer                                :: ilat,jlat,iorb,jorb,ispin,jspin,is,js,ibath
  integer                                :: kp,k1,k2,k3,k4
  integer                                :: ialfa,ibeta,indx
  real(8)                                :: sg1,sg2,sg3,sg4
  complex(8)                             :: htmp,htmpup,htmpdw
  logical                                :: Jcondition
  integer                                :: Nfoo,Nfoo2
  real(8),dimension(:,:,:,:),allocatable :: diag_hybr ![Nlat,Nspin,Norb,Nbath]
  real(8),dimension(:,:,:,:),allocatable :: bath_diag ![Nlat,Nspin,Norb/1,Nbath]

contains



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine vca_buildh_main(isector,Hmat)
    integer                               :: isector
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Htmp_up,Htmp_dw,Hrdx
    integer,dimension(Ns)                 :: ibup,ibdw
    integer,dimension(2*Ns_Ud)            :: Indices    ![2-2*Norb] 
    integer,dimension(Ns_Ud,Ns_Orb)       :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath] 
    integer,dimension(Nlat,Norb)          :: Nup,Ndw    !
    !
    nup=zero
    ndw=zero
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hstatus)stop "buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"vca_buildh_main","Hmat")
    !
    if(spH0d%status)call sp_delete_matrix(spH0d)
    !
    !Get diagonal hybridization, bath energy
    if(Nbath>0)then
      include "VCA_HAMILTONIAN/diag_hybr_bath.f90"         
    endif        
    !    
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,Dim)
    else
       call sp_init_matrix(spH0d,Dim)
    endif
#else
    call sp_init_matrix(spH0d,Dim)
#endif
    call sp_init_matrix(spH0dws(1),DimDw)
    call sp_init_matrix(spH0ups(1),DimUp)
    !
    !
    !-----------------------------------------------!
    !LOCAL CLUSTER HAMILTONIAN TERMS
    include "VCA_HAMILTONIAN/stored/H_cluster_local.f90"
    !
    !UP CLUSTER TERMS
    include "VCA_HAMILTONIAN/stored/H_cluster_up.f90"
    !
    !DW CLUSTER TERMS
    include "VCA_HAMILTONIAN/stored/H_cluster_dw.f90"
    !-----------------------------------------------!


    if(present(Hmat))then
       Hmat = zero
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=zero
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=zero
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat)
       else
          call sp_dump_matrix(spH0d,Hmat)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat)
#endif
       call sp_dump_matrix(spH0ups(1),Htmp_up)
       call sp_dump_matrix(spH0dws(1),Htmp_dw)
       Hmat = Hmat + kronecker_product(Htmp_dw,one*eye(DimUp)) ! kron(Htmp_dw,eye(DimUp))
       Hmat = Hmat + kronecker_product(one*eye(DimDw),Htmp_up) ! kron(eye(DimDw),Htmp_up)
       !
       deallocate(Htmp_up,Htmp_dw)
    endif
    !
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    return

    !if(present(Hmat))then
       !call assert_shape(Hmat,[dim,dim],"buildH_c","Hmat")
       !Hmat=zero
       !call sp_dump_matrix(spH0,Hmat)
    !endif
  end subroutine vca_buildh_main


  !####################################################################
  !        TODO: ORBS
  !####################################################################






  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  !+------------------------------------------------------------------+
  subroutine spMatVec_main(Nloc,v,Hv)
    integer                         :: Nloc
    complex(8),dimension(Nloc)      :: v
    complex(8),dimension(Nloc)      :: Hv
    complex(8)                      :: val
    integer                         :: i,iup,idw,j,jup,jdw,jj
    !
    !
    Hv=zero
    !
    do i = 1,Nloc
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       enddo
    enddo
    !
    !DW:
    do iup=1,DimUp
       !
       do idw=1,DimDw
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dws(1)%row(idw)%Size
             jup = iup
             jdw = spH0dws(1)%row(idw)%cols(jj)
             val = spH0dws(1)%row(idw)%vals(jj)
             j     = jup +  (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !UP:
    do idw=1,DimDw
       !
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j =  jup + (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
  end subroutine spMatVec_main


  !####################################################################
  !        TODO: ORBS
  !####################################################################




  !####################################################################
  !       MPI MATVEC
  !####################################################################

#ifdef _MPI
  subroutine spMatVec_mpi_main(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    !
    integer                             :: N
    complex(8),dimension(:),allocatable :: vt,Hvt
    complex(8)                          :: val
    integer                             :: i,iup,idw,j,jup,jdw,jj
    !local MPI
    integer                             :: irank
    !
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=zero
    do i=1,Nloc                 !==spH0%Nrow
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(i)
       end do
    end do
    !
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do idw=1,MpiQdw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          hxv_up: do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j   = jup + (idw-1)*DimUp
             Hv(i) = Hv(i) + val*v(j)
          end do hxv_up
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    !
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    allocate(vt(mpiQup*DimDw)) ;vt=zero
    allocate(Hvt(mpiQup*DimDw));Hvt=zero
    call vector_transpose_MPI(DimUp,MpiQdw,v,DimDw,MpiQup,vt)
    Hvt=zero    
    do idw=1,MpiQup             !<= Transposed order:  column-wise DW <--> UP  
       do iup=1,DimDw           !<= Transposed order:  column-wise DW <--> UP
          i = iup + (idw-1)*DimDw
          hxv_dw: do jj=1,spH0dws(1)%row(iup)%Size
             jup = spH0dws(1)%row(iup)%cols(jj)
             jdw = idw             
             j   = jup + (jdw-1)*DimDw
             val = spH0dws(1)%row(iup)%vals(jj)
             Hvt(i) = Hvt(i) + val*vt(j)
          end do hxv_dw
       enddo
    end do
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=zero
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv = Hv + Vt
  end subroutine spMatVec_mpi_main


  !####################################################################
  !        TODO: ORBS
  !####################################################################

#endif

END MODULE VCA_HAMILTONIAN_SPARSE_HxV
