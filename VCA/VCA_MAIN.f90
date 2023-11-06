module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_OBSERVABLES
  USE VCA_GREENS_FUNCTIONS
  USE VCA_BATH_SETUP
  USE VCA_EIGENSPACE
  USE VCA_HAMILTONIAN
  USE VCA_OMEGA
  !
  USE SF_LINALG,  only: eigh,diag,kron,eye
  USE SF_IOTOOLS, only: str,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape,sort_array
  implicit none
  private


  !>INIT VCA SOLVER
  interface vca_init_solver
     module procedure :: vca_init_solver_serial
#ifdef _MPI
     module procedure :: vca_init_solver_mpi
#endif
  end interface vca_init_solver
  !>
  public :: vca_init_solver

  !>VCA SOLVER
  interface vca_solve
     module procedure :: vca_solve_serial
#ifdef _MPI
     module procedure :: vca_solve_mpi
#endif
  end interface vca_solve
  !>
  public :: vca_solve




contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver_serial(bath)
    real(8),intent(inout),optional :: bath(:)
    logical,save                   :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    if(present(bath))then
       if(.not.check_bath_dimension(bath))stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       bath=0d0
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
       call vca_get_bath(vca_bath,bath)
    else
      write(LOGfile,"(A)") "Bath not present, setting Nbath to 0"
      Nbath=0
    endif
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    if(isetup)call setup_global
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
    isetup=.false.
    !
  end subroutine vca_init_solver_serial

#ifdef _MPI
  subroutine vca_init_solver_mpi(MpiComm,bath)
    integer                                     :: MpiComm
    real(8),dimension(:),intent(inout),optional :: bath
    !complex(8),intent(in)                      :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                                     :: check 
    logical,save                                :: isetup=.true.
    integer                                     :: i
    !
    call vca_set_mpicomm(MpiComm)
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    if(present(bath))then
       if(.not.check_bath_dimension(bath))stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       bath=0d0
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
       call vca_get_bath(vca_bath,bath)
    else
      write(LOGfile,"(A)") "Bath not present, setting Nbath to 0"
      Nbath=0
    endif
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    !Init bath:
    !call set_hloc(Hloc)
    !
    !
    if(isetup)call setup_global
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
    isetup=.false.
    !
    call vca_del_MpiComm()
    !
  end subroutine vca_init_solver_mpi
#endif




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_solve_serial(Hloc,Hk,bath,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    complex(8),intent(in),dimension(:,:,:,:,:,:)   :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: Hk ![Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nktot]
    real(8),intent(inout),dimension(:),optional    :: bath
    !
    real(8),optional                               :: Uloc_ii(Nlat,Norb)
    real(8),optional                               :: Ust_ii(Nlat)
    real(8),optional                               :: Jh_ii(Nlat)
    real(8),optional                               :: Jp_ii(Nlat)
    real(8),optional                               :: Jx_ii(Nlat)
    !
    integer                                        :: Lk,Nsites
    integer                                        :: ilat,jlat
    integer                                        :: iorb,jorb
    integer                                        :: ispin

    integer                                        :: i,j
    real(8)                                        :: Tr,omega_integral,omegaprime
    integer                                        :: unit

    !
    write(LOGfile,"(A)")"SOLVING VCA"
    !
    if(rank(Hloc) .ne. 6) STOP "STOP: wrong cluster matrix dimensions"
    if(rank(Hk)   .ne. 7) STOP "STOP: wrong lattice matrix dimensions"
    !
    !INITIALIZE THE INTERNAL INTERACTION COEFFICIENTS
    !
    if(allocated(Uloc_per_site))deallocate(Uloc_per_site)
    if(allocated(Ust_per_site))deallocate(Ust_per_site)
    if(allocated(Jh_per_site))deallocate(Jh_per_site)
    if(allocated(Jx_per_site))deallocate(Jx_per_site)
    if(allocated(Jp_per_site))deallocate(Jp_per_site)
    !
    allocate(Uloc_per_site(Nlat,Norb));Uloc_per_site=zero
    allocate(Ust_per_site(Nlat));Ust_per_site=zero
    allocate(Jx_per_site(Nlat));Jx_per_site=zero
    allocate(Jp_per_site(Nlat));Jp_per_site=zero
    allocate(Jh_per_site(Nlat));Jh_per_site=zero
    !
    do ilat=1,Nlat
      do iorb=1,Norb
        if(present(Uloc_ii))then
          Uloc_per_site(ilat,iorb)=Uloc_ii(ilat,iorb)
        else
          Uloc_per_site(ilat,iorb)=Uloc(iorb)
        endif
      enddo
      !
      if(present(Ust_ii))then
        Ust_per_site(ilat)=Ust_ii(ilat)
      else
        Ust_per_site(ilat)=Ust
      endif
      !
      if(present(Jh_ii))then
        Jh_per_site(ilat)=Jh_ii(ilat)
      else
        Jh_per_site(ilat)=Jh
      endif
      !
      if(present(Jx_ii))then
        Jx_per_site(ilat)=Jx_ii(ilat)
      else
        Jx_per_site(ilat)=Jx
      endif
      !
      if(present(Jp_ii))then
        Jp_per_site(ilat)=Jp_ii(ilat)
      else
        Jp_per_site(ilat)=Jp
      endif
    enddo
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath,vca_bath)
       call vca_write_bath(vca_bath,LOGfile)
       call vca_save_bath(vca_bath,used=.true.)
    endif
    !
    !
    select case(vca_sparse_H)
      case (.true.)
        spHtimesV_p => spMatVec_main
      case (.false.)
        spHtimesV_p => directMatVec_main
      case default
        stop "vca_solve_single ERROR: vca_sparse_H undefined"
    end select
    !
    !
    !GENERATE THE CLUSTER HAMILTONIAN AND THE HOPPING MATRIX FOR THE LATTICE
    !
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(impHk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Hk,7)))
    impHloc=zero
    impHk=zero
    call vca_set_Hcluster(Hloc)
    call vca_set_Hk(Hk)
    !
    !GET CLUSTER GREEN'S FUNCTION AND GROUND STATE ENERGY
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    if(print_observables)then
      call observables_cluster()    !obtain impurity observables as thermal averages.
      call get_custom_observables()    !obtain impurity observables as thermal averages.
    endif
    !
    !CALCULATE THE VARIATIONAL GRAND POTENTIAL
    !
    omegaprime=0.d0
    omega_integral=0.d0
    !
    if(finiteT)then
      do i=1,state_list%size
        omegaprime=omegaprime+exp(-beta*(es_return_energy(state_list,i)-state_list%emin))
      enddo
      omegaprime=state_list%emin-(1.d0/beta)*log(omegaprime)
      omega_integral=frequency_integration_finite_t()
    else
      omegaprime=state_list%emin
      omega_integral=frequency_integration()
    endif
    !
    sft_potential = omegaprime-omega_integral
    !
    write(LOGfile,"(A,10f18.12,A)")"EGS PER SITE",state_list%emin/NLAT
    write(LOGfile,"(A,10f18.12,A)")"OMEGA POTENTIAL PER SITE=",(state_list%emin-omega_integral)/NLAT
    open(free_unit(unit),file="SFT_potential.vca",position='append')
    write(unit,*)sft_potential,omegaprime,-omega_integral
    close(unit)
    !
    !CLEAN UP
    !
    call es_delete_espace(state_list)
    nullify(spHtimesV_p)
    !
    !
    deallocate(impHloc)
    deallocate(impHk)
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
  end subroutine vca_solve_serial


#ifdef _MPI

  subroutine vca_solve_mpi(MpiComm,Hloc,Hk,bath,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    complex(8),intent(in),dimension(:,:,:,:,:,:)   :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: Hk   ![Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nktot]
    real(8),intent(inout),dimension(:),optional    :: bath
    !
    real(8),optional                               :: Uloc_ii(Nlat,Norb)
    real(8),optional                               :: Ust_ii(Nlat)
    real(8),optional                               :: Jh_ii(Nlat)
    real(8),optional                               :: Jp_ii(Nlat)
    real(8),optional                               :: Jx_ii(Nlat)
    !
    integer                                        :: Lk,Nsites
    integer                                        :: ilat,jlat
    integer                                        :: iorb,jorb
    integer                                        :: ispin

    integer                                        :: i,j
    real(8)                                        :: Tr,omega_integral,omegaprime
    integer                                        :: unit
    integer                                        :: MpiComm
    logical                                        :: check
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    !
    call vca_set_MpiComm(MpiComm)
    !
    if(rank(Hloc) .ne. 6) STOP "STOP: wrong cluster matrix dimensions"
    if(rank(Hk)   .ne. 7) STOP "STOP: wrong lattice matrix dimensions"
    !    
    !INITIALIZE THE INTERNAL INTERACTION COEFFICIENTS
    !
    if(allocated(Uloc_per_site))deallocate(Uloc_per_site)
    if(allocated(Ust_per_site))deallocate(Ust_per_site)
    if(allocated(Jh_per_site))deallocate(Jh_per_site)
    if(allocated(Jx_per_site))deallocate(Jx_per_site)
    if(allocated(Jp_per_site))deallocate(Jp_per_site)
    !
    allocate(Uloc_per_site(Nlat,Norb));Uloc_per_site=zero
    allocate(Ust_per_site(Nlat));Ust_per_site=zero
    allocate(Jx_per_site(Nlat));Jx_per_site=zero
    allocate(Jp_per_site(Nlat));Jp_per_site=zero
    allocate(Jh_per_site(Nlat));Jh_per_site=zero
    !
    do ilat=1,Nlat
      do iorb=1,Norb
        if(present(Uloc_ii))then
          Uloc_per_site(ilat,iorb)=Uloc_ii(ilat,iorb)
        else
          Uloc_per_site(ilat,iorb)=Uloc(iorb)
        endif
      enddo
      !
      if(present(Ust_ii))then
        Ust_per_site(ilat)=Ust_ii(ilat)
      else
        Ust_per_site(ilat)=Ust
      endif
      !
      if(present(Jh_ii))then
        Jh_per_site(ilat)=Jh_ii(ilat)
      else
        Jh_per_site(ilat)=Jh
      endif
      !
      if(present(Jx_ii))then
        Jx_per_site(ilat)=Jx_ii(ilat)
      else
        Jx_per_site(ilat)=Jx
      endif
      !
      if(present(Jp_ii))then
        Jp_per_site(ilat)=Jp_ii(ilat)
      else
        Jp_per_site(ilat)=Jp
      endif
    enddo
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath,vca_bath)
       call vca_write_bath(vca_bath,LOGfile)
       call vca_save_bath(vca_bath,used=.true.)
    endif
    !
    select case(vca_sparse_H)
      case (.true.)
        spHtimesV_p => spMatVec_main
      case (.false.)
        spHtimesV_p => directMatVec_main
      case default
        stop "vca_solve_single ERROR: vca_sparse_H undefined"
    end select
    !
    !
    !GENERATE THE CLUSTER HAMILTONIAN AND THE HOPPING MATRIX FOR THE LATTICE
    !
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(impHk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Hk,7)))
    impHloc=zero
    impHk=zero
    call vca_set_Hcluster(Hloc)
    call vca_set_Hk(Hk)

    !
    !GET CLUSTER GREEN'S FUNCTION AND GROUND STATE ENERGY
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    if(print_observables)then
      call observables_cluster()       !obtain impurity observables as thermal averages.
      call get_custom_observables()    !obtain impurity observables as thermal averages.
    endif
    !call save_gfprime("gfprime",use_formatted=.true.)
    !call read_gfprime("gfprime",use_formatted=.true.)
    !call reconstruct_g()
    !
    !CALCULATE THE VARIATIONAL GRAND POTENTIAL
    !
    omegaprime=0.d0
    omega_integral=0.d0
    !
    if(finiteT)then
      do i=1,state_list%size
        omegaprime=omegaprime+exp(-beta*(es_return_energy(state_list,i)-state_list%emin))
      enddo
      omegaprime=state_list%emin-(1.d0/beta)*log(omegaprime)
      if(MPIMASTER) omega_integral=frequency_integration_finite_t()
    else
      omegaprime=state_list%emin
      if(MPIMASTER) omega_integral=frequency_integration()
    endif
    !
    if(MPIMASTER) sft_potential = omegaprime-omega_integral
    call Bcast_MPI(MpiComm,sft_potential)
    !
    if(MPIMASTER)then
      write(LOGfile,"(A,10f18.12,A)")"EGS PER SITE",omegaprime/NLAT
      write(LOGfile,"(A,10f18.12,A)")"OMEGA POTENTIAL PER SITE=",(omegaprime-omega_integral)/NLAT
      open(free_unit(unit),file="SFT_potential.vca",position='append')
      write(unit,*)sft_potential,omegaprime,-omega_integral
      close(unit)
    endif
    !
    call barrier_mpi(MpiComm)
    !
    !CLEAN UP
    !
    call es_delete_espace(state_list)
    nullify(spHtimesV_p)
    !
    deallocate(impHloc)
    deallocate(impHk)
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    call vca_del_MpiComm()
    !
  end subroutine vca_solve_mpi
#endif


end module VCA_MAIN








! Nexc     = Qcluster%Nexc
! Nexc_sys = Ncopies*Nexc
! call allocate_Qmatrix(Qsystem)
! !
! print*,"Nexc_sys=",Nexc_sys
! allocate(Mmat(Nexc_sys,Nexc_sys))
! allocate(Lvec(Nexc_sys));Lvec=zero
! !
! Nc = vca_get_cluster_dimension(present(bath))
! !
! Mmat=zero
! !
! do icopy=1,Ncopies
!    do jcopy=1,Ncopies
!       i1 = 1 + (icopy-1)*Nc
!       i2 = icopy*Nc
!       j1 = 1 + (jcopy-1)*Nc
!       j2 = jcopy*Nc
!       !
!       ii1 = 1 + (icopy-1)*Nexc
!       ii2 = icopy*Nexc
!       jj1 = 1 + (jcopy-1)*Nexc
!       jj2 = jcopy*Nexc
!       !
!       Mmat(ii1:ii2,jj1:jj2)  = matmul(Qcluster%cdg,  matmul(Vmat(i1:i2,j1:j2),Qcluster%c))
!    enddo
! enddo
! !
! Mmat = Mmat + kron(eye(Ncopies),diag( Qcluster%poles ))
! !
! call eigh(Mmat,Lvec)!>from here on M-->S
! !
! Qsystem%poles = Lvec
! call sort_array(Qcluster%poles)
! Lvec=zero
! do iexc=1,Nexc
!    do icopy=1,Ncopies
!       i = icopy + (iexc-1)*Ncopies
!       Lvec(i) = Qcluster%poles(iexc)
!    enddo
! enddo
! !
! Tr=0d0
! do i=1,Nexc_sys
!    !
!    arg1 = beta*Lvec(i)
!    arg2 = beta*Qsystem%poles(i)
!    if(arg1 < -20d0 .OR. arg2 < -20d0)then
!       Tr  = Tr - beta*( Lvec(i) - Qsystem%poles(i) )
!    else
!       Tr  = Tr + log( (1d0 + exp(-beta*Lvec(i)))/(1d0 + exp(-beta*Qsystem%poles(i))) )
!    endif
! enddo
! !
