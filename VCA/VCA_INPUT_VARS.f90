MODULE VCA_INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  implicit none

  !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile


  !INPUT VARIABLES (to be exported in MAIN module)
  !input variables
  !=========================================================
  integer              :: Norb         !# of lattice orbitals per site
  integer              :: Nspin        !# spin degeneracy (max 2)
  integer              :: Nlat         !# size of cluster
  integer              :: Nbath        !# of bath sites (per orbital or not depending on bath_type)
  real(8),dimension(5) :: Uloc         !local interactions
  real(8)              :: Ust          !intra-orbitals interactions
  real(8)              :: Jh           !J_Hund: Hunds' coupling constant 
  real(8)              :: Jx           !J_X: coupling constant for the spin-eXchange interaction term
  real(8)              :: Jp           !J_P: coupling constant for the Pair-hopping interaction term 
  real(8)              :: xmu          !chemical potential
  real(8)              :: beta         !inverse temperature
  real(8)              :: eps          !broadening
  real(8)              :: bandwidth    !noninteracting original bandwidth
  real(8)              :: wini,wfin    !
  logical              :: chiflag      !
  logical              :: Jhflag       !spin-exchange and pair-hopping flag.
  logical              :: HFmode       !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  logical              :: HFshift      !flag for half-filling bath local energy
  real(8)              :: cutoff       !cutoff for spectral summation
  real(8)              :: gs_threshold !Energy threshold for ground state degeneracy loop up
  real(8)              :: sb_field     !symmetry breaking field
  logical              :: print_Sigma  !flag to print impurity Green`s functions
  logical              :: print_impG   !flag to print impurity Green`s functions
  logical              :: print_impG0  !flag to print impurity Green`s functions
  logical              :: print_observables  !flag to calculate and print observables
  logical              :: vca_twin    !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  logical              :: vca_sectors          !flag to reduce sector scan for the spectrum to specific sectors +/- vca_sectors_shift
  integer              :: vca_sectors_shift    !shift to the vca_sectors scan
  ! character(len=7)     :: bath_type           !flag to set bath type: normal (1bath/imp), hybrid(1bath)  !FIXME: MAYBE ADD HYBRID
  !
  !
  character(len=12)    :: lanc_method         !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only), DVDSON (no MPI)
  real(8)              :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer              :: lanc_niter          !Max number of Lanczos iterations
  integer              :: lanc_niter_spectrum !Max number of iterations in the spectrum annealing
  integer              :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer              :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer              :: lanc_nstates_step   !Number of states added at each step to determine the optimal spectrum size at finite T
  integer              :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
  real(8)              :: lanc_spectrum_threshold  !Threshold for the spectrum annealing error.
  logical              :: vca_sparse_H         !flag to select storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE
  logical              :: vca_total_ud         !flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw
  logical              :: vca_solve_offdiag_gf !flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal 
  logical              :: vca_gf_symmetric     !flag to select the calculation of the off-diagonal impurity GF. if T the 2-channel method is used 
  !

  real(8)              :: nread        !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)              :: nerr         !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)              :: ndelta       !initial chemical potential step
  integer              :: niter        !
  integer              :: verbose      !


  !Some parameters for function dimension:
  !=========================================================
  integer              :: Lmats
  integer              :: Lreal
  integer              :: Ltau


  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)   :: Hfile,HLOCfile
  integer,save         :: LOGfile





contains




  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine vca_read_input(INPUTunit,comm)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    integer,optional :: comm
    logical          :: master=.true.
    integer          :: i,rank=0
#ifdef _MPI
    if(present(comm))then
       master=get_Master_MPI(comm)
       rank  =get_Rank_MPI(comm)
    endif
#endif
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of orbitals per cluster site.")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy")
    call parse_input_variable(Nlat,"NLAT",INPUTunit,default=1,comment="Number of cluster copies tiling the system")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=0,comment="Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)")
    !
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2.d0,0.d0,0.d0,0.d0,0.d0],comment="Values of the local interaction per orbital")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(jhflag,"JHFLAG",INPUTunit,default=.false.,comment="Flag to include full rotational invariant terms: SE, PH.")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="Spin-Exchance coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="Pair-Hopping coupling")
    !
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 it is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(hfshift,"HFSHIFT",INPUTunit,default=.true.,comment="Half-filling shift of the bath local energy: if TRUE, zero is the half-filling chemical potential.")
    !
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=5000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of Real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000,comment="Number of Imaginary time points.")
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(chiflag,"CHIFLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(print_Sigma,"PRINT_SIGMA",INPUTunit,default=.true.,comment="Flag to print impurity Self-energy")
    call parse_input_variable(print_impG,"PRINT_IMPG",INPUTunit,default=.true.,comment="Flag to print impurity interacting Greens function")
    call parse_input_variable(print_impG0,"PRINT_IMPG0",INPUTunit,default=.true.,comment="Flag to print impurity non-interacting Greens function")
    call parse_input_variable(print_observables,"PRINT_OBSERVABLES",INPUTunit,default=.true.,comment="Flag to calculate and print observables")
    ! call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal',comment="flag to set bath type: normal (1bath/imp), hybrid(1bath)")
    !
    !
    call parse_input_variable(vca_twin,"VCA_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    !
    call parse_input_variable(vca_sectors,"VCA_SECTORS",INPUTunit,default=.false.,comment="flag to reduce sector scan for the spectrum to specific sectors +/- vca_sectors_shift.")
    call parse_input_variable(vca_sectors_shift,"VCA_SECTORS_SHIFT",INPUTunit,1,comment="shift to vca_sectors")
    call parse_input_variable(vca_sparse_H,"VCA_SPARSE_H",INPUTunit,default=.true.,comment="flag to select storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    call parse_input_variable(vca_total_ud,"VCA_TOTAL_UD",INPUTunit,default=.true.,comment="flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw")
    call parse_input_variable(vca_solve_offdiag_gf,"VCA_SOLVE_OFFDIAG_GF",INPUTunit,default=.false.,comment="flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal")
    call parse_input_variable(vca_gf_symmetric,"VCA_GF_SYMMETRIC",INPUTunit,default=.false.,comment="flag to select the calculation of the off-diagonal impurity GF. if T the 2-channel method is used")
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only)") 
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=6,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=3,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=5,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_niter_spectrum,"LANC_NITER_SPECTRUM",INPUTunit,default=10,comment="Number of iterations in spectrum annealing.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=0.000000000001d0,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=256,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    call parse_input_variable(lanc_spectrum_threshold,"LANC_SPECTRUM_THRESHOLD",INPUTunit,default=1d-5,comment="Threshold for the spectrum annealing error.")
    !
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    !
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")
    call parse_input_variable(verbose,"VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 3:all.")
    !
#ifdef _MPI
    if(present(comm))then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".vca")
          do i=1,get_Size_MPI(comm)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
      call print_input()
      call save_input_file(INPUTunit)
      call scifor_version()
      call code_version(revision)
    endif
  end subroutine vca_read_input

  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop



END MODULE VCA_INPUT_VARS
