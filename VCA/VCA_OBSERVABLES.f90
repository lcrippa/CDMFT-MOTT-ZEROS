MODULE VCA_OBSERVABLES
  USE VCA_INPUT_VARS
  USE VCA_BATH_FUNCTIONS
  USE VCA_VARS_GLOBAL
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_IO, only:vca_gf_cluster
  USE VCA_GF_SHARED, only:max_exc
  USE VCA_AUX_FUNX
  USE VCA_EIGENSPACE
  USE VCA_HAMILTONIAN
  !
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE SF_INTEGRATE
  !
  implicit none
  private
  !
  interface add_custom_observable
     module procedure :: add_custom_observable_local
     module procedure :: add_custom_observable_kdep
  end interface add_custom_observable
  
  public :: observables_cluster
  public :: init_custom_observables
  public :: add_custom_observable
  public :: get_custom_observables
  public :: clear_custom_observables

  logical,save                                             :: iolegend=.true.
  real(8),dimension(:,:),allocatable                       :: dens,dens_up,dens_dw
  real(8),dimension(:,:),allocatable                       :: docc
  real(8),dimension(:,:),allocatable                       :: magz
  real(8),dimension(:,:,:),allocatable                     :: sz2
  real(8),dimension(:,:,:),allocatable                     :: zimp,simp
  real(8)                                                  :: Egs
  real(8)                                                  :: s2tot
  real(8)                                                  :: Ei
  real(8)                                                  :: integrationR
  !
  integer                                                  :: r,m,k,k1,k2,Nkpts
  integer                                                  :: iup,idw
  integer                                                  :: jup,jdw
  integer                                                  :: mup,mdw
  real(8)                                                  :: sgn,sgn1,sgn2,sg1,sg2
  real(8)                                                  :: gs_weight
  !
  real(8)                                                  :: peso
  real(8)                                                  :: norm
  !
  integer                                                  :: i,j,ii
  integer                                                  :: isector,jsector
  integer                                                  :: idim,idimUP,idimDW
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)          :: sij
  !
  complex(8),dimension(:),pointer                          :: state_cvec
  logical                                                  :: Jcondition



contains 

  !+---------------------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties for the cluster
  !+---------------------------------------------------------------------------------+

  subroutine observables_cluster()
    integer :: iorb
    write(LOGfile,"(A)")"Get cluster observables:"
    call vca_lanc_observables
  end subroutine observables_cluster

  subroutine init_custom_observables(N,Nktot_)
    integer                      :: N,Nktot_
    !
    if(MpiMaster)then
      custom_o%N_filled=0
      custom_o%N_asked=N
      allocate(custom_o%item(N))
      custom_o%init=.true.
      Nkpts=Nktot_
    endif
    !
  end subroutine init_custom_observables
    
  subroutine add_custom_observable_local(o_name,sij)
    integer                               :: i
    complex(8),dimension(:,:,:,:,:,:)     :: sij
    character(len=*)                      :: o_name
    !
    if(MpiMaster .and. custom_o%init)then
      if(custom_o%N_filled .gt. custom_o%N_asked)then
        STOP "add_custom_observable: too many observables given"
        call clear_custom_observables
      endif
      !
      custom_o%N_filled=custom_o%N_filled+1
      custom_o%item(custom_o%N_filled)%o_name=o_name
      custom_o%item(custom_o%N_filled)%o_value=0.d0
      !
      allocate(custom_o%item(custom_o%N_filled)%sij(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts))
      do i=1,size(custom_o%item(custom_o%N_filled)%sij,7)
        custom_o%item(custom_o%N_filled)%sij(:,:,:,:,:,:,i)=sij
      enddo
    else
      STOP "add_custom_observable: custom observables not initialized"
    endif
  end subroutine add_custom_observable_local


  subroutine add_custom_observable_kdep(o_name,sijk)
    integer                               :: i
    complex(8),dimension(:,:,:,:,:,:,:)   :: sijk
    character(len=*)                      :: o_name
    !
    if(MpiMaster .and. custom_o%init)then
      if(custom_o%N_filled .gt. custom_o%N_asked)then
        STOP "add_custom_observable: too many observables given"
        call clear_custom_observables
      endif
      !
      custom_o%N_filled=custom_o%N_filled+1
      custom_o%item(custom_o%N_filled)%o_name=o_name
      custom_o%item(custom_o%N_filled)%o_value=0.d0
      !
      allocate(custom_o%item(custom_o%N_filled)%sij(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts))
      custom_o%item(custom_o%N_filled)%sij=sijk
    else
      STOP "add_custom_observable: custom observables not initialized"
    endif
  end subroutine add_custom_observable_kdep


  subroutine get_custom_observables()
    integer            :: i
    !
    if(MpiMaster .and. custom_o%init)then
      if(custom_o%N_filled .eq. 0)then
        write(Logfile,*)"WARNING! Custom observables initialized but none given."
        RETURN
      endif
      !
      write(LOGfile,*)"Calculating custom observables"
      !
      allocate(sij(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts))
      sij=zero
      !
      do i=1,custom_o%N_filled
        sij=custom_o%item(i)%sij
        if(finiteT) then
          custom_o%item(i)%o_value=calculate_observable_integral_finite_t()
        else
          custom_o%item(i)%o_value=calculate_observable_integral_zero_t()
        endif
        write(LOGfile,"(A,10f18.12,A)")reg(custom_o%item(i)%o_name)//" = ",custom_o%item(i)%o_value
      enddo
      call write_custom_legend()
      call write_custom_observables()
      deallocate(sij)
    endif
    !
  end subroutine get_custom_observables
  

  subroutine clear_custom_observables()
    integer                       :: i
    if(MpiMaster .and. custom_o%init)then 
      do i=1,custom_o%N_filled
        deallocate(custom_o%item(i)%sij)
        custom_o%item(i)%o_name=""
        custom_o%item(i)%o_value=0.d0
      enddo
      custom_o%N_asked=0
      custom_o%N_filled=0
      custom_o%init=.false.
    endif
  end subroutine clear_custom_observables

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties for the cluster
  !+-------------------------------------------------------------------+
  subroutine vca_lanc_observables()
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer                             :: i,j,ii
    integer                             :: istate
    integer                             :: isector,jsector
    integer                             :: idim,jdim,ilat
    integer                             :: isz,jsz
    integer                             :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                             :: numstates
    integer                             :: r,m,k
    real(8)                             :: sgn,sgn1,sgn2
    real(8)                             :: weight
    real(8)                             :: Ei
    real(8)                             :: peso
    real(8)                             :: norm
    real(8),dimension(Nlat,Norb)        :: nup,ndw,Sz,nt
    complex(8),dimension(:),pointer     :: state_cvec
    type(sector_map)                    :: Hi(2*Ns_Ud)
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb))
    allocate(docc(Nlat,Norb))
    allocate(magz(Nlat,Norb),sz2(Nlat,Norb,Norb))
    allocate(simp(Nlat,Norb,Nspin),zimp(Nlat,Norb,Nspin))
    !
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    simp    = 0.d0
    zimp    = 0.d0
    nup     = 0.d0
    ndw     = 0.d0
    Egs     = state_list%emin
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       idim    = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !norm=sqrt(dot_product(state_cvec,state_cvec))
       !
       !if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then 
          call build_sector(isector,Hi)
          do i=1,idim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = HI(ii)%map(Indices(ii))
                mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb)
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             IbUp = Breorder(Nups)
             IbDw = Breorder(Ndws)
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !Get operators:
             do ilat=1,Nlat
                do iorb=1,Norb
                   nup(ilat,iorb)= ibup(imp_state_index(ilat,iorb))
                   ndw(ilat,iorb)= ibdw(imp_state_index(ilat,iorb))
                   sz(ilat,iorb) = (nup(ilat,iorb) - ndw(ilat,iorb))/2d0
                   nt(ilat,iorb) =  nup(ilat,iorb) + ndw(ilat,iorb)
                enddo
             enddo
             !
             !Evaluate averages of observables:
             do ilat=1,Nlat
                do iorb=1,Norb
                   dens(ilat,iorb)     = dens(ilat,iorb)      +  nt(ilat,iorb)*gs_weight
                   dens_up(ilat,iorb)  = dens_up(ilat,iorb)   +  nup(ilat,iorb)*gs_weight
                   dens_dw(ilat,iorb)  = dens_dw(ilat,iorb)   +  ndw(ilat,iorb)*gs_weight
                   docc(ilat,iorb)     = docc(ilat,iorb)      +  nup(ilat,iorb)*ndw(ilat,iorb)*gs_weight
                   magz(ilat,iorb)     = magz(ilat,iorb)      +  (nup(ilat,iorb)-ndw(ilat,iorb))*gs_weight
                   sz2(ilat,iorb,iorb) = sz2(ilat,iorb,iorb)  +  (sz(ilat,iorb)*sz(ilat,iorb))*gs_weight
                   do jorb=iorb+1,Norb
                      sz2(ilat,iorb,jorb) = sz2(ilat,iorb,jorb)  +  (sz(ilat,iorb)*sz(ilat,jorb))*gs_weight
                      sz2(ilat,jorb,iorb) = sz2(ilat,jorb,iorb)  +  (sz(ilat,jorb)*sz(ilat,iorb))*gs_weight
                   enddo
                enddo
                s2tot = s2tot  + (sum(sz))**2*gs_weight
             enddo
          enddo
          call delete_sector(isector,HI)
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    if(MPIMASTER)then
       call get_szr
       if(iolegend)call write_legend
       call write_observables()
       !
       do ilat=1,Nlat
          write(LOGfile,"(A,10f18.12,f18.12,A)")"dens site "//str(ilat)//" "//reg(file_suffix)//"=",(dens(ilat,iorb),iorb=1,Norb),sum(dens(ilat,:))
       enddo
       write(LOGfile,"(A,10f18.12,f18.12,A)")"dens  avg  "//reg(file_suffix)//" =",(sum(dens(:,iorb))/Nlat,iorb=1,Norb),sum(dens)/Nlat
       !
       write(LOGfile,"(A,10f18.12,A)")"docc       "//reg(file_suffix)//" =",(sum(docc(:,iorb))/Nlat,iorb=1,Norb)
       if(Nspin==2)write(LOGfile,"(A,10f18.12,A)") "mag "//reg(file_suffix)//"=",(sum(magz(:,iorb))/Nlat,iorb=1,Norb)
       !
       imp_dens_up=dens_up
       imp_dens_dw=dens_dw
       imp_dens   =dens
       imp_docc   =docc
    endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,dens_up)
       call Bcast_MPI(MpiComm,dens_dw)
       call Bcast_MPI(MpiComm,dens)
       call Bcast_MPI(MpiComm,docc)
    endif
#endif
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2)
    deallocate(simp,zimp)
  end subroutine vca_lanc_observables


!+---------------------------------------------------------------------------------+
!PURPOSE  : Evaluate and print out custom observable
!+---------------------------------------------------------------------------------+
  
  !T=0
  function calculate_observable_integral_zero_t() result(out_2)
    integer                                                   :: inf
    real(8)                                                   :: out_2,spin_multiplicity
    !
    out_2=0.d0
    spin_multiplicity=3.d0-Nspin
    !
    call quad(sum_observable_kmesh,a=0.0d0,inf=1,verbose=(VERBOSE>=3),result=out_2,strict=.false.)
    !
    out_2=spin_multiplicity*out_2/(pi)
    return
  end function calculate_observable_integral_zero_t

  !T finite

  function calculate_observable_integral_finite_t() result(out_2)
    integer                         :: inf,Nmax,ii
    real(8)                         :: out_2,spin_multiplicity,omegamax,integralpart
    !
    !1) Find the real omegamax
    nmax=int(2*(abs(max_exc)+bandwidth)*beta/pi)
    if (mod(nmax,2)==0)then
      nmax=nmax/2    
    else
      nmax=(nmax+1)/2
    endif
    integrationR=2*(nmax+1)*pi/beta
    !2) Evaluate discrete sum
    !
    out_2=0.d0
    do ii=0,Nmax
      out_2=out_2+dreal(sum_observable_kmesh_complex(xi*(2*ii+1)*pi/beta))
    enddo
    !
    out_2=2.d0*(1/beta)*out_2
    !
    !3) Evaluate integral part
    integralpart=0.d0
    call quad(integral_contour,a=-pi,b=pi,verbose=(VERBOSE>=3),key=6,result=integralpart,strict=.false.)
    !
    !4) Sum all
    out_2=out_2+integralpart
    !5) Spin trick
    spin_multiplicity=3.d0-Nspin 
    out_2=spin_multiplicity*out_2
    return
  end function calculate_observable_integral_finite_t


 function integral_contour(zeta) result(f)
    real(8)                 :: zeta,f
    complex(8)              :: w,fermi
    !
    w=integrationR*exp(xi*zeta)
    if(dreal((w-XMU)*beta)>= 100)then
      fermi=0.d0
    else
      fermi=(1/(exp(beta*(w-XMU))+1))
    endif
    !
    f=dreal((1.d0/pi)*w*fermi*sum_observable_kmesh_complex(w))
 end function integral_contour

  !+-------------------------------------------------------------------+
  !PURPOSE  : sum on k-vectors
  !+-------------------------------------------------------------------+


  function sum_observable_kmesh(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    real(8)                                                  :: omega
    real(8)                                                  :: out_1
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    complex(8),allocatable,dimension(:,:,:,:,:,:)            :: gfprime
    complex(8),allocatable,dimension(:,:)                    :: gfprime_lso
    complex(8),allocatable,dimension(:,:)                    :: Gk_lso
    complex(8),allocatable,dimension(:,:,:,:,:,:)            :: deltamat 
    !
    out_1=0.d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    if(allocated(gfprime))deallocate(gfprime)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(Gk_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(deltamat(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    gfprime_lso=zero
    Gk_lso=zero
    !
    !    
    call vca_gf_cluster(xi*omega,gfprime)
    gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
    call inv(gfprime_lso)
    !
    if(Nbath>0)then
      deltamat=delta_bath_freq(xi*omega,vca_bath)
    endif
    !
    do ii=1,size(impHk,7)
       Gk_lso=vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc-deltamat,Nlat,Nspin,Norb)    
       Gk_lso=gfprime_lso-Gk_lso
       call inv(Gk_lso)
       out_1=out_1+DREAL(trace(matmul(vca_nnn2lso_reshape(sij(:,:,:,:,:,:,ii),Nlat,Nspin,Norb),Gk_lso))-&
         trace(vca_nnn2lso_reshape(sij(:,:,:,:,:,:,ii),Nlat,Nspin,Norb))/(xi*omega-1.1d0))
    enddo
    out_1=out_1/(size(impHk,7))
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    return
    !
  end function sum_observable_kmesh

  function sum_observable_kmesh_complex(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    complex(8)                                               :: omega
    complex(8)                                               :: out_1
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    complex(8),allocatable,dimension(:,:,:,:,:,:)            :: gfprime
    complex(8),allocatable,dimension(:,:)                    :: gfprime_lso
    complex(8),allocatable,dimension(:,:)                    :: Gk_lso
    complex(8),allocatable,dimension(:,:,:,:,:,:)            :: deltamat 
    !
    out_1=0.d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    if(allocated(gfprime))deallocate(gfprime)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(Gk_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(deltamat(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    gfprime_lso=zero
    Gk_lso=zero
    !
    !    
    call vca_gf_cluster(omega,gfprime)
    gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
    call inv(gfprime_lso)
    !
    if(Nbath>0)then
      deltamat=delta_bath_freq(omega,vca_bath)
    endif
    !
    do ii=1,size(impHk,7)
       Gk_lso=vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc-deltamat,Nlat,Nspin,Norb)    
       Gk_lso=gfprime_lso-Gk_lso
       call inv(Gk_lso)
       out_1=out_1+DREAL(trace(matmul(vca_nnn2lso_reshape(sij(:,:,:,:,:,:,ii),Nlat,Nspin,Norb),Gk_lso)))  !!FIXME: CHECK
    enddo
    out_1=out_1/(size(impHk,7))
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    return
    !
  end function sum_observable_kmesh_complex


  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ilat,ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             simp(ilat,iorb,ispin) = dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1)) - &
                  wm1*(dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,2))-&
                  dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
             zimp(ilat,iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1))/wm1 ))
          enddo
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.vca")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"egs",&
         ((reg(txtfy(5*Norb+1+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((6+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy((7+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.vca")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend


  subroutine write_custom_legend()
    integer :: unit,i
    unit = free_unit()
    open(unit,file="custom_observables_info.vca")
    write(unit,"(A1,90(A10,6X))")"#",(reg(txtfy(i))//reg(custom_o%item(i)%o_name),i=1,custom_o%N_filled)
    close(unit)
  end subroutine write_custom_legend

  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: ilat,iorb,jorb,ispin
    do ilat=1,Nlat
      unit = free_unit()     
      open(unit,file="parameters_last"//reg(file_suffix)//"_site"//str(ilat,3)//".vca")
      write(unit,"(90F15.9)")xmu,beta,(uloc_per_site(ilat,iorb),iorb=1,Norb),Ust_per_site(ilat),Jh_per_site(ilat),Jx_per_site(ilat),Jp_per_site(ilat)
      close(unit)
    enddo
    !
    do ilat=1,Nlat
       unit = free_unit()
       open(unit,file="observables_last"//reg(file_suffix)//"_site"//str(ilat,3)//".vca")
       write(unit,"(90(F15.9,1X))")&
            (dens(ilat,iorb),iorb=1,Norb),&
            (docc(ilat,iorb),iorb=1,Norb),&
            (dens_up(ilat,iorb),iorb=1,Norb),&
            (dens_dw(ilat,iorb),iorb=1,Norb),&
            (magz(ilat,iorb),iorb=1,Norb),&
            egs,&
            ((sz2(ilat,iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
    enddo
  end subroutine write_observables

  subroutine write_custom_observables()
    integer :: i
    integer :: unit
    unit = free_unit()
    open(unit,file="custom_observables_all.vca",position='append')
    write(unit,"(90(F15.9,1X))")&
         (custom_o%item(i)%o_value,i=1,custom_o%N_filled)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="custom_observables_last.vca")
    write(unit,"(90(F15.9,1X))")&
         (custom_o%item(i)%o_value,i=1,custom_o%N_filled)
    close(unit)
    !
  end subroutine write_custom_observables

END MODULE VCA_OBSERVABLES









