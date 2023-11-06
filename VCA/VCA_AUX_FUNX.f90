MODULE VCA_AUX_FUNX
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_BATH_SETUP
  USE SF_TIMER
  USE SF_LINALG
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg
  implicit none
  private



  interface vca_lso2nnn_reshape
     module procedure :: d_nlso2nnn_scalar
     module procedure :: c_nlso2nnn_scalar
  end interface vca_lso2nnn_reshape


  interface vca_nnn2lso_reshape
     module procedure :: d_nnn2nlso_scalar
     module procedure :: c_nnn2nlso_scalar
  end interface vca_nnn2lso_reshape

  interface vca_so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface vca_so2nn_reshape

  interface vca_nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface vca_nn2so_reshape

  interface vca_set_Hcluster
     module procedure :: set_Hcluster_nnn
     module procedure :: set_Hcluster_lso
  end interface vca_set_Hcluster

  interface vca_set_Hk
     module procedure :: set_Hk_nnn
     module procedure :: set_Hk_lso
  end interface vca_set_Hk

  interface vca_print_Hcluster
     module procedure :: print_Hcluster_nnn
     module procedure :: print_Hcluster_lso
  end interface vca_print_Hcluster

#if __GNUC__ > 6
  interface read(unformatted)
    procedure read_unformatted
  end interface read(unformatted)
  
  interface write(unformatted)
    procedure write_unformatted
  end interface write(unformatted)

 interface read(formatted)
    procedure read_formatted
  end interface read(formatted)
  
  interface write(formatted)
    procedure write_formatted
  end interface write(formatted)
#endif

  public :: vca_get_cluster_dimension
  !
  public :: vca_set_Hcluster
  public :: vca_set_Hk
  public :: vca_print_Hcluster
  !
  public :: vca_lso2nnn_reshape
  public :: vca_nnn2lso_reshape
  public :: vca_so2nn_reshape
  public :: vca_nn2so_reshape
  !
  public :: save_gfprime
  public :: read_gfprime
  !
  public :: vca_search_variable
  public :: search_chemical_potential
  !
  public :: print_embedded_H_lso



contains

  !##################################################################
  !                   DIMENSION PROCEDURES
  !##################################################################
  function vca_get_cluster_dimension(with_bath) result(Ncluster)
    logical,optional :: with_bath
    logical          :: bool
    integer          :: Ns
    integer          :: Ncluster
    !
    bool = .false. ; if(present(with_bath))bool = with_bath
    !
    !Count how many levels are there in the cluster:
    Ns = Nlat*Norb
    if(bool)then
       ! select case(bath_type)
       ! case default
       Ns = (Nbath+1)*Nlat*Norb !Norb per site plus Nbath per orb per site
       ! case ('hybrid')
       !    Ns = Nbath+Nlat*Norb     !Norb per site plus shared Nbath sites FIXME: MAYBE ADD HYBRID
       ! end select
    endif
    !
    !Count the spin:
    Ncluster  = Nspin*Ns
    !
  end function vca_get_cluster_dimension




  !##################################################################
  !                   HCLUSTER ROUTINES
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_Hcluster_nnn(hloc,file)
    complex(8),dimension(:,:,:,:,:,:) :: hloc
    character(len=*),optional         :: file
    integer                           :: ilat,jlat
    integer                           :: iorb,jorb
    integer                           :: ispin,jspin
    integer                           :: unit
    character(len=32)                 :: fmt
    !
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"print_Hcluster_nnn","Hloc")
    !
    unit=LOGfile;
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc to file :"//reg(file)
    endif
    write(fmt,"(A,I0,A)")"(",Nlat*Nspin*Norb,"F9.2)"
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,fmt)(((dreal(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)),jlat=1,Nlat),jspin=1,Nspin),jorb=1,Norb)
          enddo
       enddo
    enddo
    write(unit,*)""
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,fmt)(((dimag(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)),jlat=1,Nlat),jspin=1,Nspin),jorb=1,Norb)
          enddo
       enddo
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hcluster_nnn
  !
  subroutine print_Hcluster_lso(hloc,file)
    complex(8),dimension(:,:) :: hloc
    character(len=*),optional :: file
    integer                   :: unit,is,js
    character(len=32)         :: fmt
    !
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"print_Hcluster_lso","Hloc")
    !
    unit=LOGfile;
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc to file :"//reg(file)
    endif
    write(fmt,"(A,I0,A)")"(",Nlat*Nspin*Norb,"A)"
    do is=1,Nlat*Nspin*Norb
       write(unit,fmt)(str(Hloc(is,js),3)//" ",js=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hcluster_lso

  subroutine print_embedded_H_lso(hloc,file)
    complex(8),dimension(:,:) :: hloc
    character(len=*),optional :: file
    integer                   :: unit,is,js
    character(len=32)         :: fmt
    !
    !
    unit=LOGfile;
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc to file :"//reg(file)
    endif
    write(fmt,"(A,I0,A)")"(",Nlat*Nspin*Norb*(Nbath+1),"A)"
    write(unit,"(A)")"REAL"
    do is=1,Nlat*Nspin*Norb*(Nbath+1)
       write(unit,fmt)(str(DIMAG(hloc(is,js)),5)//" ",js=1,Nlat*Nspin*Norb*(Nbath+1))
    enddo
    write(unit,"(A)")"IMAG"
    do is=1,Nlat*Nspin*Norb*(Nbath+1)
       write(unit,fmt)(str(DREAL(hloc(is,js)),5)//" ",js=1,Nlat*Nspin*Norb*(Nbath+1))
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_embedded_H_lso


  !+------------------------------------------------------------------+
  !PURPOSE  : allocate code-internal hopping matrices
  !+------------------------------------------------------------------+
  subroutine set_Hcluster_nnn(hloc)
    complex(8),dimension(:,:,:,:,:,:) :: Hloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"set_Hcluster_nnn","Hloc")
    !
    impHloc = Hloc
    !
    write(LOGfile,"(A)")"Set Hcluster: done"
    if(verbose>2)call vca_print_Hcluster(impHloc)
  end subroutine set_Hcluster_nnn

  subroutine set_Hcluster_lso(hloc)
    complex(8),dimension(:,:) :: Hloc ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"set_Hcluster_lso","Hloc")
    !
    impHloc = vca_lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
    !
    write(LOGfile,"(A)")"Set Hcluster: done"
    if(verbose>2)call vca_print_Hcluster(impHloc)
  end subroutine set_Hcluster_lso


  subroutine set_Hk_nnn(hloc)
    complex(8),dimension(:,:,:,:,:,:,:) :: Hloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Nk]
    !
    impHk = Hloc
    !
    write(LOGfile,"(A)")"Set Hk: done"
    !if(verbose>2)call vca_print_Hcluster(impHk)
  end subroutine set_Hk_nnn

  subroutine set_Hk_lso(hloc)
    complex(8),dimension(:,:,:) :: Hloc ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    integer                     :: i
    !
    do i=1,size(Hloc,3)
        impHk(:,:,:,:,:,:,i) = vca_lso2nnn_reshape(Hloc(:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    write(LOGfile,"(A)")"Set Hk: done"
    !if(verbose>2)call vca_print_Hcluster(impHk)
  end subroutine set_Hk_lso


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from/to the
  ! [Nlso][Nlso] and
  ! [Nlat][Nlat][Norb][Norb][Nspin][Nspin]
  ! shapes.
  ! - dble & cmplx
  ! 0-Reshape a scalar array, dble & cmplx
  ! 1-Reshape a function array (:,:,:,:,:,:,1:L)
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn_scalar(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn_scalar
  !
  function d_nnn2nlso_scalar(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso_scalar
  !
  function c_nlso2nnn_scalar(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn_scalar
  !
  function c_nnn2nlso_scalar(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso_scalar

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn

 function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso

#if __GNUC__ > 6

  !##################################################################
  !##################################################################
  ! ROUTINES TO READ AND WRITE CLUSTER GREEN FUNCTION
  ! unformatted and formatted I/O
  !##################################################################
  !##################################################################


!+-------------------------------------------------------------------+
!PURPOSE  : write overload for GFmatrix type (formatted)
!+-------------------------------------------------------------------+

subroutine write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
    class(GFmatrix), intent(in)         :: dtv
    integer, intent(in)                 :: unit
    integer, intent(out)                :: iostat
    character(*), intent(in)            :: iotype
    integer, intent(in)                 :: v_list(:)
    integer                             :: Nexc,iexc,Ichan,ilat,jlat,iorb,ispin,istate
    integer                             :: Nchan,Nstates
    character(*), intent(inout)         :: iomsg
    !
    !
    Nstates = size(dtv%state)
    write (unit, *,IOSTAT=iostat, IOMSG=iomsg) Nstates
    do istate=1,Nstates
      Nchan = size(dtv%state(istate)%channel)
      write (unit, *,IOSTAT=iostat, IOMSG=iomsg) Nchan
      do ichan=1,Nchan
        write (unit, *,IOSTAT=iostat, IOMSG=iomsg) size(dtv%state(istate)%channel(ichan)%poles)
        write (unit, *,IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%poles
        write (unit, *,IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%weight
      enddo
      write (unit, *,IOSTAT=iostat, IOMSG=iomsg) "\n"
    enddo
    !
end subroutine write_formatted

!+-------------------------------------------------------------------+
!PURPOSE  : read overload for GFmatrix type (formatted)
!+-------------------------------------------------------------------+

subroutine read_formatted(dtv, unit,iotype, v_list, iostat, iomsg)
    class(GFmatrix), intent(inout)                :: dtv
    integer, intent(in)                           :: unit
    integer, intent(out)                          :: iostat
    character(*), intent(in)                      :: iotype
    integer, intent(in)                           :: v_list(:)
    character(*), intent(inout)                   :: iomsg
    logical                                       :: alloc
    integer                                       :: ichan,Nchan,Nlanc,istate,Nstates
    !
    read (unit,*,IOSTAT=iostat, IOMSG=iomsg) Nstates
    call GFmatrix_allocate(dtv,Nstate=Nstates)
    do istate=1,Nstates
      read (unit,*,IOSTAT=iostat, IOMSG=iomsg) Nchan
      call GFmatrix_allocate(dtv,istate=istate,Nchan=Nchan)
      do ichan=1,Nchan
        read (unit,*, IOSTAT=iostat, IOMSG=iomsg) Nlanc
        call GFmatrix_allocate(dtv,istate=istate,ichan=ichan,Nexc=Nlanc)
        read (unit, *, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%poles
        read (unit, *, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%weight
      enddo
    enddo
    !
end subroutine read_formatted



!+-------------------------------------------------------------------+
!PURPOSE  : write overload for GFmatrix type (unformatted)
!+-------------------------------------------------------------------+

subroutine write_unformatted(dtv, unit, iostat, iomsg)
    class(GFmatrix), intent(in)         :: dtv
    integer, intent(in)                 :: unit
    integer, intent(out)                :: iostat
    integer                             :: Nexc,iexc,Ichan,ilat,jlat,iorb,ispin,istate
    integer                             :: Nchan,Nstates
    character(*), intent(inout)         :: iomsg
    !
    !
    Nstates = size(dtv%state)
    write (unit, IOSTAT=iostat, IOMSG=iomsg) Nstates
    do istate=1,Nstates
      Nchan = size(dtv%state(istate)%channel)
      write (unit, IOSTAT=iostat, IOMSG=iomsg) Nchan
      do ichan=1,Nchan
        write (unit, IOSTAT=iostat, IOMSG=iomsg) size(dtv%state(istate)%channel(ichan)%poles), dtv%state(istate)%channel(ichan)%poles, dtv%state(istate)%channel(ichan)%weight
      enddo
    enddo
    !
end subroutine write_unformatted

!+-------------------------------------------------------------------+
!PURPOSE  : read overload for GFmatrix type (unformatted)
!+-------------------------------------------------------------------+

subroutine read_unformatted(dtv, unit, iostat, iomsg)
    class(GFmatrix), intent(inout)                :: dtv
    integer, intent(in)                           :: unit
    integer, intent(out)                          :: iostat
    character(*), intent(inout)                   :: iomsg
    logical                                       :: alloc
    integer                                       :: ichan,Nchan,Nlanc,istate,Nstates
    !
    read (unit, IOSTAT=iostat, IOMSG=iomsg) Nstates
    call GFmatrix_allocate(dtv,Nstate=Nstates)
    do istate=1,Nstates
      read (unit, IOSTAT=iostat, IOMSG=iomsg) Nchan
      call GFmatrix_allocate(dtv,istate=istate,Nchan=Nchan)
      do ichan=1,Nchan
        read (unit, IOSTAT=iostat, IOMSG=iomsg) Nlanc
        call GFmatrix_allocate(dtv,istate=istate,ichan=ichan,Nexc=Nlanc)
        read (unit, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%poles
        read (unit, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%weight
      enddo
    enddo
    !
end subroutine read_unformatted

#endif

!+-------------------------------------------------------------------+
!PURPOSE  : Save cluster GF to file
!+-------------------------------------------------------------------+

subroutine save_gfprime(file,used,use_formatted)

  character(len=*),optional :: file
  character(len=256)        :: file_
  logical,optional          :: used
  logical                   :: used_
  logical,optional          :: use_formatted
  logical                   :: use_formatted_
  character(len=16)         :: extension
  integer                   :: unit_,Nchannel,Nexc,ichan,iexc,ilat,jlat,ispin,iorb,jorb
  !
#if __GNUC__ > 6
  if(.not.allocated(impGmatrix))stop "vca_gf_cluster ERROR: impGmatrix not allocated!"
  used_=.false.;if(present(used))used_=used
  use_formatted_=.false.;if(present(use_formatted))use_formatted_=use_formatted
  extension=".restart";if(used_)extension=".used"
  file_=str(str(file)//str(file_suffix)//str(extension))
  unit_=free_unit()
  !
  if(use_formatted_)then
    open(unit_,file=str(file_),access='sequential')
  else
    open(unit_,file=str(file_),form='unformatted',access='sequential')
  endif
  !
  do ilat=1,Nlat
    do jlat=1,Nlat
      do ispin=1,Nspin
        do iorb=1,Norb
          do jorb=1,Norb
            if(use_formatted_)then
              write(unit_,*)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
            else
              write(unit_)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  close(unit_)
#else
  print*,"Rear/write overloading requires Gfortran 6+"
#endif
end subroutine save_gfprime

!+-------------------------------------------------------------------+
!PURPOSE  : Read cluster GF from file
!+-------------------------------------------------------------------+

subroutine read_gfprime(file,used,use_formatted)
  character(len=*),optional :: file
  character(len=256)        :: file_
  logical,optional          :: used
  logical                   :: used_
  logical,optional          :: use_formatted
  logical                   :: use_formatted_
  character(len=16)         :: extension
  integer                   :: unit_,Nchannel,Nexc,ichan,iexc,ilat,jlat,ispin,iorb,jorb
  !
#if __GNUC__ > 6
  if(allocated(impGmatrix))deallocate(impGmatrix)
  allocate(impGmatrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
  used_=.false.;if(present(used))used_=used
  use_formatted_=.false.;if(present(use_formatted))use_formatted_=use_formatted
  extension=".restart";if(used_)extension=".used"
  file_=str(str(file)//str(file_suffix)//str(extension))
  unit_=free_unit()
  !
  if(use_formatted_)then
    open(unit_,file=str(file_),access='sequential')
  else
    open(unit_,file=str(file_),form='unformatted',access='sequential')
  endif
  !
  rewind(unit_)
  !
  do ilat=1,Nlat
    do jlat=1,Nlat
      do ispin=1,Nspin
        do iorb=1,Norb
          do jorb=1,Norb
            if(use_formatted_)then
              read(unit_,*)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
            else
              read(unit_)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  close(unit_)
#else
  print*,"Rear/write overloading requires Gfortran 6+"
#endif
end subroutine read_gfprime


  !##################################################################
  !##################################################################
  ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! can be used to fix any other *var so that  *ntmp == nread
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine vca_search_variable(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8),save          :: chich
    real(8),save          :: nold
    real(8),save          :: var_new
    real(8),save          :: var_old
    real(8)               :: var_sign
    !
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer               :: unit
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    !  
    if(count==1)then
       chich = ndelta        !~0.2
       inquire(file="var_compressibility.restart",EXIST=bool)
       if(bool)then
          open(free_unit(unit),file="var_compressibility.restart")
          read(unit,*)chich
          close(unit)
       endif
       var_old = var
    endif
    !
    ndiff=ntmp-nread
    !
    !Get 'charge compressibility"
    if(count>1)chich = (ntmp-nold)/(var-var_old)
    !
    !Add here controls on chich: not to be too small....
    !
    !update chemical potential
    var_new = var - ndiff/chich
    !
    !
    !re-define variables:
    nold    = ntmp
    var_old = var
    var     = var_new
    !
    !Print information
    write(LOGfile,"(A9,F16.9,A,F15.9)")  "n    = ",ntmp,"| instead of",nread
    write(LOGfile,"(A9,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nerr
    var_sign = (var-var_old)/abs(var-var_old)
    if(var_sign>0d0)then
       write(LOGfile,"(A9,ES16.9,A4)")"shift = ",ndiff/chich," ==>"
    else
       write(LOGfile,"(A9,ES16.9,A4)")"shift = ",ndiff/chich," <=="
    end if
    write(LOGfile,"(A9,F16.9)")"var  = ",var
    !
    !Save info about search variable iteration:
    open(free_unit(unit),file="search_variable_iteration_info"//reg(file_suffix)//".ed",position="append")
    if(count==1)write(unit,*)"#var,ntmp,ndiff"
    write(unit,*)var,ntmp,ndiff
    close(unit)
    !
    !If density is not converged set convergence to .false.
    if(abs(ndiff)>nerr)converged=.false.
    !
    write(LOGfile,"(A18,I5)")"Search var count= ",count
    write(LOGfile,"(A19,L2)")"Converged       = ",converged
    print*,""
    !
    open(free_unit(unit),file="var_compressibility.used")
    write(unit,*)chich
    close(unit)
    !
  end subroutine vca_search_variable





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine search_chemical_potential(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2
    logical,save          :: ireduce=.true.
    integer               :: unit
    !
    ndiff=ntmp-nread
    nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    if(count>2)then
       do i=1,2
          nindex_old(i+1)=nindex_old(i)
       enddo
    endif
    nindex_old(1)=nindex
    !
    if(ndiff >= nth)then
       nindex=-1
    elseif(ndiff <= -nth)then
       nindex=1
    else
       nindex=0
    endif
    !
    ndelta_old=ndelta
    bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
    !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
    if(bool)then
       ndelta=ndelta_old*nratio !decreasing the step
    else
       ndelta=ndelta_old
    endif
    !
    if(ndelta_old<1.d-9)then
       ndelta_old=0.d0
       nindex=0
    endif
    !update chemical potential
    var=var+dble(nindex)*ndelta
    !xmu=xmu+dble(nindex)*ndelta
    !
    !Print information
    write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
    if(nindex>0)then
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," ==>"
    elseif(nindex<0)then
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," <=="
    else
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," == "
    endif
    write(LOGfile,"(A,f15.9)")"var  = ",var
    write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
    unit=free_unit()
    open(unit,file="search_mu_iteration"//reg(file_suffix)//".ed",position="append")
    write(unit,*)var,ntmp,ndiff
    close(unit)
    !
    !check convergence within actual threshold
    !if reduce is activetd
    !if density is in the actual threshold
    !if DMFT is converged
    !if threshold is larger than nerror (i.e. this is not last loop)
    bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
    if(bool)then
       nth_magnitude_old=nth_magnitude        !save old threshold magnitude
       nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
       nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
       count=0                                !reset the counter
       converged=.false.                      !reset convergence
       ndelta=ndelta_old*nratio                  !reduce the delta step
       !
    endif
    !
    !if density is not converged set convergence to .false.
    if(abs(ntmp-nread)>nth)converged=.false.
    !
    !check convergence for this threshold
    !!---if smallest threshold-- NO MORE
    !if reduce is active (you reduced the treshold at least once)
    !if # iterations > max number
    !if not yet converged
    !set threshold back to the previous larger one.
    !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
    bool=ireduce.AND.(count>niter).AND.(.not.converged)
    if(bool)then
       ireduce=.false.
       nth=10.d0**(nth_magnitude_old)
    endif
    !
    write(LOGfile,"(A,I5)")"count= ",count
    write(LOGfile,"(A,L2)")"Converged=",converged
    print*,""
    !
  end subroutine search_chemical_potential















END MODULE VCA_AUX_FUNX




