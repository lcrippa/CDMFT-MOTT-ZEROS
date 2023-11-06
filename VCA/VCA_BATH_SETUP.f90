MODULE VCA_BATH_SETUP
   USE SF_CONSTANTS, only: zero
   USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
   USE SF_LINALG, only: eye,inv
   USE SF_MISC, only: assert_shape
   USE VCA_INPUT_VARS
   USE VCA_VARS_GLOBAL
   implicit none


   private



   !##################################################################
   !
   !     USER BATH ROUTINES:
   !
   !##################################################################
   public :: vca_get_bath_dimension




   !##################################################################
   !
   !     VCA BATH ROUTINES:
   !
   !##################################################################
   !VCA BATH procedures:
   public :: vca_allocate_bath               !INTERNAL (for effective_bath)
   public :: vca_deallocate_bath             !INTERNAL (for effective_bath)
   public :: vca_init_bath                   !INTERNAL (for effective_bath)
   public :: vca_write_bath                  !INTERNAL (for effective_bath)
   public :: vca_save_bath                   !INTERNAL (for effective_bath)
   public :: vca_get_bath                    !INTERNAL (for effective_bath)
   public :: vca_set_bath                    !INTERNAL (for effective_bath)
   public :: check_bath_dimension            !INTERNAL (for effective_bath)
   public :: set_bath_component              !INTERNAL (for effective_bath)


contains




   !##################################################################
   !
   !     USER BATH ROUTINES:
   !
   !##################################################################
   !+-------------------------------------------------------------------+
   !PURPOSE  : Inquire the correct bath size to allocate the 
   ! the bath array in the calling program.
   !+-------------------------------------------------------------------+
   function vca_get_bath_dimension(ispin_) result(bath_size)
      integer,optional               :: ispin_
      integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo
      !
      ! select case(bath_type)
      ! case default
      !e:[Nlat][Norb][Nspin][Nbath] + v:[Nlat][Norb][Nspin][Nbath]
      bath_size = Nlat*Norb*Nbath + Nlat*Norb*Nbath
      if(.not.present(ispin_))bath_size=Nspin*bath_size
      ! case('hybrid')
      !    !e:[1][1][Nspin][Nbath] + v:[Nlat][Norb][Nspin][Nbath]
      !    bath_size = Nbath + Nlat*Norb*Nbath
      !    if(.not.present(ispin_))bath_size=Nspin*bath_size
      ! end select
   end function vca_get_bath_dimension










   !##################################################################
   !
   !     VCA BATH ROUTINES:
   !
   !##################################################################
   !+-------------------------------------------------------------------+
   !PURPOSE  : Allocate the VCA bath
   !+-------------------------------------------------------------------+
   subroutine vca_allocate_bath(vca_bath_)
      type(effective_bath) :: vca_bath_
      !
      if(vca_bath_%status)call vca_deallocate_bath(vca_bath_)
      !
      if(Nbath==0)stop "Allocate_vca_bath ERROR: Nbath==0"
      !
      ! select case(bath_type)
      ! case default!normal_normal
      !
      allocate(vca_bath_%e(Nlat,Nspin,Norb,Nbath))  !local energies of the bath per site,orb
      allocate(vca_bath_%v(Nlat,Nspin,Norb,Nbath))  !same-spin hybridization per site,orb
      !
      ! case('hybrid')                            !hybrid_normal
      !    !
      !    allocate(vca_bath_%e(1,1,Nspin,Nbath))        !local energies of the bath stand-alone
      !    allocate(vca_bath_%v(Nlat,Nspin,Norb,Nbath))  !same-spin hybridization, connects site,orb to bath sites
      !    !
      !    ! case('normal_hybrid')
      !    !    allocate(vca_bath_%e(Nlat,1,Nspin,Nbath))        !local energies of the bath stand-alone
      !    !    allocate(vca_bath_%v(Nlat,Nspin,Norb,Nbath))  !same-spin hybridization, connects site,orb to bath sites
      !    !
      ! end select
      vca_bath_%e=0d0
      vca_bath_%v=0d0
      vca_bath_%status=.true.
   end subroutine vca_allocate_bath





   !+-------------------------------------------------------------------+
   !PURPOSE  : Deallocate the VCA bath
   !+-------------------------------------------------------------------+
   subroutine vca_deallocate_bath(vca_bath_)
      type(effective_bath) :: vca_bath_
      if(allocated(vca_bath_%e))   deallocate(vca_bath_%e)
      if(allocated(vca_bath_%v))   deallocate(vca_bath_%v)
      vca_bath_%status=.false.
   end subroutine vca_deallocate_bath





   !+------------------------------------------------------------------+
   !PURPOSE  : Initialize the VCA loop, builindg H parameters and/or 
   !reading previous (converged) solution
   !+------------------------------------------------------------------+
   subroutine vca_init_bath(vca_bath_)
      type(effective_bath) :: vca_bath_
      integer              :: i,unit,flen,Nh
      integer              :: io,jo,ispin,iorb,jorb,jspin,ilat
      logical              :: IOfile
      real(8)              :: de,hwband
      character(len=21)    :: space
      real,dimension(Nbath):: ran
      !
      hwband=1d0
      !
      if(.not.vca_bath_%status)stop "vca_init_bath error: bath not allocated"
      !
      if(Nbath==0)stop "VCA_allocate_bath ERROR: Nbath==0"
      !
      ! !Get energies:
      ! if(Nbath==1)then
      !    vca_bath_%e(:,:,:,Nbath)= 0d0
      ! else
      !    vca_bath_%e(:,:,:,1)    =-hwband
      !    vca_bath_%e(:,:,:,Nbath)= hwband
      ! endif
      ! Nh=Nbath/2
      ! if(mod(Nbath,2)==0.and.Nbath>=4)then
      !    de=hwband/max(Nh-1,1)
      !    vca_bath_%e(:,:,:,Nh)  = -1.d-3
      !    vca_bath_%e(:,:,:,Nh+1)=  1.d-3
      !    do i=2,Nh-1
      !       vca_bath_%e(:,:,:,i)        =-hwband + (i-1)*de
      !       vca_bath_%e(:,:,:,Nbath-i+1)= hwband - (i-1)*de
      !    enddo
      ! elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
      !    de=hwband/Nh
      !    vca_bath_%e(:,:,:,Nh+1)= 0.0d0
      !    do i=2,Nh
      !       vca_bath_%e(:,:,:,i)        =-hwband + (i-1)*de
      !       vca_bath_%e(:,:,:,Nbath-i+1)= hwband - (i-1)*de
      !    enddo
      ! endif
      ! ! call random_number(ran)
      ! ! forall(ilat=1:size(vca_bath_%e,1),iorb=1:size(vca_bath_%e,2),ispin=1:Nspin)&
      ! !      vca_bath_%e(ilat,ispin,iorb,:) = ran !impHloc(ilat,ilat,ispin,iorb,ispin,iorb) + ran/10d0
      ! !
      ! !Get spin-keep yhbridizations
      ! do i=1,Nbath
      !    vca_bath_%v(:,:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))
      ! enddo
      !
      !
      vca_bath_%e = 1d0
      vca_bath_%v = 0d0
      !
      !
      !Read from file if exist:
      !
      inquire(file="bath.restart",exist=IOfile)
      if(IOfile)then
         write(LOGfile,"(A)")"Reading bath from file bath.restart"
         unit = free_unit()
         !
         open(unit,file="bath.restart")
         !
         read(unit,*)          !read the header:
         ! select case(bath_type)
         ! case default
         do ilat=1,Nlat
            do i=1,Nbath
               read(unit,*)((&
                  vca_bath_%e(ilat,ispin,iorb,i),&
                  vca_bath_%v(ilat,ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
            enddo
         enddo
         !
         ! case ('hybrid')
         !    do ilat=1,Nlat
         !       do i=1,Nbath
         !          read(unit,*)(&
         !               vca_bath_%e(1,1,ispin,i),&
         !               (&
         !               vca_bath_%v(ilat,ispin,iorb,i),&
         !               iorb=1,Norb),ispin=1,Nspin)
         !       enddo
         !    enddo
         ! end select
         close(unit)
      endif
   end subroutine vca_init_bath




   !+-------------------------------------------------------------------+
   !PURPOSE  : write out the bath to a given unit with 
   ! the following column formatting: 
   ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
   !+-------------------------------------------------------------------+
   subroutine vca_write_bath(vca_bath_,unit)
      type(effective_bath) :: vca_bath_
      integer,optional     :: unit
      integer              :: unit_
      integer              :: i
      character(len=7)     :: sitestring
      integer              :: io,jo,ispin,iorb,ilat
      complex(8)           :: hybr_aux
      complex(8)           :: hrep_aux(Nspin*Norb,Nspin*Norb)
      !
      unit_=LOGfile;if(present(unit))unit_=unit
      !
      if(.not.vca_bath_%status)stop "write_vca_bath error: bath not allocated"
      !
      ! select case(bath_type)
      ! case default
      if(unit_==LOGfile) then
         write(unit_, "(A7)") ""
         write(unit_, "(A7)", advance="no") ""
         write(unit_,"(90(A12,1X))")((&
            "Ek_l"//str(iorb)//"_s"//str(ispin),"Vk_l"//str(iorb)//"_s"//str(ispin),&
            iorb=1,Norb),ispin=1,Nspin)
         do ilat=1,Nlat
            sitestring="Site "//str(ilat,2)
            do i=1,Nbath
               write(unit_, "(A7)", advance="no") sitestring
               write(unit_,"(90(F12.4,1X))")((&
                  vca_bath_%e(ilat,ispin,iorb,i),&
                  vca_bath_%v(ilat,ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
               sitestring=""
            enddo
         enddo
         write(unit_, "(A7)") ""
      else
         write(unit_,"(90(A21,1X))")((&
            "#Ek_l"//str(iorb)//"_s"//str(ispin),"Vk_l"//str(iorb)//"_s"//str(ispin),&
            iorb=1,Norb),ispin=1,Nspin)
         do ilat=1,Nlat
            do i=1,Nbath
               write(unit_,"(90(F21.12,1X))")((&
                  vca_bath_%e(ilat,ispin,iorb,i),&
                  vca_bath_%v(ilat,ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
            enddo
         enddo

      endif
      !
      ! case('hybrid')
      !    !
      !    write(unit_,"(90(A21,1X))")(&
      !         "#Ek_s"//reg(txtfy(ispin)),("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
      !         iorb=1,Norb),ispin=1,Nspin)
      !    do ilat=1,Nlat
      !       do i=1,Nbath
      !          write(unit_,"(90(F21.12,1X))")(&
      !               vca_bath_%e(ilat,1,ispin,i),&
      !               (vca_bath_%v(ilat,ispin,iorb,i),&
      !               iorb=1,Norb),ispin=1,Nspin)
      !       enddo
      !    enddo
      ! end select
   end subroutine vca_write_bath






   !+-------------------------------------------------------------------+
   !PURPOSE  : save the bath to a given file using the write bath
   ! procedure and formatting: 
   !+-------------------------------------------------------------------+
   subroutine vca_save_bath(vca_bath_,file,used)
      type(effective_bath)      :: vca_bath_
      character(len=*),optional :: file
      character(len=256)        :: file_
      logical,optional          :: used
      logical                   :: used_
      character(len=16)         :: extension
      integer                   :: unit_
      !
      if(.not.vca_bath_%status)stop "save_vca_bath error: bath is not allocated"
      !
      used_=.false.;if(present(used))used_=used
      extension=".restart";if(used_)extension=".used"
      !
      file_="bath"//reg(extension)
      if(present(file))file_=reg(file)
      unit_=free_unit()
      open(unit_,file=reg(file_))
      call vca_write_bath(vca_bath_,unit_)
      close(unit_)
   end subroutine vca_save_bath




   !+-------------------------------------------------------------------+
   !PURPOSE  : set the bath components from a given user provided 
   ! bath-array 
   !+-------------------------------------------------------------------+
   subroutine vca_set_bath(bath_,vca_bath_)
      real(8),dimension(:)   :: bath_
      type(effective_bath)   :: vca_bath_
      integer                :: stride,io,jo,i
      integer                :: ispin,iorb,jorb,jspin,ibath,ilat
      logical                :: check
      complex(8)             :: hrep_aux(Nspin*Norb,Nspin*Norb)
      complex(8)             :: U(Nspin*Norb,Nspin*Norb)
      complex(8)             :: Udag(Nspin*Norb,Nspin*Norb)
      real(8)                :: element_R,element_I,eps_k,lambda_k
      !
      if(.not.vca_bath_%status)stop "set_vca_bath error: bath not allocated"
      !
      check = check_bath_dimension(bath_)
      if(.not.check)stop "set_vca_bath error: wrong bath dimensions"
      !
      ! select case(bath_type)
      ! case default
      stride = 0
      do ilat=1,Nlat
         do iorb=1,Norb
            do ispin=1,Nspin
               do i=1,Nbath
                  io = stride + i + (iorb-1)*Nbath + (ilat-1)*Nbath*Norb + (ispin-1)*Nbath*Norb*Nlat
                  vca_bath_%e(ilat,ispin,iorb,i) = bath_(io)
               enddo
            enddo
         enddo
      enddo
      stride = Nlat*Nspin*Norb*Nbath
      do ilat=1,Nlat
         do iorb=1,Norb
            do ispin=1,Nspin
               do i=1,Nbath
                  io = stride + i + (iorb-1)*Nbath + (ilat-1)*Nbath*Norb + (ispin-1)*Nbath*Norb*Nlat
                  vca_bath_%v(ilat,ispin,iorb,i) = bath_(io)
               enddo
            enddo
         enddo
      enddo
      !
      if(hfmode)then
         do ilat=1,Nlat
            do iorb=1,Norb
               do ispin=1,Nspin
                  do i=1,Nbath
                     if( .not. hfshift) vca_bath_%e(ilat,ispin,iorb,i) = vca_bath_%e(ilat,ispin,iorb,i)-(0.5d0*Uloc_per_site(ilat,iorb)+0.5d0*Ust_per_site(ilat)*(Norb-1)+0.5d0*(Ust_per_site(ilat)-Jh_per_site(ilat))*(Norb-1))
                  enddo
               enddo
            enddo
         enddo
      endif
      !
      ! case ('hybrid')
      !    stride = 0
      !    do ispin=1,Nspin
      !       do i=1,Nbath
      !          io = stride + i + (ispin-1)*Nbath
      !          vca_bath_%e(1,1,ispin,i) = bath_(io)
      !       enddo
      !    enddo
      !    stride = Nspin*Nbath
      !    do ilat=1,Nlat
      !       do iorb=1,Norb
      !          do ispin=1,Nspin
      !             do i=1,Nbath
      !                io = stride + i + index_stride_los(ilat,ispin,iorb)
      !                vca_bath_%v(ilat,ispin,iorb,i) = bath_(io)
      !             enddo
      !          enddo
      !       enddo
      !    enddo
      ! end select
   end subroutine vca_set_bath


   !+-------------------------------------------------------------------+
   !PURPOSE  : let the user modify bath components
   !+-------------------------------------------------------------------+


   subroutine set_bath_component(bath_,ilat,ispin,iorb,e_component,v_component)
      real(8),dimension(:),allocatable        :: bath_
      integer                                 :: stride,io,jo,i
      integer                                 :: ispin,iorb,jorb,jspin,ibath,ilat
      logical                                 :: check, is_e
      complex(8)                              :: hrep_aux(Nspin*Norb,Nspin*Norb)
      complex(8)                              :: U(Nspin*Norb,Nspin*Norb)
      complex(8)                              :: Udag(Nspin*Norb,Nspin*Norb)
      real(8),dimension(Nbath),optional       :: e_component,v_component
      !
      !
      !
      ! select case(bath_type)
      ! case default
      if(present(e_component))then
         stride = 0
         do i=1,Nbath
            io = stride + i + (iorb-1)*Nbath + (ilat-1)*Nbath*Norb + (ispin-1)*Nbath*Norb*Nlat
            bath_(io) = e_component(i)
         enddo
      endif
      if(present(v_component))then
         stride = Nlat*Nspin*Norb*Nbath
         do i=1,Nbath
            io = stride + i + (iorb-1)*Nbath + (ilat-1)*Nbath*Norb + (ispin-1)*Nbath*Norb*Nlat
            bath_(io) = v_component(i)
         enddo
      endif
      !
      ! case ('hybrid')
      !    stride = 0
      !    do ispin=1,Nspin
      !       do i=1,Nbath
      !          io = stride + i + (ispin-1)*Nbath
      !          vca_bath_%e(1,1,ispin,i) = bath_(io)
      !       enddo
      !    enddo
      !    stride = Nspin*Nbath
      !    do ilat=1,Nlat
      !       do iorb=1,Norb
      !          do ispin=1,Nspin
      !             do i=1,Nbath
      !                io = stride + i + index_stride_los(ilat,ispin,iorb)
      !                vca_bath_%v(ilat,ispin,iorb,i) = bath_(io)
      !             enddo
      !          enddo
      !       enddo
      !    enddo
      ! end select
   end subroutine set_bath_component




   !+-------------------------------------------------------------------+
   !PURPOSE  : copy the bath components back to a 1-dim array 
   !+-------------------------------------------------------------------+
   subroutine vca_get_bath(vca_bath_,bath_)
      type(effective_bath)   :: vca_bath_
      real(8),dimension(:)   :: bath_
      complex(8)             :: hrep_aux(Nspin*Norb,Nspin*Norb)
      complex(8)             :: U(Nspin*Norb,Nspin*Norb)
      complex(8)             :: Udag(Nspin*Norb,Nspin*Norb)
      integer                :: stride,io,jo,i,ilat
      integer                :: ispin,iorb,jorb,jspin,ibath
      logical                :: check
      !
      if(.not.vca_bath_%status)stop "get_vca_bath error: bath not allocated"
      !
      check=check_bath_dimension(bath_)
      if(.not.check)stop "get_vca_bath error: wrong bath dimensions"
      !
      ! select case(bath_type)
      ! case default
      stride = 0
      do ilat=1,Nlat
         do iorb=1,Norb
            do ispin=1,Nspin
               do i=1,Nbath
                  io = stride + i + (iorb-1)*Nbath + (ilat-1)*Nbath*Norb + (ispin-1)*Nbath*Norb*Nlat
                  bath_(io) = vca_bath_%e(ilat,ispin,iorb,i)
               enddo
            enddo
         enddo
      enddo
      stride = Nlat*Nspin*Norb*Nbath
      do ilat=1,Nlat
         do iorb=1,Norb
            do ispin=1,Nspin
               do i=1,Nbath
                  io = stride + i + (iorb-1)*Nbath + (ilat-1)*Nbath*Norb + (ispin-1)*Nbath*Norb*Nlat
                  bath_(io) = vca_bath_%v(ilat,ispin,iorb,i)
               enddo
            enddo
         enddo
      enddo
      !
      ! case ('hybrid')
      !    stride = 0
      !    do ispin=1,Nspin
      !       do i=1,Nbath
      !          io = stride + i + (ispin-1)*Nbath
      !          bath_(io) = vca_bath_%e(1,1,ispin,i)
      !       enddo
      !    enddo
      !    stride = Nspin*Nbath
      !    do ilat=1,Nlat
      !       do iorb=1,Norb
      !          do ispin=1,Nspin
      !             do i=1,Nbath
      !                io = stride + i + index_stride_los(ilat,ispin,iorb)
      !                bath_(io) = vca_bath_%v(ilat,ispin,iorb,i)
      !             enddo
      !          enddo
      !       enddo
      !    enddo
      !
      ! end select
   end subroutine vca_get_bath







   !##################################################################
   !
   !     USER BATH CHECKS:
   !
   !##################################################################
   !+-------------------------------------------------------------------+
   !PURPOSE  : Check if the dimension of the bath array are consistent
   !+-------------------------------------------------------------------+
   function check_bath_dimension(bath_) result(bool)
      real(8),dimension(:)           :: bath_
      integer                        :: Ntrue
      logical                        :: bool
      Ntrue = vca_get_bath_dimension()
      bool  = ( size(bath_) == Ntrue )
   end function check_bath_dimension






END MODULE VCA_BATH_SETUP
