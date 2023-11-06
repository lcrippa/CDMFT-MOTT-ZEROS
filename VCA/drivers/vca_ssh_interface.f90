program vca_ssh
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  implicit none
  integer                                         :: Nlso,Ndimer,Nk
  real(8)                                         :: vhop,whop,eps0,energy_offset,wmixing,eps_shift,eps_left
  logical                                         :: converged
  real(8)                                         :: wband
  !
  logical                                         :: wloop,wmin
  integer                                         :: nloop,iloop
  real(8)                                         :: vhop_var,whop_var,dummy_omega
  real(8),dimension(:),allocatable                :: vhop_array,whop_array,omega_array
  integer,dimension(1)                            :: min_loc
  !
  real(8),dimension(:,:),allocatable              :: Tsys,Tref,Vmat,Htb,Mmat,dens
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8)                                      :: iw
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: h_cluster
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: h_k
  character(len=16)                               :: finput
  !
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master
  !
  real(8),dimension(:,:),allocatable              :: U_loc_per_site
  real(8)                                         :: U2
  real(8)                                         :: interface_coupling
  integer                                         :: interface_offset
  integer                                         :: ilat,iorb
  !

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(vhop,"vhop",finput,default=0.25d0,comment="intra-dimer hopping")
  call parse_input_variable(whop,"whop",finput,default=0.25d0,comment="inter-dimer hopping")
  call parse_input_variable(eps0,"eps0",finput,default=0.00d0,comment="difference in on-site energy")
  call parse_input_variable(eps_left,"eps_left",finput,default=0.00d0,comment="energy shift on left side of the chain")
  call parse_input_variable(eps_shift,"eps_shift",finput,default=0.00d0,comment="shift of on-site energy for all sites")
  call parse_input_variable(Ndimer,"Ndimer",finput,default=1,comment="number of dimers")
  call parse_input_variable(Nk,"Nk",finput,default=10,comment="Number of k point for BZ integration")
  call parse_input_variable(wloop,"wloop",finput,default=.false.,comment="T: includes loop over ts")
  call parse_input_variable(wmin,"wmin",finput,default=.false.,comment="T: includes global minimization")
  call parse_input_variable(nloop,"NLOOP",finput,default=5)
  call parse_input_variable(interface_coupling,"INTERFACE_COUPLING",finput,default=1.0d0)
  call parse_input_variable(interface_offset,"INTERFACE_OFFSET",finput,default=0)
  call parse_input_variable(U2,"U2",finput,default=0.0d0)
  !
  call vca_read_input(trim(finput),comm)
  !
  !Add DMFTtools CTRL Variables:
  !
  call add_ctrl_var(Nlat,"NLAT")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !if (modulo(Ndimer,2) /= 0 ) then
  ! print *, "Ndimer has to be even!"
  ! stop
  !endif
  Nlat=Ndimer*2
  Nlso=Nlat*Nspin*Norb
  !
  vhop_var=vhop
  whop_var=whop
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  !Set U_loc for each site
  allocate(U_loc_per_site(Nlat,Norb))
  U_loc_per_site = zero
  do ilat=1,Nlat
    do iorb = 1, Norb
      if (ilat<= Nlat/2+interface_offset) then
         U_loc_per_site(ilat,iorb) = Uloc(iorb)
      else
         U_loc_per_site(ilat,iorb) = U2
      endif
    enddo
  enddo
  print *, U_loc_per_site
  if (abs(interface_offset) > 1)then
   print *, "interface_offset larger than 1 not implemented!"
   stop
  endif
  !
  !INIT SOLVER:
  !
  
  call vca_init_solver(comm)
  !print_observables=.false.
  !
   dummy_omega = solve_vca(vhop,whop)
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  
  call finalize_MPI()

!+------------------------------------------------------------------+
!FUNCTIONS
!+------------------------------------------------------------------+

contains


  function solve_vca(vij,wij) result(Omega)
    real(8)                      :: vij,wij
    real(8)                      :: Omega
    integer ::ik
    !
    vhop_var=vij
    whop_var=wij
    !print*,""
    !print*,"----- INIT -----"
    !
    !print*,"VHOP_VAR = ",vhop_var
    h_cluster = lso2nnn(Hloc_model(Nlso,vhop_var,whop_var,eps0),Nlat,Nspin,Norb)
    call generate_hk()
    print *, h_cluster
    call vca_solve(comm,h_cluster,h_k,Uloc_ii = U_loc_per_site)
    call vca_get_sft_potential(omega)
    !
    !print*,"------ DONE ------"
    !print*,""
    !
  end function solve_vca

   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d BHZ model
   !+------------------------------------------------------------------+


   function Hloc_model(N,vhop_,whop_,eps0_) result (H0)
      integer                                               :: N,ispin,idimer
      real(8)                                               :: vhop_,whop_,eps0_
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: H0
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do idimer=1,Ndimer
          hopping_matrix(2*idimer-1,2*idimer  ,ispin,ispin,:,:) = -vhop_
          hopping_matrix(2*idimer  ,2*idimer-1,ispin,ispin,:,:) = -vhop_
          !
          if(idimer < Ndimer) then
             hopping_matrix(2*idimer  ,2*idimer+1,ispin,ispin,:,:) = -whop_
             hopping_matrix(2*idimer+1,2*idimer  ,ispin,ispin,:,:) = -whop_
          endif
          !
          if(idimer > Ndimer) then
             hopping_matrix(2*idimer-1,2*idimer-2,ispin,ispin,:,:) = -whop_
             hopping_matrix(2*idimer-2,2*idimer-1,ispin,ispin,:,:) = -whop_
          endif
          ! adjust the coupling at the interface
          if (idimer == Ndimer/2+interface_offset) then
            hopping_matrix(2*idimer -interface_offset ,2*idimer+1-interface_offset,ispin,ispin,:,:) = -whop_*interface_coupling
            hopping_matrix(2*idimer+1-interface_offset,2*idimer  -interface_offset,ispin,ispin,:,:) = -whop_*interface_coupling
          endif
          !set on-site energy terms for the interacting half
          if (idimer <= Ndimer/2+interface_offset) then
            hopping_matrix(2*idimer -1 ,2*idimer-1,ispin,ispin,:,:) = eps_shift + eps0_ -Uloc(1)/2 +eps_left
            hopping_matrix(2*idimer,2*idimer,ispin,ispin,:,:) = eps_shift + eps0_ -Uloc(1)/2 +eps_left
          else
            hopping_matrix(2*idimer -1 ,2*idimer-1,ispin,ispin,:,:) =eps_shift -eps0_ -U2/2
            hopping_matrix(2*idimer,2*idimer,ispin,ispin,:,:) =eps_shift -eps0_ -U2/2
          endif
         enddo

      enddo
      !
      H0=nnn2lso(hopping_matrix,Nlat,Nspin,Norb)
      !
   end function hloc_model


   function hk_model(kpoint,N) result(Hk)
      integer                                                      :: Nispin,ispin,N
      real(8),dimension(:)                                         :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
      complex(8),dimension(N,N)                                    :: hk
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
            hopping_matrix(1,Nlat,ispin,ispin,:,:) = hopping_matrix(1,Nlat,ispin,ispin,:,:) + whop*exp(-xi*kpoint(1)*Ndimer)
            hopping_matrix(Nlat,1,ispin,ispin,:,:) = hopping_matrix(Nlat,1,ispin,ispin,:,:) + whop*exp( xi*kpoint(1)*Ndimer)
      enddo
      !
      Hk=nnn2lso(hopping_matrix,Nlat,Nspin,Norb)+hloc_model(N,vhop,whop,eps0)
      !
   end function hk_model

   subroutine generate_hk()
      integer                                     :: ik
      real(8),dimension(Nk,1)                     :: kgrid
      real(8),dimension(2)                        :: e1,e2,bk1,bk2
      real(8)                                     :: bklen
      complex(8),allocatable                      :: Hk_lso(:,:,:)
      !
      e1 = [1d0, 0d0]
      call TB_set_ei(eix=e1)
      bklen=2d0*pi
      bk1=bklen*[1d0, 0d0]
      call TB_set_bk(bkx=bk1)
      !
      call TB_build_kgrid([Nk],kgrid)
      kgrid(:,1)=kgrid(:,1)/Ndimer
      !
      if(allocated(hk_lso))deallocate(hk_lso)
      !
      allocate(Hk_lso(Nlso,Nlso,Nk))

      if(allocated(h_k))deallocate(h_k)
      allocate(h_k(Nlat,Nlat ,NSpin,NSpin,Norb,Norb,Nk))
      hk_lso=zero
      !
      !call TB_build_model(Hk_lso,hk_model,Nlso,kgrid)
      !
      do ik=1,Nk
        h_k(:,:,:,:,:,:,ik)=lso2nnn(hk_lso(:,:,ik),Nlat,Nspin,Norb)
      enddo
      !
   end subroutine generate_hk

   !+------------------------------------------------------------------+
   !PURPOSE  : Auxilliary reshape functions
   !+------------------------------------------------------------------+

   function lso2nnn(Hlso,Nlat_,Nspin_,Norb_) result(Hnnn)
      integer                                                     :: ilat,jlat,Nlat_,Nspin_,Norb_
      complex(8),dimension(Nlat_*Nspin_*Norb_,Nlat_*Nspin_*Norb_) :: Hlso
      complex(8),dimension(Nlat_,Nlat_,Nspin_,Nspin_,Norb_,Norb_) :: Hnnn
      integer                                                     :: iorb,jorb
      integer                                                     :: ispin,jspin
      integer                                                     :: is,js
      Hnnn=zero
      do ilat=1,Nlat_
         do jlat=1,Nlat_
            do ispin=1,Nspin_
               do jspin=1,Nspin_
                  do iorb=1,Norb_
                     do jorb=1,Norb_
                        is = iorb + (ilat-1)*Norb_ + (ispin-1)*Norb_*Nlat_
                        js = jorb + (jlat-1)*Norb_ + (jspin-1)*Norb_*Nlat_
                        Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnn


   function nnn2lso(Hnnn,Nlat_,Nspin_,Norb_) result(Hlso)
      integer                                                     :: ilat,jlat,Nlat_,Nspin_,Norb_
      complex(8),dimension(Nlat_,Nlat_,Nspin_,Nspin_,Norb_,Norb_) :: Hnnn
      complex(8),dimension(Nlat_*Nspin_*Norb_,Nlat_*Nspin_*Norb_) :: Hlso
      integer                                                     :: iorb,jorb
      integer                                                     :: ispin,jspin
      integer                                                     :: is,js
      Hlso=zero
      do ilat=1,Nlat_
         do jlat=1,Nlat_
            do ispin=1,Nspin_
               do jspin=1,Nspin_
                  do iorb=1,Norb_
                     do jorb=1,Norb_
                        is = iorb + (ilat-1)*Norb_ + (ispin-1)*Norb_*Nlat_
                        js = jorb + (jlat-1)*Norb_ + (jspin-1)*Norb_*Nlat_
                        Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso

end program vca_ssh



