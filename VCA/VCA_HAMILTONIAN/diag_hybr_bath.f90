  
  if(allocated(diag_hybr))deallocate(diag_hybr)
  if(allocated(bath_diag))deallocate(bath_diag)
  !select case (bath_type)
  !case default
     Nfoo = size(vca_bath%e,3)
     Nfoo2= size(vca_bath%e,1)  
     allocate(diag_hybr(Nlat,Nspin,Norb,Nbath));diag_hybr=0d0
     allocate(bath_diag(Nfoo2,Nspin,Nfoo,Nbath));bath_diag=0d0
     do ibath=1,Nbath
       do ispin=1,Nspin
         do ilat=1,Nlat             
           do iorb=1,Norb
             diag_hybr(ilat,ispin,iorb,ibath)=vca_bath%v(ilat,ispin,iorb,ibath)
           enddo
         enddo
         do ilat=1,Nfoo2  
           do iorb=1,Nfoo
             bath_diag(ilat,ispin,iorb,ibath)=vca_bath%e(ilat,ispin,iorb,ibath)
           enddo
         enddo
       enddo
     enddo
!  case ("replica")
!     Nfoo = Norb
!     allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
!     allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0
!     do ibath=1,Nbath
!        do ispin=1,Nspin
!           do iorb=1,Norb
!              diag_hybr(ispin,iorb,ibath)=vca_bath%vr(ibath)
!              bath_diag(ispin,iorb,ibath)=vca_bath%h(ispin,ispin,iorb,iorb,ibath)
!           enddo
!        enddo
!     enddo
  !end select
