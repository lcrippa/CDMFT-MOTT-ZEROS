  subroutine vca_get_gimp_matsubara_full(Gmats)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats = impGmats
  end subroutine vca_get_gimp_matsubara_full

  subroutine vca_get_gimp_matsubara_ij(Gmats,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats = impGmats(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_gimp_matsubara_ij
