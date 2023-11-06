!NORMAL, MATSUBARA SELF-ENEGRGY
  subroutine vca_get_sigma_matsubara_full(Smats)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats = impSmats
  end subroutine vca_get_sigma_matsubara_full

  subroutine vca_get_sigma_matsubara_ij(Smats,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats = impSmats(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_sigma_matsubara_ij
