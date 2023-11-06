!NORMAL, REAL SELF-ENERGY

  subroutine vca_get_sigma_real_full(Sreal)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal = impSreal
  end subroutine vca_get_sigma_real_full

  subroutine vca_get_sigma_real_ij(Sreal,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal = impSreal(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_sigma_real_ij


