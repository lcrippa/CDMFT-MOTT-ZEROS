!NORMAL, REAL GREEN'S FUNCTION

  subroutine vca_get_gimp_real_full(Greal)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal = impGreal
  end subroutine vca_get_gimp_real_full

  subroutine vca_get_gimp_real_ij(Greal,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal = impGreal(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_gimp_real_ij

