  subroutine vca_get_dens(dens)
    real(8),dimension(Nlat,Norb) :: dens
    dens = imp_dens
  end subroutine vca_get_dens
