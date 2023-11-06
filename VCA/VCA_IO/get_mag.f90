  subroutine vca_get_mag(mag)
    real(8),dimension(Nlat,Norb) :: mag
    mag = imp_dens_up - imp_dens_dw
  end subroutine vca_get_mag
