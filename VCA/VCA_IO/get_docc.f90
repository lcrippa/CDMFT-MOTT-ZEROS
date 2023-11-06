  subroutine vca_get_docc(docc)
    real(8),dimension(Nlat,Norb) :: docc
    docc = imp_docc
  end subroutine vca_get_docc
