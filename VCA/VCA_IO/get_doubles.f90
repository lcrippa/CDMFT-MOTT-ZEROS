subroutine vca_get_doubles_(docc)
  real(8),dimension(4) :: docc
  docc = [vca_Dust,vca_Dund,vca_Dse,vca_Dph]
end subroutine vca_get_doubles_

subroutine vca_get_dust_(docc)
  real(8) :: docc
  docc = vca_Dust
end subroutine vca_get_dust_

subroutine vca_get_dund_(docc)
  real(8) :: docc
  docc = vca_Dund
end subroutine vca_get_dund_

subroutine vca_get_dse_(docc)
  real(8) :: docc
  docc = vca_Dse
end subroutine vca_get_dse_

subroutine vca_get_dph_(docc)
  real(8) :: docc
  docc = vca_Dph
end subroutine vca_get_dph_

subroutine vca_get_doubles_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "vca_get_doubles error: required N_sites > evaluated N_sites"
     yii=ddii(:,:)
  endif
end subroutine vca_get_doubles_lattice

subroutine vca_get_dust_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "vca_get_dust error: required N_sites > evaluated N_sites"
     yii=ddii(:,1)
  endif
end subroutine vca_get_dust_lattice

subroutine vca_get_dund_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "vca_get_dund error: required N_sites > evaluated N_sites"
     yii=ddii(:,2)
  endif
end subroutine vca_get_dund_lattice

subroutine vca_get_dse_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "vca_get_dse error: required N_sites > evaluated N_sites"
     yii=ddii(:,3)
  endif
end subroutine vca_get_dse_lattice

subroutine vca_get_dph_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "vca_get_dph error: required N_sites > evaluated N_sites"
     yii=ddii(:,4)
  endif
end subroutine vca_get_dph_lattice
