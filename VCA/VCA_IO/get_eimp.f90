subroutine vca_get_eimp_(eimp)
  real(8),dimension(4) :: eimp
  eimp = [vca_Epot,vca_Eint,vca_Ehartree,vca_Eknot]
end subroutine vca_get_eimp_

subroutine vca_get_epot_(eimp)
  real(8) :: eimp
  eimp = vca_Epot
end subroutine vca_get_epot_

subroutine vca_get_eint_(eimp)
  real(8) :: eimp
  eimp = vca_Eint
end subroutine vca_get_eint_

subroutine vca_get_ehartree_(eimp)
  real(8) :: eimp
  eimp = vca_Ehartree
end subroutine vca_get_ehartree_

subroutine vca_get_eknot_(eimp)
  real(8) :: eimp
  eimp = vca_Eknot
end subroutine vca_get_eknot_




subroutine vca_get_eimp_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "vca_get_eimp error: required N_sites > evaluated N_sites"
     yii=eii
  endif
end subroutine vca_get_eimp_lattice

subroutine vca_get_epot_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "vca_get_epot error: required N_sites > evaluated N_sites"
     yii=eii(:,1)
  endif
end subroutine vca_get_epot_lattice

subroutine vca_get_eint_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "vca_get_eint error: required N_sites > evaluated N_sites"
     yii=eii(:,2)
  endif
end subroutine vca_get_eint_lattice

subroutine vca_get_ehartree_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "vca_get_ehartree error: required N_sites > evaluated N_sites"
     yii=eii(:,3)
  endif
end subroutine vca_get_ehartree_lattice

subroutine vca_get_eknot_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "vca_get_knot error: required N_sites > evaluated N_sites"
     yii=eii(:,4)
  endif
end subroutine vca_get_eknot_lattice
