MODULE VCA_GREENS_FUNCTIONS
  !
  USE VCA_GF_SHARED
  USE VCA_GF_NORMAL
  !
  implicit none
  private 


  public :: build_gf_cluster



contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Build the GF of the cluster 
  !+------------------------------------------------------------------+
  subroutine build_gf_cluster()
    !
    call allocate_grids
    !    
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    !
    write(LOGfile,"(A)")"Get cluster Greens functions:"
    call build_gf_normal()
    call build_sigma_normal()
    !
    if(print_Sigma)call vca_print_impSigma()
    if(print_impG)call vca_print_impG()
    if(print_impG0)call vca_print_impG0()
    !
    call deallocate_grids
    !
  end subroutine build_gf_cluster


end MODULE VCA_GREENS_FUNCTIONS













