MODULE VCA

  USE VCA_INPUT_VARS
  

  USE VCA_AUX_FUNX, only:                        &
       !vca_lso2nnn_reshape                     ,&
       !vca_nnn2lso_reshape                     ,&
       !vca_so2nn_reshape                       ,&
       !vca_nn2so_reshape                       ,&
       !vca_set_Hcluster                        ,&
       search_chemical_potential


  USE VCA_IO, only:                              &
       vca_get_sigma_matsubara                  ,&
       vca_get_gimp_matsubara                   ,&
       vca_get_sigma_realaxis                   ,&
       vca_get_gimp_realaxis                    ,&
       vca_get_dens                             ,&
       vca_get_mag                              ,&
       vca_get_docc                             ,&
       vca_get_sft_potential                    ,&
       vca_read_impsigma                        ,&
       vca_gf_cluster


  USE VCA_BATH_SETUP, only:                     &
       vca_get_bath_dimension                  ,&
       set_bath_component            


  USE VCA_MAIN, only:                           &
       vca_init_solver                         ,&
       vca_solve

  USE VCA_OBSERVABLES, only:                    &
       init_custom_observables                , &
       clear_custom_observables               , &
       add_custom_observable


END MODULE VCA
