!NORMAL, MATSUBARA SELF-ENEGRGY
subroutine ed_get_sigma_matsubara_1(Smats)
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:,:,:) = impSmats(:,:,:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_1

subroutine ed_get_sigma_matsubara_2(Smats)
  complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats),intent(inout) :: Smats
  integer  :: io,jo,ilat,jlat,iorb,jorb,ispin,jspin
  do ilat=1,Nlat
    do jlat=1,Nlat
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  io = index_stride_lso(ilat,ispin,iorb)
                  jo = index_stride_lso(jlat,jspin,jorb)
                  Smats(io,jo,:) = impSmats(ilat,jlat,ispin,jspin,iorb,jorb,:)
               enddo
            enddo
         enddo
      enddo
    enddo
  enddo
end subroutine ed_get_sigma_matsubara_2

subroutine ed_get_sigma_matsubara_3(Smats,ilat,jlat,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lmats),intent(inout) :: Smats
  integer                          :: ilat,jlat,iorb,jorb,ispin,jspin
  Smats(:) = impSmats(ilat,jlat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_sigma_matsubara_3



