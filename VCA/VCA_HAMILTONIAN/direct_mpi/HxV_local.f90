  do i=1,Nloc
     iup = iup_index(i+mpiIshift,DimUp)
     idw = idw_index(i+mpiIshift,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     ibup = bdecomp(mup,Ns)
     ibdw = bdecomp(mdw,Ns)
     !
     do ilat=1,Nlat
        do iorb=1,Norb
           nup(ilat,iorb)=dble(ibup(imp_state_index(ilat,iorb)))
           ndw(ilat,iorb)=dble(ibdw(imp_state_index(ilat,iorb)))
        enddo
     enddo
     !> HxV_imp: Diagonal Elements, i.e. local part
     htmp = zero
     do ilat=1,Nlat
       do iorb=1,Norb
          htmp = htmp + impHloc(ilat,ilat,1,1,iorb,iorb)*Nup(ilat,iorb)
          htmp = htmp + impHloc(ilat,ilat,Nspin,Nspin,iorb,iorb)*Ndw(ilat,iorb)
          htmp = htmp - xmu*( Nup(ilat,iorb)+Ndw(ilat,iorb) )
       enddo
     enddo
     !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
     !
     !density-density interaction: same orbital, opposite spins:
     ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
      do ilat=1,Nlat
         do iorb=1,Norb
            htmp = htmp + Uloc_per_site(ilat,iorb)*nup(ilat,iorb)*ndw(ilat,iorb)
         enddo
      enddo
      if(Norb>1)then
         !density-density interaction: different orbitals, opposite spins:
         ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
         ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
         do ilat=1,Nlat
            do iorb=1,Norb
               do jorb=iorb+1,Norb
                  htmp = htmp + Ust_per_site(ilat)*(nup(ilat,iorb)*ndw(ilat,jorb) + nup(ilat,jorb)*ndw(ilat,iorb))
               enddo
            enddo
         enddo
         !density-density interaction: different orbitals, parallel spins
         ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
         ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
         do ilat=1,Nlat
            do iorb=1,Norb
               do jorb=iorb+1,Norb
                  htmp = htmp + (Ust_per_site(ilat)-Jh_per_site(ilat))*(nup(ilat,iorb)*nup(ilat,jorb) + ndw(ilat,iorb)*ndw(ilat,jorb))
               enddo
            enddo
         enddo
      endif
      !if using the Hartree-shifted chemical potential: mu=0 for half-filling
      !sum up the contributions of hartree terms:
      if(hfmode)then
         do ilat=1,Nlat
            do iorb=1,Norb
               htmp = htmp - 0.5d0*Uloc_per_site(ilat,iorb)*(nup(ilat,iorb)+ndw(ilat,iorb)) + 0.25d0*Uloc_per_site(ilat,iorb)
            enddo
         enddo
         if(Norb>1)then
            do ilat=1,Nlat
               do iorb=1,Norb
                  do jorb=iorb+1,Norb
                     htmp=htmp-0.5d0*Ust_per_site(ilat)*(nup(ilat,iorb)+ndw(ilat,iorb)+nup(ilat,jorb)+ndw(ilat,jorb))+0.25d0*Ust_per_site(ilat)
                     htmp=htmp-0.5d0*(Ust_per_site(ilat)-Jh_per_site(ilat))*(nup(ilat,iorb)+ndw(ilat,iorb)+nup(ilat,jorb)+ndw(ilat,jorb))+0.25d0*(Ust_per_site(ilat)-Jh_per_site(ilat))
                  enddo
               enddo
            enddo
         endif
      endif
     !
     !
     !> H_Bath: local bath energy contribution.
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     if(Nbath>0)then
       do ilat=1,size(bath_diag,1)
         do iorb=1,size(bath_diag,3)
           do kp=1,Nbath
               ialfa = getBathStride(ilat,iorb,kp)
               htmp =htmp + bath_diag(ilat,1    ,iorb,kp)*ibup(ialfa) !UP
               htmp =htmp + bath_diag(ilat,Nspin,iorb,kp)*ibdw(ialfa) !DW
               htmp = htmp - xmu*dble(ibup(ialfa)+ibdw(ialfa))
            enddo
         enddo
       enddo
     endif
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
  enddo


