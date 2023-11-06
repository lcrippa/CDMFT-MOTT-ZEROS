  do jdw=1,MpiQup
     do jup=1,DimDw
        mdw  = Hs(2)%map(jup)
        ibdw  = bdecomp(mdw,Ns)
        !
        j    = jup + (jdw-1)*DimDw
        !
        !
        !> H_imp: Off-diagonal elements, i.e. non-local part. 
        !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
         do ilat=1,Nlat
             do jlat=1,Nlat
              do iorb=1,Norb
                 do jorb=1,Norb
                       is = imp_state_index(ilat,iorb) ! 1 is actually not used
                       js = imp_state_index(jlat,jorb)
                    Jcondition = &
                         (impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                         (ibdw(js)==1) .AND. (ibdw(is)==0)
                    if (Jcondition) then
                       call c(js,mdw,k1,sg1)
                       call cdg(is,k1,k2,sg2)
                       iup = binary_search(Hs(2)%map,k2)
                       idw = jdw             
                       i   = iup + (idw-1)*DimDw
                       htmp = impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2
                       !
                       Hvt(i) = Hvt(i) + htmp*vt(j)
                       !
                    endif
                 enddo
              enddo
            enddo
          enddo
        !>H_hyb: hopping terms for a given spin (imp <--> bath)
        if(Nbath>0)then
          do ilat=1,Nlat
            do iorb=1,Norb
               do kp=1,Nbath
                  ialfa=getBathStride(ilat,iorb,kp)
                  is = imp_state_index(ilat,iorb) ! 1 is actually not used
                  !
                  if( (diag_hybr(ilat,Nspin,iorb,kp)/=0d0) .AND. &
                       (ibdw(is)==1) .AND. (ibdw(ialfa)==0) )then
                     call c(is,mdw,k1,sg1)
                     call cdg(ialfa,k1,k2,sg2)
                     iup = binary_search(Hs(2)%map,k2)
                     idw = jdw             
                     i   = iup + (idw-1)*DimDw
                     htmp=diag_hybr(ilat,Nspin,iorb,kp)*sg1*sg2
                     !
                     hvt(i) = hvt(i) + htmp*vt(j)
                     !
                  endif
                  if( (diag_hybr(ilat,Nspin,iorb,kp)/=0d0) .AND. &
                       (ibdw(is)==0) .AND. (ibdw(ialfa)==1) )then
                     call c(ialfa,mdw,k1,sg1)
                     call cdg(is,k1,k2,sg2)
                     iup = binary_search(Hs(2)%map,k2)
                     idw = jdw             
                     i   = iup + (idw-1)*DimDw
                     htmp=diag_hybr(ilat,Nspin,iorb,kp)*sg1*sg2
                     !
                     hvt(i) = hvt(i) + htmp*vt(j)
                     !
                  endif
               enddo
            enddo
          enddo
        endif
     enddo
  enddo

