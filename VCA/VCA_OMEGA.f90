MODULE VCA_OMEGA
  !
  USE VCA_GF_SHARED
  USE VCA_BATH_FUNCTIONS
  USE VCA_GF_NORMAL
  USE SF_INTEGRATE
  USE SF_LINALG
  USE SF_IOTOOLS, only: splot
  !
  implicit none
  private 

  public                                            :: sum_kmesh
  public                                            :: frequency_integration
  public                                            :: frequency_integration_finite_t

  public                                            :: reconstruct_g
  public                                            :: test_g

  !public                                            :: frequency_integration_sample

  !public                                            :: test_ksum

  complex(8),allocatable,dimension(:,:,:,:,:,:)     :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: gftest ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  real(8)                                           :: integrationR



contains

 
  !+------------------------------------------------------------------+
  !PURPOSE  : Function that does the k-sum in the RBZ of the VCA variational function(al)
  !+------------------------------------------------------------------+
  function sum_kmesh(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    complex(8)                                               :: omega
    real(8)                                                  :: out_1
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    complex(8),allocatable,dimension(:,:,:,:,:,:)            :: deltamat 
    !
    out_1=0.d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    if(allocated(gfprime))deallocate(gfprime)
    if(allocated(deltamat))deallocate(deltamat)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(deltamat(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    deltamat=zero
    !
    !    
    call vca_gf_cluster(omega,gfprime)
    if(Nbath>0)then
      deltamat=delta_bath_freq(omega,vca_bath)
    endif
    !
    do ii=1,size(impHk,7)
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(deltamat+impHloc-impHk(:,:,:,:,:,:,ii),Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
       out_1=out_1+log(abs(det(tmp_mat)))
    enddo
    out_1=out_1/size(impHk,7)
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    return
    !
  end function sum_kmesh


  function sum_kmesh_complex(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    complex(8)                                               :: omega
    complex(8)                                               :: out_1
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    complex(8),allocatable,dimension(:,:,:,:,:,:)            :: deltamat 
    !
    out_1=0.d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    if(allocated(gfprime))deallocate(gfprime)
    if(allocated(deltamat))deallocate(deltamat)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(deltamat(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    deltamat=zero
    !
    !    
    call vca_gf_cluster(omega,gfprime)
    if(Nbath>0)then
      deltamat=delta_bath_freq(omega,vca_bath)
    endif
    !
    do ii=1,size(impHk,7)
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(deltamat+impHloc-impHk(:,:,:,:,:,:,ii),Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
       out_1=out_1+log(det(tmp_mat))
    enddo
    out_1=out_1/size(impHk,7)
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    return
    !
  end function sum_kmesh_complex


  !+------------------------------------------------------------------+
  !PURPOSE  : Do the frequency sum at zero T
  !+------------------------------------------------------------------+

  function frequency_integration() result(out_2)
    integer                         :: inf
    real(8)                         :: out_2,spin_multiplicity
    !
    out_2=0.d0
    spin_multiplicity=3.d0-Nspin 
    !
    !
    write(LOGfile,"(A)")"Calculating Omega"
    call quad(imaginary_axis,a=0.0d0,inf=1,verbose=(verbose>=3),result=out_2,strict=.false.)
    !
    out_2=spin_multiplicity*out_2/pi 
    return
  end function frequency_integration


 function imaginary_axis(zeta) result(f)
    real(8)                 :: zeta,f
    complex(8)              :: w
    !
    w=xi*zeta
    !
    f=sum_kmesh(w)
 end function imaginary_axis


  !+------------------------------------------------------------------+
  !PURPOSE  : Do the frequency sum at finite T
  !+------------------------------------------------------------------+

  function frequency_integration_finite_t() result(out_2)
    integer                         :: inf,Nmax,ii
    real(8)                         :: out_2,spin_multiplicity,omegamax,integralpart
    !
    !1) Find the real omegamax
    nmax=int(2*(abs(max_exc)+bandwidth)*beta/pi)
    if (mod(nmax,2)==0)then
      nmax=nmax/2    
    else
      nmax=(nmax+1)/2
    endif
    integrationR=2*(nmax+1)*pi/beta
    print*,"NMAX=",nmax
    print*,"INTEGRATION R=",integrationR
    !2) Evaluate discrete sum
    !
    out_2=0.d0
    do ii=0,Nmax
      out_2=out_2+dreal(sum_kmesh_complex(xi*(2*ii+1)*pi/beta))
    enddo
    !
    out_2=2.d0*(1/beta)*out_2
    print*,"SUM PART = ",out_2
    !
    !3) Evaluate integral part
    integralpart=0.d0
    call quad(integral_contour,a=0.0d0,b=pi,verbose=(verbose>=3),key=6,result=integralpart,strict=.false.)
    !
    print*,"INTEGRAL PART = ",integralpart
    !4) Sum all
    out_2=out_2+integralpart
    !5) Spin trick
    spin_multiplicity=3.d0-Nspin 
    out_2=spin_multiplicity*out_2
    return
  end function frequency_integration_finite_t


 function integral_contour(zeta) result(f)
    real(8)                 :: zeta,f
    complex(8)              :: w,fermi
    !
    w=integrationR*exp(xi*zeta)
    if(dreal((w-XMU)*beta)>= 100)then
      fermi=0.d0
    else
      fermi=(1/(exp(beta*(w-XMU))+1))
    endif
    !
    f=dreal((1/pi)*w*fermi*sum_kmesh_complex(w))
    !print*,zeta,f,fermi,sum_kmesh(w)
 end function integral_contour



 !+------------------------------------------------------------------+
  !PURPOSE  : DEBUG
  !+------------------------------------------------------------------+

  subroutine reconstruct_g
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    character(len=64)                                        :: suffix
    integer                                                  :: ilat,jlat,iorb,jorb,ispin,ifreq
    !
    allocate(gftest(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    gftest=0
    !
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    !
    do ifreq=1,size(wm)
        call vca_gf_cluster(xi*wm(ifreq),gftest(:,:,:,:,:,:,ifreq))
    enddo
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                call splot("Gtest"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,gftest(ilat,jlat,ispin,ispin,iorb,jorb,:))
             enddo
            enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(gftest))deallocate(gftest)
    !
  end subroutine reconstruct_g


  subroutine test_g
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    character(len=64)                                        :: suffix
    integer                                                  :: ilat,jlat,iorb,ispin,ifreq
    !
    allocate(gftest(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    gftest=0
    !
    !
    call vca_gf_cluster(xi*0.5d0,gftest(:,:,:,:,:,:,1))
    if(allocated(gftest))deallocate(gftest)
    !
  end subroutine test_g

  !subroutine test_ksum
  !  integer                                                  :: ii,jj,kk
  !  real(8)                                                  :: omega
  !  real(8)                                                  :: out_1
  !  complex(8),allocatable,dimension(:,:)                    :: tmp_mat
  !  !
  !  out_1=0.d0
  !  omega=1.0d0
  !  !
  !  !
  !  if(allocated(tmp_mat))deallocate(tmp_mat)
  !  if(allocated(gfprime))deallocate(gfprime)
  !  !
  !  allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
  !  allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
  !  !
  !  tmp_mat=zero
  !  gfprime=zero
  !  !    
  !  call vca_gf_cluster(xi*omega,gfprime)
  !  do ii=1,size(impHk,7)
  !     print*,"ITERATION",ii
  !     !print*,"OMEGA",xi*omega
  !     !print*,"AAAAAAAAA",matmul(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
  !     !print*,"BBBBBBBBB",vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb)
  !     tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
  !     out_1=out_1+log(abs(REAL(det(tmp_mat))))
  !     print*,"CCCCCCCCC",(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb))
  !  enddo
  !  out_1=out_1/size(impHk,7) !*(pi**Ndim)
  !  !
  !  deallocate(tmp_mat)
  !  deallocate(gfprime)
  ! return
  ! !
  !end subroutine test_ksum

  !function frequency_integration_sample() result(out_2)
  !    integer                          :: N,i
  !    real(8),dimension(:),allocatable :: x,func
  !    real(8)                          :: out_2,a,b,spin_multiplicity
  !    !
  !    a=0.0001d0
  !    b=99999.d0
  !    N=1000
  !   out_2=0.d0
  !    spin_multiplicity=2.d0
  !    allocate(x(N),func(N))
  !    x = linspace(a,b,N)
  !    do i=1,N
  !       func(i) = sum_kmesh(xi*x(i))
  !    enddo
  !    call quad(func,a,b,Ninterp=3,key=6,epsabs=0d0,epsrel=1d-4,verbose=.true.,result=out_2)
  !    out_2=spin_multiplicity*out_2/pi 
  !    deallocate(x,func)
  !    return
  !end function frequency_integration_sample

end MODULE VCA_OMEGA













