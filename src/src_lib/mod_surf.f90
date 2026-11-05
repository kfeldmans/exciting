module modsurf

!  use mod_atoms, only: natmtot
implicit none

integer              :: naxismax
integer,allocatable  :: naxis(:) !ias
integer,allocatable  :: axisy(:,:),axisz(:,:) !iax,ias
real(8),allocatable  :: axisx_sph(:,:) !iax,ias
complex(8),allocatable :: ylm_mat(:,:,:) !lm,iax,ias
complex(8),allocatable :: ylm_tmat(:,:,:) !lm,iax,ias
real(8),allocatable   :: raxis(:,:,:)! 3 ,iax,ias
complex(8),allocatable :: interp_exp(:,:,:) !ngrid(1),iax,ais
contains

subroutine generate_surf_grid(lmax)
  use modmpi, only: mpiglobal
  use invert, only: zinvert_lapack
  use m_linalg, only: zlsp
  use modinput
  use mod_Gvector, only: ngrtot, ngrid, ngvec, cfunir
  use mod_atoms, only: nspecies, natoms, atposc, idxas,natmtot
  use mod_muffin_tin, only: rmt
  use constants, only: zzero, pi
  implicit none
  
  real(8) :: a1(3),a2(3),a3(3),r0(3),rv(3),rv1(3)
  real(8) :: xmu
  real(8) :: ratom(3), disc
  integer, intent(In) :: lmax
  integer :: lmmax
  integer :: is,ia,ias,yi,zi,p1,p2,p3
  integer :: ir,ir2,iq
  real(8) :: irf, xx
  integer :: nax,iax
  integer :: tmpsize
  real(8),allocatable  :: axisx_sph_tmp(:,:), tp_tmp(:,:,:) ! 2,iax,ias
  real(8),allocatable ::  raxis_tmp(:,:,:) !3,iax,ias
  integer,allocatable  :: axisy_tmp(:,:),axisz_tmp(:,:)
  complex(8),allocatable  :: A(:,:),Ainv(:,:),AA(:,:)
  complex(8), allocatable :: interp_tmp(:,:), zftmat(:,:)
  complex(8) :: ii,zt1
  integer :: lmax_ylm, lmax_tmp
  integer,allocatable  :: lmax_ias(:)
  complex(8),allocatable :: ylm_mat_large(:,:,:),ylm_tmat_large(:,:,:)
  ii=cmplx(0d0,1d0,8)

  lmax_ylm=20
  lmmax=(lmax+1)**2


  !write(*,*)"Režģa izmērs",ngrid
  !write(*,*)"x ass garums",sqrt(sum(input%structure%crystal%basevect(:, 1)**2))


  a1=input%structure%crystal%basevect(:, 1)/ngrid(1)
  a2=input%structure%crystal%basevect(:, 2)/ngrid(2)
  a3=input%structure%crystal%basevect(:, 3)/ngrid(3)
tmpsize=ngrid(2)*ngrid(3)*2
  Allocate(axisy_tmp(tmpsize,natmtot),axisz_tmp(tmpsize,natmtot),axisx_sph_tmp(tmpsize,natmtot))
  allocate(naxis(natmtot))
  allocate(tp_tmp(2,tmpsize,natmtot))
  allocate(raxis_tmp(3,tmpsize,natmtot))
  ia=1
  is=1
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      nax=0
      do p1=-1,1
        do p2=-1,1
          do p3=-1,1
            ratom=atposc (:, ia, is) + p1*input%structure%crystal%basevect(:, 1) + p2*input%structure%crystal%basevect(:, 2) + p3*input%structure%crystal%basevect(:, 3)
            do yi=0,ngrid(2)-1
              do zi=0,ngrid(3)-1

                r0 = yi*a2 + zi*a3
                rv=r0-ratom

                disc= sum(a1*rv)**2 - sum(a1**2)*(sum(rv**2)-rmt(is)**2)
              
                if (disc.gt.1d-9) then
                  xmu= ( -sum(a1*rv) - sqrt(disc) )/ sum(a1**2)
                  if ((xmu.ge.0d0).and.(xmu.lt.(ngrid(1)))) then
                    nax = nax+1
                    axisx_sph_tmp(nax,ias)=xmu
                    axisy_tmp(nax,ias)=yi
                    axisz_tmp(nax,ias)=zi
                    rv1=r0 + a1*xmu - ratom
                    raxis_tmp(:,nax,ias)=rv1
                    !write(*,*)"rmt1",sqrt(sum((rv1)**2)) 
                    tp_tmp(1,nax,ias) = atan2(sqrt(rv1(1)**2+rv1(2)**2),rv1(3))
                    tp_tmp(2,nax,ias) = atan2(rv1(2),rv1(1))
 
                  endif

                  xmu= ( -sum(a1*rv) + sqrt(disc) )/ sum(a1**2)
                  if ((xmu.ge.0d0).and.(xmu.lt.(ngrid(1)))) then
                    nax = nax+1
                    axisx_sph_tmp(nax,ias)=xmu
                    axisy_tmp(nax,ias)=yi
                    axisz_tmp(nax,ias)=zi
                    rv1=r0 + a1*xmu - ratom
                    raxis_tmp(:,nax,ias)=rv1
                    !write(*,*)"rmt2",sqrt(sum((rv1)**2)) 
                    tp_tmp(1,nax,ias) = atan2(sqrt(rv1(1)**2+rv1(2)**2),rv1(3))
                    tp_tmp(2,nax,ias) = atan2(rv1(2),rv1(1))
                  endif
                  
                endif
              enddo !zi
            enddo !yi
          enddo !p3
        enddo !p2
      enddo !p1
      naxis(ias)=nax
    enddo !ia
  enddo !is
  naxismax=maxval(naxis(:))
! write(*,*)"maxasis",naxismax
! write(*,*)"naxis",naxis
! write(*,*)"natmtot", natmtot

Allocate(axisy(naxismax,natmtot))
allocate(axisz(naxismax,natmtot))
allocate(axisx_sph(naxismax,natmtot))
allocate(ylm_mat(lmmax,naxismax,natmtot))
allocate(ylm_tmat(lmmax,naxismax,natmtot))
Allocate(raxis(3,naxismax,natmtot))

allocate(lmax_ias(natmtot))
do ias=1, natmtot
  lmax_tmp=int(sqrt(dble(naxis(ias)/2d0)))-1 !! 2d0 or 1.6d0 ????
  if (lmax_tmp.ge.lmax_ylm)then
    lmax_ias(ias)=lmax_ylm
  else if(lmax_tmp.le.lmax) then
    lmax_ias(ias)=lmax
  else
    lmax_ias(ias)=lmax_tmp
  endif
enddo

Allocate(ylm_mat_large((lmax_ylm+1)**2,naxismax,natmtot))
Allocate(ylm_tmat_large((lmax_ylm+1)**2,naxismax,natmtot))



ylm_mat=zzero
axisy=0
axisz=0
axisx_sph=0d0
do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    ratom=atposc (:, ia, is) 
    axisy(1:naxis(ias),ias)=axisy_tmp(1:naxis(ias),ias)
    axisz(1:naxis(ias),ias)=axisz_tmp(1:naxis(ias),ias)
    axisx_sph(1:naxis(ias),ias)=axisx_sph_tmp(1:naxis(ias),ias)
    
    do iax=1, naxis(ias)
      raxis(:,iax,ias)=raxis_tmp(:,iax,ias)+ratom
      Call genylm (lmax_ias(ias), tp_tmp(:, iax,ias), ylm_mat_large(1:(lmax_ias(ias)+1)**2,iax,ias))
    enddo
    
  enddo
enddo

do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    call inv_svd(naxis(ias),(lmax_ias(ias)+1)**2,transpose(ylm_mat_large(1:(lmax_ias(ias)+1)**2,1:naxis(ias),ias)),ylm_tmat_large(1:(lmax_ias(ias)+1)**2,1:naxis(ias),ias))
    ylm_mat(1:lmmax,1:naxis(ias),ias)=ylm_mat_large(1:lmmax,1:naxis(ias),ias)
    ylm_tmat(1:lmmax,1:naxis(ias),ias)=ylm_tmat_large(1:lmmax,1:naxis(ias),ias)
  enddo!ia
enddo!is

deallocate(ylm_mat_large)
deallocate(ylm_tmat_large)

allocate(zftmat(ngrid(1),ngrid(1)))
  do iax=0,ngrid(1)-1
    do iq=0,ngrid(1)-1 
      zftmat(iq+1,iax+1)=dcmplx(cos(2*iax*iq*pi/dble(ngrid(1))), -sin(2*iax*iq*pi/dble(ngrid(1))) )
    enddo
  enddo
  zt1=1/dble(ngrid(1))

allocate(interp_tmp(ngrid(1),naxismax))
allocate(interp_exp(ngrid(1),naxismax,natmtot))
interp_tmp=zzero

do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    
    do iax=1,naxis(ias)
      xx=axisx_sph(iax,ias)
      do ir=1, ceiling(dble(ngrid(1))/2d0)
        ir2=ir-1
        !write(*,*)ir,",",ir2
        irf=2d0*pi*dble(ir2)/dble(ngrid(1))
        interp_tmp(ir,iax)=exp(ii*irf*xx)
      enddo

      if (mod(ngrid(1),2).eq.0)then
        ir=ngrid(1)/2+1
        ir2=ngrid(1)/2
        !write(*,*)"vidus",ir,ir2
        irf=2d0*pi*dble(ir2)/dble(ngrid(1))
        interp_tmp(ir,iax)=cos (irf*xx)
      endif

      do ir=ceiling((dble(ngrid(1))-1)/2)+2,ngrid(1)
        ir2=ngrid(1)-ir+1
        !write(*,*)ir,",",ir2
        irf=2d0*pi*dble(-ir2)/dble(ngrid(1))
        interp_tmp(ir,iax)=exp(ii*irf*xx)
      enddo
    enddo !iax
    Call zgemm ('N', 'N', ngrid(1), naxismax, ngrid(1), &
               & zt1, zftmat, ngrid(1), interp_tmp(1,1), ngrid(1), zzero, &
               & interp_exp(1,1,ias), ngrid(1))

  enddo !ia
enddo !is

deallocate(interp_tmp,zftmat)




deallocate(axisy_tmp,axisz_tmp,axisx_sph_tmp,tp_tmp,raxis_tmp)
! open(11,file='surf.dat',status='replace')
! open(12,file='surf2.dat',status='replace')
!   do is=1, nspecies
!     do ia=1, natoms(is)
!       ias=idxas(ia,is)
!       do iax=1, naxis(ias)
!         rv(:)=axisx_sph(iax,ias)*a1 + axisy(iax,ias)*a2 + axisz(iax,ias)*a3
!         write(11,*)rv(1),rv(2),rv(3)
!         write(12,*)raxis(1,iax,ias),raxis(2,iax,ias),raxis(3,iax,ias)
!       enddo
!     enddo
!   enddo
! close(11)
! close(12)
if(mpiglobal%is_root) then
open(11,file='surf_info.dat',status='replace')
write(11,*)"lmmax:",lmmax
write(11,*)"is, ia, number_of_ponts_on_MTsurface,lmax_ylm,lmmax, naxis/lmmax"
do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    write(11,*)is,ia,naxis(ias),lmax_ias(ias),(lmax_ias(ias)+1)**2,real(naxis(ias))/real((lmax_ias(ias)+1)**2)
  enddo
enddo
close(11)
endif

deallocate(lmax_ias)
end subroutine 

subroutine surf_pot(lmax,zvclir,igfft,qvec,vlm)
  use modinput
  use m_linalg, only:zlsp
  use mod_Gvector, only: ngrtot, ngrid, ngvec, cfunir,vgc
  use mod_atoms, only: nspecies, natoms, atposc, idxas,natmtot
  use mod_muffin_tin, only: rmt
  use constants, only: zzero, zone, pi
  implicit none
  integer  , intent(in) :: lmax
  complex(8) ,intent(inout) :: zvclir(:)
  integer, intent(in) :: igfft(:)
  real(8),intent(in) :: qvec(3)
  complex(8) ,intent(out) :: vlm(:,:)!lm,ias
  complex(8) :: zvaxft(ngrid(1))
  integer :: ig1,ig2,iax,is,ia,ias,iy,iz 
  complex(8),allocatable :: vmu(:,:),zvaxft_all(:,:,:),interp2(:,:,:)
  complex(8),allocatable :: zftmat(:,:)
  complex(8) :: vmu1(naxismax)
  real(8) :: xx, irf,ir2,ta,tb,rv(3)
  complex(8) ::zt1, ii
  integer :: ir, lmmaxvr,iq
 
  complex(8) :: phase
  complex(8), allocatable :: zftvec(:)
  complex(8), external :: zdotu

  ii=cmplx(0d0,1d0,8)
  


call timesec(ta)
  lmmaxvr=size(vlm,1)

  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      
      do iax=1,naxis(ias)
        iz=axisz(iax,ias)
        iy=axisy(iax,ias)
        ig1=iz*ngrid(2)*ngrid(1) + iy*ngrid(1) + 1

        rv(:)=raxis(:,iax,ias)
        phase=exp(ii*sum(qvec*rv))
        vmu1(iax)=phase*zdotu(ngrid(1),zvclir(ig1),1,interp_exp(1,iax,ias),1)
      enddo !iax

      Call zgemv ('N', lmmaxvr, naxis(ias), zone, ylm_tmat(1,1,ias), lmmaxvr, vmu1, 1, zzero, vlm(1,ias), 1)

    enddo !ia
  enddo !is

call timesec(tb)
!write(*,*)"lm_komplekts:",tb-ta
! open(1, file = 'vlm_new.dat', status = 'replace')
! do ias=1, natmtot
! do ir=1,(lmax+1)**2
!   write(1,*)dble(vlm(ir,ias)),",",aimag(vlm(ir,ias)),",", ias, ir
! enddo
! enddo
! close(1)
! stop
end subroutine


subroutine inv_svd(m,n,A_in,Ainv)
  implicit none
  integer  , intent(in) :: m, n
  complex(8) ,intent(in) :: A_in(m,n)
  complex(8) ,intent(out) :: Ainv(n,m)
  complex(8) :: A(m,n),U(m,m),VT(n,n)
  real(8)    :: S(n),Sinv(n,m)
  real(8)    :: rwork(5*m) !can be decreased

  integer :: info ,work_size, lm
  complex(8),allocatable  :: work(:) 
  real(8),allocatable     :: Smat(:,:)
  complex(8),allocatable  :: E(:,:)
  
  A=A_in
  work_size=5*m
  allocate(work(work_size))
  call zgesvd('A', 'A', m, n, A, m, S, U, m, VT, n, work, work_size, rwork, info)

  deallocate(work)
  if (info.ne.0)then
    write(*,*)"info:",info
    write(*,*)"singular value decomposition went wrong in mod_surf.f90"
    write(*,*)"lmmax=",m
    write(*,*)"number of points on the sphere =",n
    stop
  endif
  !write(*,*)"info",info
  if (.false.) then !check if decomposition was correct
    write(*,*)"matrix before SVD:"
    write(*,*)A_in(1,1:3)
    write(*,*)A_in(2,1:3)
    write(*,*)A_in(3,1:3)
    allocate(Smat(m,n))
    Smat=0d0
    do lm=1,n
      Smat(lm,lm)=S(lm)
    enddo
    A = matmul( matmul(U ,Smat) , VT )
    write(*,*)"matrix constructed back:"
    write(*,*)A(1,1:3)
    write(*,*)A(2,1:3)
    write(*,*)A(3,1:3)
    deallocate (Smat)
  endif
  Sinv=0d0

  !write(*,*)"Diagonal matrix:"
  do lm=1,n
    !write(*,*)S(lm)
    Sinv(lm,lm)=1d0/S(lm)
    if (abs(S(lm)).lt.1e-10) then
      write(*,*)"singular value decomposition went wrong in mod_surf.f90"
      write(*,*)"diagonal matrix element nr.",lm," is small:",S(lm)
      write(*,*)"lmmax=",m
      write(*,*)"number of points on the sphere =",n
      stop
      !Sinv(lm,lm)=0d0
    endif
  enddo


  Ainv = matmul( matmul(conjg(transpose(VT)),Sinv)  , conjg(transpose(U)))

  if (.false.) then !check if A^(-1)*A=E
    Allocate(E(n,n))
    E=matmul(Ainv,A)
    !write(*,*)"mazā Vienības matica:"
    !write(*,*)E(1,1:3)
    !write(*,*)E(2,1:3)
    !write(*,*)E(3,1:3)
    Write(*,*)"Diagonal of the unit matrix:"
    do lm=1,n
      write(*,*)E(lm,lm)
    enddo
    deallocate(E)
  endif

end subroutine


end module
