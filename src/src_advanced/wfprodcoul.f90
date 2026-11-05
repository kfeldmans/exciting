subroutine WFprodcoul(ist1,wf1,ist2,wf2,prod,qlm)
      use modinput
      use mod_APW_LO
      use mod_atoms
      use mod_muffin_tin
      use mod_eigenvalue_occupancy
      use constants, only : zzero, zone, fourpi
      use mod_SHT
      use mod_Gvector, only : ngrtot
      use mod_eigensystem, only : WFType
      use weinert, only : poisson_and_multipoles_mt_yukawa2
      use mod_radial
      use mod_eigensystem, only : gntyyy,gntyyyT,gntlyy,gntyyl
!      use wigner3j_symbol, only : gaunt_yyy
! !USES:
! !DESCRIPTION:
! Evaluates a product of two WFs in the real space
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
      implicit none
      integer, intent(in) :: ist1,ist2
      type (WFType) :: wf1,wf2,prod
      complex(8), intent(out) :: qlm(lmmaxvr,natmtot)

      integer :: is,ia,ias
      integer :: l1,l2,m1,m2,lm1,lm2,io1,io2,if1,if2,ilo1,ilo2,lmmaxprod,lm,ir,l,m,ilo,io,if2offset
      integer :: if1offset
      integer :: irad,jrad,irad1,irad2,ipbf,radoffset
      integer :: lmax,llow,lhi
      integer :: blkstart,chunksize,iroffset
      integer, parameter :: blksize=64
      complex(8), allocatable :: factors(:), rho(:,:), fr(:)
      complex(8) :: mtmesh(ntpll,blksize), mtmesh1(ntpll,blksize), mtrlm1(blksize,lmmaxvr)
      complex(8) :: zfmt(lmmaxvr,nrmtmax),zfmttr(nrmtmax,lmmaxvr),zpot(lmmaxvr,nrmtmax),qlm2(lmmaxvr),H3(lmmaxvr)
      real(8) :: refmttr(nrmtmax,lmmaxvr),imfmttr(nrmtmax,lmmaxvr)
      complex(8) :: zt,zt2
      real(8) :: ta,tb
      complex(8),allocatable :: H(:,:,:), F(:,:), HH(:,:,:), FF(:,:), H2(:,:)
      real(8),allocatable :: reFF(:,:),imFF(:,:)
      real(8),external :: oldgaunt,oldwigner3j
 
      if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
      if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
 
! call timesec(ta)
qlm=0d0 !(:,ias)
      lmax= input%groundstate%lmaxvr

      allocate(H(lmmaxvr,maxradial,maxradial))
      allocate(HH(lmmaxvr,maxradial2,maxradial))
      allocate(H2(lmmaxvr,maxradial2))
      allocate(F(maxpbf,lmmaxvr))
      allocate(FF(maxpbf2,lmmaxvr))
      allocate(reFF(maxpbf2,lmmaxvr))
      allocate(imFF(maxpbf2,lmmaxvr))

      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)

if (.true.) then 
          H=0d0
          F=0d0
          
          if2=0
          do irad2=1,nradial(is)
            l2=lrad(irad2,is)
            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              if2=if2+1
              if1=0

              do irad1=1,nradial(is)
                l1=lrad(irad1,is)
                do m1=-l1,l1
                  lm1=idxlm(l1,m1)
                  if1=if1+1
                  zt=wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))

                  do L=abs(l1-l2),min(lmax,l1+l2),2
                    M=m1-m2
                    if ((M.le.L).and.(M.ge.-L)) then
                      LM=idxlm(L,M)
                      H(LM,irad1,irad2)=H(LM,irad1,irad2)+gntyyl(L,lm2,lm1)*zt !wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))                
!                      H(LM,irad1,irad2)=H(LM,irad1,irad2)+gntyyy(LM,lm2,lm1)*zt !wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))                
!                      H(LM,irad1,irad2)=H(LM,irad1,irad2)+gntyyy(LM,lm2,lm1)*wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))                
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo

          do irad2=1,nradial(is)
            l2=lrad(irad2,is)
            do irad1=1,nradial(is)
              l1=lrad(irad1,is)
              do L=abs(l1-l2),min(lmax,l1+l2),2
                do M=-L,L
                  LM=idxlm(L,M)
                  F(1:npbf(L,ias),LM)=F(1:npbf(L,ias),LM)+H(LM,irad1,irad2)*uproducts(1:npbf(L,ias),L,irad1,irad2,ias) 
                enddo
              enddo
            enddo
          enddo

!-----------------

          qlm(:,ias)=0d0
          do l=0,lmax
            do m=-l,l
              lm=idxlm(l,m)
              zt=0d0
              do ipbf=1,npbf(l,ias)
                qlm(lm,ias)=qlm(lm,ias)+pbfmultipoles(ipbf,l,ias)*F(ipbf,lm)
              enddo
            enddo
          enddo

!-----------------



!          HH=0d0
          FF=0d0
!          if2=0
          if2offset=0
          do irad2=1,nradial(is)
!            H2=0d0
            l2=lrad(irad2,is)

            irad1=0
            do l1=0,lmax
              do ipbf=1,npbf(l1,ias) !nradial2(ias)
                H3=0d0
                irad1=irad1+1

                do m1=-l1,l1
                  lm1=idxlm(l1,m1)
                  zt2=F(ipbf,lm1)

                  if2=if2offset
                  do m2=-l2,l2
                    lm2=idxlm(l2,m2)
                    if2=if2+1
                    zt=zt2*wf2%mt(if2,ist2,ias)

                    do L=abs(l1-l2),min(lmax,l1+l2),2

                      M=m1+m2
                      if ((M.le.L).and.(M.ge.-L)) then
                        LM=idxlm(L,M)
!                        H3(LM)=H3(LM)+gntlyy(L,lm1,lm2)*zt
                        H3(LM)=H3(LM)+gntlyy(L,lm2,lm1)*zt
                      endif 
                    enddo
                  enddo
                enddo

                do L=abs(l1-l2),min(lmax,l1+l2),2
                  do M=-L,L
                    LM=idxlm(L,M)
                    FF(1:npbf2(L,ias),LM)=FF(1:npbf2(L,ias),LM)+H3(LM)*uproducts2(1:npbf2(L,ias),L,irad1,irad2,ias) 
                  enddo
                enddo
                
              enddo
            enddo
            if2offset=if2

!            do irad1=1,nradial2(ias)
!              l1=lrad2(irad1,ias)
!              do L=abs(l1-l2),min(lmax,l1+l2),2
!                do M=-L,L
!                  LM=idxlm(L,M)
!                  FF(1:npbf2(L,ias),LM)=FF(1:npbf2(L,ias),LM)+H2(LM,irad1)*uproducts2(1:npbf2(L,ias),L,irad1,irad2,ias) 
!                enddo
!              enddo
!            enddo

          enddo

!          do irad2=1,nradial(is)
!            l2=lrad(irad2,is)
!            radoffset=0
!            do l1=0,lmax
!              do L=abs(l1-l2),min(lmax,l1+l2),2
!                do M=-L,L
!                  LM=idxlm(L,M)
!                  irad1=radoffset
!                  do ipbf=1,npbf(l1,ias)
!                    irad1=irad1+1 
!                    FF(1:npbf2(L,ias),LM)=FF(1:npbf2(L,ias),LM)+HH(LM,irad1,irad2)*uproducts2(1:npbf2(L,ias),L,irad1,irad2,ias) 
!                  enddo
!                enddo
!              enddo
!              radoffset=radoffset+npbf(l1,ias)
!            enddo
!          enddo

!          do irad2=1,nradial(is)
!            l2=lrad(irad2,is)
!            do irad1=1,nradial2(ias)
!              l1=lrad2(irad1,ias)
!              do L=abs(l1-l2),min(lmax,l1+l2),2
!                do M=-L,L
!                  LM=idxlm(L,M)
!                  FF(1:npbf2(L,ias),LM)=FF(1:npbf2(L,ias),LM)+HH(LM,irad1,irad2)*uproducts2(1:npbf2(L,ias),L,irad1,irad2,ias) 
!!                  call zaxpy(npbf2(L,ias), HH(LM,irad1,irad2), uproducts2(1,L,irad1,irad2,ias),1,FF(1,LM),1)
!!                 call daxpy(npbf2(L,ias),  dble(HH(LM,irad1,irad2)), uproducts2(1,L,irad1,irad2,ias),1,reFF(1,LM),1)
!!                 call daxpy(npbf2(L,ias), dimag(HH(LM,irad1,irad2)), uproducts2(1,L,irad1,irad2,ias),1,imFF(1,LM),1)
!                enddo
!              enddo
!            enddo
!          enddo

!-----------------
 
          refmttr=0d0
          imfmttr=0d0

          do L=0,lmax
            do M=-L,L
              LM=idxlm(L,M)
!              zt=0d0
              do ipbf=1,npbf2(L,ias)
!                call daxpy(nrmt(is), reFF(ipbf,LM),productbasis2(1,ipbf,L,ias),1,refmttr(1,LM),1)
!                call daxpy(nrmt(is), imFF(ipbf,LM),productbasis2(1,ipbf,L,ias),1,imfmttr(1,LM),1)
                call daxpy(nrmt(is), dble(FF(ipbf,LM)),productbasis2(1,ipbf,L,ias),1,refmttr(1,LM),1)
                call daxpy(nrmt(is),dimag(FF(ipbf,LM)),productbasis2(1,ipbf,L,ias),1,imfmttr(1,LM),1)
!                zfmttr(1:nrmt(is),lm)=zfmttr(1:nrmt(is),lm)+pbfpotential(1:nrmt(is),ipbf,l,ias)*F(ipbf,lm) 
              enddo
            enddo
          enddo

          do ir= 1,nrmt(is)
!            zpot(1:lmmaxvr,ir)= zfmttr(ir,1:lmmaxvr)
            prod%mtrlm(1:lmmaxvr,ir,ias,1)=dcmplx(refmttr(ir,1:lmmaxvr),imfmttr(ir,1:lmmaxvr))
!            zpot(1:lmmaxvr,ir)=dcmplx(refmttr(ir,1:lmmaxvr),imfmttr(ir,1:lmmaxvr))
          enddo


!-----------------


          
else
          chunksize=blksize
          do iroffset=1,nrmt(is),blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset

!-----------
! expand WFs in spherical harmonics
!-----------
            mtrlm1=0d0
            ias=idxas(ia,is)

            if1=0
            do irad = 1, nradial(is)
              l = lrad (irad,is)
              do m = - l, l
                if1=if1+1
                lm = idxlm (l, m)
                mtrlm1(1:chunksize,lm)=mtrlm1(1:chunksize,lm)+wf1%mt(if1,ist1,ias)*radial(iroffset:iroffset+chunksize-1,irad,ias)
              enddo
            enddo


            Call zgemm ('N', 'T', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
                       & mtrlm1, blksize, zzero, mtmesh, ntpll)
            do ir=1,chunksize
                 mtmesh(1:ntpll,ir)=mtmesh(1:ntpll,ir)*conjg(wf2%mtmesh(1:ntpll,iroffset+ir-1,ias,ist2))
            enddo
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, zfmt(1,iroffset),lmmaxvr) 
          enddo
          zpot=0d0
          call poisson_and_multipoles_mt_yukawa2( lmax, nrmt(is), spr(1:nrmt(is),is), zfmt, zpot, qlm(:,ias), is)



if(.false.)then
write(*,*) 'debug'
 open(11,file='mt_test.dat',status='replace')
 do ir=1, nrmt(1)
  write(11,*) spr(ir,is),dble(zfmt(1,ir)),dble(zpot(1,ir))  
 enddo
 close(11)
 stop
endif

          
          do iroffset=1,nrmt(is),blksize
            chunksize=blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset
            Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll, zpot(1,iroffset), lmmaxvr, zzero, mtmesh, ntpll)                        
            do ir=1,chunksize
              mtmesh(1:ntpll,ir)=mtmesh(1:ntpll,ir)*wf2%mtmesh(1:ntpll,iroffset+ir-1,ias,ist2)
            enddo
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, prod%mtrlm(1,iroffset,ias,1),lmmaxvr)
          enddo          
          
if(.false.)then
write(*,*) 'debug'
 open(11,file='mt_test.dat',status='replace')
 do ir=1, nrmt(1)
  write(11,*) spr(ir,1),dble(zpot(1,ir)), dble(prod%mtrlm(1,ir,1,1))
 enddo
 close(11)
 stop
endif
endif
 
        enddo
      enddo
      



! call timesec(tb)
!write(*,*) tb-ta
 
end subroutine WFprodcoul
