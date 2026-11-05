subroutine WFprodcoul2(ist1,wf1,ist2,wf2,prod,qlm)
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
      integer :: blkstart,chunksize,iroffset,LMoffset
      integer, parameter :: blksize=64
      complex(8), allocatable :: factors(:), rho(:,:), fr(:)
      complex(8) :: mtmesh(ntpll,blksize), mtmesh1(ntpll,blksize), mtrlm1(blksize,lmmaxvr)
      complex(8) :: zfmt(lmmaxvr,nrmtmax),zfmttr(nrmtmax,lmmaxvr),zpot(lmmaxvr,nrmtmax),qlm2(lmmaxvr)!,H3(lmmaxvr)
      real(8) :: refmttr(nrmtmax,lmmaxvr),imfmttr(nrmtmax,lmmaxvr)
      complex(8) :: zt,zt2
      real(8) :: ta,tb
      complex(8),allocatable :: H(:,:,:), HH(:,:,:)! , FF(:,:)! , F(:,:), H2(:,:)
!      complex(8) :: F(maxpbf,lmmaxvr),H2(maxpbf,lmmaxvr)
      complex(8) :: F(maxpbf,lmmaxvr),H2(maxpbf,lmmaxvr), FF(maxpbf2,lmmaxvr)
      real(8),allocatable :: reFF(:,:),imFF(:,:)
      real(8),external :: oldgaunt,oldwigner3j
 
      if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
      if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
 
! call timesec(ta)
qlm=0d0 !(:,ias)
      lmax= input%groundstate%lmaxvr

      allocate(H(lmmaxvr,maxradial,maxradial))
!      allocate(HH(lmmaxvr,maxradial2,maxradial))
!      allocate(H2(maxpbf,lmmaxvr))
!      allocate(H2(maxradial2,lmmaxvr))
!      allocate(F(maxpbf,lmmaxvr))
!      allocate(FF(maxpbf2,lmmaxvr))
!      allocate(reFF(maxpbf2,lmmaxvr))
!      allocate(imFF(maxpbf2,lmmaxvr))

      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)

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



          FF=0d0
          if2offset=0
          do irad2=1,nradial(is)
            l2=lrad(irad2,is)

            radoffset=0
            do l1=0,lmax

              do L=abs(l1-l2),min(lmax,l1+l2),2
                LMoffset=idxlm(L,-L)+L
                do LM=idxlm(L,-L),idxlm(L,-L)+2*L
                  H2(:,LM)=0d0
                enddo


                if2=if2offset
                do m2=-l2,l2
                  lm2=idxlm(l2,m2)
                  if2=if2+1

                  do m1=-l1,l1
                    lm1=idxlm(l1,m1)



                    zt=gntlyy(L,lm2,lm1)*wf2%mt(if2,ist2,ias)

                    M=m1+m2
                    if ((M.le.L).and.(M.ge.-L)) then
                      LM=idxlm(L,M)
!                      call zaxpy(npbf(l1,ias),zt,F(1,lm1),1,H2(1,LM),1)
                      do ipbf=1,npbf(l1,ias) !nradial2(ias)
                        H2(ipbf,LM)=H2(ipbf,LM)+zt*F(ipbf,lm1)
                      enddo
                      

                    endif 

                enddo

                              
              enddo

                do LM=idxlm(L,-L),idxlm(L,-L)+2*L
                  do jrad=1,npbf(l1,ias)
                    FF(1:npbf2(0,ias),LM)=FF(1:npbf2(0,ias),LM)+H2(jrad,LM)*uproducts2(1:npbf2(0,ias),0,radoffset+jrad,irad2,ias)
                  enddo
                enddo

        enddo
              radoffset=radoffset+npbf(l1,ias) !irad1
            enddo
            if2offset=if2

          enddo


!-----------------
 
          refmttr=0d0
          imfmttr=0d0

          do L=0,lmax
            do M=-L,L
              LM=idxlm(L,M)
!              do ipbf=1,npbf2(L,ias)
!                call daxpy(nrmt(is), dble(FF(ipbf,LM)),productbasis2(1,ipbf,L,ias),1,refmttr(1,LM),1)
!                call daxpy(nrmt(is),dimag(FF(ipbf,LM)),productbasis2(1,ipbf,L,ias),1,imfmttr(1,LM),1)
              do ipbf=1,npbf2(0,ias)
                call daxpy(nrmt(is), dble(FF(ipbf,LM)),productbasis2(1,ipbf,0,ias),1,refmttr(1,LM),1)
                call daxpy(nrmt(is),dimag(FF(ipbf,LM)),productbasis2(1,ipbf,0,ias),1,imfmttr(1,LM),1)
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

          
 
        enddo
      enddo
      
deallocate(H)
!deallocate(FF)
!deallocate(H2,F)


! call timesec(tb)
!write(*,*) tb-ta
 
end subroutine WFprodcoul2
