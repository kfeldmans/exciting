subroutine WFprodcoul3(ia,is,ist1,wf1,ist2,wf2,FF,qlm)
!subroutine WFprodcoul3(ist1,wf1,ist2,wf2,FF,qlm)
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
      use mod_eigensystem, only : gntyyy,gntyyyT,gntlyy,gntyyl,gntyly
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
      integer, intent(in) :: ist1,ist2 !,ia,is
      type (WFType) :: wf1,wf2
      complex(8), intent(out) :: qlm(lmmaxvr) !,natmtot)
      complex(8), intent(out) :: FF(maxpbf2,lmmaxvr) !,natmtot)
    
      integer :: is,ia
      integer :: ias
      integer :: l1,l2,m1,m2,lm1,lm2,io1,io2,ilo1,ilo2,lmmaxprod,lm,ir,l,m,ilo,io
      integer :: irad,jrad,irad1,irad2,ipbf,radoffset, padpbf2
      integer :: if1,if2,ifoffset1,ifoffset2,ifl2
      integer :: radoffset1, radoffset2
      integer :: lmax
      integer :: blkstart,chunksize,iroffset,LMoffset
      complex(8) :: zt,zt2
      real(8) :: ta,tb
!      complex(8),allocatable :: T(:,:,:)
      complex(8),allocatable :: T2(:,:)
      complex(8) :: F(maxpbf,lmmaxvr),U(maxpbf2,lmmaxvr), H(lmmaxvr,maxradial) ,G(maxpbf2)
      real(8) :: reG(maxpbf2),imG(maxpbf2)
 
 
! call timesec(ta)
!qlm=0d0 !(:,ias)
      
      lmax= input%groundstate%lmaxvr


!      allocate(T(maxpbf2,2*lmax+1,2*lmax+1))

          ias=idxas(ia,is)
          padpbf2=npbf2(0,ias)
          if ((padpbf2/4)*4.lt.padpbf2) padpbf2=(padpbf2/4+1)*4
          allocate(T2(padpbf2,(2*lmax+1)*(2*lmax+1)))

          F=0d0

          if2=1
          do irad2=1,nradial(is)
            H=0d0
            l2=lrad(irad2,is)
            ifoffset2=if2+l2

            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              zt=conjg(wf2%mtordered(ifoffset2+m2,ist2,ias))
              
              if1=1
              do irad1=1,nradial(is)
                l1=lrad(irad1,is)
                ifoffset1=if1+l1
!                lm1=idxlm(l1,-l1)+l1
                lm1=l1*(l1+1)+1

                do L=abs(l1-l2),min(lmax,l1+l2),2
!                  LM=idxlm(L,-L)+L-m2
!                  LM=L*(L+1)+1-m2
                  do m1=max(-l1,m2-L),min(l1,m2+L)
!                    H(LM+m1,irad1)=H(LM+m1,irad1)+gntyly(lm1+m1,L,lm2)*wf1%mtordered(ifoffset1+m1,ist1,ias)*zt
!                    H(L*(L+1)+1-m2+m1,irad1)=H(L*(L+1)+1-m2+m1,irad1)+gntyly(lm1+m1,L,lm2)*wf1%mtordered(ifoffset1+m1,ist1,ias)*zt
                    H(L*(L+1)+1-m2+m1,irad1)=H(L*(L+1)+1-m2+m1,irad1)+gntyly(lm1+m1,L,lm2)*wf1%mtordered(ifoffset1+m1,ist1,ias)*zt
                  enddo
                enddo
                if1=if1+2*l1+1
              enddo
!              H(:,1:nradial(is),irad2)=H(:,1:nradial(is),irad2)+G(:,1:nradial(is))*conjg(wf2%mtordered(ifoffset2+m2,ist2,ias))
            enddo
            if2=if2+2*l2+1

            do irad1=1,nradial(is)
              l1=lrad(irad1,is)
              do L=abs(l1-l2),min(lmax,l1+l2),2
!                LM=idxlm(L,-L)+L
!               LM=L*(L+1)+1
                do M=-L,L
                  F(1:npbf(L,ias),L*(L+1)+1+M)=F(1:npbf(L,ias),L*(L+1)+1+M)+H(L*(L+1)+1+M,irad1)*uproducts(1:npbf(L,ias),L,irad1,irad2,ias) 
                enddo
              enddo
            enddo
          enddo

!-----------------

          qlm=0d0
          do l=0,lmax
            do m=-l,l
              lm=idxlm(l,m)
              zt=0d0
              do ipbf=1,npbf(l,ias)
                qlm(lm)=qlm(lm)+pbfmultipoles(ipbf,l,ias)*F(ipbf,lm)
              enddo
            enddo
          enddo

!-----------------


          ifl2=1
          radoffset2=0
          do l2=0,lmax 
            lm2=idxlm(l2,-l2)+l2
            radoffset1=0
            do l1=0,lmax
              lm1=idxlm(l1,-l1)+l1
!              T(:,1:2*l2+1,1:2*l1+1)=0d0
              T2(:,1:(2*l1+1)*(2*l2+1))=0d0

              if2=ifl2
              do irad2=1,nlradial(l2,is)
                ifoffset2=if2+l2

                do m1=-l1,l1
                  reG=0d0
                  imG=0d0
                  do irad1=1,npbf(l1,ias)
                    ReG(1:padpbf2)=ReG(1:padpbf2)+dble(F(irad1,lm1+m1))*uproducts3(1:padpbf2,radoffset1+irad1,radoffset2+irad2,ias)
                    ImG(1:padpbf2)=ImG(1:padpbf2)+dimag(F(irad1,lm1+m1))*uproducts3(1:padpbf2,radoffset1+irad1,radoffset2+irad2,ias)
                  enddo
                  G(1:padpbf2)=dcmplx(ReG(1:padpbf2),ImG(1:padpbf2))
                  do m2=-l2,l2
!                    T(1:padpbf2,l2+m2+1,l1+m1+1)=T(1:padpbf2,l2+m2+1,l1+m1+1)+G(1:padpbf2)*wf2%mtordered(ifoffset2+m2,ist2,ias)
                    T2(1:padpbf2,l2+1+(l1+m1)*(2*l2+1)+m2)=T2(1:padpbf2,l2+1+(l1+m1)*(2*l2+1)+m2)+G(1:padpbf2)*wf2%mtordered(ifoffset2+m2,ist2,ias)
!                    T2(1:npbf2(0,ias),l2+1+(l1+m1)*(2*l2+1)+m2)=T2(1:npbf2(0,ias),l2+1+(l1+m1)*(2*l2+1)+m2)+G(1:npbf2(0,ias))*wf2%mtordered(ifoffset2+m2,ist2,ias)
                  enddo
                enddo
                if2=if2+2*l2+1
              enddo


              radoffset1=radoffset1+npbf(l1,ias)

              do L=abs(l1-l2),min(lmax,l1+l2),2
                do m1=-l1,l1
                  LM=L*(L+1)+m1+1
                  do m2=max(-l2,-m1-L),min(l2,-m1+L)
!                    FF(1:npbf2(0,ias),LM+m2)=FF(1:npbf2(0,ias),LM+m2)+T(1:npbf2(0,ias),l2+m2+1,l1+m1+1)*gntlyy(lm2+m2,lm1+m1,L)
!                    FF(1:npbf2(0,ias),LM+m2)=FF(1:npbf2(0,ias),LM+m2)+T2(1:npbf2(0,ias),l2+1+(l1+m1)*(2*l2+1)+m2)*gntlyy(lm2+m2,lm1+m1,L)
                    FF(1:padpbf2,LM+m2)=FF(1:padpbf2,LM+m2)+T2(1:padpbf2,l2+1+(l1+m1)*(2*l2+1)+m2)*gntlyy(lm2+m2,lm1+m1,L)
                  enddo
                enddo
              enddo

            enddo

            radoffset2=radoffset2+nlradial(l2,is)

            ifl2=if2

          enddo

      
!deallocate(T)
deallocate(T2)


! call timesec(tb)
!write(*,*) tb-ta
 
end subroutine WFprodcoul3
