subroutine expandinpb(FF,prod)
      use modinput
      use mod_APW_LO
      use mod_atoms
      use mod_muffin_tin
      use mod_eigenvalue_occupancy
      use constants, only : zzero, zone, fourpi
      use mod_SHT
      use mod_Gvector, only : ngrtot
      use mod_eigensystem, only : WFType
      use mod_radial
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
      complex(8), intent(in) :: FF(maxpbf2,lmmaxvr,natmtot)
      type (WFType) :: prod

      integer :: is,ia,ias
      integer :: lmax,lm,l,m,ir
      integer :: ipbf
      real(8) :: refmttr(nrmtmax,lmmaxvr),imfmttr(nrmtmax,lmmaxvr)
      real(8) :: ta,tb
 
      if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
      if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
 
! call timesec(ta)

      lmax= input%groundstate%lmaxvr

      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)

 
          refmttr=0d0
          imfmttr=0d0

          do L=0,lmax
            do M=-L,L
              LM=idxlm(L,M)
              do ipbf=1,npbf2(0,ias)
                call daxpy(nrmt(is), dble(FF(ipbf,LM,ias)),productbasis2(1,ipbf,0,ias),1,refmttr(1,LM),1)
                call daxpy(nrmt(is),dimag(FF(ipbf,LM,ias)),productbasis2(1,ipbf,0,ias),1,imfmttr(1,LM),1)
              enddo
            enddo
          enddo

          do ir= 1,nrmt(is)
            prod%mtrlm(1:lmmaxvr,ir,ias,1)=dcmplx(refmttr(ir,1:lmmaxvr),imfmttr(ir,1:lmmaxvr))
          enddo

        enddo
      enddo
      
! call timesec(tb)
!write(*,*) tb-ta
 
end subroutine expandinpb
