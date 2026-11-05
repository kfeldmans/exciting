!
!
!
!BOP
! !ROUTINE: WFRelease
! !INTERFACE:
!
!
subroutine orderWF(wf)
! !USES:
! !DESCRIPTION:
! Generates wave functions in the WFType representation from first- and second-variational eigenvectors.
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
!     use mod_kpoint, only : vkl
!     use mod_eigenvalue_occupancy, only : nstfv,nstsv
!     use mod_gkvector, only : ngk,vgkl,gkc,tpgkc,sfacgk,ngkmax
!     use mod_APW_LO, only : apwordmax
!     use mod_atoms, only : natmtot
!     use mod_muffin_tin, only : lmmaxapw
!     use modinput, only : input
 Use modinput
 Use mod_eigensystem
 use mod_kpoint
 use mod_eigenvalue_occupancy
 use mod_gkvector
 use mod_APW_LO
 use mod_atoms
 use mod_muffin_tin

 use modgw, only : kqset, Gkqset
 use mod_Gvector, only : ngrid, ngrtot, igfft
 Use mod_lattice, only : omega
 use constants, only : zzero, zone

!      use modmpi


implicit none
type (WFType) :: wf

integer :: ia,is,ias,if3,io,l,m,lm,wfsize,l1,l3,m1,m3,lm1,lm3,j1,j3,lmax,ilo
integer :: ifun, irad, iapw, napw(nspecies)
real (8) :: t1
integer :: ifg,igk,j
integer, allocatable :: tmpmap(:)

lmax=input%groundstate%lmaxapw
wfsize=wf%maxaa+wf%maxnlo
write(*,*) wfsize,wf%maxaa,wf%maxnlo
if (.not.allocated(wf%mtordered)) Allocate(wf%mtordered(wfsize,nstsv,natmtot))
allocate(tmpmap(wfsize))

Do is = 1, nspecies

  napw(is)=0
  do l=0,lmax
    napw(is)=napw(is)+apword(l,is)
  enddo


  irad=0
  ifun=1
  do l=0,lmax
    do io = 1, apword (l, is)
      irad=irad+1
      tmpmap(irad)=ifun
      ifun=ifun+2*l+1
    enddo
  enddo
  Do ilo = 1, nlorb (is)
    l = lorbl (ilo, is)
    irad=irad+1
    tmpmap(irad)=ifun
    ifun=ifun+2*l+1
  End Do

  Do ia = 1, natoms (is)
    ias = idxas (ia, is)

       iapw=0
       ifun=1
       do l=0,lmax
         do io = 1, apword (l, is)
           iapw=iapw+1
           wf%mtordered(ifun:ifun+2*l,:,ias)=wf%mt(tmpmap(iapw):tmpmap(iapw)+2*l,:,ias)
           ifun=ifun+2*l+1
         enddo
         Do ilo = 1, nlorb (is)
           if (l.eq.lorbl (ilo, is)) then
             wf%mtordered(ifun:ifun+2*l,:,ias)=wf%mt(tmpmap(napw(is)+ilo):tmpmap(napw(is)+ilo)+2*l,:,ias)
             ifun=ifun+2*l+1
           endif
         
         End Do
       enddo







  End Do
End Do

end subroutine orderWF

