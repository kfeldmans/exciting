
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> Overlap and Hamiltonian variables
Module mod_eigensystem
  implicit none
  !> order of overlap and Hamiltonian matrices for each k-point
  Integer, pointer :: nmat_ptr(:, :)
  Integer, Allocatable, target :: nmat(:, :)
  !> maximum nmat over all k-points
  integer, pointer :: nmatmax_ptr
  Integer, target :: nmatmax
  !> size of packed matrices
  Integer, Allocatable :: npmat(:, :)
  !> index to the position of the local-orbitals in the H and O matrices
  Integer, Allocatable :: idxlo(:, :, :)
  !> APW-local-orbital overlap integrals
  Real(8), Allocatable :: oalo(:, :, :)
  !> local-orbital-local-orbital overlap integrals
  Real(8), Allocatable :: ololo(:, :, :)
  !> APW-APW Hamiltonian integrals
  Real(8), Allocatable :: haa(:, :, :, :, :, :)
  Complex(8), Allocatable :: haaij(:, :, :)
  integer :: haaijSize
  !> local-orbital-APW Hamiltonian integrals
  Real(8), Allocatable :: hloa(:, :, :, :, :)
  Complex(8), Allocatable :: haloij(:, :, :)
  integer, allocatable :: haloijSize(:)
  !> local-orbital-local-orbital Hamiltonian integrals
  Real(8), Allocatable :: hlolo(:, :, :, :)
  Complex(8), Allocatable :: hloloij(:, :, :)
  ! The stuff for the linearized Koelling-Harmon
  !> APW-APW Hamiltonian integrals
  Real(8), Allocatable :: h1aa(:, :, :, :)
  !> local-orbital-APW Hamiltonian integrals
  Real(8), Allocatable :: h1loa(:, :, :)
  !> local-orbital-local-orbital Hamiltonian integrals
  Real(8), Allocatable :: h1lolo(:, :, :)
  logical :: h1on
  !> complex Gaunt coefficient array
  Complex(8), Allocatable :: gntyry(:, :, :), gntryy(:, :, :), gntnonz(:)
  real(8), allocatable :: gntyyy(:,:,:)
  !> list of non-zero Gaunt coefficients
  Integer, Allocatable :: gntnonzlm1(:), gntnonzlm2(:), gntnonzlm3(:), gntnonzlindex(:), gntnonzl2index(:, :)
  !> another compact list of Gaunt coefficients
  complex(8), allocatable :: listgnt(:, :, :)
  integer, allocatable :: indgnt(:, :, :)
  !gntnonzlist (:, :, :)
  !> tseqit is .true. if the first-variational secular equation is to be solved iteratively
  Logical :: tseqit
  !> number of secular equation iterations per self-consistent loop
  Integer :: nseqit
  !> iterative solver step length
  Real(8) :: tauseq
  !> ARPACK seed vector
  complex(8), allocatable :: arpackseed(:, :)
  !> pace
  complex(8), allocatable :: pace(:, :, :)

! singular components for the Davidson algorithm
  integer :: nsingular
  complex(8),allocatable :: singular(:,:,:)
  real(8),allocatable :: evalsingular(:,:)

  !> Matrix-elements for muffin-tin functions
  Type MTHamiltonianType
    complex(8), pointer :: aa(:, :, :), alo(:, :, :), loa(:, :, :), lolo(:, :, :)
  end type MTHamiltonianType

  !>  Muffin-tin Hamiltonian components
  Type MTHamiltonianList
    Type(MTHamiltonianType) :: main
    Type(MTHamiltonianType) :: spinless
    Type(MTHamiltonianType) :: alpha
    Type(MTHamiltonianType) :: beta
    Type(MTHamiltonianType) :: ab
    Type(MTHamiltonianType) :: ba
    integer :: maxnlo
    integer :: maxaa
    integer, allocatable :: losize(:)
  contains
    procedure :: release => MTRelease    !! Release all memory
  end type MTHamiltonianList

  !> Muffin-tin Hamlitonian for SCF
  Type(MTHamiltonianList) :: mt_hscf

! Wave functions in different representations
  Type WFType
    complex(8), allocatable :: mt(:,:,:)      !! coefficients of MT functions
    complex(8), allocatable :: mtrlm(:,:,:,:) !! expansion in spherical harmonics
    complex(8), allocatable :: mtmesh(:,:,:,:)!! values on a real-space mesh
    complex(8), allocatable :: gk(:,:)        !! coefficients of PWs
    complex(8), allocatable :: ir(:,:)        !! WFs on an FFT grid
    integer :: maxnlo                     !! maximum number of LOs over all species
    integer :: maxaa                      !! maximum number of augmenting functions for (L)APW
    integer,allocatable :: losize(:)      !! number of LOs for each species
  end type WFType


  !> Relativity settings
  integer :: level_nr, level_zora, level_iora
  parameter(level_nr=0, level_zora=1, level_iora=2)

Contains

!
!
!
!BOP
! !ROUTINE: MTNullify
! !INTERFACE:
!
!
  subroutine MTNullify(mt_h)
! !USES:
!      Use modinput
!      Use modmain
! !DESCRIPTION:
! Initialises all parts of the muffin-tin Hamiltonian or overlap.
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!   Modified Dezember 2020 (Ronaldo)
!EOP
!BOC
    Implicit None
    Type(MTHamiltonianList), Intent(Inout) :: mt_h

    nullify (mt_h%main%aa)
    nullify (mt_h%main%alo)
    nullify (mt_h%main%loa)
    nullify (mt_h%main%lolo)

    nullify (mt_h%spinless%aa)
    nullify (mt_h%spinless%alo)
    nullify (mt_h%spinless%loa)
    nullify (mt_h%spinless%lolo)

    nullify (mt_h%alpha%aa)
    nullify (mt_h%alpha%alo)
    nullify (mt_h%alpha%loa)
    nullify (mt_h%alpha%lolo)

    nullify (mt_h%beta%aa)
    nullify (mt_h%beta%alo)
    nullify (mt_h%beta%loa)
    nullify (mt_h%beta%lolo)

    nullify (mt_h%ab%aa)
    nullify (mt_h%ab%alo)
    nullify (mt_h%ab%loa)
    nullify (mt_h%ab%lolo)

    nullify (mt_h%ba%aa)
    nullify (mt_h%ba%alo)
    nullify (mt_h%ba%loa)
    nullify (mt_h%ba%lolo)

  end subroutine MTNullify

!
!
!
!BOP
! !FUNCTION: MaxAPWs
! !INTERFACE:
!
!
  integer Function MaxAPWs()
! !USES:
    Use modinput
    Use mod_APW_LO
    Use mod_atoms
    Use mod_muffin_tin
! !DESCRIPTION:
! Finds the maximum number of APWs across all muffin-tins.
!
! !REVISION HISTORY:
!   Created January 2021 (Andris)
!EOP
!BOC
    Implicit None
    Type(MTHamiltonianList) :: mt_h
    integer ::  io, if1, l, m, lm, maxsize, is, ias

    maxsize = 0
    Do is = 1, nspecies
      if1 = 0
      Do l = 0, input%groundstate%lmaxmat
        Do m = -l, l
          lm = idxlm(l, m)
          Do io = 1, apword(l, is)
            if1 = if1 + 1
          End Do
        End Do
      End Do
      maxsize = max(if1, maxsize)
    End do

    MaxAPWs = maxsize

  end function MaxAPWs

!
!
!
!BOP
! !ROUTINE: MTInitAll
! !INTERFACE:
!
!
  subroutine MTInitAll(mt_h)
! !USES:
    Use modinput
    Use mod_APW_LO
    Use mod_atoms
    Use mod_muffin_tin
! !DESCRIPTION:
! Initialises all parts of the muffin-tin Hamiltonian or overlap.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
    Implicit None
    Type(MTHamiltonianList) :: mt_h
    integer ::  io, ilo, if1, l, m, lm, l1, m1, lm1, l3, m3, lm3, is, ias
!      integer, external :: MaxAPWs

    mt_h%maxaa = MaxAPWs()

    if (allocated(mt_h%losize)) deallocate (mt_h%losize)
    allocate (mt_h%losize(nspecies))
    mt_h%maxnlo = 0
    Do is = 1, nspecies
      ias = idxas(1, is)
      ilo = nlorb(is)
      if (ilo .gt. 0) then
        l1 = lorbl(ilo, is)
        lm1 = idxlm(l1, l1)
        l3 = lorbl(1, is)
        lm3 = idxlm(l3, -l3)
        mt_h%losize(is) = idxlo(lm1, ilo, ias) - idxlo(lm3, 1, ias) + 1
        if (mt_h%maxnlo .lt. mt_h%losize(is)) mt_h%maxnlo = mt_h%losize(is)
      else
        mt_h%losize(is) = 0
      end if
    End do

    call MTInit(mt_h%spinless, mt_h%maxaa, mt_h%maxnlo)
!      call MTRedirect(mt_h%main,mt_h%alpha)

  end subroutine MTInitAll

!
!
!
!BOP
! !ROUTINE: MTInit
! !INTERFACE:
!
!
  subroutine MTInit(mt_block, maxaa, maxnlo)
! !USES:
    Use mod_atoms
!      Use modmain
! !DESCRIPTION:
! Initialises a part of the muffin-tin Hamiltonian or overlap.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
    Implicit None
    Type(MTHamiltonianType) :: mt_block
    integer, intent(in) :: maxaa, maxnlo

    if (.not. associated(mt_block%aa)) then
      allocate (mt_block%aa(maxaa, maxaa, natmtot))
      mt_block%aa = 0d0
      if (maxnlo .gt. 0) then
        allocate (mt_block%loa(maxnlo, maxaa, natmtot))
        allocate (mt_block%alo(maxaa, maxnlo, natmtot))
        allocate (mt_block%lolo(maxnlo, maxnlo, natmtot))
        mt_block%loa = 0d0
        mt_block%alo = 0d0
        mt_block%lolo = 0d0
      end if
    else
      mt_block%aa = 0d0
      if (maxnlo .gt. 0) then
        mt_block%loa = 0d0
        mt_block%alo = 0d0
        mt_block%lolo = 0d0
      end if

    end if

  end subroutine MTInit

!
!
!
!BOP
! !ROUTINE: MTRedirect
! !INTERFACE:
!
!
  subroutine MTRedirect(mt_blockA, mt_blockB)
! !USES:
    Use mod_atoms
!      Use modmain
! !DESCRIPTION:
! Redirects pointers of mt\_blockA to mt\_blockB of muffin-tin Hamiltonians or overlaps.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
    Implicit None
    Type(MTHamiltonianType), intent(in)  :: mt_blockB
    Type(MTHamiltonianType), intent(out) :: mt_blockA

    if (associated(mt_blockB%aa)) mt_blockA%aa => mt_blockB%aa
    if (associated(mt_blockB%alo)) mt_blockA%alo => mt_blockB%alo
    if (associated(mt_blockB%loa)) mt_blockA%loa => mt_blockB%loa
    if (associated(mt_blockB%lolo)) mt_blockA%lolo => mt_blockB%lolo

  end subroutine MTRedirect
!
!
!
!BOP
! !ROUTINE: MTCopy
! !INTERFACE:
!
!
  subroutine MTCopy(mt_blockA, mt_blockB)
! !USES:
    use mod_atoms
    use modmpi, only: terminate_mpi_env, mpiglobal

!      Use modmain
! !DESCRIPTION:
! Copies the contents of mt\_blockA to mt\_blockB of muffin-tin Hamiltonians or overlaps.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
    Implicit None
    Type(MTHamiltonianType) :: mt_blockA, mt_blockB
    integer :: maxaa, maxnlo

    ! Check if mt_blockA%aa is associated to some target
    if (.not. associated(mt_blockA%aa)) then
      call terminate_mpi_env(mpiglobal, &
        & 'Error (MTCopy): mt_blockA%aa is not allocated')
    end if
    ! Check if mt_blockB%aa is associated to some target
    if (.not. associated(mt_blockB%aa)) then
      call terminate_mpi_env(mpiglobal, &
        & 'Error (MTCopy): mt_blockB%aa is not allocated')
    end if

    ! Now we can safely make mt_blockB%aa points to the same target as mt_blockA%aa
    mt_blockB%aa = mt_blockA%aa

    ! We will repeat this procedure below

    ! This first check is to verify if there are any local orbitals
    if (associated(mt_blockA%lolo)) then
      if (.not. associated(mt_blockA%alo)) then
        call terminate_mpi_env(mpiglobal, &
          & 'Error (MTCopy): mt_blockA%alo is not allocated')
      end if
      if (.not. associated(mt_blockB%alo)) then
        call terminate_mpi_env(mpiglobal, &
          & 'Error (MTCopy): mt_blockB%alo is not allocated')
      end if
      mt_blockB%alo = mt_blockA%alo

      if (.not. associated(mt_blockA%loa)) then
        call terminate_mpi_env(mpiglobal, &
          & 'Error (MTCopy): mt_blockA%loa is not allocated')
      end if
      if (.not. associated(mt_blockB%loa)) then
        call terminate_mpi_env(mpiglobal, &
          & 'Error (MTCopy): mt_blockB%loa is not allocated')
      end if
      mt_blockB%loa = mt_blockA%loa
      if (.not. associated(mt_blockB%lolo)) then
        call terminate_mpi_env(mpiglobal, &
          & 'Error (MTCopy): mt_blockB%lolo is not allocated')
      end if
      mt_blockB%lolo = mt_blockA%lolo
    end if
  end subroutine MTCopy


  !> Release all memory used by the muffin-tin Hamiltonian or overlap
  subroutine MTRelease(this)
    implicit none
    class(MTHamiltonianList), Intent(Inout) :: this

    if (associated(this%spinless%aa)) deallocate (this%spinless%aa)
    if (associated(this%spinless%alo)) deallocate (this%spinless%alo)
    if (associated(this%spinless%loa)) deallocate (this%spinless%loa)
    if (associated(this%spinless%lolo)) deallocate (this%spinless%lolo)

    if (associated(this%alpha%aa)) deallocate (this%alpha%aa)
    if (associated(this%alpha%alo)) deallocate (this%alpha%alo)
    if (associated(this%alpha%loa)) deallocate (this%alpha%loa)
    if (associated(this%alpha%lolo)) deallocate (this%alpha%lolo)

    if (associated(this%beta%aa)) deallocate (this%beta%aa)
    if (associated(this%beta%alo)) deallocate (this%beta%alo)
    if (associated(this%beta%loa)) deallocate (this%beta%loa)
    if (associated(this%beta%lolo)) deallocate (this%beta%lolo)

    if (associated(this%ab%aa)) deallocate (this%ab%aa)
    if (associated(this%ab%alo)) deallocate (this%ab%alo)
    if (associated(this%ab%loa)) deallocate (this%ab%loa)
    if (associated(this%ab%lolo)) deallocate (this%ab%lolo)

    if (associated(this%ba%aa)) deallocate (this%ba%aa)
    if (associated(this%ba%alo)) deallocate (this%ba%alo)
    if (associated(this%ba%loa)) deallocate (this%ba%loa)
    if (associated(this%ba%lolo)) deallocate (this%ba%lolo)

    if (allocated(this%losize)) deallocate (this%losize)

  end subroutine MTRelease

!
!
!
!BOP
! !ROUTINE: WFInit
! !INTERFACE:
!
!
  subroutine WFInit(wf)
! !USES:
! !DESCRIPTION:
! Initialised an WFType datastructure
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
    Use modinput
    Use mod_APW_LO
    Use mod_atoms
    Use mod_muffin_tin

    implicit none
    type(WFType) :: wf

    integer :: haaijsize
    integer :: if1,l,m,lm,is,ias,io,ilo,l1,l3,lm1,lm3


    haaijSize=0
    Do is = 1, nspecies
      if1=0
      Do l = 0, input%groundstate%lmaxmat
        Do m = - l, l
          lm = idxlm (l, m)
          Do io = 1, apword (l, is)
            if1=if1+1
          End Do
        End Do
      End Do
      if (if1.gt.haaijSize) haaijSize=if1
    Enddo
    wf%maxaa=haaijSize


    if (allocated(wf%losize)) deallocate(wf%losize)
    allocate(wf%losize(nspecies))
    wf%maxnlo=0
    Do is = 1, nspecies
      ias=idxas (1, is)
      ilo=nlorb (is)
      if (ilo.gt.0) then
        l1 = lorbl (ilo, is)
        lm1 = idxlm (l1, l1)
        l3 = lorbl (1, is)
        lm3 = idxlm (l3, -l3)
        wf%losize(is)=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
        if (wf%maxnlo.lt.wf%losize(is)) wf%maxnlo=wf%losize(is)
      else
        wf%losize(is)=0
      endif
    Enddo

!     nullify(wf%mt)
!     nullify(wf%gk)
!     nullify(wf%ir)
!     nullify(wf%mtrlm)
!     nullify(wf%mtmesh)

  end subroutine WFInit


!
!
!
!BOP
! !ROUTINE: WFRelease
! !INTERFACE:
!
!
  subroutine WFRelease(wf)
! !USES:
! !DESCRIPTION:
! Release all memory used by wave functions in the WFType representation
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!EOP
!BOC
    implicit none
    type(WFType) :: wf

    if (allocated(wf%mt)) then
      deallocate(wf%mt)
    endif
!     nullify(wf%mt)
    if (allocated(wf%gk)) then
      deallocate(wf%gk)
    endif
!     nullify(wf%gk)
    if (allocated(wf%ir)) then
      deallocate(wf%ir)
    endif
!     nullify(wf%ir)
    if (allocated(wf%mtrlm)) then
      deallocate(wf%mtrlm)
    endif
!     nullify(wf%mtrlm)
    if (allocated(wf%mtmesh)) then
      deallocate(wf%mtmesh)
    endif
!     nullify(wf%mtmesh)

  end subroutine WFRelease

!
!
!
!BOP
! !ROUTINE: genWFinMT
! !INTERFACE:
!
!
  subroutine genWFinMT(wf)
    use modinput
    use mod_APW_LO
    use mod_atoms
    use mod_muffin_tin
    use mod_eigenvalue_occupancy
! !USES:
! !DESCRIPTION:
! Evaluates WF in terms of spherical harmonics
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
    implicit none
    type (WFType) :: wf
    integer :: j,is,ia,ias,ir,l,m,lm,if1,ilo,io

    if (.not.allocated(wf%mtrlm)) allocate(wf%mtrlm(lmmaxvr,nrmtmax,natmtot,nstsv))
    wf%mtrlm=0d0
    do j=1,nstsv
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          if1=0
! APW part
          do l=0,input%groundstate%lmaxvr
            do io = 1, apword (l, is)
              do m=-l,l
                lm=idxlm(l,m)
                if1=if1+1
                wf%mtrlm(lm,1:nrmt(is),ias,j)=wf%mtrlm(lm,1:nrmt(is),ias,j)+wf%mt(if1,j,ias)*apwfr(1:nrmt(is),1,io,l,ias)
              enddo
            enddo
          enddo

! local-orbital functions
          Do ilo = 1, nlorb (is)
            l = lorbl (ilo, is)
            Do m = - l, l
              if1=if1+1
              lm = idxlm (l, m)
              wf%mtrlm(lm,1:nrmt(is),ias,j)=wf%mtrlm(lm,1:nrmt(is),ias,j)+wf%mt(if1,j,ias)*lofr(1:nrmt(is),1,ilo,ias)
            End Do
          End Do

        enddo
      enddo
    enddo




  end subroutine genWFinMT



!
!
!
!BOP
! !ROUTINE: genWFonMesh
! !INTERFACE:
!
!
  subroutine genWFonMesh(wf)
    use modinput
    use mod_APW_LO
    use mod_atoms
    use mod_muffin_tin
    use mod_eigenvalue_occupancy
    use constants, only : zzero, zone
    use mod_SHT
! !USES:
! !DESCRIPTION:
! Evaluates WF on a real-space mesh
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
    implicit none
    type (WFType) :: wf
    integer :: j,is,ia,ias

    if (.not.allocated(wf%mtmesh)) allocate(wf%mtmesh(ntpll,nrmtmax,natmtot,nstsv))

    do j=1,nstsv
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          Call zgemm ('N', 'N', ntpll, nrmt(is), lmmaxvr, &
            & zone, zbshthf, ntpll, wf%mtrlm(1,1,ias,j), lmmaxvr, zzero, &
            & wf%mtmesh(1,1,ias,j), ntpll)
        enddo
      enddo
    enddo
  end subroutine genWFonMesh

!
!
!
!BOP
! !ROUTINE: genWFonMeshOne
! !INTERFACE:
!
!
  subroutine genWFonMeshOne(wf)
    use modinput
    use mod_APW_LO
    use mod_atoms
    use mod_muffin_tin
    use mod_eigenvalue_occupancy
    use constants, only : zzero, zone
    use mod_SHT
! !USES:
! !DESCRIPTION:
! Evaluates WF on a real-space mesh
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
    implicit none
    type (WFType) :: wf
    integer :: is,ia,ias

    if (.not.allocated(wf%mtmesh)) allocate(wf%mtmesh(ntpll,nrmtmax,natmtot,1)) ! was lmmaxvr instead of lmmaxhf/ntpll

    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        Call zgemm ('N', 'N', ntpll, nrmt(is), lmmaxvr, &
          & zone, zbshthf, ntpll, wf%mtrlm(1,1,ias,1), lmmaxvr, zzero, &
          & wf%mtmesh(1,1,ias,1), ntpll) ! Genshtmat3

        !  Call zgemm ('N', 'N', lmmaxhf, nrmt(is), lmmaxvr, &
        ! & zone, zbshthf, lmmaxhf, wf%mtrlm(1,1,ias,1), lmmaxvr, zzero, &
        ! & wf%mtmesh(1,1,ias,1), lmmaxhf) ! Genshtmat2

        ! Call zgemm ('N', 'N', lmmaxvr, nrmt(is), lmmaxvr, &
        ! & zone, zbshtvr, lmmaxvr, wf%mtrlm(1,1,ias,1), lmmaxvr, zzero, &
        ! & wf%mtmesh(1,1,ias,1), lmmaxvr) ! Genshtmat
      enddo
    enddo
  end subroutine genWFonMeshOne

!
!
!
!BOP
! !ROUTINE: WFprod
! !INTERFACE:
!
!
  subroutine WFprod(ist1,wf1,ist2,wf2,prod)
    use modinput
    use mod_APW_LO
    use mod_atoms
    use mod_muffin_tin
    use mod_eigenvalue_occupancy
    use constants, only : zzero, zone
    use mod_SHT
    use mod_Gvector, only : ngrtot
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
    integer :: is,ia,ias
    integer :: l1,l3,m1,m3,lm1,lm3,lm2,io1,io2,if1,if3,if1old,if3old,ilo1,ilo2,lmmaxprod
    complex(8), allocatable :: factors(:), rho(:,:), fr(:)
    complex(8) :: zt
    real(8) :: ta,tb




    if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
    if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))


    call timesec(ta)

    lmmaxprod= lmmaxvr
    prod%mtrlm=0d0
    allocate(factors(lmmaxvr))
    allocate(rho(nrmtmax,lmmaxvr))
    allocate(fr(nrmtmax))

!do if1=1,wf1%maxaa
!
! write(*,*)
!
!enddo


    Do is = 1, nspecies
!        n = lmmaxvr * nrmt (is)
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        rho=0d0

! APW-APW part
!if (.true.) then
        if1=0
        Do l1 = 0, input%groundstate%lmaxapw
          Do io1 = 1, apword (l1, is)
            if3=0
            if1old=if1

            Do l3 = 0, input%groundstate%lmaxapw
              Do io2 = 1, apword (l3, is)
                fr(:)=apwfr (:, 1, io1, l1, ias) * apwfr (:, 1, io2, l3, ias)
                if1=if1old
                if3old=if3
                factors=0d0
                Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  if1=if1+1
                  if3=if3old
                  Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    if3=if3+1
!                      i=1
!                      do while(indgnt(i,lm3,lm1).ne.0)
!                        lm2=indgnt(i,lm3,lm1)
                    zt=conjg(wf1%mt(if1,ist1,ias))*wf2%mt(if3,ist2,ias)
                    do lm2=1,lmmaxprod
                      !                       if (dble(gntyyy(lm1,lm2,lm3)*conjg(gntyyy(lm1,lm2,lm3))).gt.1d-10) then
                      factors(lm2) = factors(lm2) + zt*gntyyy(lm2,lm1,lm3)
!                        endif
                    enddo
!                        if (lm2.le.lmmaxvr) factorsnew(lm2,1)=factorsnew(lm2,1)+dble(mt_dm%main%ff(if1,if3,ias)*conjg(listgnt(i,lm3,lm1)))
!                        i=i+1
!                      enddo
                  enddo
                enddo

                do lm2=1,lmmaxprod
                  if (factors(lm2).ne.0d0) then
!                        rho(1,lm2) = rho(1,lm2) + fr(1)*factors(lm2)
                    rho(1:nrmt(is),lm2) = rho(1:nrmt(is),lm2) + fr(1:nrmt(is))*factors(lm2)
                  endif

                enddo

              End Do
            End Do


          End Do
        End Do

        if (wf1%losize(is).gt.0) then
!APW-LO part

          if1=0
          Do l1 = 0, input%groundstate%lmaxapw
            Do io1 = 1, apword (l1, is)
              if3=0
              if1old=if1
              Do ilo2 = 1, nlorb (is)
                l3 = lorbl (ilo2, is)
                fr(:)=apwfr (:, 1, io1, l1, ias) * lofr (:, 1, ilo2, ias)
                if1=if1old
                if3old=if3
                factors=0d0
                Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  if1=if1+1
                  if3=if3old
                  Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    if3=if3+1

!                      i=1
!                      do lm2=1,lmmaxvr
!                        factorsnew(lm2,1)=factorsnew(lm2,1)+2d0*conjg(gntryy(lm2,lm1,lm3))*dble(mt_dm%main%ff(if1,maxaa+if3,ias))
!                      enddo
                    zt=conjg(wf1%mt(if1,ist1,ias))*wf2%mt(wf2%maxaa+if3,ist2,ias)
                    do lm2=1,lmmaxprod
!                        if (dble(gntyyy(lm1,lm2,lm3)*conjg(gntyyy(lm1,lm2,lm3))).gt.1d-10) then
                      factors(lm2) = factors(lm2) + zt*gntyyy(lm2,lm1,lm3)
!                        endif
                    enddo
                    zt=conjg(wf1%mt(wf1%maxaa+if3,ist1,ias))*wf2%mt(if1,ist2,ias)
                    do lm2=1,lmmaxprod
!                        if (dble(gntyyy(lm1,lm2,lm3)*conjg(gntyyy(lm1,lm2,lm3))).gt.1d-10) then
                      factors(lm2) = factors(lm2) + zt*gntyyy(lm2,lm3,lm1)
!                        endif
                    enddo

!                      do while(indgnt(i,lm3,lm1).ne.0)
!                        lm2=indgnt(i,lm3,lm1)
!                        if (lm2.le.lmmaxvr) factorsnew(lm2,1)=factorsnew(lm2,1)+2d0*dble(mt_dm%main%ff(if1,maxaa+if3,ias)*conjg(listgnt(i,lm3,lm1)))
!                        i=i+1
!                      enddo


                  enddo
                enddo

                do lm2=1,lmmaxprod
                  if (factors(lm2).ne.0d0) then
                    rho(1:nrmt(is),lm2) = rho(1:nrmt(is),lm2) + fr(1:nrmt(is))*factors(lm2)
!                        rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+frnew(1:nrmt(is),1)*factorsnew(lm2,1)
                  endif

                enddo
!                  do lm2=1,lmmaxvr
!                    if (factorsnew(lm2,1).ne.0d0) then
!                      rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+fr(1:nrmt(is),1)*factorsnew(lm2,1)
!                    endif
!                  enddo

              End Do
            End Do
          End Do


!LO-LO part
          if1=0
          Do ilo1 = 1, nlorb (is)
            l1 = lorbl (ilo1, is)
            if3=0
            if1old=if1
            Do ilo2 = 1, nlorb (is)
              l3 = lorbl (ilo2, is)
              fr(:)=lofr (:, 1, ilo1, ias) * lofr (:, 1, ilo2, ias)
              if1=if1old
              if3old=if3
              factors=0d0
              Do m1 = - l1, l1
                lm1 = idxlm (l1, m1)
                if1=if1+1
                if3=if3old
                Do m3 = - l3, l3
                  lm3 = idxlm (l3, m3)
                  if3=if3+1
!                      i=1
                  zt=conjg(wf1%mt(wf1%maxaa+if1,ist1,ias))*wf2%mt(wf2%maxaa+if3,ist2,ias)
                  do lm2=1,lmmaxprod
!                        if (dble(gntyyy(lm1,lm2,lm3)*conjg(gntyyy(lm1,lm2,lm3))).gt.1d-10) then
                    factors(lm2) = factors(lm2) + zt*gntyyy(lm2,lm1,lm3)
!                        endif
                  enddo
!                      do while(indgnt(i,lm3,lm1).ne.0)
!                        lm2=indgnt(i,lm3,lm1)
!                        if (lm2.le.lmmaxvr) factorsnew(lm2,1)=factorsnew(lm2,1)+dble(mt_dm%main%ff(maxaa+if1,maxaa+if3,ias)*conjg(listgnt(i,lm3,lm1)))
!                        i=i+1
!                      enddo

                enddo
              enddo

              do lm2=1,lmmaxprod
                if (factors(lm2).ne.0d0) then
                  rho(1:nrmt(is),lm2) = rho(1:nrmt(is),lm2) + fr(1:nrmt(is))*factors(lm2)
!                        rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+frnew(1:nrmt(is),1)*factorsnew(lm2,1)
                endif
              enddo
!                do lm2=1,lmmaxvr
!                  if (factorsnew(lm2,1).ne.0d0) then
!                    rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+frnew(1:nrmt(is),1)*factorsnew(lm2,1)
!                  endif
!                enddo

            End Do
          End Do

        endif




        do lm2=1,lmmaxvr
          prod%mtrlm(lm2,1:nrmt(is),ias,1)=rho(1:nrmt(is),lm2)
        enddo

      End Do
    End Do


    deallocate(fr,factors)
    deallocate(rho)


    call timesec(tb)
!write(*,*) tb-ta

    prod%ir(:,1)=conjg(wf1%ir(:,ist1))*wf2%ir(:,ist2)
  end subroutine WFprod


!
!
!
!BOP
! !ROUTINE: WFprodrs
! !INTERFACE:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THE ORIGINAL exciting-neon WFprodrs subroutine BEGINS HERE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine WFprodrs(ist1,wf1,ist2,wf2,prod)
!      use modinput
!      use mod_APW_LO
!      use mod_atoms
!      use mod_muffin_tin
!      use mod_eigenvalue_occupancy
!      use constants, only : zzero, zone
!      use mod_SHT
!      use mod_Gvector, only : ngrtot
! ! !USES:
! ! !DESCRIPTION:
! ! Evaluates a product of two WFs in the real space
! !
! ! !REVISION HISTORY:
! !   Created 2021 (Andris)
! !EOP
! !BOC
!      implicit none
!      integer, intent(in) :: ist1,ist2
!      type (WFType) :: wf1,wf2,prod
!      integer :: is,ia,ias
!      integer :: l1,l3,m1,m3,lm1,lm3,lm2,io1,io2,if1,if3,if1old,if3old,ilo1,ilo2,lmmaxprod
!      complex(8), allocatable :: factors(:), rho(:,:), fr(:)
!      complex(8) :: zt
!      real(8) :: ta,tb

!      if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
!      if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
! call timesec(ta)
!      if (.not.allocated(prod%mtmesh)) allocate(prod%mtmesh(ntpll,nrmtmax,natmtot,1))
!      prod%mtmesh(:,:,:,1)=conjg(wf1%mtmesh(:,:,:,ist1))*wf2%mtmesh(:,:,:,ist2)

!      do is=1,nspecies
!        do ia=1,natoms(is)
!          ias=idxas(ia,is)
!         !  Call zgemm ('N', 'N', lmmaxvr, nrmt(is), lmmaxvr, zone, zfshtvr, lmmaxvr, prod%mtmesh(1,1,ias,1), lmmaxvr, zzero, prod%mtrlm(1,1,ias,1) , lmmaxvr) ! Genshtmat
!         !  Call zgemm ('N', 'N', lmmaxvr, nrmt(is), lmmaxhf, zone, zfshthf, lmmaxhf, prod%mtmesh(1,1,ias,1), lmmaxhf, zzero, prod%mtrlm(1,1,ias,1) , lmmaxvr) ! Genshtmat2
!          Call zgemm ('N', 'N', lmmaxvr, nrmt(is), ntpll, zone, zfshthf, lmmaxvr, prod%mtmesh(1,1,ias,1), ntpll, zzero, prod%mtrlm(1,1,ias,1) , lmmaxvr) ! Genshtmat3
!        enddo
!      enddo

!      deallocate(prod%mtmesh)
! call timesec(tb)
! !write(*,*) tb-ta

!      !prod%ir(:,1)=conjg(wf1%ir(:,ist1))*wf2%ir(:,ist2)
!      end subroutine WFprodrs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THE ORIGINAL exciting-neon WFprodrs subroutine ENDS HERE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VERSION OF WFprodrs subroutine that works faster STARTS HERE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine WFprodrs(ist1,wf1,ist2,wf2,prod)
    use modinput
    use mod_APW_LO
    use mod_atoms
    use mod_muffin_tin
    use mod_eigenvalue_occupancy
    use constants, only : zzero, zone
    use mod_SHT
    use mod_Gvector, only : ngrtot
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
    integer :: is,ia,ias
    integer :: l1,l3,m1,m3,lm1,lm3,lm2,io1,io2,if1,if3,if1old,if3old,ilo1,ilo2,lmmaxprod,lm,ir
    integer :: blkstart,chunksize,iroffset
    integer, parameter :: blksize=64
    complex(8), allocatable :: factors(:), rho(:,:), fr(:), mtmesh(:,:), angmesh(:)
    real(8), allocatable :: REzfshthf(:,:), IMzfshthf(:,:), REmtmesh(:,:), IMmtmesh(:,:),TMPzfshthf(:,:),TMPmtmesh(:,:)
    real(8) :: REresult(lmmaxvr,blksize),IMresult(lmmaxvr,blksize),TMPresult(lmmaxvr,blksize)
    complex(8) :: zt
    real(8) :: ta,tb

    if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
    if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))

    call timesec(ta)

    allocate(mtmesh(ntpll,blksize))
    allocate(angmesh(ntpll))

! needed for option 3 (below)
    if (.true.) then
      allocate(REmtmesh(ntpll,blksize))
      allocate(IMmtmesh(ntpll,blksize))
      allocate(TMPmtmesh(ntpll,blksize))

      allocate(REzfshthf(lmmaxvr,ntpll))
      allocate(IMzfshthf(lmmaxvr,ntpll))
      allocate(TMPzfshthf(lmmaxvr,ntpll))

      REzfshthf=dble(zfshthf)
      IMzfshthf=dimag(zfshthf)
      TMPzfshthf=REzfshthf+IMzfshthf
    endif

    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)

! Pick your poison

! Option 1: zgemv to transform to spherical harmonics shell by shell.
! Bad because of zgemv.
        if (.false.) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ir,lm,angmesh)
!$OMP DO
          do ir=1,nrmt(is)
            do lm=1,ntpll
              angmesh(lm)=conjg(wf1%mtmesh(lm,ir,ias,ist1))*wf2%mtmesh(lm,ir,ias,ist2)
            enddo
            call zgemv('N', lmmaxvr, ntpll, zone, zfshthf, lmmaxvr, angmesh, 1, zzero, prod%mtrlm(1,ir,ias,1), 1)
          enddo
!$OMP END DO
!$OMP END PARALLEL
        endif
! Option 2: blocking allows to maintain the memory access locality and to replace zgemv with zgemm
! Slightly longer code because of blocking, but it's faster than option 1
        if (.false.) then
          chunksize=blksize
          do iroffset=1,nrmt(is),blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset
            do ir=iroffset,iroffset+chunksize-1
              do lm=1,ntpll
                mtmesh(lm,ir-iroffset+1)=conjg(wf1%mtmesh(lm,ir,ias,ist1))*wf2%mtmesh(lm,ir,ias,ist2)
              enddo
            enddo
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, prod%mtrlm(1,iroffset,ias,1) , lmmaxvr) ! Genshtmat3
          enddo

        endif
! Option 3: keep the blocking introduced above, but replace zgemm with 3 dgemms.
! 3 dgemms are faster than 1 zgemm, but additional memory transfers have extra cost.
! Option 3 beats option 2 on a single core, but does it on all cores?
! Again, even longer code.
        if (.true.) then
          chunksize=blksize
          do iroffset=1,nrmt(is),blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset
            do ir=iroffset,iroffset+chunksize-1
              do lm=1,ntpll
                zt=conjg(wf1%mtmesh(lm,ir,ias,ist1))*wf2%mtmesh(lm,ir,ias,ist2)
                REmtmesh(lm,ir-iroffset+1)=dble(zt)
                IMmtmesh(lm,ir-iroffset+1)=dimag(zt)
                TMPmtmesh(lm,ir-iroffset+1)=REmtmesh(lm,ir-iroffset+1)+IMmtmesh(lm,ir-iroffset+1)
              enddo
            enddo

            Call dgemm ('N', 'N', lmmaxvr, chunksize, ntpll, 1d0, REzfshthf, lmmaxvr, REmtmesh, ntpll, 0d0, REresult , lmmaxvr)
            Call dgemm ('N', 'N', lmmaxvr, chunksize, ntpll, 1d0, IMzfshthf, lmmaxvr, IMmtmesh, ntpll, 0d0, TMPresult , lmmaxvr)
            Call dgemm ('N', 'N', lmmaxvr, chunksize, ntpll, 1d0, TMPzfshthf, lmmaxvr, TMPmtmesh, ntpll, 0d0, IMresult , lmmaxvr)

            IMresult=IMresult-REresult
            IMresult=IMresult-TMPresult
            REresult=REresult-TMPresult

            prod%mtrlm(1:lmmaxvr,iroffset:iroffset+chunksize-1,ias,1)=dcmplx(REresult(1:lmmaxvr,1:chunksize),IMresult(1:lmmaxvr,1:chunksize))

          enddo
        endif

      enddo
    enddo
    deallocate(mtmesh,angmesh)
    deallocate(REzfshthf,IMzfshthf,REmtmesh,IMmtmesh,TMPmtmesh,TMPzfshthf)

    call timesec(tb)
    !write(*,*) tb-ta

    !prod%ir(:,1)=conjg(wf1%ir(:,ist1))*wf2%ir(:,ist2)
  end subroutine WFprodrs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VERSION OF WFprodrs subroutine that works faster ENDS HERE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prodshrs(zfun1,zfun2,zres)
    use modinput
    use mod_APW_LO
    use mod_atoms
    use mod_muffin_tin
    use mod_eigenvalue_occupancy
    use constants, only : zzero, zone
    use mod_SHT
    use mod_Gvector, only : ngrtot
    ! !USES:
    ! !DESCRIPTION:
    ! Calculates a product of two complex functions given in a muffin tin.
    ! The first and the second function are given in spherical harmonics and on a spherical grid, respectively.
    ! The output is produced in spherical harmonics.
    ! !REVISION HISTORY:
    !   Created 2024 (Andris)
    !EOP
    !BOC
    implicit none
    complex(8), intent(in) :: zfun1(lmmaxvr,nrmtmax,natmtot)
    complex(8), intent(in) :: zfun2(ntpll,nrmtmax,natmtot)
    complex(8), intent(out) :: zres(lmmaxvr,nrmtmax,natmtot)



!     integer, intent(in) :: ist1,ist2
!     type (WFType) :: wf1,wf2,prod
    integer :: is,ia,ias
    integer :: ir,lm
    integer :: blkstart,chunksize,iroffset
    integer, parameter :: blksize=64
    complex(8) :: mtmesh(ntpll,blksize)

    real(8) :: Ash(lmmaxvr,blksize),Bsh(lmmaxvr,blksize),Csh(lmmaxvr,blksize)
    real(8) :: Asc(ntpll,blksize),Bsc(ntpll,blksize),Csc(ntpll,blksize)
    real(8) :: Afsht(lmmaxvr,ntpll),Bfsht(lmmaxvr,ntpll),Cfsht(lmmaxvr,ntpll)
    real(8) :: Absht(ntpll,lmmaxvr),Bbsht(ntpll,lmmaxvr),Cbsht(ntpll,lmmaxvr)
    real(8) :: ref2, imf2

!,prodmesh
!      real(8), allocatable :: REzfshthf(:,:), IMzfshthf(:,:), REmtmesh(:,:), IMmtmesh(:,:),TMPzfshthf(:,:),TMPmtmesh(:,:)
!      real(8) :: REresult(lmmaxvr,blksize),IMresult(lmmaxvr,blksize),TMPresult(lmmaxvr,blksize)
!      complex(8) :: zt
!      real(8) :: ta,tb


!      allocate(mtmesh(ntpll,blksize))

    Afsht=dble(zfshthf)
    Bfsht=dimag(zfshthf)
    Cfsht=Afsht+Bfsht

    Absht=dble(zbshthf)
    Bbsht=dimag(zbshthf)
    Cbsht=Absht+Bbsht


    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! transform function zfun1 into the real space
        chunksize=blksize
        do iroffset=1,nrmt(is),blksize
          if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset

          if (.false.) then
            Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll, zfun1(1,iroffset,ias), lmmaxvr, zzero, mtmesh, ntpll)
            mtmesh(:,1:chunksize)=mtmesh(:,1:chunksize)*zfun2(:,iroffset:iroffset+chunksize-1,ias)
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, zres(1,iroffset,ias) , lmmaxvr)
          else
            Ash(1:lmmaxvr,1:chunksize)=dble(zfun1(1:lmmaxvr,iroffset:iroffset+chunksize-1,ias))
            Bsh(1:lmmaxvr,1:chunksize)=dimag(zfun1(1:lmmaxvr,iroffset:iroffset+chunksize-1,ias))
            Csh(1:lmmaxvr,1:chunksize)=Ash(1:lmmaxvr,1:chunksize)+Bsh(1:lmmaxvr,1:chunksize)

            Call dgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, Absht, ntpll, Ash, lmmaxvr, zzero, Asc, ntpll)
            Call dgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, Bbsht, ntpll, Bsh, lmmaxvr, zzero, Csc, ntpll)
            Call dgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, Cbsht, ntpll, Csh, lmmaxvr, zzero, Bsc, ntpll)

! the imaginary part of zfun1 in SH
            Bsc=Bsc-Asc-Csc
! the real part of zfun2 in SH
            Asc=Asc-Csc
! back up copy of the real part
            Csc=Asc

            do ir=1,chunksize
              do lm=1,ntpll
                ref2=dble(zfun2(lm,iroffset+ir-1,ias))
                imf2=dimag(zfun2(lm,iroffset+ir-1,ias))
                Asc(lm,ir)=ref2*Asc(lm,ir)-imf2*Bsc(lm,ir)
                Bsc(lm,ir)=ref2*Bsc(lm,ir)+imf2*Csc(lm,ir)
                Csc(lm,ir)=Asc(lm,ir)+Bsc(lm,ir)
              enddo
            enddo
!                Asc(1:ntpll,1:chunksize)=dble(zfun2(1:ntpll,iroffset:iroffset+chunksize-1,ias))*Asc(1:ntpll,1:chunksize)-dimag(zfun2(1:ntpll,iroffset:iroffset+chunksize-1,ias))*Bsc(1:ntpll,1:chunksize)
!                Bsc(1:ntpll,1:chunksize)=dble(zfun2(1:ntpll,iroffset:iroffset+chunksize-1,ias))*Bsc(1:ntpll,1:chunksize)+dimag(zfun2(1:ntpll,iroffset:iroffset+chunksize-1,ias))*Csc(1:ntpll,1:chunksize)
!                Csc=Asc+Bsc

            Call dgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, Afsht, lmmaxvr, Asc, ntpll, zzero, Ash, lmmaxvr)
            Call dgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, Bfsht, lmmaxvr, Bsc, ntpll, zzero, Csh, lmmaxvr)
            Call dgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, Cfsht, lmmaxvr, Csc, ntpll, zzero, Bsh, lmmaxvr)

            Bsh=Bsh-Ash-Csh
            Ash=Ash-Csh

            zres(1:lmmaxvr,iroffset:iroffset+chunksize-1,ias)=dcmplx(Ash(1:lmmaxvr,1:chunksize),Bsh(1:lmmaxvr,1:chunksize))

          endif

        enddo

! calculate the product zfun1*zfun2
! transform back to spherical harmonics
      enddo
    enddo

  end subroutine prodshrs

end module mod_eigensystem
