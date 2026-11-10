! Coefficient-space Fock exchange inside the muffin-tin spheres
! (WP3 of research-log/handoff-plan-claude-work.md). Since step 3.5 FockExchange and
! calcACE work only with the projections W(f,ias,n) built from these tables; calc_vxnl
! stops if an assumption of fock_mt_assumptions_ok is not met.
!
! Sphere functions f of species is are numbered in the order of wf%mtordered
! (orderwf.f90): the radial functions irad = 1..nradial(is) in the order of mod_radial's
! radial(:,irad,ias) / lrad(irad,is) (by l; for each l first the APW orders
! io = 1..apword(l,is), then the local orbitals with lorbl(ilo,is) = l in increasing ilo),
! each with m = -l..l. Function f belongs to irad = f_irad(f,is), l = f_l, m = f_m,
! lm = f_lm = idxlm(l,m); fstart(irad,is) is its f for m = -l.
!
! Maps from the loop orders of calcACE (calcace.f90) to f:
!   apw2f(if3,is): APW loop  l = 0..lmaxmat, m = -l..l, io = 1..apword(l,is)
!   lo2f(ilom,is): LO loop   ilo = 1..nlorb(is), m = -l..l  (counter per atom; calcACE's
!                  loindex runs on over atoms, so loindex = (LOs of earlier atoms) + ilom)
!
! Radial quadrature (Kristaps 11-10, gate 3.1): the spline weights of fderiv, so that
! the coefficient path reproduces the old fderiv integrals (calcACE, zfmtinp) and the
! normalisation of the radial functions (genapwfr, genlofr). fderiv(-1) is linear in
! its data (clamped spline, src_lib/spline.f90), so int_0^R f dr = sum_ir wspl(ir,is) f(ir)
! with wspl(ir,is) = fderiv(-1) of the unit vector e_ir, taken at R. init_fock_mt prints
! the check against fderiv for u_irad^2 r^2.
!
! Tables, built once per run next to mod_radial's tables (same guard as WP2), with wspl:
!   uwr(ir,irad,ias)      = wspl(ir,is) r^2 u_irad(r)
!   Rb(irad,L,jrad,ias)   = int u_irad (r/R_alpha)^L u_jrad r^2 dr,  L = 0..lmaxvr
!   Tpb2(irad,p,ias)      = int u_irad productbasis2_p r^2 dr,      p = 1..npbf2(0,ias)
!   gnt_mt                = <Y_lm | Y_LM | Y_l'm'> = int conj(Y_lm) Y_LM Y_l'm' dOmega
!                           (complex Y_lm; module gaunt, math/gaunt.f90, gaunt_yyy),
!                           sparse in LM: for lm1 = lm (l <= lmaxapw) and lm3 = l'm'
!                           (l' <= lmaxapw), i = 1..gnt_mt%num(lm1,lm3) gives
!                           LM = gnt_mt%lm2(i,lm1,lm3) (L <= lmaxvr, M = m - m') and the
!                           value gnt_mt%val(i,lm1,lm3). l_of_lm(LM) = L.
!   The lm index of module gaunt is the sequential l, m = -l..l index; init_fock_mt
!   asserts that it equals idxlm.
!
! project_mt(is, ias, zf, w): w(f) = int u_irad(f)(r) zf(lm(f), r) r^2 dr for one atom.

module mod_fock_mt
  use mod_radial, only: radial, nradial, lrad, maxradial, productbasis2, npbf2, maxpbf2
  use gaunt, only: non_zero_gaunt_real, non_zero_gaunt_yyy
  implicit none
  private

  !> .true. when the assumptions hold and the tables below are built
  logical, public :: fock_mt_on = .false.
  integer, public :: nfmax = 0, lmaxpot = -1
  integer, allocatable, public :: nfun(:)                                  ! (is)
  integer, allocatable, public :: f_irad(:,:), f_l(:,:), f_m(:,:), f_lm(:,:) ! (f,is)
  integer, allocatable, public :: fstart(:,:)                              ! (irad,is)
  integer, allocatable, public :: napwf(:), apw2f(:,:)                     ! (is), (if3,is)
  integer, allocatable, public :: nlof(:), lo2f(:,:)                       ! (is), (ilom,is)
  integer, allocatable, public :: l_of_lm(:)                               ! (lm)
  real(8), allocatable, public :: wspl(:,:)                                ! (ir,is)
  real(8), allocatable, public :: uwr(:,:,:)                               ! (ir,irad,ias)
  real(8), allocatable, public :: Rb(:,:,:,:)                              ! (irad,0:L,jrad,ias)
  real(8), allocatable, public :: Tpb2(:,:,:)                              ! (irad,p,ias)
  type(non_zero_gaunt_real), public :: gnt_mt

  public :: fock_mt_assumptions_ok, init_fock_mt, release_fock_mt, project_mt

contains

  !> Checks the assumptions under which the coefficient-space path is valid
  !> (handoff plan, section 0). With report, prints each one that is not met.
  logical function fock_mt_assumptions_ok(report) result(ok)
    use modinput, only: input
    use mod_atoms, only: nspecies, spr
    use mod_muffin_tin, only: nrmt, nrcmt, rcmt
    logical, intent(in) :: report
    integer :: is

    ok = .true.
    if (input%groundstate%lmaxmat /= input%groundstate%lmaxapw .or. &
        input%groundstate%lmaxvr /= input%groundstate%lmaxapw) call fail('lmaxmat = lmaxapw = lmaxvr')
    if (input%groundstate%lradstep /= 1) then
      call fail('lradstep = 1')
    else
      do is = 1, nspecies
        if (nrcmt(is) /= nrmt(is)) then
          call fail('nrcmt = nrmt')
        else if (any(rcmt(1:nrmt(is),is) /= spr(1:nrmt(is),is))) then
          call fail('rcmt = spr')
        end if
      end do
    end if
    if (trim(input%groundstate%hybrid%erfcapprox) /= 'none') call fail('erfcapprox = none')
    if (trim(input%groundstate%hybrid%singularity) /= 'exc0d') call fail('singularity = exc0d')
    if (.not. input%groundstate%hybrid%rpseudo) call fail('rpseudo = true')
    if (.not. input%groundstate%hybrid%rsurf) call fail('rsurf = true')
    if (associated(input%groundstate%spin)) call fail('no spin')
    if (input%groundstate%hybrid%updateRadial) call fail('updateRadial = false')

  contains

    subroutine fail(what)
      character(*), intent(in) :: what
      ok = .false.
      if (report) write(*,'("Info(mod_fock_mt): assumption not met: ",a)') what
    end subroutine fail

  end function fock_mt_assumptions_ok


  !> Builds the index maps and tables. Call after mod_radial's init_radial and
  !> init_productbasis_v3 (it reads radial, lrad, nradial, productbasis2, npbf2).
  !> If an assumption is not met, nothing is built and fock_mt_on stays .false.
  subroutine init_fock_mt
    use modinput, only: input
    use mod_APW_LO, only: apword, apwordmax, nlorb, lorbl
    use mod_atoms, only: nspecies, natoms, natmtot, idxas, spr
    use mod_muffin_tin, only: nrmt, nrmtmax, rmt, idxlm, lmmaxvr
    integer :: is, ia, ias, irad, jrad, l, m, io, ilo, f, if3, ilom, nr, p, nf, lmax, ll, lm, ir
    integer, allocatable :: irad_apw(:,:), irad_lo(:)
    real(8), allocatable :: r2(:), rl(:), f1(:), gr(:), cf(:,:)
    real(8) :: ta, tb, mb, dw(2)

    call timesec(ta)
    call release_fock_mt
    if (.not. fock_mt_assumptions_ok(.true.)) then
      write(*,'("Info(mod_fock_mt): coefficient-space path off (calc_vxnl stops)")')
      return
    end if
    lmax = input%groundstate%lmaxapw
    lmaxpot = input%groundstate%lmaxvr

    ! the lm order of module gaunt (sequential) must be idxlm
    lm = 0
    do l = 0, lmax
      do m = -l, l
        lm = lm + 1
        if (idxlm(l,m) /= lm) call stop_msg('idxlm is not the sequential (l,m) index')
      end do
    end do
    allocate(l_of_lm((lmaxpot+1)**2))
    do l = 0, lmaxpot
      l_of_lm(l*l+1:(l+1)**2) = l
    end do
    if (lmmaxvr /= (lmaxpot+1)**2) call stop_msg('lmmaxvr /= (lmaxvr+1)**2')

    ! --- sphere functions f, in the order of wf%mtordered
    allocate(nfun(nspecies), fstart(maxradial,nspecies))
    fstart = 0
    nfmax = 0
    do is = 1, nspecies
      nf = 0
      do irad = 1, nradial(is)
        fstart(irad,is) = nf + 1
        nf = nf + 2*lrad(irad,is) + 1
      end do
      nfun(is) = nf
      nfmax = max(nfmax, nf)
    end do
    allocate(f_irad(nfmax,nspecies), f_l(nfmax,nspecies), f_m(nfmax,nspecies), f_lm(nfmax,nspecies))
    f_irad = 0; f_l = -1; f_m = 0; f_lm = 0
    do is = 1, nspecies
      do irad = 1, nradial(is)
        l = lrad(irad,is)
        do m = -l, l
          f = fstart(irad,is) + l + m
          f_irad(f,is) = irad
          f_l(f,is) = l
          f_m(f,is) = m
          f_lm(f,is) = idxlm(l,m)
        end do
      end do
    end do

    ! --- maps from calcACE's loop orders to f; the irad order is re-derived with the
    !     loops of init_radial / orderWF and checked against lrad
    allocate(irad_apw(apwordmax,0:lmax), irad_lo(max(1,maxval(nlorb(1:nspecies)))))
    allocate(napwf(nspecies), nlof(nspecies), apw2f(nfmax,nspecies), lo2f(nfmax,nspecies))
    apw2f = 0; lo2f = 0
    do is = 1, nspecies
      irad_apw = 0; irad_lo = 0
      irad = 0
      do l = 0, lmax
        do io = 1, apword(l,is)
          irad = irad + 1
          irad_apw(io,l) = irad
          if (lrad(irad,is) /= l) call stop_msg('APW radial order differs from lrad')
        end do
        do ilo = 1, nlorb(is)
          if (lorbl(ilo,is) == l) then
            irad = irad + 1
            irad_lo(ilo) = irad
            if (lrad(irad,is) /= l) call stop_msg('LO radial order differs from lrad')
          end if
        end do
      end do
      if (irad /= nradial(is)) call stop_msg('number of radial functions differs from nradial')
      do ilo = 1, nlorb(is)
        if (irad_lo(ilo) == 0) call stop_msg('local orbital with l > lmaxapw')
      end do
      if3 = 0
      do l = 0, lmax
        do m = -l, l
          do io = 1, apword(l,is)
            if3 = if3 + 1
            apw2f(if3,is) = fstart(irad_apw(io,l),is) + l + m
          end do
        end do
      end do
      napwf(is) = if3
      ilom = 0
      do ilo = 1, nlorb(is)
        l = lorbl(ilo,is)
        do m = -l, l
          ilom = ilom + 1
          lo2f(ilom,is) = fstart(irad_lo(ilo),is) + l + m
        end do
      end do
      nlof(is) = ilom
      if (napwf(is) + nlof(is) /= nfun(is)) call stop_msg('APW + LO count differs from nfun')
    end do
    deallocate(irad_apw, irad_lo)

    ! --- spline weights of fderiv(-1): wspl(ir,is) = integral of the unit vector e_ir
    allocate(wspl(nrmtmax,nspecies), gr(nrmtmax), cf(3,nrmtmax))
    allocate(r2(nrmtmax), rl(nrmtmax), f1(nrmtmax))
    wspl = 0d0
    do is = 1, nspecies
      nr = nrmt(is)
      do ir = 1, nr
        f1(1:nr) = 0d0
        f1(ir) = 1d0
        call fderiv(-1, nr, spr(:,is), f1, gr, cf)
        wspl(ir,is) = gr(nr)
      end do
    end do
    ! check: sum wspl u^2 r^2 vs fderiv for every radial function of the first atom
    dw = 0d0
    do is = 1, nspecies
      nr = nrmt(is)
      ias = idxas(1,is)
      do irad = 1, nradial(is)
        f1(1:nr) = radial(1:nr,irad,ias)**2*spr(1:nr,is)**2
        call fderiv(-1, nr, spr(:,is), f1, gr, cf)
        dw(1) = max(dw(1), abs(dot_product(wspl(1:nr,is), f1(1:nr)) - gr(nr)))
        dw(2) = max(dw(2), abs(gr(nr)))
      end do
    end do
    write(*,'("Info(mod_fock_mt): spline weights vs fderiv, max rel. diff = ",es10.2)') dw(1)/dw(2)

    ! --- radial tables
    allocate(uwr(nrmtmax,maxradial,natmtot), Rb(maxradial,0:lmaxpot,maxradial,natmtot))
    allocate(Tpb2(maxradial,max(1,maxpbf2),natmtot))
    uwr = 0d0; Rb = 0d0; Tpb2 = 0d0
    do is = 1, nspecies
      nr = nrmt(is)
      r2(1:nr) = spr(1:nr,is)**2
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        do irad = 1, nradial(is)
          uwr(1:nr,irad,ias) = wspl(1:nr,is)*r2(1:nr)*radial(1:nr,irad,ias)
        end do
        ! (r/R)^L as in match_bound_mt2 (weinert.f90): repeated multiplication by r/R
        rl(1:nr) = 1d0
        do ll = 0, lmaxpot
          do jrad = 1, nradial(is)
            do irad = 1, nradial(is)
              f1(1:nr) = radial(1:nr,irad,ias)*rl(1:nr)*radial(1:nr,jrad,ias)*r2(1:nr)
              Rb(irad,ll,jrad,ias) = dot_product(wspl(1:nr,is), f1(1:nr))
            end do
          end do
          rl(1:nr) = rl(1:nr)*(spr(1:nr,is)/rmt(is))
        end do
        do p = 1, npbf2(0,ias)
          do irad = 1, nradial(is)
            f1(1:nr) = radial(1:nr,irad,ias)*productbasis2(1:nr,p,0,ias)*r2(1:nr)
            Tpb2(irad,p,ias) = dot_product(wspl(1:nr,is), f1(1:nr))
          end do
        end do
      end do
    end do
    deallocate(r2, rl, f1, gr, cf)

    ! --- Gaunt table <Y_lm | Y_LM | Y_l'm'>
    gnt_mt = non_zero_gaunt_yyy(lmax, lmaxpot, lmax)

    fock_mt_on = .true.
    call timesec(tb)
    mb = 8d0*(size(wspl) + size(uwr) + size(Rb) + size(Tpb2) + size(gnt_mt%val))/1d6 &
       + 4d0*(size(gnt_mt%lm2) + size(gnt_mt%num))/1d6
    write(*,'("Info(mod_fock_mt): nfun per species ",20i5)') nfun
    write(*,'("Info(mod_fock_mt): APW, LO functions per species ",20(2i5,2x))') &
      (napwf(is), nlof(is), is = 1, nspecies)
    write(*,'("Info(mod_fock_mt): Gaunt non-zeros ",i10,", tables ",f10.2," MB")') sum(gnt_mt%num), mb
    write(*,'("fock_mt_init : ",f12.4)') tb - ta

  contains

    subroutine stop_msg(what)
      character(*), intent(in) :: what
      write(*,*)
      write(*,'("Error(init_fock_mt): ",a)') what
      write(*,*)
      stop
    end subroutine stop_msg

  end subroutine init_fock_mt


  subroutine release_fock_mt
    if (allocated(nfun)) deallocate(nfun)
    if (allocated(f_irad)) deallocate(f_irad)
    if (allocated(f_l)) deallocate(f_l)
    if (allocated(f_m)) deallocate(f_m)
    if (allocated(f_lm)) deallocate(f_lm)
    if (allocated(fstart)) deallocate(fstart)
    if (allocated(napwf)) deallocate(napwf)
    if (allocated(apw2f)) deallocate(apw2f)
    if (allocated(nlof)) deallocate(nlof)
    if (allocated(lo2f)) deallocate(lo2f)
    if (allocated(l_of_lm)) deallocate(l_of_lm)
    if (allocated(wspl)) deallocate(wspl)
    if (allocated(uwr)) deallocate(uwr)
    if (allocated(Rb)) deallocate(Rb)
    if (allocated(Tpb2)) deallocate(Tpb2)
    if (allocated(gnt_mt%num)) deallocate(gnt_mt%num)
    if (allocated(gnt_mt%lm2)) deallocate(gnt_mt%lm2)
    if (allocated(gnt_mt%val)) deallocate(gnt_mt%val)
    nfmax = 0
    lmaxpot = -1
    fock_mt_on = .false.
  end subroutine release_fock_mt


  !> w(f) = int_0^R u_irad(f)(r) zf(lm(f), r) r^2 dr for the sphere functions of atom ias
  !> of species is, with the spline weights wspl. zf(lm, ir) on the muffin-tin grid
  !> (nrcmt = nrmt asserted). Only w(1:nfun(is)) is meaningful; the rest is set to zero.
  !> Reads module tables only, so it may be called inside parallel regions.
  subroutine project_mt(is, ias, zf, w)
    use mod_muffin_tin, only: nrmt
    integer, intent(in) :: is, ias
    complex(8), intent(in) :: zf(:,:)
    complex(8), intent(out) :: w(:)
    integer :: irad, l, lm0, f0, ir

    w = (0d0, 0d0)
    do irad = 1, nradial(is)
      l = lrad(irad,is)
      lm0 = l*l + 1
      f0 = fstart(irad,is)
      do ir = 1, nrmt(is)
        w(f0:f0+2*l) = w(f0:f0+2*l) + zf(lm0:lm0+2*l, ir)*uwr(ir,irad,ias)
      end do
    end do
  end subroutine project_mt

end module mod_fock_mt
