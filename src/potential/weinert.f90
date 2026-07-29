!> This module provides functionalities to compute the electrostatic potential
!> from a given charge density distribution, i.e., solving Poisson's equation,
!> using Weinert's method [1].
!>
!> [1]: [M. Weinert. Solution of Poisson's equation: Beyond Ewald-type methods.
!> *J. Math. Phys.* **22**, 2433(1981)](https://doi.org/10.1063/1.524800)
module weinert
  use precision, only: dp

  implicit none
  private

  public :: surface_ir, multipoles_ir, poisson_ir
  public :: poisson_and_multipoles_mt, match_bound_mt
  public :: poisson_and_multipoles_mt_yukawa, multipoles_ir_yukawa, pseudocharge_gspace_yukawa, poisson_ir_yukawa
  public :: poisson_mt_yukawa, pseudocharge_rspace_matrix, pseudocharge_rspace_new
  public :: multipoles_ir2, multipoles_ir3, multipoles_ir4, multipoles_ir5, pseudocharge_gspace2, pseudocharge_gspace3, poisson_ir2, surface_ir2, surface_ir3
contains

  !> This subroutine evaluates an interstitial function given by the Fourier series
  !> \[ f({\bf r}) = \sum_{\bf G} \hat{f}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
  !> on the surface of the muffin-tin spheres.
  !>
  !> The function evaluated on the surface of the muffin-tin sphere \(\alpha\) is given by
  !> \[ f^{{\rm SF},\alpha}(\theta,\phi) = f({\bf \tau}_\alpha + R_\alpha \hat{\bf r}(\theta,\phi)) =
  !> \sum_{l,m} f^{{\rm SF},\alpha}_{lm} \, Y_{lm}(\theta,\phi) \;,\]
  !> with
  !> \[ f^{{\rm SF},\alpha}_{lm} = 4\pi \, {\rm i}^l \sum_{\bf G} {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}
  !> \hat{f}({\bf G+p}) \, j_l(|{\bf G+p}| R_\alpha) \, Y_{lm}^\ast({\bf \widehat{G+p}}) \;,\]
  !> where \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) which is centered at \({\bf \tau}_\alpha\)
  !> and \(j_l(x)\) are the spherical Bessel functions.
  subroutine surface_ir( lmax, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, fsf)
    use mod_atoms, only: nspecies, natoms, idxas
    use mod_muffin_tin, only: rmt
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function
    complex(dp), intent(in) :: fig(:)
    !> function values on the muffin-tin sphere surfaces given by \(f^{{\rm SF},\alpha}_{lm}\)
    complex(dp), intent(out) :: fsf(:,:)

    integer :: is, ia, ias

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas( ia, is)
        call surface_ir_single_mt( lmax, rmt(is), ngp, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
          fig, fsf(:,ias))
      end do
    end do
  end subroutine surface_ir
  !> Same as [[surface_ir(subroutine)]] but for a single muffin-tin sphere.
  subroutine surface_ir_single_mt( lmax, rmt, ngp, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, fsf)
    use constants, only: zzero, fourpi
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> muffin-tin radius \(R_\alpha\)
    real(dp), intent(in) :: rmt
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:)
    !> bounds for integer components of \({\bf G}\)
    integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function
    complex(dp), intent(in) :: fig(:)
    !> function values on the muffin-tin sphere surfaces given by \(f^{{\rm SF},\alpha}_{lm}\)
    complex(dp), intent(out) :: fsf(:)

    integer :: l, m, lm, igp, ifg, ig(3)
    complex(dp) :: z1, z2, zil

    real(dp), allocatable :: rl3(:)

    fsf = zzero

    allocate( rl3(0:lmax))

    rl3(0) = rmt**3
    do l = 1, lmax
      rl3(l) = rl3(l-1)*rmt
    end do
!$omp parallel default(shared) private(igp,ig,ifg,z1,z2,zil,l,m,lm) reduction(+:fsf)
!$omp do
    do igp = 1, ngp
      ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
      ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
      z1 = fourpi*sfacgp(igp)*fig(ifg)
      zil = cmplx( 1._dp, 0._dp, dp)
      do l = 0, lmax
        z2 = z1*zil*jlgpr( l, igp)
        do m = -l, l
          lm = l*(l+1) + m + 1
          fsf(lm) = fsf(lm) + z2*conjg( ylmgp( lm, igp))
        end do
        zil = cmplx( -aimag(zil), dble(zil), dp)
      end do
    end do
!$omp end do
!$omp end parallel

    deallocate( rl3)
  end subroutine surface_ir_single_mt

  !> This subroutine calculates the multipole moments corresponding to the extension of
  !> an interstitial function given by the Fourier series
  !> \[ f({\bf r}) = \sum_{\bf G} \hat{f}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
  !> into the interior of the muffin-tin spheres.
  !>
  !> The multipole moments corresponding to the extension of the function inside the the muffin-tin sphere \(\alpha\)
  !> are given by
  !> \[ q^{{\rm IR},\alpha}_{lm} = 4\pi \, {\rm i}^l \, R_\alpha^{l+3} \sum_{\bf G} {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}
  !> \hat{f}({\bf G+p}) \, \frac{j_{l+1}(|{\bf G+p}| R_\alpha)}{|{\bf G+p}| R_\alpha} \, Y_{lm}^\ast({\bf \widehat{G+p}}) \;,\]
  !> where \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) which is centered at \({\bf \tau}_\alpha\)
  !> and \(j_l(x)\) are the spherical Bessel functions.
  subroutine multipoles_ir( lmax, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, qlm)
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas
    use mod_muffin_tin, only: rmt
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: fig(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qlm(:,:)

    integer :: is, ia, ias

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas( ia, is)
        call multipoles_ir_single_mt( lmax, rmt(is), ngp, gpc, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
          fig, qlm(:,ias), epslat=input%structure%epslat)
      end do
    end do
  end subroutine multipoles_ir
  !> Same as [[multipoles_ir(subroutine)]] but for a single muffin-tin sphere.
  subroutine multipoles_ir_single_mt( lmax, rmt, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, qlm, epslat)
    use constants, only: zzero, fourpi, y00
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> muffin-tin radius \(R_\alpha\)
    real(dp), intent(in) :: rmt
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:)
    !> bounds for integer components of \({\bf G}\)
    integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: fig(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qlm(:)
    !> threshold below which \(|{\bf G+p}|\) is considered zero
    real(dp), optional, intent(in) :: epslat

    integer :: i, l, m, lm, igp, ngpf, ngpz, ifg, ig(3)
    real(dp) :: eps
    complex(dp) :: z1, z2, zil, figzero

    integer, allocatable :: igp_finite(:), igp_zero(:)
    real(dp), allocatable :: rl3(:)

    eps = 1.d-12
    if( present( epslat)) eps = epslat

    qlm = zzero
    figzero = zzero

    allocate( rl3(0:lmax))

    rl3(0) = rmt**3
    do l = 1, lmax
      rl3(l) = rl3( l-1)*rmt
    end do

    igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > eps, i=1, ngp)])
    ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg,z1,z2,zil,l,m,lm) reduction(+:qlm)
!$omp do
    do i = 1, ngpf
      igp = igp_finite(i)
      ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
      ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
      z1 = fourpi*sfacgp(igp)*fig(ifg)/(gpc(igp)*rmt)
      zil = cmplx( 1._dp, 0._dp, dp)
      do l = 0, lmax
        z2 = z1*zil*rl3(l)*jlgpr( l+1, igp)
        do m = -l, l
          lm = l*(l+1) + m + 1
          qlm(lm) = qlm(lm) + z2*conjg( ylmgp( lm, igp))
        end do
        zil = cmplx( -aimag(zil), dble(zil), dp)
      end do
    end do
!$omp end do
!$omp end parallel
    igp_zero = pack( [(i, i=1, ngp)], [(gpc(i) <= eps, i=1, ngp)])
    ngpz = size( igp_zero)
    do i = 1, ngpz
      igp = igp_zero(i)
      ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
      ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
      qlm(1) = qlm(1) + fourpi/3._dp*rl3(0)*y00*fig(ifg)
    end do

    deallocate( rl3)
  end subroutine multipoles_ir_single_mt

  !> This subroutine solves Poisson's equation for a given complex charge density
  !> contained in an isolated muffin-tin sphere using the Green's function approach
  !> and calculates the multipole moments of the charge distribution.
  !>
  !> The general electrostatic potential arising from an arbitrary charge distribution
  !> \[ n({\bf r}) = n({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) =
  !> \sum_{l,m} n^\alpha_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
  !> inside the muffin-tin sphere \(\alpha\) is given by
  !> \[ v_{\rm sph}({\bf r}) = v_{\rm sph}({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) =
  !> \sum_{l,m} v_{\rm sph}[n^\alpha_{lm}](r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
  !> with
  !> \[ v_{\rm sph}[n^\alpha_{lm}](r) = \frac{4\pi}{2l+1} \left[
  !> \frac{1}{r^{l+1}} \int_0^r s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s +
  !> r^l \int_r^{R_\alpha} \frac{n^\alpha_{lm}(s)}{s^{l-1}} \, {\rm d}s -
  !> \frac{r^l}{R_\alpha^{2l+1}} \int_0^{R_\alpha} s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s \right] \;.\]
  !>
  !> In addition, the multipole moments of the charge distribution in the sphere are calculated as
  !> \[ q^{{\rm MT},\alpha}_{lm} = \int_0^{R_\alpha} s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s \,. \]
  !>
  !> @note \(v_{\rm sph}({\rm r})\) as defined above goes to zero on the muffin-tin surface.@endnote
  subroutine poisson_and_multipoles_mt( lmax, nr, r, zrhomt, zvclmt, qlm)
    use constants, only: fourpi
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> number of radial grid points
    integer, intent(in) :: nr
    !> radial grid
    real(dp), intent(in) :: r(:)
    !> complex charge distribution \(n^\alpha_{lm}(r)\)
    complex(dp), intent(in) :: zrhomt(:,:)
    !> complex electrostatic potential \(v_{\rm sph}[n^\alpha_{lm}](r)\)
    complex(dp), intent(out) :: zvclmt(:,:)
    !> multipole moments of the charge distribution \(q^{{\rm MT},\alpha}_{lm}\)
    complex(dp), intent(out) :: qlm(:)

    integer :: l, m, lm, ir
    real(dp) :: t1, t2, t3, t4, t5

    real(dp), allocatable :: ri(:), rl(:), ril1(:), cf(:,:)
    real(dp), allocatable :: fr(:,:), gr(:,:)

    allocate( ri(nr), rl(nr), ril1(nr), cf(3,nr))
    allocate( fr(nr,4), gr(nr,4))

    ! initialise r^l and r^(-l-1)
    do ir = 1, nr
      rl(ir) = 1._dp
      ri(ir) = 1._dp/r(ir)
      ril1(ir) = ri(ir)
    end do
    lm = 0
    do l = 0, lmax
      t1 = fourpi/dble( 2*l+1)
      do m = -l, l
        lm = lm + 1
        do ir = 1, nr
          t2 = rl(ir)*r(ir)*r(ir)    ! r^(2+l)
          t3 = ril1(ir)*r(ir)*r(ir)  ! r^(1-l)
          t4 = dble( zrhomt(lm, ir))
          t5 = aimag( zrhomt(lm, ir))
          fr(ir,1) = t2*t4
          fr(ir,2) = t2*t5
          fr(ir,3) = t3*t4
          fr(ir,4) = t3*t5
        end do
        call fderiv( -1, nr, r, fr(:,1), gr(:,1), cf)
        call fderiv( -1, nr, r, fr(:,2), gr(:,2), cf)
        call fderiv( -1, nr, r, fr(:,3), gr(:,3), cf)
        call fderiv( -1, nr, r, fr(:,4), gr(:,4), cf)
        qlm(lm) = cmplx( gr(nr,1), gr(nr,2), 8)
        do ir = 1, nr
          t2 = ril1(ir)*gr(ir,1) + rl(ir)*(gr(nr,3) - gr(ir,3) - gr(nr,1)*ril1(nr)/rl(nr))
          t3 = ril1(ir)*gr(ir,2) + rl(ir)*(gr(nr,4) - gr(ir,4) - gr(nr,2)*ril1(nr)/rl(nr))
          zvclmt( lm, ir) = t1*cmplx( t2, t3, 8)
        end do
      end do
      ! update r^l and r^(-l-1)
      if( l < lmax) then
        do ir = 1, nr
          rl(ir) = rl(ir)*r(ir)
          ril1(ir) = ril1(ir)*ri(ir)
        end do
      end if
    end do

    deallocate( ri, rl, ril1, cf, fr, gr)
    return
  end subroutine poisson_and_multipoles_mt

  !> This subroutine solves Poisson's equation for a complex charge density given in the
  !> interstitial region by
  !> \[ n({\bf r}) = \sum_{\bf G} \hat{n}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
  !> and having multipole moments \(q^\alpha_{lm}\) inside the muffin-tin spheres.
  !>
  !> This is done by constructing a pseudodensity with the given multipole moments and having a
  !> rapidly converging Fourier series. The Fourier components of the pseudodensity are given by
  !> \[ \hat{n}^{\rm ps}({\bf G+p}) = \hat{n}({\bf G+p}) + \sum_\alpha \hat{n}^{{\rm ps},\alpha}({\bf G+p}) \;.\]
  !> See [[pseudodensity_ir_single_mt(subroutine)]] for further details.
  !>
  !> From the pseudodensity the electrostatic potential in the interstitial region is obtained
  !> by solving Poisson's equation in reciprocal space, i.e.,
  !> \[ V({\bf r}) = \sum_{\bf G} \hat{V}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \;, \]
  !> with
  !> \[ \hat{V}({\bf G+p}) = 4\pi \frac{\hat{n}^{\rm ps}({\bf G+p})}{|{\bf G+p}|^2} \;.\]
  subroutine poisson_ir( lmax, npsden, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, zrhoig, qlm, zvclig,&
    & cutoff_in,hybrid_in,rpseudo_in,rpseudomat)
    use modinput
    use mod_lattice, only: omega
    use mod_atoms, only: nspecies, natoms, idxas
    Use mod_kpoint, only: nkptnr
    use mod_muffin_tin, only: rmt
    use constants, only: zzero, fourpi, twopi
    use mod_Gvector, only: ngrid
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> pseudodensity expansion order
    integer, intent(in) :: npsden
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> on entry: Fourier components \(\hat{n}({\bf G+p})\) of the interstitial charge density on the FFT grid
    !> on exit: Fourier components \(\hat{n}^{\rm ps}({\bf G+p})\) of the pseudodensity on the FFT grid
    complex(dp), intent(inout) :: zrhoig(:)
    !> muffin-tin multipole moments \(q^\alpha_{lm}\) of the charge density
    complex(dp), intent(in) :: qlm(:,:)
    !> option for using coulomb cutoff for solving Poisson's equation
    logical, optional, intent(in) :: cutoff_in
    logical, optional, intent(in) :: rpseudo_in
    logical, optional, intent(in) :: hybrid_in
    !> Fourier components \(\hat{V}({\bf G+p})\) of the interstitial electrostatic potential on the FFT grid
    complex(dp), intent(out) :: zvclig(:)
    complex(dp),optional, intent(in) :: rpseudomat(:,:,:)
    integer :: i, is, ia, ias, igp, ngpf, ifg, ig(3),ii
    real(8) :: r_c
    integer, allocatable :: igp_finite(:)

    logical :: cutoff, hybrid, rpseudo

    if (present(hybrid_in)) then
      hybrid=hybrid_in
    else
      hybrid=.false.
    endif

    if (present(cutoff_in)) then
      cutoff=cutoff_in
    else
      cutoff=.false.
    endif

    if (present(rpseudo_in)) then
      rpseudo=rpseudo_in
    else
      rpseudo=.false.
    endif

    if (rpseudo) then

      call pseudocharge_rspace_new(lmax,qlm,rpseudomat,zrhoig)
    else

      ! add Fourier components of pseudodensity from multipole moments
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call pseudodensity_ir_single_mt( lmax, rmt(is), omega, npsden, ngp, gpc, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
            zrhoig, qlm(:,ias), epslat=input%structure%epslat)
        end do
      end do

      if(hybrid) then
        write(*,*)"pseudo_ref.dat"
        open(11,file='pseudo_ref.dat',status='replace')
        do ii=1, ngp
          i=igfft(ii)
          write(11,*)dble(zrhoig(i)),imag(zrhoig(i))
        enddo
        close(11)
        stop
      endif



    endif

    ! solve Poisson's equation in reciprocal space
!      zvclig = zzero
    igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > input%structure%epslat, i=1, ngp)])
    ngpf = size( igp_finite)



    if (cutoff) then
      r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0
      zvclig = zrhoig*twopi*r_c**2
!$omp parallel default(shared) private(i,igp,ig,ifg)
!$omp do
      ! cutoff correction for > epslat
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        zvclig(ifg) = fourpi*zrhoig(ifg)*(1d0-cos(gpc(igp)*r_c))/(gpc(igp)**2)
      end do
!$omp end do
!$omp end parallel
    else
      zvclig = zzero
      ! without cutoff correction
!$omp parallel default(shared) private(i,igp,ig,ifg)
!$omp do
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        zvclig(ifg) = fourpi*zrhoig(ifg)/(gpc(igp)**2)
      end do
!$omp end do
!$omp end parallel
    end if
  end subroutine poisson_ir

  !> This subroutine computes the Fourier components of a quickly converging pseudodensity with multipole moments \(q_{lm}\)
  !> in muffin-tin \(\alpha\) and adds them to the input argument `zrhoig`.
  !>
  !> The Fourier components to be added are given by
  !> \[ \hat{n}^{{\rm ps},\alpha}({\bf G+p}) = \frac{4\pi}{\Omega} {\rm e}^{-{\rm i}({\bf G+p}) \cdot {\bf \tau}_\alpha}
  !> \sum_{l,m} \left( -\frac{\rm i}{R_\alpha} \right)^l \frac{(2l + 2N + 3)!!}{(2l + 1)!!}
  !> \frac{j_{l+N+1}(|{\bf G+p}|R_\alpha)}{(|{\bf G+p}|R_\alpha)^{N+1}} \, q^\alpha_{lm} \, Y_{lm}(\widehat{\bf G+p}) \;,\]
  !> where \(\Omega\) is the unit cell volume, \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\)
  !> which is centered at \({\bf \tau}_\alpha\) and \(j_l(x)\) are the spherical Bessel functions.
  !> \(N\) is the expansion order of the pseudodensity.
  subroutine pseudodensity_ir_single_mt( lmax, rmt, omega, npsden, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, zrhoig, qlm, epslat)
    use constants, only: zzero, fourpi, zi, y00
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> muffin-tin radius \(R_\alpha\)
    real(dp), intent(in) :: rmt
    !> unit cell volume \(\Omega\)
    real(dp), intent(in) :: omega
    !> pseudodensity expansion order
    integer, intent(in) :: npsden
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:)
    !> bounds for integer components of \({\bf G}\)
    integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{n}^{\rm ps}({\bf G+p})\) of the pseudodensity on the FFT grid
    complex(dp), intent(inout) :: zrhoig(:)
    !> muffin-tin multipole moments \(q^\alpha_{lm}\) of the charge density
    complex(dp), intent(in) :: qlm(:)
    !> threshold below which \(|{\bf G+p}|\) is considered zero
    real(dp), optional, intent(in) :: epslat

    integer :: i, l, m, lm, igp, ngpf, ngpz, ifg, ig(3)
    real(dp) :: eps, fpo, t1, t2
    complex(dp) :: z1, z2

    integer, allocatable :: igp_finite(:), igp_zero(:)
    complex(dp), allocatable :: zrp(:)

    real(dp), external :: factnm

    eps = 1.d-12
    if( present( epslat)) eps = epslat

    fpo = fourpi/omega
    allocate( zrp( (lmax+1)**2))

    ! add Fourier components of pseudodensity from multipole moments
    t1 = 1.d0
    do l = 0, lmax
      t2 = t1*factnm( 2*l + 2*npsden + 3, 2)/factnm( 2*l+1, 2)
      t1 = t1/rmt
      do m = -l, l
        lm = l*(l+1) + m + 1
        zrp(lm) = conjg( zi**l)*t2*qlm(lm)
      end do
    end do

    igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > eps, i=1, ngp)])
    ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg,l,m,lm,z1,z2,t1,t2)
!$omp do
    do i = 1, ngpf
      igp = igp_finite(i)
      ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
      ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
      z1 = fpo*conjg( sfacgp(igp))
      t1 = (gpc(igp)*rmt)**(npsden + 1)
      do l = 0, lmax
        t2 = jlgpr( l+npsden+1, igp)/t1
        z2 = zzero
        do m = -l, l
          lm = l*(l+1) + m + 1
          z2 = z2 + zrp(lm)*ylmgp( lm, igp)
        end do
        zrhoig(ifg) = zrhoig(ifg) + z1*t2*z2
      end do
    end do
!$omp end do
!$omp end parallel
    igp_zero = pack( [(i, i=1, ngp)], [(gpc(i) <= eps, i=1, ngp)])
    ngpz = size( igp_zero)
    do i = 1, ngpz
      igp = igp_zero(i)
      ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
      ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
      zrhoig(ifg) = zrhoig(ifg) + fpo/factnm( 2*npsden + 3, 2)*zrp(1)*y00
    end do
    deallocate( zrp)
  end subroutine pseudodensity_ir_single_mt

  !> Given a complex muffin-tin function
  !> \[ f({\bf r}) = f({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) =
  !> \sum_{l,m} f^\alpha_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
  !> and an interstitial function given on the muffin-tin sphere surface by
  !> \(f^{{\rm SF},\alpha}_{lm}\) according to subroutine [[surface_ir(subroutine)]],
  !> this subroutine matches the muffin-tin function to the interstial function by
  !> adding the homogeneous function
  !> \[ f^{\rm hom}({\bf r}) = f^{\rm hom}({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) =
  !> \sum_{l,m} \left[ f^{{\rm SF},\alpha}_{lm} - f^\alpha_{lm}(R_\alpha) \right] \,
  !> \left( \frac{r_\alpha}{R_\alpha} \right)^l \, Y_{lm}(\hat{\bf r}_\alpha) \]
  !> to the muffin-tin function.
  subroutine match_bound_mt( lmax, nr, r, rmt, firsf, fmt, yukawa, zbessi)
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> number of radial grid points
    integer, intent(in) :: nr
    !> radial grid
    real(dp), intent(in) :: r(:)
    !> muffin-tin radius \(R_\alpha\)
    real(dp), intent(in) :: rmt
    !> interstitial function on muffin-tin sphere surface \(f^{{\rm SF},\alpha}_{lm}\)
    complex(dp), intent(in) :: firsf(:)
    !> muffin-tin function \(f^\alpha_{lm}(r)\)
    complex(dp), intent(inout) :: fmt(:,:)

    logical, optional, intent(in) :: yukawa
    complex(dp),optional, intent(in) :: zbessi(:,0:) !nrmtmax, 0:input%groundstate%lmaxvr+input%groundstate%npsden+1
    integer :: l, m, lm, ir
    complex(dp) :: df,zt1

    real(dp), allocatable :: rr(:,:)


    if ((present(yukawa)).and.yukawa) then

      lm = 0
      do l = 0, lmax
        do m = -l, l
          lm = lm + 1
          zt1 = firsf(lm) - fmt (lm, nr)
          Do ir = 1, nr
            fmt (lm, ir) = fmt (lm, ir) + zt1 &
              & * zbessi(ir,l)/zbessi(nr,l)


          End Do !ir
        enddo
      enddo



    else


      allocate( rr(nr,2), source=1._dp)
      do ir = 1, nr
        rr(ir,2) = r(ir)/rmt
      end do

      lm = 0
      do l = 0, lmax
        do m = -l, l
          lm = lm + 1
          df = firsf(lm) - fmt(lm,nr)
          do ir = 1, nr
            fmt(lm,ir) = fmt(lm,ir) + df*rr(ir,1)
          end do
        end do
        do ir = 1, nr
          rr(ir,1) = rr(ir,1)*rr(ir,2)
        end do
      end do

      deallocate( rr)
    endif
  end subroutine match_bound_mt









  subroutine poisson_and_multipoles_mt_yukawa( lmax, nr, r, zrhomt, zvclmt, qlm, is, yukawa_in ,zlambda, il, kl)
    use modinteg
    use constants, only: fourpi
!> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
!> number of radial grid points
    integer, intent(in) :: nr
!> radial grid
    real(dp), intent(in) :: r(:)
!> complex charge distribution \(n^\alpha_{lm}(r)\)
    complex(dp), intent(in) :: zrhomt(:,:)
!> complex electrostatic potential \(v_{\rm sph}[n^\alpha_{lm}](r)\)
    complex(dp), intent(out) :: zvclmt(:,:)
!> multipole moments of the charge distribution \(q^{{\rm MT},\alpha}_{lm}\)
    complex(dp), intent(out) :: qlm(:)
    integer, intent(in) :: is
    logical, optional, intent(In) :: yukawa_in
    Complex (8),optional,Intent (In) :: il(:,0:), kl(:,0:),zlambda


    integer :: l, m, lm


    complex(dp) :: zt1,zt2
    logical :: yukawa
!external functions
    Real (8) :: factnm
    External factnm

    if(present(yukawa_in))then
      yukawa=yukawa_in
      call poisson_mt_yukawa( lmax, nr, r, zrhomt(:, :), zvclmt, is, &
        & yukawa_in=yukawa,zlambda=zlambda, il=il(:,:), kl=kl(:,:))
    else
      yukawa=.false.
      call poisson_mt_yukawa( lmax, nr, r, zrhomt(:, :), zvclmt, is)

    endif



    if(yukawa)then
      lm=0
      Do l = 0, lmax
        zt1 = factnm (2*l+1, 2) / (zlambda ** l)
        zt2 = 1d0 / ( kl(nr,l)* fourpi * zlambda)
        Do m = - l, l
          lm = lm + 1
          qlm (lm) = zt1 * zt2 * zvclmt (lm,nr)
        End Do
      enddo
    else
      lm=0
      Do l = 0, lmax
        Do m = - l, l
          lm = lm + 1
          qlm (lm) = (2*l+1) * r(nr)**(l+1) *zvclmt (lm,nr) / fourpi
        End Do
      enddo
    endif

!write(*,*)"mans"
!do ir=1, nr
    !write(*,*)dble(il(ir, 0)),",",imag(il(ir, 0))
    !write(*,*)dble(zrhomt(1, ir)),",",imag(zrhomt(1, ir))
    !write(*,*)dble(zvclmt (1, ir)),",", imag(zvclmt (1, ir))
!enddo
!stop
!write(*,*)qlm(1:3)


  end subroutine poisson_and_multipoles_mt_yukawa


  subroutine multipoles_ir_yukawa( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qi, zlambda,zilmt)
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil,zzero
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qi(:,:)
    complex(dp), intent(in) :: zlambda
    complex(dp), intent(in) :: zilmt(0:,:)

    integer :: is, ia, ias, ig, ifg, l, lm, m
    complex(dp) :: zt1,zt2,tbessi,zt3
    real(dp) :: t1,tbessj


!external functions
    Real (8) :: factnm
    External factnm

    qi (:, :) = zzero


    do is = 1, nspecies
      call msbesselic1(rmt(is)*zlambda,tbessi)
      do ia = 1, natoms(is)
        ias = idxas( ia, is)

        Do ig = 1, ngvec
          ifg = igfft (ig)

          If (gpc(ig) .Gt. input%structure%epslat) Then
            zt1 = zvclir (ifg) * sfacgp (ig, ias)/((gpc(ig) ** 2)+(zlambda ** 2))

            lm = 0
            Do l = 0, input%groundstate%lmaxvr
              if (l .eq. 0) then
                call sbessel1 (gpc(ig)*rmt(is), tbessj)
                zt2 = zt1 *  fourpi * zil(l) * (rmt(is) ** 2) * factnm (2*l+1, 2) /(zlambda ** (l))&
                  &* ((zlambda * jlgpr (l, ig, is) *  tbessi )&
                  &- (gpc(ig) *tbessj* zilmt(l,is)))
                !zt2 = zt1 *  fourpi * (rmt(is) ** 2) &
                !&* ((zlambda * jlgpr (0, ig, is) * (zilmt(1,is)+ (zilmt(0,is)/(zlambda*rmt(is)))))&
                !&- (gpc(ig) * ((jlgpr(0, ig, is)/(gpc(ig)*rmt(is)))-jlgpr(1, ig, is)) * zilmt(0,is)))
              else
                zt2 = zt1 *  fourpi * zil(l) * (rmt(is) ** 2) * factnm (2*l+1, 2) /(zlambda ** (l))&
                  &* ((zlambda * jlgpr (l, ig, is) * zilmt(l-1,is))&
                  &- (gpc(ig) * jlgpr (l-1, ig, is) * zilmt(l,is)))
              endif

              Do m = - l, l
                lm = lm + 1
                qi (lm, ias) = qi (lm, ias) + zt2 * conjg (ylmgp(lm, ig))
              End Do
            End Do
          Else
            zt3 = fourpi * y00 * (rmt (is) ** 2) * zilmt(1,is) / zlambda
            qi (1,ias) = qi (1,ias) + zt3 * zvclir (ifg)
          End If
        enddo
      end do
    end do

  end subroutine multipoles_ir_yukawa


  subroutine pseudocharge_gspace_yukawa( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qlm, zlambda, zilmt)
    use mod_lattice, only: omega
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(inout) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(in) :: qlm(:,:)
    complex(dp), intent(in) :: zlambda

    complex(dp), intent(in) :: zilmt(0:,:)

    complex(dp) :: zrp((lmax+1)**2)


    integer :: is, ia, ias, ig, ifg, l, lm, m,npsd
    complex(dp) :: zt1,zsum1, zsum2
    real(dp) :: t1,z1,fpo
!external functions
    Real (8) :: factnm
    External factnm


    fpo = fourpi/omega
    npsd=input%groundstate%npsden

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas( ia, is)
        lm=0
        Do l = 0, lmax
          t1 = 1d0 / (factnm (2*l+1, 2))
          Do m = - l, l
            lm = lm + 1
            zrp (lm) = qlm(lm, ias) * t1
            !write(*,*)"qlm",qlm(lm, ias),t1
          End Do
        End Do



        Do ig = 1, ngvec
          ifg = igfft (ig)
          If (gpc(ig) .Gt. input%structure%epslat) Then
            z1 = gpc (ig)
            zt1 = fpo * conjg (sfacgp(ig, ias))/(z1 ** (npsd+1))
            lm = 0
            Do l = 0, lmax
              zsum1 = (jlgpr (npsd+l+1, ig, is)*(zlambda ** (l+npsd+1))/ zilmt(npsd+l+1,is))* conjg(zil(l))
              lm = lm + 1
              zsum2 = zrp (lm) * ylmgp (lm, ig)
              Do m = - l + 1, l
                lm = lm + 1
                zsum2 = zsum2 + zrp (lm) * ylmgp (lm, ig)
              End Do
              zvclir (ifg) = zvclir (ifg) +  zt1 * zsum1 * zsum2
            End Do

          Else
            !Sie divi labojumi no oriģinālā szpotcoul:
            zt1 = (fpo * y00 * (rmt(is) ** (npsd+1))&
              &* (zlambda ** (npsd+1))) / ((factnm (2*npsd+3, 2))* zilmt(npsd+1,is))

            zvclir (ifg) = zvclir (ifg) + zt1 * qlm(1, ias)

          End If
        End Do
      End Do
    End Do

!do ig=1,20
!  ifg = igfft (ig)
!  write(*,*)zvclir (ifg)
!enddo
  end subroutine pseudocharge_gspace_yukawa

  subroutine poisson_ir_yukawa( lmax, ngvec, gpc, igfft, zrhoir, zlambda,zvclir,cutoff)
    use mod_lattice, only: omega
    Use mod_kpoint, only: nkptnr
    use modinput
    use constants, only: zzero
    ! use mod_atoms, only: nspecies, natoms, idxas
    ! use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zrhoir(:)

    complex(dp), intent(in) :: zlambda
    complex(dp), intent(out) :: zvclir(:)
    logical, intent(in) :: cutoff


    integer :: ig, ifg
    real(dp) :: r_c
!external functions
    Real (8) :: factnm
    External factnm

    zvclir=zzero
    r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0
    Do ig = 1, ngvec
      ifg = igfft (ig)

      if (cutoff) then

        If (gpc(ig) .Gt. input%structure%epslat) Then
          zvclir (ifg) = fourpi * zrhoir (ifg) *&
            (1d0 - exp(-zlambda*r_c)* ( zlambda * sin(gpc(ig)*r_c) / gpc(ig) + cos(gpc(ig)*r_c) ) )/&
            ((gpc(ig) ** 2) + (zlambda ** 2))
        Else
          zvclir (ifg) = fourpi * zrhoir(ifg)* (1d0 - exp(-zlambda*r_c)*(zlambda*r_c + 1))/(zlambda ** 2)
        End If
        !write(*,*)ig,zvclir (ifg)
      else

        zvclir (ifg) = fourpi * zrhoir (ifg) / ((gpc(ig) ** 2) + (zlambda ** 2))

        !write(*,*)ig,zvclir (ifg)
      endif

    End Do
!stop
  end subroutine poisson_ir_yukawa




  subroutine poisson_mt_yukawa( lmax, nr, r, zrhomt, zvclmt,is, yukawa_in,zlambda, il, kl)
    use modinteg
    use constants, only: fourpi,zzero
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> number of radial grid points
    integer, intent(in) :: nr
    !> radial grid
    real(dp), intent(in) :: r(:)
    !> complex charge distribution \(n^\alpha_{lm}(r)\)
    complex(dp), intent(in) :: zrhomt(:,:)
    !> complex electrostatic potential \(v_{\rm sph}[n^\alpha_{lm}](r)\)
    complex(dp), intent(out) :: zvclmt(:,:)
    integer, intent(in) :: is

    Complex (8),optional,Intent (In) :: il(:,0:), kl(:,0:),zlambda
    logical, optional, intent(In) :: yukawa_in

    integer :: l, m, lm, ir

    complex(dp) :: f1(nr),f2(nr),g1(nr),g2(nr),zt1,zrho_tmp(nr)
    real(dp)  :: t1
    complex(dp) ,allocatable :: rl(:),ril1(:),ri(:),r2(:),tr1(:),tr2(:)
    logical :: yukawa

    if(present(yukawa_in))then
      yukawa=yukawa_in
    else
      yukawa=.false.
    endif

    zvclmt=zzero

    if(yukawa) then
      zt1=fourpi * zlambda
      allocate(r2(nr))
      r2=r(:nr)*r(:nr)
      lm = 0
      Do l = 0, lmax
        Do m= -l, l
          lm = lm + 1
          zrho_tmp = zrhomt(lm, :nr)
          f1 = il(:nr,l) * r2 * zrho_tmp
          call integ_cf (nr, is, f1, g1, mt_integw)
          f1 = kl(:nr,l) * g1
          f2 = kl(:nr,l) * r2 * zrho_tmp
          call integ_cf (nr, is, f2, g2, mt_integw)
          f2= il(:nr,l) * (g2(nr)-g2)
          zvclmt (lm, :nr)=zt1 * (f1+f2)
        Enddo
      enddo
      deallocate(r2)
    else !not Yukawa
      allocate(rl(nr),ril1(nr),ri(nr),r2(nr),tr1(nr),tr2(nr))

      rl(:) = 1d0          !r^l
      ri(:) = 1d0/r(:nr)
      ril1(:) = 1d0/r(:nr) ! r^(-l-1)
      r2=r(:nr)*r(:nr)

      lm = 0
      Do l = 0, lmax
        t1 = fourpi/(2*l+1)
        tr1=rl*r2
        tr2=ril1*r2
        Do m= -l, l
          lm = lm + 1
          zrho_tmp = zrhomt(lm, :nr)
          f1=zrho_tmp*tr1
          call integ_cf (nr, is, f1, g1, mt_integw)
          f1=g1*ril1 !/ r(:nr)**(l+1)

          f2=zrho_tmp*tr2 !/r(:nr)**(l-1)
          call integ_cf (nr, is, f2, g2, mt_integw)
          f2= rl * (g2(nr)-g2)

          zvclmt (lm, :nr)=(f1+f2)*t1
        Enddo
        ! update r^l and r^(-l-1)
        if( l < lmax) then
          rl = rl*r(:nr)
          ril1 = ril1*ri
        end if
      enddo
      deallocate(rl,ril1,ri,r2,tr1,tr2)







      ! lm = 0
      ! Do l = 0, lmax
      !   Do m= -l, l
      !     lm = lm + 1
      !     f1=zrhomt(lm, :nr)*r(:nr)**(l+2)
      !     call integ_cf (nr, is, f1, g1, mt_integw)
      !     f1=g1 / r(:nr)**(l+1)

      !     f2=zrhomt(lm, :nr)/r(:nr)**(l-1)
      !     call integ_cf (nr, is, f2, g2, mt_integw)
      !     f2= r(:nr)**(l)* (g2(nr)-g2)

      !     !f3=-r(:nr)**(l)*g1(nr)/r(nr)**(2*l+1)

      !     zvclmt (lm, :nr)=fourpi *(f1+f2)/(2*l+1)
      !   Enddo
      ! enddo

    endif! if yukawa or not



  end subroutine poisson_mt_yukawa


  subroutine pseudocharge_rspace_matrix(lmax,npsd,kvec,mat,yukawa_in,zlambda,zilmt,zbessi)

    use modrspace
    use modinput
    use modinteg
    use mod_atoms, only: nspecies, natoms, idxas, natmtot,atposc,spr
    use mod_Gvector, only: ngrid,ngrtot
    use mod_muffin_tin, only: rmt, nrmtmax,nrmt
    use constants, only: zzero, fourpi, y00
    implicit none

    integer, intent(in) :: lmax
    integer, intent(in) :: npsd
    Real(8), Intent (In) :: kvec(3)
    Complex (8), Intent (Out) ::  mat(:,:,:) !(lm,irg,ias)
    logical, optional, intent(in) :: yukawa_in
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)
    complex(dp),optional, intent(in) :: zbessi(:,0:,:)


    complex(dp) :: phase
    logical :: yukawa

    complex (8) :: zf1(nrmtmax),zf2(nrmtmax)
    complex (8) :: zt1, zt2
    complex (8) :: rinteg(0:lmax,nspecies) !(l,is)
    complex (8) :: zldep(0:lmax)
    Real (8) :: rr,t1,t2,rv_abs,ratom(3), rv(3) , ldep(0:lmax)

    integer :: is, ia, ias, lm, l,ig,m,igr
    integer :: ir


    !external functions
    Real (8) :: factnm
    External factnm

    mat(:,:,:)=zzero
    !write(*,*)"rspace metode ar matricu"



    if(present(yukawa_in))then
      yukawa=yukawa_in
    else
      yukawa=.false.
    endif

    if (yukawa)then
      do is=1, nspecies
        Do l = 0, lmax
          zf1=zzero
          do ir=1, nrmt(is)
            rr=spr(ir,is)/rmt(is)
            zf1(ir)=(spr(ir,is)/rmt(is))**l*(1d0-(spr(ir,is)/rmt(is))**2)**npsd*zbessi(ir,l,is)*spr(ir,is)**2
          enddo
          call integ_cv (nrmt(is), is, zf1(1:nrmt(is)),rinteg(l,is), mt_integw)
          !rinteg(l,is)=zf2(nrmt(is))
        enddo
      enddo
    endif


    if (yukawa)then
      Do l = 0, lmax
        zldep(l)=zlambda**l/factnm(2*l+1, 2)
      enddo
    else
      Do l = 0, lmax
        ldep(l)=factnm( 2*l + 2*npsd + 3, 2)/(2**(npsd)*factnm (npsd, 1)*factnm( 2*l+1, 2))
      enddo
    endif


    do is=1, nspecies
      do ia=1, natoms(is)
        ias=idxas(ia,is)
        ratom=atposc (:, ia, is)
        do igr=1, rgrid_nmtpoints(ias)
          ig=rgrid_mt_map(igr,ias)
          rv=rgrid_mt_rv(:,ias,igr)
          phase=exp(-cmplx(0,1,8)*sum(kvec*(rv+ratom)))
          rv_abs=rgrid_mt_rabs(ias,igr)

          if(yukawa) then
            lm=0
            Do l = 0, lmax
              zt1=zldep(l)*(rv_abs/rmt(is))**l*(1d0-(rv_abs/rmt(is))**2)**npsd/rinteg(l,is)
              if (rv_abs .Gt. input%structure%epslat) then
                Do m = - l, l
                  lm = lm + 1
                  mat(lm,igr,ias) = mat(lm,igr,ias) + zt1*rgrid_zylm(lm,ias,igr)*phase
                enddo !m
              else
                lm=1
                mat(lm,igr,ias) = mat(lm,igr,ias) + zt1*y00*phase
                exit ! exit l loop just to get the l=0 contribution
              endif
            enddo ! l



          else !if not Yukawa



            lm=0
            Do l = 0,lmax
              t1=ldep(l)*(rv_abs/rmt(is))**l*(1d0-rv_abs**2/rmt(is)**2)**npsd/(rmt(is)**(l+3))
              if (rv_abs .Gt. input%structure%epslat) then
                Do m = - l, l
                  lm = lm + 1
                  mat(lm,igr,ias) = mat(lm,igr,ias)+t1*rgrid_zylm(lm,ias,igr)*phase
                enddo ! m
              else
                lm=1
                mat(lm,igr,ias) = mat(lm,igr,ias) + t1*y00*phase
                exit   ! exit l loop just to get the l=0 contribution
              endif
            enddo ! l

          endif !if yukawa or not

        enddo !igr

      enddo !ia
    enddo !is

    ! ! do is=1, nspecies
    ! !   do ia=1, natoms(is)
    ! !     ias=idxas(ia,is)
    ! !     do igr=1, rgrid_nmtpoints(ias)
    ! !       ig=rgrid_mt_map(ias,igr)
    ! !       lm=0
    ! !       Do l = 0,lmax
    ! !         Do m = - l, l
    ! !           lm = lm+1
    ! !           zvclir(ig)=zvclir(ig) + qlm(lm,ias)*mat(lm,ias,igr)
    ! !         enddo ! m
    ! !       enddo
    ! !     enddo !igr
    ! !   enddo !ia
    ! ! enddo !is

    !   call zfftifc( 3, ngrid, -1, zvclir)


  end subroutine pseudocharge_rspace_matrix

  subroutine pseudocharge_rspace_new(lmax,qlm,rpseudomat,zvclir)
    use modrspace, only: rgrid_mt_map, rgrid_nmtpoints
    use mod_atoms, only: nspecies, natoms, natmtot,idxas
    use mod_Gvector, only: ngrtot, ngrid
    implicit none

    integer, intent(in) :: lmax
    complex(8), intent(in) :: qlm((lmax+1)**2,natmtot)
    complex(8), intent(in) :: rpseudomat(:,:,:)
    Complex (8), Intent (InOut) :: zvclir(ngrtot)
    complex (8) :: zt1
    integer :: is,ia,ias,l,m,lm,ig,igr,lmmax
    real(8) :: ta,tb,tc,td

    ! external functions
    Complex (8) :: zdotc
    External :: zdotc

    call timesec(ta)
    call zfftifc( 3, ngrid, 1, zvclir)
    call timesec(tb)

    lmmax=(lmax+1)**2
! !$OMP PARALLEL DEFAULT(NONE) PRIVATE(is,ia,ias,igr,ig,zt1) SHARED(idxas,nspecies,natoms,rgrid_mt_map,qlm,rpseudomat,rgrid_nmtpoints,zvclir)
    do is=1, nspecies
      do ia=1, natoms(is)
        ias=idxas(ia,is)

! !$OMP DO SCHEDULE(DYNAMIC)
        do igr=1, rgrid_nmtpoints(ias)
          ig=rgrid_mt_map(igr,ias)
          zt1=sum(qlm(:,ias)*rpseudomat(:,igr,ias))
          !zt1=zdotc(lmmax, qlm(:,ias), 1, rpseudomat(:,igr,ias), 1)
          zvclir(ig)=zvclir(ig)+zt1
        enddo !igr
! !$OMP END DO NOWAIT
      enddo !ia
    enddo !is
! !$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do is=1, nspecies
    !   do ia=1, natoms(is)
    !     ias=idxas(ia,is)
    !     do igr=1, rgrid_nmtpoints(ias)
    !       ig=rgrid_mt_map(igr,ias)
    !       lm=0
    !       Do l = 0,lmax
    !         Do m = - l, l
    !           lm = lm+1
    !           zvclir(ig)=zvclir(ig) + qlm(lm,ias)*rpseudomat(lm,igr,ias)
    !         enddo ! m
    !       enddo
    !     enddo !igr
    !   enddo !ia
    ! enddo !is


    call timesec(tc)
    call zfftifc( 3, ngrid, -1, zvclir)
    call timesec(td)
    !write(*,*)"fft 2x:",td-tc+tb-ta
    !write(*,*)"rekinasana :",tc-tb

  end subroutine pseudocharge_rspace_new

! This version will need reordering of array indices in jlgpr, ylmgp, sfacgp
  subroutine multipoles_ir5( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qi, zlambda,zilmt)
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas,natmtot
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt, idxlm
    use constants, only: fourpi,y00,zil,zzero
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qi(:,:)
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)

    integer :: is, ia, ias, ig, ifg, l, lm, m
    complex(dp) :: zt1,zt2,tbessi,zt3
    complex(4) :: ct1
    real(dp) :: t1,tbessj
    logical :: yukawa

    integer :: firstnonzeroG,lmmaxvr
    Complex (8) :: sf(natmtot)
    Complex (8) :: zlm ((lmax+1)**2), zlm2((lmax+1)**2)
    real(8):: rmtl (0:lmax+3, nspecies)

    integer, parameter :: blksize=64
    Complex (8) :: zmat1((lmax+1)**2,blksize) ,zmat2((lmax+1)**2,blksize)
    Complex (4) :: cmat1((lmax+1)**2,blksize) ,cmat2((lmax+1)**2,blksize)
    Complex (4) :: cmat3(blksize,(lmax+1)**2) ,cmat4(blksize,(lmax+1)**2)
    complex (4) :: cqi((lmax+1)**2,natmtot), csf(blksize,natmtot)
!    real(4) :: sjl(0:lmax+1,blksize)
    real(4) :: sjl(blksize,0:lmax+1)
    integer :: igoffset, chunksize, sfld

    !external functions
    Real (8) :: factnm
    External factnm

    lmmaxvr=(lmax+1)**2
    sfld=size(sfacgp,dim=1)

    if (present(zlambda)) then
      yukawa=.true.
    else
      yukawa=.false.
    endif




    cqi (:, :) = 0e0

    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    ! compute (R_mt)^l
    Do is = 1, nspecies
      rmtl (0, is) = 1.d0
      Do l = 1, input%groundstate%lmaxvr + 3
        rmtl (l, is) = rmtl (l-1, is) * rmt (is)
      End Do
    End Do

    chunksize=blksize
!x$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,zt1,t1,lm,m,l,is,ia,ias,zlm,zlm2,sf) SHARED(gpc,sfacgp,zvclir,rmt,ngvec,ylmgp,jlgpr,zil,nspecies,natoms,idxas,firstnonzeroG,lmax,lmmaxvr) REDUCTION(+:qi)
!x$OMP DO

! Is the alignment worth the hassle?
! The unaligned version goes as follows:
!          Do igoffset = 1, ngvec,blksize
!            if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1
!            Do ig=igoffset,igoffset+chunksize-1
!             t1 = 1.d0 / gpc(ig)
!             zt1 = zvclir (ig) * t1
!             Do lm=1,lmmaxvr
!               zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
!             End Do !lm
!            End Do
!            ...............
! If it makes any difference (for better or worse), it should be visible in large systems.
! For now, sticking with the aligned code.
    Do igoffset = 1, ngvec,blksize
      if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1

      if (gpc(igoffset).lt. input%structure%epslat) Then
        Do ig=igoffset,igoffset+chunksize-1
          If (.not.(gpc(ig) .lt. input%structure%epslat)) Then
            t1 = 1.d0 / gpc(ig)
          Else
            t1 = 0d0
          End If
          zt1 = zvclir (ig) * t1
          Do lm=1,lmmaxvr
!                zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
            cmat3(ig-igoffset+1,lm)=cmplx(zt1*conjg (ylmgp(lm, ig)),kind=4)
          End Do !lm
        End Do
      else
        Do ig=igoffset,igoffset+chunksize-1
          t1 = 1.d0 / gpc(ig)
!              zt1 = zvclir (ig) * t1
          ct1=cmplx(zvclir (ig) * t1, kind=4)
          Do lm=1,lmmaxvr
!                zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
!               cmat1(lm,ig-igoffset+1)=cmplx(zt1*conjg (ylmgp(lm, ig)),kind=4)
            cmat3(ig-igoffset+1,lm)=ct1*conjg (cmplx(ylmgp(lm, ig),kind=4))
          End Do !lm
        End Do
      endif
      csf(1:chunksize,1:natmtot)=cmplx(sfacgp(igoffset:igoffset+chunksize-1,1:natmtot),kind=4)
! End of the differing part

      Do is = 1, nspecies

        Do ig=igoffset,igoffset+chunksize-1
          Do l = 0, lmax
            sjl(ig-igoffset+1,l)=real(jlgpr (l+1, ig, is),kind=4)
          End Do
        End Do
!                 lm = 0
        Do l = 0, lmax
          Do lm=idxlm(l,-l),idxlm(l,l)
            Do ig=1,chunksize
              cmat4(ig,lm)= sjl(ig,l) * cmat3(ig,lm)
            End Do ! ig
          End Do
!                   lm=lm+1+2*l
        End Do !l
        Call cgemm('T','N', lmmaxvr, natoms(is), chunksize, (1.0e0,0.0), cmat4, blksize, csf(1,idxas (1, is)), blksize, (1.0e0,0.0), cqi(1,idxas (1, is)), lmmaxvr)

      End Do !is

    End Do
!x$OMP END DO
!x$OMP END PARALLEL
    qi=cmplx(cqi,kind=8)

    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        lm = 0
        do l = 0, lmax
          qi (lm+1:lm+1+2*l,ias) = qi (lm+1:lm+1+2*l,ias) * fourpi * zil (l) * rmtl (l+3, is) / rmt(is) ! šo va riznest ārā
          lm=lm+1+2*l
        enddo
      End Do
    End Do

    if (gpc(1) .lt. input%structure%epslat) Then

      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          t1 = fourpi * y00 * rmtl (3, is) / 3.d0
          qi (1,ias) = qi (1,ias) + t1 * zvclir (1)
        End Do
      End Do
    endif





  end subroutine multipoles_ir5

  subroutine multipoles_ir4( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qi, zlambda,zilmt)
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas,natmtot
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil,zzero
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qi(:,:)
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)

    integer :: is, ia, ias, ig, ifg, l, lm, m
    complex(dp) :: zt1,zt2,tbessi,zt3
    complex(4) :: ct1
    real(dp) :: t1,tbessj
    logical :: yukawa

    integer :: firstnonzeroG,lmmaxvr
    Complex (8) :: sf(natmtot)
    Complex (8) :: zlm ((lmax+1)**2), zlm2((lmax+1)**2)
    real(8):: rmtl (0:lmax+3, nspecies)

    integer, parameter :: blksize=64
    Complex (8) :: zmat1((lmax+1)**2,blksize) ,zmat2((lmax+1)**2,blksize)
    Complex (4) :: cmat1((lmax+1)**2,blksize) ,cmat2((lmax+1)**2,blksize)
    complex (4) :: cqi((lmax+1)**2,natmtot), csf(blksize,natmtot)
    real(4) :: sjl(0:lmax+1,blksize)
    integer :: igoffset, chunksize, sfld

    !external functions
    Real (8) :: factnm
    External factnm

    lmmaxvr=(lmax+1)**2
    sfld=size(sfacgp,dim=1)

    if (present(zlambda)) then
      yukawa=.true.
    else
      yukawa=.false.
    endif




    cqi (:, :) = 0e0

    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    ! compute (R_mt)^l
    Do is = 1, nspecies
      rmtl (0, is) = 1.d0
      Do l = 1, input%groundstate%lmaxvr + 3
        rmtl (l, is) = rmtl (l-1, is) * rmt (is)
      End Do
    End Do

    chunksize=blksize
!x$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,zt1,t1,lm,m,l,is,ia,ias,zlm,zlm2,sf) SHARED(gpc,sfacgp,zvclir,rmt,ngvec,ylmgp,jlgpr,zil,nspecies,natoms,idxas,firstnonzeroG,lmax,lmmaxvr) REDUCTION(+:qi)
!x$OMP DO

! Is the alignment worth the hassle?
! The unaligned version goes as follows:
!          Do igoffset = 1, ngvec,blksize
!            if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1
!            Do ig=igoffset,igoffset+chunksize-1
!             t1 = 1.d0 / gpc(ig)
!             zt1 = zvclir (ig) * t1
!             Do lm=1,lmmaxvr
!               zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
!             End Do !lm
!            End Do
!            ...............
! If it makes any difference (for better or worse), it should be visible in large systems.
! For now, sticking with the aligned code.
    Do igoffset = 1, ngvec,blksize
      if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1

      if (gpc(igoffset).lt. input%structure%epslat) Then
        Do ig=igoffset,igoffset+chunksize-1
          If (.not.(gpc(ig) .lt. input%structure%epslat)) Then
            t1 = 1.d0 / gpc(ig)
          Else
            t1 = 0d0
          End If
          zt1 = zvclir (ig) * t1
          Do lm=1,lmmaxvr
            cmat1(lm,ig-igoffset+1)=cmplx(zt1*conjg (ylmgp(lm, ig)),kind=4)
          End Do !lm
        End Do
      else
        Do ig=igoffset,igoffset+chunksize-1
          t1 = 1.d0 / gpc(ig)
          ct1=cmplx(zvclir (ig) * t1, kind=4)
          Do lm=1,lmmaxvr
            cmat1(lm,ig-igoffset+1)=ct1*conjg (cmplx(ylmgp(lm, ig),kind=4))
          End Do !lm
        End Do
      endif
      csf(1:chunksize,1:natmtot)=cmplx(sfacgp(igoffset:igoffset+chunksize-1,1:natmtot),kind=4)
! End of the differing part

      Do is = 1, nspecies
        Do ig=igoffset,igoffset+chunksize-1
          sjl(0:lmax,ig-igoffset+1)=real(jlgpr (1:lmax+1, ig, is),kind=4)
        End Do
        Do ig=1,chunksize
          lm = 0
          Do l = 0, lmax
            cmat2(lm+1:lm+1+2*l,ig)= sjl(l,ig) * cmat1(lm+1:lm+1+2*l,ig)
            lm=lm+1+2*l
          End Do !l
        End Do ! ig
        Call cgemm('N','N', lmmaxvr, natoms(is), chunksize, (1.0e0,0.0), cmat2, lmmaxvr, csf(1,idxas (1, is)), blksize, (1.0e0,0.0), cqi(1,idxas (1, is)), lmmaxvr)

      End Do !is

    End Do
!x$OMP END DO
!x$OMP END PARALLEL
    qi=cmplx(cqi,kind=8)

    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        lm = 0
        do l = 0, lmax
          qi (lm+1:lm+1+2*l,ias) = qi (lm+1:lm+1+2*l,ias) * fourpi * zil (l) * rmtl (l+3, is) / rmt(is) ! šo va riznest ārā
          lm=lm+1+2*l
        enddo
      End Do
    End Do

    if (gpc(1) .lt. input%structure%epslat) Then

      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          t1 = fourpi * y00 * rmtl (3, is) / 3.d0
          qi (1,ias) = qi (1,ias) + t1 * zvclir (1)
        End Do
      End Do
    endif

  end subroutine multipoles_ir4

  subroutine multipoles_ir3( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qi, zlambda,zilmt)
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas,natmtot
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil,zzero
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qi(:,:)
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)

    integer :: is, ia, ias, ig, ifg, l, lm, m
    complex(dp) :: zt1,zt2,tbessi,zt3
    real(dp) :: t1,tbessj
    logical :: yukawa

    integer :: firstnonzeroG,lmmaxvr
    Complex (8) :: sf(natmtot)
    Complex (8) :: zlm ((lmax+1)**2), zlm2((lmax+1)**2)
    real(8):: rmtl (0:lmax+3, nspecies)

    integer, parameter :: blksize=64
    Complex (8) :: zmat1((lmax+1)**2,blksize) ,zmat2((lmax+1)**2,blksize)
    integer :: igoffset, chunksize, sfld

    !external functions
    Real (8) :: factnm
    External factnm

    lmmaxvr=(lmax+1)**2
    sfld=size(sfacgp,dim=1)

    if (present(zlambda)) then
      yukawa=.true.
    else
      yukawa=.false.
    endif




    qi (:, :) = zzero

    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    if (.not.yukawa) then
      ! compute (R_mt)^l
      Do is = 1, nspecies
        rmtl (0, is) = 1.d0
        Do l = 1, input%groundstate%lmaxvr + 3
          rmtl (l, is) = rmtl (l-1, is) * rmt (is)
        End Do
      End Do

      chunksize=blksize
!x$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,zt1,t1,lm,m,l,is,ia,ias,zlm,zlm2,sf) SHARED(gpc,sfacgp,zvclir,rmt,ngvec,ylmgp,jlgpr,zil,nspecies,natoms,idxas,firstnonzeroG,lmax,lmmaxvr) REDUCTION(+:qi)
!x$OMP DO

! Is the alignment worth the hassle?
! The unaligned version goes as follows:
!          Do igoffset = 1, ngvec,blksize
!            if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1
!            Do ig=igoffset,igoffset+chunksize-1
!             t1 = 1.d0 / gpc(ig)
!             zt1 = zvclir (ig) * t1
!             Do lm=1,lmmaxvr
!               zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
!             End Do !lm
!            End Do
!            ...............
! If it makes any difference (for better or worse), it should be visible in large systems.
! For now, sticking with the aligned code.
      Do igoffset = 1, ngvec,blksize
        if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1

        if (gpc(igoffset).lt. input%structure%epslat) Then
          Do ig=igoffset,igoffset+chunksize-1
            If (.not.(gpc(ig) .lt. input%structure%epslat)) Then
              t1 = 1.d0 / gpc(ig)
            Else
              t1 = 0d0
            End If
            zt1 = zvclir (ig) * t1
            Do lm=1,lmmaxvr
              zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
            End Do !lm
          End Do
        else
          Do ig=igoffset,igoffset+chunksize-1
            t1 = 1.d0 / gpc(ig)
            zt1 = zvclir (ig) * t1
            Do lm=1,lmmaxvr
              zmat1(lm,ig-igoffset+1)=zt1*conjg (ylmgp(lm, ig))
            End Do !lm
          End Do
        endif
! End of the differing part

        Do is = 1, nspecies
          Do ig=igoffset,igoffset+chunksize-1
            lm = 0
            Do l = 0, lmax
              zmat2(lm+1:lm+1+2*l,ig-igoffset+1)=jlgpr (l+1, ig, is) * zmat1(lm+1:lm+1+2*l,ig-igoffset+1)
              lm=lm+1+2*l
            End Do !l
          End Do ! ig
          Call zgemm('N','N', lmmaxvr, natoms(is), chunksize, (1.0D0,0.0), zmat2, lmmaxvr, sfacgp(igoffset,idxas (1, is)), sfld, (1.0D0,0.0), qi(1,idxas (1, is)), lmmaxvr)
        End Do !is

      End Do
!x$OMP END DO
!x$OMP END PARALLEL

      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm = 0
          do l = 0, lmax
            qi (lm+1:lm+1+2*l,ias) = qi (lm+1:lm+1+2*l,ias) * fourpi * zil (l) * rmtl (l+3, is) / rmt(is)
            lm=lm+1+2*l
          enddo
        End Do
      End Do

      if (gpc(1) .lt. input%structure%epslat) Then

        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = fourpi * y00 * rmtl (3, is) / 3.d0
            qi (1,ias) = qi (1,ias) + t1 * zvclir (1)
          End Do
        End Do
      endif




    else ! yukawa

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,zt1,t1,lm,m,l,is,ia,ias,zlm,zlm2,sf) SHARED(gpc,zlambda,sfacgp,zvclir,rmt,ngvec,ylmgp,jlgpr,zil,zilmt,nspecies,natoms,idxas,firstnonzeroG,lmax,lmmaxvr) REDUCTION(+:qi)
!$OMP DO
      Do ig = firstnonzeroG, ngvec
        sf(:)=sfacgp(ig,:)
        zt1 = 1.d0 / (gpc(ig)**2 + zlambda**2)
        zt1 = zvclir (ig) * zt1
        zlm(:)=zt1*conjg (ylmgp(:, ig))

        Do is = 1, nspecies
          zlm2(1)=(zlambda * jlgpr (0, ig, is) * (zilmt(1,is) + (zilmt(0,is)/(zlambda*rmt(is)))))&
            &- (gpc(ig) * ((jlgpr(0, ig, is)/(gpc(ig)*rmt(is)))-jlgpr(1, ig, is)) * zilmt(0,is))
          zlm2(1) = zlm2(1) * zlm(1)

          lm = 1
          Do l = 1, lmax
            zlm2(lm+1:lm+1+2*l) = (zlambda * jlgpr (l, ig, is) * zilmt(l-1,is))-(gpc(ig) * jlgpr (l-1, ig, is) * zilmt(l,is))
            zlm2(lm+1:lm+1+2*l)=zlm2(lm+1:lm+1+2*l)*zlm(lm+1:lm+1+2*l)
            lm=lm+1+2*l
          End Do !l

          Do ia = 1, natoms(is)
            ias = idxas (ia, is)
            qi (1:lmmaxvr,ias) = qi (1:lmmaxvr,ias) + sf (ias) * zlm2(1:lmmaxvr)
          End Do !ia
        End Do !is
      End Do !ig
!$OMP END DO
!$OMP END PARALLEL

      Do is = 1, nspecies
        t1=fourpi*rmt(is)**2
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm = 0
          do l = 0, lmax
            qi (lm+1:lm+1+2*l,ias) = qi (lm+1:lm+1+2*l,ias) * t1 * zil (l) * zlambda**(-l) * factnm (2*l+1, 2)
            lm=lm+1+2*l
          enddo
        End Do
      End Do

      if (gpc(1) .lt. input%structure%epslat) Then


        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)

            zt1 = fourpi * y00 * (rmt (is) ** 2) * zilmt(1,is) / zlambda
            qi (1,ias) = qi (1,ias) + zt1 * zvclir (1)

          End Do
        End Do
      endif

    endif! yukawa or not

  end subroutine multipoles_ir3

  subroutine multipoles_ir2( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qi, zlambda,zilmt)
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas,natmtot
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil,zzero
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(out) :: qi(:,:)
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)

    integer :: is, ia, ias, ig, ifg, l, lm, m
    complex(dp) :: zt1,zt2,tbessi,zt3
    real(dp) :: t1,tbessj
    logical :: yukawa

    integer :: firstnonzeroG,lmmaxvr
    Complex (8) :: sf(natmtot)
    Complex (8) :: zlm ((lmax+1)**2), zlm2((lmax+1)**2)
    real(8):: rmtl (0:lmax+3, nspecies)
    !external functions
    Real (8) :: factnm
    External factnm

    lmmaxvr=(lmax+1)**2

    if (present(zlambda)) then
      yukawa=.true.
    else
      yukawa=.false.
    endif




    qi (:, :) = zzero

    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    if (.not.yukawa) then
      ! compute (R_mt)^l
      Do is = 1, nspecies
        rmtl (0, is) = 1.d0
        Do l = 1, input%groundstate%lmaxvr + 3
          rmtl (l, is) = rmtl (l-1, is) * rmt (is)
        End Do
      End Do

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,zt1,t1,lm,m,l,is,ia,ias,zlm,zlm2,sf) SHARED(gpc,sfacgp,zvclir,rmt,ngvec,ylmgp,jlgpr,zil,nspecies,natoms,idxas,firstnonzeroG,lmax,lmmaxvr) REDUCTION(+:qi)
!$OMP DO
      Do ig = firstnonzeroG, ngvec
        sf(:)=sfacgp(ig,:)
        t1 = 1.d0 / gpc(ig)
        zt1 = zvclir (ig) * t1
        Do lm=1,lmmaxvr
          zlm(lm)=zt1*conjg (ylmgp(lm, ig))
        End Do !lm

        Do is = 1, nspecies
          lm = 0
          Do l = 0, lmax
            zlm2(lm+1:lm+1+2*l) = jlgpr (l+1, ig, is) * zlm(lm+1:lm+1+2*l)
            lm=lm+1+2*l
          End Do !l

          Do ia = 1, natoms(is)
            ias = idxas (ia, is)
            Do lm=1,lmmaxvr
              qi (lm,ias) = qi (lm,ias) + sf (ias) * zlm2(lm)
            End Do !lm
          End Do !ia
        End Do !is
      End Do !ig
!$OMP END DO
!$OMP END PARALLEL

      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm = 0
          do l = 0, lmax
            qi (lm+1:lm+1+2*l,ias) = qi (lm+1:lm+1+2*l,ias) * fourpi * zil (l) * rmtl (l+3, is) / rmt(is) ! šo va riznest ārā
            lm=lm+1+2*l
          enddo
        End Do
      End Do

      if (gpc(1) .lt. input%structure%epslat) Then

        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = fourpi * y00 * rmtl (3, is) / 3.d0
            qi (1,ias) = qi (1,ias) + t1 * zvclir (1)
          End Do
        End Do
      endif




    else ! yukawa

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,zt1,t1,lm,m,l,is,ia,ias,zlm,zlm2,sf) SHARED(gpc,zlambda,sfacgp,zvclir,rmt,ngvec,ylmgp,jlgpr,zil,zilmt,nspecies,natoms,idxas,firstnonzeroG,lmax,lmmaxvr) REDUCTION(+:qi)
!$OMP DO
      Do ig = firstnonzeroG, ngvec
        sf(:)=sfacgp(ig,:)
        zt1 = 1.d0 / (gpc(ig)**2 + zlambda**2)
        zt1 = zvclir (ig) * zt1
        zlm(:)=zt1*conjg (ylmgp(:, ig))

        Do is = 1, nspecies
          zlm2(1)=(zlambda * jlgpr (0, ig, is) * (zilmt(1,is) + (zilmt(0,is)/(zlambda*rmt(is)))))&
            &- (gpc(ig) * ((jlgpr(0, ig, is)/(gpc(ig)*rmt(is)))-jlgpr(1, ig, is)) * zilmt(0,is))
          zlm2(1) = zlm2(1) * zlm(1)

          lm = 1
          Do l = 1, lmax
            zlm2(lm+1:lm+1+2*l) = (zlambda * jlgpr (l, ig, is) * zilmt(l-1,is))-(gpc(ig) * jlgpr (l-1, ig, is) * zilmt(l,is))
            zlm2(lm+1:lm+1+2*l)=zlm2(lm+1:lm+1+2*l)*zlm(lm+1:lm+1+2*l)
            lm=lm+1+2*l
          End Do !l

          Do ia = 1, natoms(is)
            ias = idxas (ia, is)
            qi (1:lmmaxvr,ias) = qi (1:lmmaxvr,ias) + sf (ias) * zlm2(1:lmmaxvr)
          End Do !ia
        End Do !is
      End Do !ig
!$OMP END DO
!$OMP END PARALLEL

      Do is = 1, nspecies
        t1=fourpi*rmt(is)**2
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm = 0
          do l = 0, lmax
            qi (lm+1:lm+1+2*l,ias) = qi (lm+1:lm+1+2*l,ias) * t1 * zil (l) * zlambda**(-l) * factnm (2*l+1, 2)
            lm=lm+1+2*l
          enddo
        End Do
      End Do

      if (gpc(1) .lt. input%structure%epslat) Then


        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)

            zt1 = fourpi * y00 * (rmt (is) ** 2) * zilmt(1,is) / zlambda
            qi (1,ias) = qi (1,ias) + zt1 * zvclir (1)

          End Do
        End Do
      endif

    endif! yukawa or not

  end subroutine multipoles_ir2


  subroutine pseudocharge_gspace2( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qlm, zlambda, zilmt)
    use mod_lattice, only: omega
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas, natmtot
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(inout) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(in) :: qlm(:,:)
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)

    complex(dp) :: zrp((lmax+1)**2)


    integer :: is, ia, ias, ig, ifg, l, lm, m,npsd
    complex(dp) :: zt1,zsum1, zsum2
    real(dp) :: t1,t2, t3,z1,fpo

    logical :: yukawa
    integer :: lmmaxvr, firstnonzeroG
    complex(dp) :: zrp2((lmax+1)**2,natmtot), sf(natmtot)
    Complex (8) :: zlm ((lmax+1)**2), zlm2((lmax+1)**2), zlm3((lmax+1)**2)
    !external functions
    Real (8) :: factnm, one_over_rmtl(0:lmax,nspecies)
    External factnm
    lmmaxvr=(lmax+1)**2
    if (present(zlambda))then
      yukawa=.true.
    else
      yukawa=.false.
    endif

    fpo = fourpi/omega
    npsd=input%groundstate%npsden

    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    if(.not.yukawa) then
      ! compute (R_mt)^l
      Do is = 1, nspecies
        Do l = 0, lmax
          one_over_rmtl(l, is) = rmt(is)**(-l)
        End Do
      End Do



      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm = 0
          Do l = 0, lmax
            t1 = factnm (2*(l+npsd)+3, 2) / factnm (2*l+1, 2)
            t1 = t1 * one_over_rmtl(l, is)
            Do m = - l, l
              lm = lm + 1
              zrp2 (lm,ias) = qlm(lm, ias) * t1 * fpo * conjg (zil(l))
            End Do !m
          End Do !l
        End Do !ia
      End Do !is





      ! add the pseudocharge and real interstitial densities in G-space
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,t1,t2,l,lm,t3,is,ia,ias,sf,zlm,zlm2,zlm3) SHARED(lmax,gpc,zvclir,rmt,ngvec,ylmgp,jlgpr,zrp2,sfacgp,firstnonzeroG,idxas,nspecies,natoms,lmmaxvr,npsd)
!$OMP DO
      Do ig = firstnonzeroG, ngvec
        sf(:)=conjg(sfacgp(ig,:))
        Do is = 1, nspecies
          t1 = gpc (ig) * rmt (is)
          t2 = 1d0/ (t1 ** (npsd+1))
          zlm(:) = t2 * ylmgp (:, ig)

          zlm2 = 0d0
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)

            zlm2(:)=zlm2(:)+zrp2(:,ias)*sf(ias)
          End Do

          lm = 0
          Do l = 0, lmax
            t3 = jlgpr (npsd+l+1, ig, is)
            zlm3(lm+1:lm+1+2*l) = t3 * zlm (lm+1:lm+1+2*l)
            lm=lm+1+2*l
          End Do
          Do lm = 1, lmmaxvr
            zvclir (ig) = zvclir (ig) + zlm2(lm) * zlm3(lm)
          End Do
        End Do !is
      End Do !ig
!$OMP END DO
!$OMP END PARALLEL



      if (firstnonzeroG.eq.2) then
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = y00 / factnm (2*npsd+3, 2)
            zvclir (1) = zvclir (1) + t1 * zrp2 (1,ias)
          End Do
        End Do
      endif




    else !if yukawa


      do is=1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm=0
          do l=0, lmax
            t1=fpo / factnm (2*l+1, 2)
            zt1=t1*zlambda**(npsd+l+1) * conjg(zil(l))
            do m=-l,l
              lm = lm + 1
              zrp2 (lm,ias)= zt1 * qlm(lm, ias) / zilmt(npsd+l+1,is)
            enddo !m
          enddo !l
        enddo!ia
      enddo!is





      ! add the pseudocharge and real interstitial densities in G-space
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,t1,l,lm,t3,is,ia,ias,sf,zlm,zlm2,zlm3) SHARED(lmax,gpc,zvclir,rmt,ngvec,ylmgp,jlgpr,zrp2,sfacgp,firstnonzeroG,idxas,nspecies,natoms,lmmaxvr,npsd)
!$OMP DO
      Do ig = firstnonzeroG, ngvec
        sf(:)=conjg(sfacgp(ig,:))
        t1 = gpc (ig) ** (-npsd-1)
        Do is = 1, nspecies


          zlm(:) = t1 * ylmgp (:, ig)

          zlm2 = 0d0
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zlm2(:)=zlm2(:)+zrp2(:,ias)*sf(ias)
          End Do

          lm = 0
          Do l = 0, lmax
            t3 = jlgpr (npsd+l+1, ig, is)
            zlm3(lm+1:lm+1+2*l) = t3 * zlm (lm+1:lm+1+2*l)
            lm=lm+1+2*l
          End Do
          Do lm = 1, lmmaxvr
            zvclir (ig) = zvclir (ig) + zlm2(lm) * zlm3(lm)
          End Do
        End Do !is
      End Do !ig
!$OMP END DO
!$OMP END PARALLEL



      if (firstnonzeroG.eq.2) then
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = rmt(is)**(npsd+1)*y00 / factnm (2*npsd+3, 2)
            zvclir (1) = zvclir (1) + t1 * zrp2 (1,ias)
          End Do
        End Do
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   do is = 1, nspecies
      !     do ia = 1, natoms(is)
      !       ias = idxas( ia, is)
      !       lm=0
      !       Do l = 0, lmax
      !         t1 = 1d0 / (factnm (2*l+1, 2))
      !         Do m = - l, l
      !             lm = lm + 1
      !             zrp (lm) = qlm(lm, ias) * t1
      !             !write(*,*)"qlm",qlm(lm, ias),t1
      !         End Do
      !       End Do

      !       Do ig = 1, ngvec
      !         ifg = ig ! igfft (ig) !becouse zvclir is sorted
      !         If (gpc(ig) .Gt. input%structure%epslat) Then
      !           z1 = gpc (ig)
      !           zt1 = fpo * conjg (sfacgp(ig, ias))/(z1 ** (npsd+1))
      !           lm = 0
      !           Do l = 0, lmax
      !               zsum1 = (jlgpr (npsd+l+1, ig, is)*(zlambda ** (l+npsd+1))/ zilmt(npsd+l+1,is))* conjg(zil(l))
      !               lm = lm + 1
      !               zsum2 = zrp (lm) * ylmgp (lm, ig)
      !               Do m = - l + 1, l
      !                 lm = lm + 1
      !                 zsum2 = zsum2 + zrp (lm) * ylmgp (lm, ig)
      !               End Do !m
      !               zvclir (ifg) = zvclir (ifg) +  zt1 * zsum1 * zsum2
      !           End Do !l

      !         Else
      !           zt1 = (fpo * y00 * (rmt(is) ** (npsd+1))&
      !             &* (zlambda ** (npsd+1))) / ((factnm (2*npsd+3, 2))* zilmt(npsd+1,is))

      !           zvclir (ifg) = zvclir (ifg) + zt1 * qlm(1, ias)
      !         End If
      !       End Do! ig
      !     End Do !ia
      !   End Do !is
    endif !if yukawa or not
    ! do ig=1,20
    !   write(*,*)zvclir (ig)
    ! enddo
    ! stop
  end subroutine pseudocharge_gspace2

  subroutine pseudocharge_gspace3( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qlm, zlambda, zilmt)
    use mod_lattice, only: omega
    use modinput
    use mod_atoms, only: nspecies, natoms, idxas, natmtot
    use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi,y00,zil
    !> maximum angular momentum \(l\)
    integer, intent(in) :: lmax
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
    !integer, intent(in) :: ivgp(:,:)
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> bounds for integer components of \({\bf G}\)
    !integer, intent(in) :: intgv(3,2)
    !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
    !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(inout) :: zvclir(:)
    !> multipole moments of the function's extension inside the muffin-tin spheres
    complex(dp), intent(in) :: qlm(:,:)
    complex(dp),optional, intent(in) :: zlambda
    complex(dp),optional, intent(in) :: zilmt(0:,:)

    complex(dp) :: zrp((lmax+1)**2)


    integer :: is, ia, ias, ig, ifg, l, lm, m,npsd
    complex(dp) :: zt1,zsum1, zsum2
    real(dp) :: t1,t2, t3,z1,fpo

    logical :: yukawa
    integer :: lmmaxvr, firstnonzeroG
    complex(dp) :: zrp2((lmax+1)**2,natmtot), sf(natmtot)
    Complex (8) :: zlm ((lmax+1)**2), zlm2((lmax+1)**2), zlm3((lmax+1)**2)
    !external functions
    Real (8) :: one_over_rmtl(0:lmax,nspecies)
    Real (8), external :: factnm
    integer, parameter :: blksize=64
    Real (8) :: prefg(blksize)
    Real (8) :: prefsp(nspecies)
    Complex (8) :: zmat1((lmax+1)**2,blksize)
    integer :: igoffset, chunksize, sfld
    Complex (8), external :: zdotu

    lmmaxvr=(lmax+1)**2
    if (present(zlambda))then
      yukawa=.true.
    else
      yukawa=.false.
    endif
    sfld=size(sfacgp,1)

    fpo = fourpi/omega
    npsd=input%groundstate%npsden

    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    if(.not.yukawa) then
      ! compute (R_mt)^l
      Do is = 1, nspecies
        Do l = 0, lmax
          one_over_rmtl(l, is) = rmt(is)**(-l)
        End Do
      End Do



      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm = 0
          Do l = 0, lmax
            t1 = factnm (2*(l+npsd)+3, 2) / factnm (2*l+1, 2)
            t1 = t1 * one_over_rmtl(l, is)
            Do m = - l, l
              lm = lm + 1
              zrp2 (lm,ias) = qlm(lm, ias) * t1 * fpo * conjg (zil(l))
            End Do !m
          End Do !l
        End Do !ia
      End Do !is

      Do is = 1, nspecies
        prefsp(is)=1d0/ (rmt(is)**(npsd+1))
      End Do



      ! add the pseudocharge and real interstitial densities in G-space
      chunksize=blksize
      Do igoffset = firstnonzeroG, ngvec,blksize
        if (igoffset+blksize-1.gt.ngvec) chunksize=ngvec-igoffset+1
        Do ig=0,chunksize-1
          prefg(ig+1) = 1d0/ ( gpc (ig+igoffset) ** (npsd+1))
        End Do
        Do is = 1, nspecies
          Call zgemm('N','C', lmmaxvr, chunksize, natoms(is), (1.0D0,0.0), zrp2(1,idxas(1,is)), lmmaxvr, sfacgp(igoffset,idxas (1, is)), sfld, (0.0D0,0.0), zmat1, lmmaxvr)
          Do ig=0,chunksize-1
            t2=prefsp(is)*prefg(ig+1)

            lm = 0
            Do l = 0, lmax
              zlm(lm+1:lm+1+2*l) = jlgpr (npsd+l+1, ig+igoffset, is) * ylmgp (lm+1:lm+1+2*l, ig+igoffset) !zlm (lm+1:lm+1+2*l)
              lm=lm+1+2*l
            End Do
            zvclir (ig+igoffset) = zvclir (ig+igoffset) + zdotu(lmmaxvr, zlm, 1, zmat1(1,ig+1), 1) * t2
          End Do !ig
        End Do !is
      End Do !ig

      if (firstnonzeroG.eq.2) then
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = y00 / factnm (2*npsd+3, 2)
            zvclir (1) = zvclir (1) + t1 * zrp2 (1,ias)
          End Do
        End Do
      endif




    else !if yukawa


      do is=1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          lm=0
          do l=0, lmax
            t1=fpo / factnm (2*l+1, 2)
            zt1=t1*zlambda**(npsd+l+1) * conjg(zil(l))
            do m=-l,l
              lm = lm + 1
              zrp2 (lm,ias)= zt1 * qlm(lm, ias) / zilmt(npsd+l+1,is)
            enddo !m
          enddo !l
        enddo!ia
      enddo!is





      ! add the pseudocharge and real interstitial densities in G-space
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,t1,l,lm,t3,is,ia,ias,sf,zlm,zlm2,zlm3) SHARED(lmax,gpc,zvclir,rmt,ngvec,ylmgp,jlgpr,zrp2,sfacgp,firstnonzeroG,idxas,nspecies,natoms,lmmaxvr,npsd)
!$OMP DO
      Do ig = firstnonzeroG, ngvec
        sf(:)=conjg(sfacgp(ig,:))
        t1 = gpc (ig) ** (-npsd-1)
        Do is = 1, nspecies


          zlm(:) = t1 * ylmgp (:, ig)

          zlm2 = 0d0
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zlm2(:)=zlm2(:)+zrp2(:,ias)*sf(ias)
          End Do

          lm = 0
          Do l = 0, lmax
            t3 = jlgpr (npsd+l+1, ig, is)
            zlm3(lm+1:lm+1+2*l) = t3 * zlm (lm+1:lm+1+2*l)
            lm=lm+1+2*l
          End Do
          Do lm = 1, lmmaxvr
            zvclir (ig) = zvclir (ig) + zlm2(lm) * zlm3(lm)
          End Do
        End Do !is
      End Do !ig
!$OMP END DO
!$OMP END PARALLEL



      if (firstnonzeroG.eq.2) then
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = rmt(is)**(npsd+1)*y00 / factnm (2*npsd+3, 2)
            zvclir (1) = zvclir (1) + t1 * zrp2 (1,ias)
          End Do
        End Do
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   do is = 1, nspecies
      !     do ia = 1, natoms(is)
      !       ias = idxas( ia, is)
      !       lm=0
      !       Do l = 0, lmax
      !         t1 = 1d0 / (factnm (2*l+1, 2))
      !         Do m = - l, l
      !             lm = lm + 1
      !             zrp (lm) = qlm(lm, ias) * t1
      !             !write(*,*)"qlm",qlm(lm, ias),t1
      !         End Do
      !       End Do

      !       Do ig = 1, ngvec
      !         ifg = ig ! igfft (ig) !becouse zvclir is sorted
      !         If (gpc(ig) .Gt. input%structure%epslat) Then
      !           z1 = gpc (ig)
      !           zt1 = fpo * conjg (sfacgp(ig, ias))/(z1 ** (npsd+1))
      !           lm = 0
      !           Do l = 0, lmax
      !               zsum1 = (jlgpr (npsd+l+1, ig, is)*(zlambda ** (l+npsd+1))/ zilmt(npsd+l+1,is))* conjg(zil(l))
      !               lm = lm + 1
      !               zsum2 = zrp (lm) * ylmgp (lm, ig)
      !               Do m = - l + 1, l
      !                 lm = lm + 1
      !                 zsum2 = zsum2 + zrp (lm) * ylmgp (lm, ig)
      !               End Do !m
      !               zvclir (ifg) = zvclir (ifg) +  zt1 * zsum1 * zsum2
      !           End Do !l

      !         Else
      !           zt1 = (fpo * y00 * (rmt(is) ** (npsd+1))&
      !             &* (zlambda ** (npsd+1))) / ((factnm (2*npsd+3, 2))* zilmt(npsd+1,is))

      !           zvclir (ifg) = zvclir (ifg) + zt1 * qlm(1, ias)
      !         End If
      !       End Do! ig
      !     End Do !ia
      !   End Do !is
    endif !if yukawa or not
    ! do ig=1,20
    !   write(*,*)zvclir (ig)
    ! enddo
    ! stop
  end subroutine pseudocharge_gspace3

  subroutine poisson_ir2(ngvec, gpc, igfft, zrhoir,zvclir,cutoff,zlambda)
    use mod_lattice, only: omega
    Use mod_kpoint, only: nkptnr
    use modinput
    use constants, only: zzero
    ! use mod_atoms, only: nspecies, natoms, idxas
    ! use mod_muffin_tin, only: rmt,nrmtmax,nrmt
    use constants, only: fourpi

    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngvec
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> map from \({\bf G}\) vector index to point in FFT grid
    integer, intent(in) :: igfft(:)
    !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
    complex(dp), intent(in) :: zrhoir(:)

    complex(dp),optional, intent(in) :: zlambda
    complex(dp), intent(out) :: zvclir(:)
    logical, intent(in) :: cutoff


    integer :: ig, ifg, firstnonzeroG
    real(dp) :: r_c
    !external functions
    Real (8) :: factnm
    External factnm
    logical :: yukawa

    if (present(zlambda)) then
      yukawa=.true.
    else
      yukawa=.false.
    endif

    zvclir=zzero
    r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0
    if (gpc(1) .lt. input%structure%epslat) Then
      firstnonzeroG=2
    else
      firstnonzeroG=1
    endif

    if (.not.yukawa)then
      If (cutoff) Then
        Do ig = firstnonzeroG, ngvec
          zvclir (ig) = fourpi * zrhoir (ig)*(1d0-cos(gpc(ig) * r_c )) / (gpc(ig)**2)
        End Do
        if (firstnonzeroG.eq.2) then
          zvclir (1) = zrhoir(1)*(fourpi*0.5d0)*r_c**2
        endif
      Else! cutof
        Do ig = firstnonzeroG, ngvec
          zvclir (ig) = fourpi * zrhoir (ig) / (gpc(ig)**2)
        End Do
        if (firstnonzeroG.eq.2) then
          zvclir (1) = 0.d0
        endif
      End If !if cutoff
    else ! if yukawa



      if (cutoff) then
        Do ig = firstnonzeroG, ngvec
          zvclir (ig) = fourpi * zrhoir (ig) *&
            (1d0 - exp(-zlambda*r_c)* ( zlambda * sin(gpc(ig)*r_c) / gpc(ig) + cos(gpc(ig)*r_c) ) )/&
            ((gpc(ig) ** 2) + (zlambda ** 2))
        End Do
        if (firstnonzeroG.eq.2) then
          zvclir (1) = fourpi * zrhoir(1)* (1d0 - exp(-zlambda*r_c)*(zlambda*r_c + 1))/(zlambda ** 2)
        endif

        ! If (gpc(ig) .Gt. input%structure%epslat) Then
        !    zvclir (ifg) = fourpi * zrhoir (ig) *&
        !       (1d0 - exp(-zlambda*r_c)* ( zlambda * sin(gpc(ig)*r_c) / gpc(ig) + cos(gpc(ig)*r_c) ) )/&
        !       ((gpc(ig) ** 2) + (zlambda ** 2))
        ! Else
        !    zvclir (ifg) = fourpi * zrhoir(ig)* (1d0 - exp(-zlambda*r_c)*(zlambda*r_c + 1))/(zlambda ** 2)
        ! End If

      else ! cutoff
        Do ig = 1, ngvec
          zvclir (ig) = fourpi * zrhoir (ig) / ((gpc(ig) ** 2) + (zlambda ** 2))
        enddo
      endif !if cutoff

    endif! if yukawa or not
  end subroutine poisson_ir2



  ! call surface_ir2( input%groundstate%lmaxvr, ngp, jlgpr, ylmgp, sfacgp, zvclir, qlmir)

  subroutine surface_ir2( lmax, ngvec,  jlgpr, ylmgp, sfacgp, vclig, vilm2)
    ! use mod_atoms, only: nspecies, natoms, idxas
    ! use mod_muffin_tin, only: rmt
    use mod_atoms, only: natmtot, nspecies,natoms, idxas
    use constants, only: fourpi,zil

    integer, intent(in) :: lmax
    integer, intent(in) :: ngvec
    real(dp), intent(in) :: jlgpr(0:,:,:)
    complex(dp), intent(in) :: ylmgp(:,:)
    complex(dp), intent(in) :: sfacgp(:,:)
    complex(dp), intent(in) :: vclig(:)

    complex(dp), intent(out) :: vilm2(:,:)

    integer :: ig, l,m, lm, is, ia, ias
    Complex (8) :: zlm ((lmax+1)**2),zlm2((lmax+1)**2)
    Complex (8) :: sf(natmtot)

    vilm2  = 0.d0
! !$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,l,lm,is,ia,ias,zlm,zlm2,sf) SHARED(vclig,ngvec,ylmgp,jlgpr,sfacgp,nspecies,natoms,idxas,lmax) REDUCTION(+:vilm2)
! !$OMP DO
    Do ig = 1, ngvec
      zlm(:) = vclig (ig) *conjg (ylmgp(:, ig))
      sf(:)=sfacgp(ig,:)
      Do is = 1, nspecies
        lm = 0
        Do l = 0, lmax
          zlm2(lm+1:lm+1+2*l) = jlgpr (l, ig, is) * zlm(lm+1:lm+1+2*l)
          lm=lm+1+2*l
        End Do

        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          vilm2 (:,ias) = vilm2 (:,ias) + sf (ias) * zlm2(:)
        End Do
      End Do
    End Do
! !$OMP END DO
! !$OMP END PARALLEL
    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        lm = 0
        Do l = 0, lmax
          Do m = - l, l
            lm = lm + 1
            vilm2 (lm,ias) = vilm2 (lm,ias) * fourpi * zil (l)
          End Do !m
        End Do!l
      enddo!ia
    enddo!ia

  end subroutine surface_ir2

  subroutine surface_ir3( lmax, ngvec,  jlgpr, ylmgp, sfacgp, vclig, vilm2)
    ! use mod_atoms, only: nspecies, natoms, idxas
    ! use mod_muffin_tin, only: rmt
    use mod_atoms, only: natmtot, nspecies,natoms, idxas
    use constants, only: fourpi,zil

    integer, intent(in) :: lmax
    integer, intent(in) :: ngvec
    real(dp), intent(in) :: jlgpr(0:,:,:)
    complex(dp), intent(in) :: ylmgp(:,:)
    complex(dp), intent(in) :: sfacgp(:,:)
    complex(dp), intent(in) :: vclig(:)

    complex(dp), intent(out) :: vilm2(:,:)

    integer :: ig, l,m, lm, is, ia, ias, lmmax, atst, sfld
    integer, parameter :: blocksize=32
    Complex (8) :: zlm ((lmax+1)**2,blocksize),zlm2((lmax+1)**2,blocksize)
    integer :: iggg, blk

!          vilm2  = 0.d0
    sfld=size(sfacgp,dim=1)
    lmmax=(lmax+1)**2
    vilm2  = 0.d0
    Do ig = 0, ngvec-1, blocksize
      blk=min(blocksize,ngvec-ig)
      do iggg=1,blk
        zlm(:,iggg) = vclig (ig+iggg) *conjg (ylmgp(:, ig+iggg))
      enddo

      Do is = 1, nspecies
        do iggg=1,blk
          lm = 0
          Do l = 0, lmax
            zlm2(lm+1:lm+1+2*l,iggg) = jlgpr (l, ig+iggg, is) * zlm(lm+1:lm+1+2*l,iggg)
            lm=lm+1+2*l
          End Do
        enddo

        atst = idxas(1, is)
        Call zgemm('N','N', lmmax, natoms(is), blk, (1.0D0,0.0), zlm2, lmmax, sfacgp(ig+1,atst), sfld, (1.0D0,0.0), vilm2(1,atst), lmmax)

      End Do
    End Do

    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        lm = 0
        Do l = 0, lmax
          Do m = - l, l
            lm = lm + 1
            vilm2 (lm,ias) = vilm2 (lm,ias) * fourpi * zil (l)
          End Do !m
        End Do!l
      enddo!ia
    enddo!ia

  end subroutine surface_ir3


end module weinert
