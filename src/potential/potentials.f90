module potentials
  use precision, only: dp

  implicit none
  private

  public :: coulomb_potential, coulomb_potential2, coulomb_potential2vFeld

contains

  !> This subroutines calcultes the Coulomb potential from a given complex charge density
  !> \(n({\bf r})\) using Weinert's method.
  !>
  !> The Coulomb potential is given by the sum of the Hartree potential and the external potential
  !> \[ V_{\rm cl}({\bf r}) = \int \frac{n({\bf r}')}{|{\bf r} - {\bf r}'|} \, {\rm d}^3 r
  !> + \sum_{\bf R} \sum_\alpha \frac{Z_\alpha}{|{\bf r} - {\bf \tau}_\alpha - {\bf R}|} \;, \]
  !> where \(Z_\alpha\) is the ionic charge of the atom \(\alpha\) at position \({\bf \tau}_\alpha\).
  !>
  !> The Coulomb potential is calculated as follows.
  !> First, Poisson's equation inside each muffin-tin
  !> sphere is solved using the Green's function method and the multipole moments of the muffin-tin
  !> charge density are obtained using the subroutine [[poisson_and_multipoles_mt(subroutine)]].
  !> Then, we add the ionic potential and monopoles to the muffin-tin potential and multipoles.
  !> In a third step, we compute the multipole moments corresponding to the extension of the interstitial
  !> charge density inside the muffin-tin spheres using the subroutine [[multipoles_ir(subroutine)]].
  !> Next, we solve Poisson's equation for the interstitial region using the subroutine [[poisson_ir(subroutine)]].
  !> This is done by constructing a pseudocharge density with the same muffin-tin multipole moments
  !> as the actual charge density and having a quickly converging Fourier series. With this pseudocharge
  !> density, we solve Poisson's equation in reciprocal space for the interstitial Coulomb potential.
  !> In the last step, we add the homogeneous solution of Poisson's equation to the muffin-tin Coulomb potential
  !> to match it with the interstitial potential on the muffin-tin sphere boundaries using the subroutines
  !> [[match_bound_mt(subroutine)]] and [[surface_ir(subroutine)]].


  subroutine coulomb_potential2( nr, r, ngp, gpc, igp0, jlgpr, ylmgp, sfacgp, zn, zrhomt, zrhoir, zvclmt, zvclir, zrho0,qvec, cutoff,&
    & hybrid_in, yukawa_in,zlambda_in,zbessi,zbessk,zilmt,rpseudo_in,rpseudomat)
    use modsurf, only: surf_pot
    use constants, only: y00,zzero
    use modinput
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_muffin_tin, only: lmmaxvr, rmt, nrmtmax
    use mod_Gvector, only: ngrtot, ngrid, ivg, intgv, igfft, ivgig
    use mod_potential_and_density, only: vmad
    use mod_convergence, only: iscl
    use weinert
!> number or radial grid points for each species
    integer, intent(in) :: nr(:)
!> radial grid for each species
    real(dp), intent(in) :: r(:,:)
!> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
!> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
!> index of shortest \({\bf G+p}\) vector
    integer, intent(in) :: igp0
!> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
!> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
!> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
!> ionic charges \(Z_\alpha\) for each species
    real(dp), intent(in) :: zn(:)
!> complex muffin-tin charge density
    complex(dp), intent(in) :: zrhomt(:,:,:)
!> complex interstitial charge density
    complex(dp), intent(In) :: zrhoir(:)
!> complex muffin-tin Coulomb potential
    complex(dp), intent(out) :: zvclmt(:,:,:) ! lm, nrmax, natoms
!> complex interstitial Coulomb potential
    complex(dp), intent(out) :: zvclir(:)
!> Fourier component of pseudocharge density for shortest \({\bf G+p}\) vector
    complex(dp), intent(out) :: zrho0
    real(8), intent(in) :: qvec(3)
!> option for using coulomb cutoff for solving Poisson's equation
    logical, optional, intent(in) :: cutoff
    logical, optional, intent(in) :: hybrid_in
    logical, optional, intent(in) :: yukawa_in
    complex(dp), optional, intent(in) :: zlambda_in
!      complex(dp), optional, intent(in) :: zbessi(nrmtmax,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies)
!      complex(dp), optional, intent(in) :: zbessk(nrmtmax,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies)
    complex(dp), optional, intent(in) :: zbessi(:,0:,:)!(nrmtmax, 0:input%groundstate%lmaxvr+input%groundstate%npsden+1, nspecies)
    complex(dp), optional, intent(in) :: zbessk(:,0:,:)
    complex(dp), optional, intent(in) :: zilmt(0:,:)
    real(dp), allocatable :: vion(:,:), vdplmt(:,:,:), vdplir(:)
    complex(dp), allocatable :: qlm(:,:), qlmir(:,:), zrhoig(:)
    logical, optional, intent(in) :: rpseudo_in
    complex(dp), optional, intent(in) :: rpseudomat(:,:,:)

    logical :: yukawa
    logical :: hybrid, rpseudo
    complex(dp) :: zlambda

    integer :: is, ia, ias, ir,j
    integer :: ig, ifg, lm

    real(dp) :: t1, ta, tb,tc,td,te,tf

    complex(dp), allocatable :: zrhoig_sort(:)
    real(dp) :: maxrhog
    integer :: ngp2

    if (present(hybrid_in)) then
      hybrid=hybrid_in
    else
      hybrid=.false.
    endif

    if (present(yukawa_in)) then
      yukawa=yukawa_in
    else
      yukawa=.false.
    endif

    if (present(zlambda_in)) then
      zlambda=zlambda_in
    endif

    if (present(rpseudo_in)) then
      rpseudo=rpseudo_in
    else
      rpseudo=.false.
    endif


    allocate( qlm( lmmaxvr, natmtot), qlmir( lmmaxvr, natmtot))
    allocate( zrhoig( ngrtot))
    allocate( vion( nrmtmax, nspecies), source=0._dp)

! solve Poisson's equation in each muffin-tin sphere
! and find muffin-tin multipoles


    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        if(yukawa) then
          call poisson_and_multipoles_mt_yukawa ( input%groundstate%lmaxvr, nr(is), r(:,is), zrhomt(:,:,ias), zvclmt(:,:,ias), qlm(:,ias), is,&
            & yukawa_in=yukawa,zlambda=zlambda, il=zbessi(1:nr(is),0:input%groundstate%lmaxvr,is), kl=zbessk(1:nr(is),0:input%groundstate%lmaxvr,is))
        else ! Coulumb
          call poisson_and_multipoles_mt_yukawa ( input%groundstate%lmaxvr, nr(is), r(:,is), zrhomt(:,:,ias), zvclmt(:,:,ias), qlm(:,ias), is)
        endif

        vmad(ias) = dble( zvclmt(1,1,ias) ) * y00
      end do
    end do

! add ionic potential and monopole
    do is = 1, nspecies
      if( zn(is) == 0.0_dp) cycle
      call potnucl( input%groundstate%ptnucl, nr(is), r(:,is), zn(is), vion(:,is))
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        zvclmt(1,1:nr(is),ias) = zvclmt(1,1:nr(is),ias) + vion(1:nr(is),is)/y00
        qlm(1,ias) = qlm(1,ias) + zn(is)*y00
      end do
    end do




! Fourier transform interstitial density to reciprocal space

    zrhoig = zrhoir
    call zfftifc( 3, ngrid, -1, zrhoig)
    allocate(zrhoig_sort(ngp))
    maxrhog=0d0
    do ig=1,ngp
      zrhoig_sort(ig)=zrhoig(igfft(ig))
      maxrhog=max(maxrhog,abs(zrhoig_sort(ig)))
    enddo

    ig=ngp
    do while (abs(zrhoig_sort(ig)).lt.1d-14*maxrhog)
      ig=ig-1
    enddo
    if (ig.eq.ngp)then
      ngp2=ig
    else
      ngp2=ig + 1
    endif
    call timesec(ta)
    if (yukawa) then
      call multipoles_ir3( input%groundstate%lmaxvr, ngp2, gpc, &
        & jlgpr, ylmgp, sfacgp, igfft, &
        zrhoig_sort, qlmir,zlambda,zilmt)
    else
      ! multipoles_ir4 - some data in single precision
      call multipoles_ir3( input%groundstate%lmaxvr, ngp2, gpc, &
        & jlgpr, ylmgp, sfacgp, igfft, &
        zrhoig_sort, qlmir)
    endif
    write(*,*) "qlmir max abs is: ", maxval(abs(qlmir))
    call timesec(tb)
!write(*,*)"multipoles:",tb-ta

! take difference of muffin-tin and interstitial multipole moments
    qlm = qlm - qlmir



    ! !solve Poisson's equation in interstitial region
    call timesec(ta)
    if(rpseudo)then

      do ig=1, ngp
        ifg = igfft (ig)
        zrhoig(ifg)=zrhoig_sort(ig)
      enddo
      call pseudocharge_rspace_new(input%groundstate%lmaxvr,qlm,rpseudomat,zrhoig)
      do ig=1, ngp
        ifg = igfft (ig)
        zrhoig_sort(ig)=zrhoig(ifg)
      enddo
      call timesec(td)
      !write(*,*)"kartošana2x :",td-tc+tb-ta
      !write(*,*)"konstrueshana :",tc-tb
    else
      if (yukawa) then
        call pseudocharge_gspace2(input%groundstate%lmaxvr, ngp, gpc, &
          & jlgpr, ylmgp, sfacgp, igfft, zrhoig_sort, qlm,zlambda,zilmt)
      else
        call pseudocharge_gspace3(input%groundstate%lmaxvr, ngp, gpc, &
          & jlgpr, ylmgp, sfacgp, igfft, zrhoig_sort, qlm)
      endif
    endif
    call timesec(tb)
!write(*,*)"pseudocharge:",tb-ta

    call timesec(ta)
    if (yukawa) then
      call poisson_ir2(ngp, gpc, igfft, zrhoig_sort, zvclir, cutoff, zlambda=zlambda)
    else
      call poisson_ir2(ngp, gpc, igfft, zrhoig_sort, zvclir, cutoff)
    endif
    call timesec(tb)
!write(*,*)"poisson_ir:", tb-ta
    ! zrhoig_sort(1:ngp)=zvclir(1:ngp)
    ! zvclir(:)=zzero
    ! do ig=1, ngp
    !  ifg = igfft (ig)
    !  zvclir(ifg)=zrhoig_sort(ig)
    ! enddo


    call timesec(ta)

    if (input%groundstate%hybrid%rsurf) then
      zrhoig_sort(1:ngp)=zvclir(1:ngp)
      zvclir(:)=zzero
      do ig=1, ngp
        ifg = igfft (ig)
        zvclir(ifg)=zrhoig_sort(ig)
      enddo
      call zfftifc( 3, ngrid, 1, zvclir)

      call surf_pot(input%groundstate%lmaxvr,zvclir,igfft,qvec,qlmir)
    else
      call surface_ir3( input%groundstate%lmaxvr, ngp, jlgpr, ylmgp, sfacgp, zvclir, qlmir)
! Fourier transform interstitial potential to real space
      zrhoig_sort(1:ngp)=zvclir(1:ngp)
      zvclir(:)=zzero
      do ig=1, ngp
        ifg = igfft (ig)
        zvclir(ifg)=zrhoig_sort(ig)
      enddo
      call zfftifc( 3, ngrid, 1, zvclir)
    endif


    call timesec(tb)
!   write(*,*)"surface_ir:",tb-ta
! stop

    ! open(1, file = 'vlm_original.dat', status = 'replace')
    ! do ias=1, natmtot
    ! do lm=1,lmmaxvr
    !   write(1,*)dble(qlmir(lm,ias)),",",aimag(qlmir(lm,ias)),",", ias, lm
    ! enddo
    ! enddo
    ! close(1)
    ! stop


    call timesec(ta)
! match muffin-tin and interstitial potential on muffin-tin boundaries
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        if (present(zbessi)) then
          call match_bound_mt( input%groundstate%lmaxvr, nr(is), r(:,is), rmt(is), qlmir(:,ias), zvclmt(:,:,ias),yukawa=yukawa, zbessi=zbessi(:,:,is) )
        else
          call match_bound_mt( input%groundstate%lmaxvr, nr(is), r(:,is), rmt(is), qlmir(:,ias), zvclmt(:,:,ias),yukawa=yukawa)
        endif
        vmad(ias) = vmad(ias) + dble( qlmir(1,ias) ) * y00 - vion(nr(is),is)
      end do
    end do
    call timesec(tb)

!write(*,*)"match_bound_mt",tb-ta



    deallocate( vion, zrhoig, qlm, qlmir,zrhoig_sort)

! add dipole correction
    if( (iscl > 0) .and. input%groundstate%dipolecorrection) then
      allocate( vdplmt( lmmaxvr, maxval(nr), natmtot), vdplir( ngrtot))
      call dipole_correction( vdplmt, vdplir)
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nr(is)
            zvclmt(1,ir,ias) = zvclmt(1,ir,ias)+vdplmt(1,ir,ias)
            zvclmt(3,ir,ias) = zvclmt(3,ir,ias)+vdplmt(3,ir,ias)
          end do
        end do
      end do
      zvclir = zvclir + vdplir
      deallocate( vdplmt, vdplir)
    end if

    if(.false.)then
      open(11,file='mt_rez.dat',status='replace')
      is=1
      lm=1
      do ir=1, nr(is)
        write(11,*)r(ir,is),",",dble(zvclmt (lm, ir, is))*y00
      enddo
      close(11)

      open(11,file='is_rez.dat',status='replace')
      Do ig = 1, ngrtot
        t1 = gpc (ig)
        write(11,*)dble(zvclir(ig))
      End Do
      close(11)
!stop
    endif



  end subroutine coulomb_potential2



  subroutine coulomb_potential2vFeld(qlmir_in,qlmir_out, nr, r, ngp, gpc, igp0, jlgpr, ylmgp, sfacgp, zn, zrhomt, zrhoir, zvclmt, zvclir, zrho0,qvec, cutoff,&
    & hybrid_in, yukawa_in,zlambda_in,zbessi,zbessk,zilmt,rpseudo_in,rpseudomat)
    use modsurf, only: surf_pot
    use constants, only: y00,zzero
    use modinput
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_muffin_tin, only: lmmaxvr, rmt, nrmtmax
    use mod_Gvector, only: ngrtot, ngrid, ivg, intgv, igfft, ivgig
    use mod_potential_and_density, only: vmad
    use mod_convergence, only: iscl
    use weinert
!> number or radial grid points for each species
    integer, intent(in) :: nr(:)
!> radial grid for each species
    real(dp), intent(in) :: r(:,:)
!> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
!> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
!> index of shortest \({\bf G+p}\) vector
    integer, intent(in) :: igp0
!> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
!> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
!> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
!> ionic charges \(Z_\alpha\) for each species
    real(dp), intent(in) :: zn(:)
!> complex muffin-tin charge density
    complex(dp), intent(in) :: zrhomt(:,:,:)
!> complex interstitial charge density
    complex(dp), intent(In) :: zrhoir(:)
!> complex muffin-tin Coulomb potential
    complex(dp), intent(out) :: zvclmt(:,:,:) ! lm, nrmax, natoms
!> complex interstitial Coulomb potential
    complex(dp), intent(out) :: zvclir(:)
!> Fourier component of pseudocharge density for shortest \({\bf G+p}\) vector
    complex(dp), intent(out) :: zrho0
    real(8), intent(in) :: qvec(3)

    !complex(dp), allocatable, intent(in) :: qlmir_in(:,:)
    complex(dp), allocatable, intent(in) :: qlmir_in(:,:)
    complex(dp), allocatable, intent(out) :: qlmir_out(:,:)
!> option for using coulomb cutoff for solving Poisson's equation
    logical, optional, intent(in) :: cutoff
    logical, optional, intent(in) :: hybrid_in
    logical, optional, intent(in) :: yukawa_in
    complex(dp), optional, intent(in) :: zlambda_in
!      complex(dp), optional, intent(in) :: zbessi(nrmtmax,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies)
!      complex(dp), optional, intent(in) :: zbessk(nrmtmax,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies)
    complex(dp), optional, intent(in) :: zbessi(:,0:,:)!(nrmtmax, 0:input%groundstate%lmaxvr+input%groundstate%npsden+1, nspecies)
    complex(dp), optional, intent(in) :: zbessk(:,0:,:)
    complex(dp), optional, intent(in) :: zilmt(0:,:)
    real(dp), allocatable :: vion(:,:), vdplmt(:,:,:), vdplir(:)
    complex(dp), allocatable :: qlm(:,:),qlmir(:,:), zrhoig(:)
    logical, optional, intent(in) :: rpseudo_in
    complex(dp), optional, intent(in) :: rpseudomat(:,:,:)

    logical :: yukawa
    logical :: hybrid, rpseudo
    complex(dp) :: zlambda

    integer :: is, ia, ias, ir,j
    integer :: ig, ifg, lm

    real(dp) :: t1, ta, tb,tc,td,te,tf

    complex(dp), allocatable :: zrhoig_sort(:)
    real(dp) :: maxrhog
    integer :: ngp2

    if (present(hybrid_in)) then
      hybrid=hybrid_in
    else
      hybrid=.false.
    endif

    if (present(yukawa_in)) then
      yukawa=yukawa_in
    else
      yukawa=.false.
    endif

    if (present(zlambda_in)) then
      zlambda=zlambda_in
    endif

    if (present(rpseudo_in)) then
      rpseudo=rpseudo_in
    else
      rpseudo=.false.
    endif


    allocate( qlm( lmmaxvr, natmtot), qlmir( lmmaxvr, natmtot))
    allocate( zrhoig( ngrtot))
    allocate( vion( nrmtmax, nspecies), source=0._dp)

! solve Poisson's equation in each muffin-tin sphere
! and find muffin-tin multipoles


    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        if(yukawa) then
          call poisson_and_multipoles_mt_yukawa ( input%groundstate%lmaxvr, nr(is), r(:,is), zrhomt(:,:,ias), zvclmt(:,:,ias), qlm(:,ias), is,&
            & yukawa_in=yukawa,zlambda=zlambda, il=zbessi(1:nr(is),0:input%groundstate%lmaxvr,is), kl=zbessk(1:nr(is),0:input%groundstate%lmaxvr,is))
        else ! Coulumb
          call poisson_and_multipoles_mt_yukawa ( input%groundstate%lmaxvr, nr(is), r(:,is), zrhomt(:,:,ias), zvclmt(:,:,ias), qlm(:,ias), is)
        endif

        vmad(ias) = dble( zvclmt(1,1,ias) ) * y00
      end do
    end do

! add ionic potential and monopole
    do is = 1, nspecies
      if( zn(is) == 0.0_dp) cycle
      call potnucl( input%groundstate%ptnucl, nr(is), r(:,is), zn(is), vion(:,is))
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        zvclmt(1,1:nr(is),ias) = zvclmt(1,1:nr(is),ias) + vion(1:nr(is),is)/y00
        qlm(1,ias) = qlm(1,ias) + zn(is)*y00
      end do
    end do




! Fourier transform interstitial density to reciprocal space

    zrhoig = zrhoir
    call zfftifc( 3, ngrid, -1, zrhoig)
    allocate(zrhoig_sort(ngp))
    maxrhog=0d0
    do ig=1,ngp
      zrhoig_sort(ig)=zrhoig(igfft(ig))
      maxrhog=max(maxrhog,abs(zrhoig_sort(ig)))
    enddo

    ig=ngp
    do while (abs(zrhoig_sort(ig)).lt.1d-14*maxrhog)
      ig=ig-1
    enddo
    if (ig.eq.ngp)then
      ngp2=ig
    else
      ngp2=ig + 1
    endif
    call timesec(ta)
#define temp_removal_1237123812038120
#ifdef temp_removal_1237123812038120
    if (yukawa) then
      call multipoles_ir3( input%groundstate%lmaxvr, ngp2, gpc, &
        & jlgpr, ylmgp, sfacgp, igfft, &
        zrhoig_sort, qlmir,zlambda,zilmt)
    else
      ! multipoles_ir4 - some data in single precision
      call multipoles_ir3( input%groundstate%lmaxvr, ngp2, gpc, &
        & jlgpr, ylmgp, sfacgp, igfft, &
        zrhoig_sort, qlmir)
    endif

#endif temp_removal_1237123812038120
    call timesec(tb)
!write(*,*)"multipoles:",tb-ta

! take difference of muffin-tin and interstitial multipole moments
    qlmir_out = qlmir
    !write(*,*) "Difference in multipole calcs: " maxval(abs(qlmir-qlmir_in))
    !if (.not. allocated(qlmir_in)) then
    ! write(*,*) "Error: qlmir_in not allocated"
    !endif
    !write(*,*) "Relative difference: ", maxval(abs(qlmir - qlmir_in))/max(maxval(abs(qlmir)),0.00000001), maxval(abs(qlmir_in)), maxval(abs(qlmir))
    qlmir = qlmir_in
    qlm = qlm - qlmir



    ! !solve Poisson's equation in interstitial region
    call timesec(ta)
    if(rpseudo)then

      do ig=1, ngp
        ifg = igfft (ig)
        zrhoig(ifg)=zrhoig_sort(ig)
      enddo
      call pseudocharge_rspace_new(input%groundstate%lmaxvr,qlm,rpseudomat,zrhoig)
      do ig=1, ngp
        ifg = igfft (ig)
        zrhoig_sort(ig)=zrhoig(ifg)
      enddo
      call timesec(td)
      !write(*,*)"kartošana2x :",td-tc+tb-ta
      !write(*,*)"konstrueshana :",tc-tb
    else
      if (yukawa) then
        call pseudocharge_gspace2(input%groundstate%lmaxvr, ngp, gpc, &
          & jlgpr, ylmgp, sfacgp, igfft, zrhoig_sort, qlm,zlambda,zilmt)
      else
        call pseudocharge_gspace3(input%groundstate%lmaxvr, ngp, gpc, &
          & jlgpr, ylmgp, sfacgp, igfft, zrhoig_sort, qlm)
      endif
    endif
    call timesec(tb)
!write(*,*)"pseudocharge:",tb-ta

    call timesec(ta)
    if (yukawa) then
      call poisson_ir2(ngp, gpc, igfft, zrhoig_sort, zvclir, cutoff, zlambda=zlambda)
    else
      call poisson_ir2(ngp, gpc, igfft, zrhoig_sort, zvclir, cutoff)
    endif
    call timesec(tb)
!write(*,*)"poisson_ir:", tb-ta
    ! zrhoig_sort(1:ngp)=zvclir(1:ngp)
    ! zvclir(:)=zzero
    ! do ig=1, ngp
    !  ifg = igfft (ig)
    !  zvclir(ifg)=zrhoig_sort(ig)
    ! enddo


    call timesec(ta)

    if (input%groundstate%hybrid%rsurf) then
      zrhoig_sort(1:ngp)=zvclir(1:ngp)
      zvclir(:)=zzero
      do ig=1, ngp
        ifg = igfft (ig)
        zvclir(ifg)=zrhoig_sort(ig)
      enddo
      call zfftifc( 3, ngrid, 1, zvclir)

      call surf_pot(input%groundstate%lmaxvr,zvclir,igfft,qvec,qlmir)
    else
      call surface_ir3( input%groundstate%lmaxvr, ngp, jlgpr, ylmgp, sfacgp, zvclir, qlmir)
! Fourier transform interstitial potential to real space
      zrhoig_sort(1:ngp)=zvclir(1:ngp)
      zvclir(:)=zzero
      do ig=1, ngp
        ifg = igfft (ig)
        zvclir(ifg)=zrhoig_sort(ig)
      enddo
      call zfftifc( 3, ngrid, 1, zvclir)
    endif


    call timesec(tb)
!   write(*,*)"surface_ir:",tb-ta
! stop

    ! open(1, file = 'vlm_original.dat', status = 'replace')
    ! do ias=1, natmtot
    ! do lm=1,lmmaxvr
    !   write(1,*)dble(qlmir(lm,ias)),",",aimag(qlmir(lm,ias)),",", ias, lm
    ! enddo
    ! enddo
    ! close(1)
    ! stop


    call timesec(ta)
! match muffin-tin and interstitial potential on muffin-tin boundaries
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        if (present(zbessi)) then
          call match_bound_mt( input%groundstate%lmaxvr, nr(is), r(:,is), rmt(is), qlmir(:,ias), zvclmt(:,:,ias),yukawa=yukawa, zbessi=zbessi(:,:,is) )
        else
          call match_bound_mt( input%groundstate%lmaxvr, nr(is), r(:,is), rmt(is), qlmir(:,ias), zvclmt(:,:,ias),yukawa=yukawa)
        endif
        vmad(ias) = vmad(ias) + dble( qlmir(1,ias) ) * y00 - vion(nr(is),is)
      end do
    end do
    call timesec(tb)

!write(*,*)"match_bound_mt",tb-ta

    deallocate( vion, zrhoig, qlm, qlmir,zrhoig_sort)

! add dipole correction
    if( (iscl > 0) .and. input%groundstate%dipolecorrection) then
      allocate( vdplmt( lmmaxvr, maxval(nr), natmtot), vdplir( ngrtot))
      call dipole_correction( vdplmt, vdplir)
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nr(is)
            zvclmt(1,ir,ias) = zvclmt(1,ir,ias)+vdplmt(1,ir,ias)
            zvclmt(3,ir,ias) = zvclmt(3,ir,ias)+vdplmt(3,ir,ias)
          end do
        end do
      end do
      zvclir = zvclir + vdplir
      deallocate( vdplmt, vdplir)
    end if

    if(.false.)then
      open(11,file='mt_rez.dat',status='replace')
      is=1
      lm=1
      do ir=1, nr(is)
        write(11,*)r(ir,is),",",dble(zvclmt (lm, ir, is))*y00
      enddo
      close(11)

      open(11,file='is_rez.dat',status='replace')
      Do ig = 1, ngrtot
        t1 = gpc (ig)
        write(11,*)dble(zvclir(ig))
      End Do
      close(11)
!stop
    endif



  end subroutine coulomb_potential2vFeld










  subroutine coulomb_potential( nr, r, ngp, gpc, igp0, jlgpr, ylmgp, sfacgp, zn, zrhomt, zrhoir, zvclmt, zvclir, zrho0, cutoff,&
    & hybrid_in, yukawa_in,zlambda_in,zbessi,zbessk,zilmt,rpseudo_in,rpseudomat)

    use constants, only: y00,zzero
    use modinput
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_muffin_tin, only: lmmaxvr, rmt, nrmtmax
    use mod_Gvector, only: ngrtot, ngrid, ivg, intgv, igfft, ivgig
    use mod_potential_and_density, only: vmad
    use mod_convergence, only: iscl
    use weinert
    !> number or radial grid points for each species
    integer, intent(in) :: nr(:)
    !> radial grid for each species
    real(dp), intent(in) :: r(:,:)
    !> total number of \({\bf G+p}\) vectors
    integer, intent(in) :: ngp
    !> lengths of \({\bf G+p}\) vectors
    real(dp), intent(in) :: gpc(:)
    !> index of shortest \({\bf G+p}\) vector
    integer, intent(in) :: igp0
    !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
    real(dp), intent(in) :: jlgpr(0:,:,:)
    !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
    complex(dp), intent(in) :: ylmgp(:,:)
    !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
    complex(dp), intent(in) :: sfacgp(:,:)
    !> ionic charges \(Z_\alpha\) for each species
    real(dp), intent(in) :: zn(:)
    !> complex muffin-tin charge density
    complex(dp), intent(in) :: zrhomt(:,:,:)
    !> complex interstitial charge density
    complex(dp), intent(In) :: zrhoir(:)
    !> complex muffin-tin Coulomb potential
    complex(dp), intent(out) :: zvclmt(:,:,:) ! lm, nrmax, natoms
    !> complex interstitial Coulomb potential
    complex(dp), intent(out) :: zvclir(:)
    !> Fourier component of pseudocharge density for shortest \({\bf G+p}\) vector
    complex(dp), intent(out) :: zrho0
    !> option for using coulomb cutoff for solving Poisson's equation
    logical, optional, intent(in) :: cutoff
    logical, optional, intent(in) :: hybrid_in
    logical, optional, intent(in) :: yukawa_in
    complex(dp), optional, intent(in) :: zlambda_in
!      complex(dp), optional, intent(in) :: zbessi(nrmtmax,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies)
!      complex(dp), optional, intent(in) :: zbessk(nrmtmax,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies)
    complex(dp), optional, intent(in) :: zbessi(:,0:,:)!(nrmtmax, 0:input%groundstate%lmaxvr+input%groundstate%npsden+1, nspecies)
    complex(dp), optional, intent(in) :: zbessk(:,0:,:)
    complex(dp), optional, intent(in) :: zilmt(0:,:)
    real(dp), allocatable :: vion(:,:), vdplmt(:,:,:), vdplir(:)
    complex(dp), allocatable :: qlm(:,:), qlmir(:,:), zrhoig(:)
    logical, optional, intent(in) :: rpseudo_in
    complex(dp), optional, intent(in) :: rpseudomat(:,:,:)

    logical :: yukawa
    logical :: hybrid, rpseudo
    complex(dp) :: zlambda

    integer :: is, ia, ias, ir,j
    integer :: ig, ifg, lm

    real(dp) :: t1


    if (present(hybrid_in)) then
      hybrid=hybrid_in
    else
      hybrid=.false.
    endif

    if (present(yukawa_in)) then
      yukawa=yukawa_in
    else
      yukawa=.false.
    endif

    if (present(zlambda_in)) then
      zlambda=zlambda_in
    endif

    if (present(rpseudo_in)) then
      rpseudo=rpseudo_in
    else
      rpseudo=.false.
    endif


    allocate( qlm( lmmaxvr, natmtot), qlmir( lmmaxvr, natmtot))
    allocate( zrhoig( ngrtot))
    allocate( vion( nrmtmax, nspecies), source=0._dp)

    ! solve Poisson's equation in each muffin-tin sphere
    ! and find muffin-tin multipoles

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)

        if(yukawa) then
          call poisson_and_multipoles_mt_yukawa ( input%groundstate%lmaxvr, nr(is), r(:,is), zrhomt(:,:,ias), zvclmt(:,:,ias), qlm(:,ias), is,&
            & yukawa_in=yukawa,zlambda=zlambda, il=zbessi(1:nr(is),0:input%groundstate%lmaxvr,is), kl=zbessk(1:nr(is),0:input%groundstate%lmaxvr,is))
          vmad(ias) = dble( zvclmt(1,1,ias) ) * y00
        else ! Coulumb
          call poisson_and_multipoles_mt ( input%groundstate%lmaxvr, nr(is), r(:,is), zrhomt(:,:,ias), zvclmt(:,:,ias), qlm(:,ias))
          vmad(ias) = dble( zvclmt(1,1,ias) ) * y00
        endif
      end do
    end do

    ! add ionic potential and monopole
    do is = 1, nspecies
      if( zn(is) == 0.0_dp) cycle
      call potnucl( input%groundstate%ptnucl, nr(is), r(:,is), zn(is), vion(:,is))
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        zvclmt(1,1:nr(is),ias) = zvclmt(1,1:nr(is),ias) + vion(1:nr(is),is)/y00
        qlm(1,ias) = qlm(1,ias) + zn(is)*y00
      end do
    end do

    ! Fourier transform interstitial density to reciprocal space
    zrhoig = zrhoir

    call zfftifc( 3, ngrid, -1, zrhoig)
    ! find multipole moments of interstitial density
    ! extended into the muffin-tin spheres





    if (yukawa) then

      !Do ig = ngp+1, ngrtot
      !  ifg = igfft (ig)
      !  zrhoig(ifg)=zzero
      !End Do

      call multipoles_ir_yukawa( input%groundstate%lmaxvr, ngp, gpc, &
        & jlgpr, ylmgp, sfacgp, igfft, &
        zrhoig, qlmir,zlambda,zilmt)
!write(*,*)"qlmir"

      !stop
    else



      call multipoles_ir( input%groundstate%lmaxvr, ngp, gpc, &
        ivg, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, &
        zrhoig, qlmir)
    endif
    ! take difference of muffin-tin and interstitial multipole moments


    qlm = qlm - qlmir



    ! solve Poisson's equation in interstitial region

    if (yukawa) then
      if (rpseudo) then
        call pseudocharge_rspace_new(input%groundstate%lmaxvr,qlm,rpseudomat,zrhoig)

      else
        call pseudocharge_gspace_yukawa(input%groundstate%lmaxvr, ngp, gpc, &
          & jlgpr, ylmgp, sfacgp, igfft, zrhoig, qlm,zlambda,zilmt)
      endif
      ! open(11,file='is_pseudo_r.dat',status='replace')
      ! Do ig = 1, ngrtot
      !     t1 = gpc (ig)
      !     ifg = igfft (ig)
      !     write(11,*)t1,",",dble(zrhoig(ifg))
      ! End Do
      ! close(11)
      !stop



      !zrho0 = zrhoig( igfft( igp0))

      call poisson_ir_yukawa(input%groundstate%lmaxvr, ngp, gpc, igfft, zrhoig, zlambda,zvclir,cutoff)


    else
      if (rpseudo) then
        call poisson_ir( input%groundstate%lmaxvr, input%groundstate%npsden, ngp, gpc, &
          ivg, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, &
          zrhoig, qlm, zvclir, cutoff,hybrid_in=hybrid,rpseudo_in=rpseudo,rpseudomat=rpseudomat)
      else
        call poisson_ir( input%groundstate%lmaxvr, input%groundstate%npsden, ngp, gpc, &
          ivg, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, &
          zrhoig, qlm, zvclir, cutoff,hybrid_in=hybrid,rpseudo_in=.false.)
      endif
      zrho0 = zrhoig( igfft( igp0))
    endif
    ! evaluate interstitial potential on muffin-tin surface
    call surface_ir( input%groundstate%lmaxvr, ngp, gpc, &
      ivg, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, &
      zvclir, qlmir)

    ! match muffin-tin and interstitial potential on muffin-tin boundaries
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        if (present(zbessi)) then
          call match_bound_mt( input%groundstate%lmaxvr, nr(is), r(:,is), rmt(is), qlmir(:,ias), zvclmt(:,:,ias),yukawa=yukawa, zbessi=zbessi(:,:,is) )
        else
          call match_bound_mt( input%groundstate%lmaxvr, nr(is), r(:,is), rmt(is), qlmir(:,ias), zvclmt(:,:,ias),yukawa=yukawa )
        endif
        vmad(ias) = vmad(ias) + dble( qlmir(1,ias) ) * y00 - vion(nr(is),is)
      end do
    end do
!if(yukawa)then
!do ig=1, 10
!  ifg = igfft (ig)
!write(*,*)"gatavs ,",zvclir(ifg)
!enddo
!stop
!endif






    ! Fourier transform interstitial potential to real space
    call zfftifc( 3, ngrid, 1, zvclir)

    deallocate( vion, zrhoig, qlm, qlmir)

    ! add dipole correction
    if( (iscl > 0) .and. input%groundstate%dipolecorrection) then
      allocate( vdplmt( lmmaxvr, maxval(nr), natmtot), vdplir( ngrtot))
      call dipole_correction( vdplmt, vdplir)
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nr(is)
            zvclmt(1,ir,ias) = zvclmt(1,ir,ias)+vdplmt(1,ir,ias)
            zvclmt(3,ir,ias) = zvclmt(3,ir,ias)+vdplmt(3,ir,ias)
          end do
        end do
      end do
      zvclir = zvclir + vdplir
      deallocate( vdplmt, vdplir)
    end if

    if(.false.)then
      open(11,file='mt_rez.dat',status='replace')
      is=1
      lm=1
      do ir=1, nr(is)
        write(11,*)r(ir,is),",",dble(zvclmt (lm, ir, is))*y00
      enddo
      close(11)

      open(11,file='is_rez.dat',status='replace')
      Do ig = 1, ngrtot
        t1 = gpc (ig)
        write(11,*)dble(zvclir(ig))
      End Do
      close(11)
      !stop
    endif



  end subroutine coulomb_potential

end module potentials
