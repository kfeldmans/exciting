module mod_Feldmans
  use modmain
  use precision, only: dp, i32
  use constants, only: pi
  use mod_sbessel, only: sbessel, sbessel_vec
  use mod_eigensystem, only: WFType
  use mod_kpointset, only: Gk_set
  use modgw, only: Gkqset, Gkset
  use mod_Basis1D
  use mod_lattice, only: omega
  implicit none
  private
  public :: hello, testing, testing_calc, set_wf_planewaves

  ! Ripped from Fortran stdlib
  real(dp), parameter :: tolerance = 4._dp*epsilon(1._dp)
  integer, parameter :: newton_iters = 100
contains

  ! Ripped from fortran standard library

  ! legegendre polynomials
  ! unspecified behaviour if n is negative
  pure elemental function legendre(n, x) result(leg)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: leg
    select case (n)
     case (0)
      leg = 1
     case (1)
      leg = x
     case default
      block
        real(dp) :: leg_down1, leg_down2
        integer :: i

        leg_down1 = x
        leg_down2 = 1

        do i = 2, n
          leg = (2*i - 1)*x*leg_down1/i - (i - 1)*leg_down2/i
          leg_down2 = leg_down1
          leg_down1 = leg
        end do
      end block
    end select
  end function legendre

  ! derivatives of legegendre polynomials
  ! unspecified behaviour if n is negative
  pure elemental function dlegendre(n, x) result(dleg)
    !! Calculates the derivatives of Legendre polynomials
    integer, intent(in) :: n !! Number of Legendre polynomials to find
    real(dp), intent(in) :: x !! Coordinate at which to calculate the polynomials
    real(dp) :: dleg

    select case (n)
     case (0)
      dleg = 0
     case (1)
      dleg = 1
     case default
      block
        real(dp) :: leg_down1, leg_down2, leg
        real(dp) :: dleg_down1, dleg_down2
        integer :: i

        leg_down1 = x
        dleg_down1 = 1

        leg_down2 = 1
        dleg_down2 = 0

        do i = 2, n
          leg = (2*i - 1)*x*leg_down1/i - (i - 1)*leg_down2/i
          dleg = dleg_down2 + (2*i - 1)*leg_down1
          leg_down2 = leg_down1
          leg_down1 = leg
          dleg_down2 = dleg_down1
          dleg_down1 = dleg
        end do
      end block
    end select
  end function dlegendre

  subroutine testing_non_zero_gaunt_count()
    integer :: l,l1,l2,m,m1,m2,count
    integer, parameter :: L_max = 10

    do l=0,L_max
      do m = -l,l
        count = 0
        do l1 = 0,L_max
          do l2=abs(l1-l),min(l+l1,L_max),2
            count = count + min(l1-m,l2)-max(-l1-m,-l2)
          end do !l2
        end do !l1
        print '(I0, " ", I0," ", I0)',l,m, count
      end do ! m
    end do !l

  end subroutine testing_non_zero_gaunt_count

  subroutine gauss_legendre(x, w, interval)
    ! Ripped from fortran standard library!
    !! Calculates grid coordinates and weights for Gauss-Legendre quadrature
    real(dp), intent(out) :: x(:) !! Quadrature coordinates
    real(dp), intent(out) :: w(:) !! Quadrature weights
    real(dp), intent(in), optional :: interval(2) !! Interval end-points, default: [0,1]

    associate (n => size(x) - 1)
      select case (n)
       case (0)
        x = 0
        w = 2
       case (1)
        x(1) = -sqrt(1._dp/3._dp)
        x(2) = -x(1)
        w = 1
       case default
        block
          integer :: i, j
          real(dp) :: leg, dleg, delta

          do i = 0, (n + 1)/2 - 1
            ! use Gauss-Chebyshev points as an initial guess
            x(i + 1) = -cos((2*i + 1)/(2._dp*n + 2._dp)*pi)
            do j = 1, newton_iters
              leg = legendre(n + 1, x(i + 1))
              dleg = dlegendre(n + 1, x(i + 1))
              delta = -leg/dleg
              x(i + 1) = x(i + 1) + delta
              if (abs(delta) <= tolerance*abs(x(i + 1))) exit
            end do
            x(n - i + 1) = -x(i + 1)

            dleg = dlegendre(n + 1, x(i + 1))
            w(i + 1) = 2._dp/((1 - x(i + 1)**2)*dleg**2)
            w(n - i + 1) = w(i + 1)
          end do

          if (mod(n, 2) == 0) then
            x(n/2 + 1) = 0

            dleg = dlegendre(n + 1, 0.0_dp)
            w(n/2 + 1) = 2._dp/(dleg**2)
          end if
        end block
      end select
    end associate

    if (present(interval)) then
      associate (a => interval(1), b => interval(2))
        x = 0.5_dp*(b - a)*x + 0.5_dp*(b + a)
        w = 0.5_dp*(b - a)*w
      end associate
    end if
  end subroutine gauss_legendre ! END Ripped from fortran standard library

  pure function linspace(min, max, N) result(grid)
    real(dp), intent(in) :: min ! Smallest value of grid
    real(dp), intent(in) :: max ! Largest value of grid
    integer, intent(in) :: N ! Number of grid points
    real(dp) :: grid(1:N) ! N eqaully spaced grid points from min to max
    integer(i32) :: i
    real(dp) :: d
    d = (max - min)/(N - 1)
    grid = [(min + d*(i - 1), i=1, N)]
  end function linspace

  subroutine hello()
    print *, "Got here2!"
  end subroutine hello

  subroutine set_wf_planewaves(wf)
    type(WFType), intent(inout) :: wf

    !wf%gk(:,:) = 0._dp
    !wf%gk(1,1) = cmplx(1._dp, 1._dp)
    !wf%gk(3,1) = 3._dp
  end subroutine set_wf_planewaves

  subroutine testing_calc(wf,q_nnLM_is)
    Type(WFType), intent(in) :: wf
    complex(dp), allocatable, intent(out) :: q_nnLM_is(:,:,:, :) ! q_nn_LM_is(n,n',LM, is) in q_nnLM_is(nstfv,nstfv,lmmaxvr, natmtot)
    real(dp), parameter :: minx = 0._dp
    real(dp) :: maxx = 1._dp
    real(dp), parameter :: minG = 0.00000001_dp
    real(dp), parameter :: maxG = 10._dp
    integer(i32), parameter :: Nx = 6000
    integer(i32), parameter :: NG = 500
    integer(i32), parameter :: n = 15 !! How many eigenvals & eigenfuncs to keep for representation
    integer(i32), parameter :: lmax = 8 !! Maximum order for spherical Bessel functions
    integer(i32), parameter :: L_max = 8 !! Maximum angular momentum to calculate integrals for
    Type(Basis1D) :: SVDb_vec(0:lmax) !! Vector of bases for each order of sBessel functions.
    real(dp) :: U_l(Nx, n, 0:lmax) !! sBessel SVD approximation bases for all orders
    real(dp) :: U_l_scaled(Nx, n, 0:lmax) !! Helper in integral calculation
    real(dp) :: A(Nx,NG)
    real(dp) :: xgrid(Nx)
    real(dp) :: weightsx(Nx)
    real(dp) :: Ggrid(NG)
    real(dp) :: jlxG(0:lmax, Nx, NG)
    real(dp), allocatable :: work(:) ! Work array for blas/lapack
    real(dp) :: koefs(n)
    real(dp) :: f(Nx)
    complex(dp), allocatable :: ylmgk(:,:) ! ylmgk(lmmaxvr,Gkqset%ngk(1,1))
    integer :: lwork
    real(dp) :: xpow(Nx)
    integer :: i,j,l,ll,l1,l2,idxg,igk,m,m1,m2,lm,ibs,ibs1,ibs2,ist,ist1,ist2,is,ias,ia, lm1,lm2,lm3
    real(dp), allocatable :: svdfac(:,:,:) ! Expansion coefficients of spherical Bessel functions in  SVD basis for angular moments, for used lengths of G vectors.
    real(dp), allocatable :: jlxg_rayl(:,:,:) !
    real(dp), allocatable :: jxg(:,:)
    complex(dp), allocatable :: B(:,:,:)
    complex(dp), allocatable :: A_is(:,:,:,:) ! A(Gkqset%ngk(1,1), lmmaxvr, nstfv,natomtot)
    complex(dp), allocatable :: B_is(:,:,:,:) ! B(n,lmmaxvr,nstfv, natomtot)
    complex(dp), allocatable :: psi_nlm_approx(:,:,:)
    complex(dp), allocatable :: psi_nlm_is_approx(:,:,:,:)
    complex(dp), allocatable :: psi_nlm(:,:,:)
    complex(dp), allocatable :: psi_nlm_is(:,:,:,:)
    real(dp) :: test_func1(Nx),test_func2(Nx)
    real(dp) :: total
    real(dp) :: omega_inv ! inverse of unit cell volume


    real(dp) :: I_Llsls(n, n, 0:lmax, 0:lmax, 0:L_max) !! Integral values between all basis elements
    complex(dp) :: psi_nlm_local(Nx)
    write(*,*) "Started testing_calc"

    omega_inv = sqrt(1/omega)

    call test_wavefunction(wf)


    write(*,*) "Testing calc started"
    write(*,*) "Found max value of radius: ", maxval(rmt(1:nspecies))
    write(*,*) "Length of only non-zero plane-wave: ", Gkqset%gkc(1,1,1)
    !maxx = maxval(rmt(1:nspecies))

    ! Generating the gridX and gridG.
    !xgrid = linspace(minx, maxx, Nx)
    call gauss_legendre(xgrid, weightsx, interval=[minx, maxx])
    Ggrid = linspace(minG, maxG, NG)
    write(*,*) "Sum of integration weights: ", sum(weightsx)

    ! Calculating the sBessel functions on the grids for all orders for generating an SVD basis.
    do idxg = 1,NG
      call sbessel_vec(lmax,maxval(rmt(1:nspecies))*xgrid * Ggrid(idxg), jlxG(:,:,idxg))
    end do

    ! Generate a basis for each l.
    do l =0, lmax
      A = jlxG(l,:,:)
      lwork = -1
      call genSVDBasis(Nx,A,weightsx,n,xgrid,work,lwork,SVDb_vec(l))
    end do

    call testing_non_zero_gaunt_count

    call gengntyyy
    write(*,*) "The Gaunt coefficient array is allocated: ", allocated(gntyyy)

    lm1 = idxlm(3,2)
    lm2 = idxlm(2,-1)
    lm3 = idxlm(3,3)
    write(*,*) "The precomputed G: ", gntyyy(lm3,lm2,lm1)
    ! The last parameter in gntyyy is for the complex conjugated sperical harmonic.
    ! We must have the sum of first two m's equal to last m.
    !write(*,*) "Gaunt coeffs: ", gntyyy(idxlm(1,1),idxlm(1,1),idxlm(2,2))
    !write(*,*) "Stop in mod_Feldmans.f90 after checking Gaunt coeffs"
    !stop


#define T413094738947394723
#ifdef T413094738947394723
    ! Checking the basis
    write(*,*) "The sum of integration weights of basis: ", sum(SVDb_vec(2)%weights)
    test_func1 = SVDb_vec(0)%basis(:,1)
    test_func2 = SVDb_vec(0)%basis(:,1)
    total = sum(test_func1*test_func2*SVDb_vec(0)%weights)
    write(*,*) "The norm of a basis function under the weights: ", total
    write(*,*) "The inner product of two basis functions under the weights: ", sum(SVDb_vec(3)%basis(:,1)*SVDb_vec(3)%basis(:,6)*SVDb_vec(3)%weights)
#endif T413094738947394723
    ! Checking the correctness of the bases and projection functions.
#define SVDProjCheck
#ifdef SVDProjCheck
    A = jlxG(3,:,:)
    koefs = xToB(SVDb_vec(3), A(:,400))
    f = BTox(SVDb_vec(3), koefs)
    write(*,*) "Total weighted abs error accross the sBessel approximation: ", sum(SVDb_vec(2)%weights*abs(A(:,400) - f))
    write(*,*) "Maximum approximation error: ", maxval(abs(A(:,400) - f))
    do i=1,1000
      koefs = xToB(SVDb_vec(2), f)
      f = BTox(SVDb_vec(2), koefs)
    end do
    write(*,*) "Maximum approximation error after 1000 repeated projections: ", maxval(abs(A(:,400) - f))

    write(*,*) "Maximum value of A(:,400): ", maxval(A(:,400))
#endif SVDProjCheck

    ! Calculating the spherical harmonics at G_ang angles, stored in ylmgk.
    allocate(ylmgk(lmmaxvr,Gkqset%ngk(1,1)))
    Do igk = 1, Gkqset%ngk(1,1)
      Call genylm(input%groundstate%lmaxvr, Gkqset%tpgkc(:, igk,1,1), ylmgk(:, igk))
    End Do

    write(*,*) "lmmaxvr: ", lmmaxvr
    write(*,*) "input%groundstate%lmaxvr: ", input%groundstate%lmaxvr

    ! Calculate the sBessel functions used in the Rayleigh expansion
    allocate(jlxg_rayl(0:lmax,Nx,Gkqset%ngk(1,1)))
    do igk=1,Gkqset%ngk(1,1)
      call sbessel_vec(lmax,xgrid * Gkqset%gkc(igk,1,1)*maxval(rmt(1:nspecies)),jlxg_rayl(:,:,igk))
    end do

    ! Calculate SVDb coefficients for each sBessel function used in Rayleigh expansion
    write(*,*) "Started finding psi_nlm through B coefficients"
    allocate(svdfac(n,Gkqset%ngk(1,1),0:lmax))
    allocate(jxg(Nx,Gkqset%ngk(1,1)))
    do l=0,lmax
      jxg = jlxg_rayl(l,:,:)
      do igk=1,Gkqset%ngk(1,1)
        svdfac(:,igk,l) = sqrt(maxval(rmt))*xToB(SVDb_vec(l), jxg(:,igk))
      end do
    end do

    ! Calculating the B coefficients
    write(*,*) "ngk(1,1): ", ngk(1,1)
    write(*,*) "Gkqset%ngk(1,1): ", Gkqset%ngk(1,1)
    allocate(B(n,lmmaxvr,nstfv))
    B = 0._dp
    ! Gkqset%sfacgk(ngk(:,:),natmmax,1,ik): ispn == 1 - index spinor,
    do igk = 1, ngk(1,1)
      do ist = 1, nstfv ! index state
        do l = 0, input%groundstate%lmaxvr
          do m= -l,l
            lm = idxlm(l,m)
            do ibs=1,n
              B(ibs,lm,ist)= B(ibs,lm,ist) + omega_inv*wf%gk(igk,ist)*fourpi* cmplx(0,1, kind=dp) ** l * (Gkqset%sfacgk(igk,1,1,1)) * (ylmgk(lm,igk))*svdfac(ibs,igk,l)
            end do !s
          end do ! m
        end do ! l
      end do ! ist
    end do ! igk


    ! Calculating the B coefficients for all atoms.
    write(*,*) "Started calclating B coefficients for all atoms"
    allocate(B_is(n,lmmaxvr,nstfv,natmtot))
    B_is = cmplx(0._dp, kind=dp)
    do is=1,nspecies
      do ia = 1,natoms(is)
        ias = idxas(ia,is)
        write(*,*) "Atom index being computed on: ", ias

        do ist = 1,nstfv
          do l = 0, input%groundstate%lmaxvr
            do m=-l,l
              lm=idxlm(l,m)
              do ibs=1,n
                do igk = 1, ngk(1,1)
                  B_is(ibs,lm,ist,ias) = B_is(ibs,lm,ist, ias) + omega_inv*wf%gk(igk,ist)*fourpi*Gkqset%sfacgk(igk,ias,1,1)* (cmplx(0,1, kind=dp) ** l)  * conjg(ylmgk(lm,igk))*svdfac(ibs,igk,l)
                end do ! igk
              end do ! ibs
            end do ! m
          end do ! l
        end do ! ist
      end do ! ia
    end do ! is
    write(*,*) "Finished calclating B coefficients for all atoms"
    write(*,*) "The maximum of all B's: ", maxval(abs(B_is))


!#define TEST_B_COEFFICIENTS_498231
#ifdef TEST_B_COEFFICIENTS_498231
    ! Finding the spherical harmonic coefficients of the wavefunction from B coefficients.
    allocate(psi_nlm_approx(Nx,lmmaxvr,nstfv))
    psi_nlm_approx = cmplx(0._dp, kind=dp)
    do ist = 1,nstfv
      do l=0,lmax
        do m=-l,l
          lm=idxlm(l,m)
          do is=1,n
            psi_nlm_approx(:,lm,ist) =  psi_nlm_approx(:,lm,ist) + B(is,lm,ist)*SVDb_vec(l)%basis(:,is)/sqrt(maxval(rmt))
          end do
        end do
      end do
    end do
    write(*,*) "Finished finding psi_nlm through B coefficients"
    ! Finding the spherical harmonic coefficients of the wavefunction directly
    write(*,*) "Started finding psi_nlm directly"
    allocate(psi_nlm(Nx,lmmaxvr,nstfv))
    psi_nlm = cmplx(0._dp, kind=dp)
    do l=0,lmax
      jxg = jlxg_rayl(l,:,:)
      do ist=1,nstfv
        do m=-l,l
          lm=idxlm(l,m)
!$OMP PARALLEL DO  PRIVATE(psi_nlm_local)
          do igk = 1, ngk(1,1)
            psi_nlm_local = omega_inv*wf%gk(igk,ist)*fourpi* cmplx(0,1, kind=dp) ** l  * conjg(ylmgk(lm,igk)) * jxg(:,igk) / Gkqset%sfacgk(igk,1,1,1)
!$OMP critical
            psi_nlm(:,lm,ist) = psi_nlm(:,lm,ist) + psi_nlm_local
!$OMP end critical
          end do
!$OMP end parallel do
        end do
      end do
    end do
    write(*,*) "Finished finding psi_nlm directly"
    write(*,*) "The difference between both approximations: ", maxval(abs(psi_nlm-psi_nlm_approx))
    write(*,*) "Maximum value of psi_nlm: ", maxval(abs(psi_nlm))
    write(*,*) "Maximum value of psi_nlm_approx: ", maxval(abs(psi_nlm_approx))
    write(*,*) "Maximum unweighted summed error of the approximation: ", maxval(sum(abs(psi_nlm_approx-psi_nlm),dim=1) )
#endif TEST_B_COEFFICIENTS_498231



!#define TEST_B_IAS_COEFFICIENTS_4313418741
#ifdef TEST_B_IAS_COEFFICIENTS_4313418741
    ! Finding the spherical harmonic coefficients of the wavefunction from B_is coefficients.
    allocate(psi_nlm_is_approx(Nx,lmmaxvr,nstfv,natmtot))
    psi_nlm_approx = cmplx(0._dp, kind=dp)
    do is = 1, nspecies
      do ia = 1,natoms(is)
        ias = idxas(ia,is)
        do ist = 1,nstfv
          do l=0,lmax
            do m=-l,l
              lm=idxlm(l,m)
              do ibs=1,n
                psi_nlm_is_approx(:,lm,ist,ias) =  psi_nlm_is_approx(:,lm,ist,ias) + B_is(ibs,lm,ist,ias)*SVDb_vec(l)%basis(:,ibs)/sqrt(maxval(rmt))
              end do
            end do
          end do
        end do
      end do
    end do

    write(*,*) "Finished finding psi_nlm_is through B coefficients"
    ! Finding the spherical harmonic coefficients of the wavefunction directly
    write(*,*) "Started finding psi_nlm_is directly"
    allocate(psi_nlm_is(Nx,lmmaxvr,nstfv,natmtot))
    psi_nlm = cmplx(0._dp, kind=dp)
    do is = 1, nspecies
      do ia = 1,natoms(is)
        ias = idxas(ia,is)
        do l=0,lmax
          jxg = jlxg_rayl(l,:,:)
          do ist=1,nstfv
            do m=-l,l
              lm=idxlm(l,m)
!$OMP PARALLEL DO  PRIVATE(psi_nlm_local)
              do igk = 1, ngk(1,1)
                psi_nlm_local = omega_inv*wf%gk(igk,ist)*fourpi* cmplx(0,1, kind=dp) ** l * conjg(ylmgk(lm,igk)) * jxg(:,igk) / Gkqset%sfacgk(igk,ias,1,1)
!$OMP critical
                psi_nlm_is(:,lm,ist,ias) = psi_nlm_is(:,lm,ist, ias) + psi_nlm_local
!$OMP end critical
              end do
!$OMP end parallel do
            end do
          end do
        end do
      end do
    end do
    write(*,*) "Finished finding psi_nlm_is directly"
    write(*,*) "The difference between both approximations: ", maxval(abs(psi_nlm_is-psi_nlm_is_approx))
    write(*,*) "Maximum value of psi_nlm_is: ", maxval(abs(psi_nlm_is))
    write(*,*) "Maximum value of psi_nlm_approx: ", maxval(abs(psi_nlm_is_approx))
    write(*,*) "Maximum unweighted summed error of the approximation: ", maxval(sum(abs(psi_nlm_is_approx-psi_nlm_is),dim=1) )
#endif TEST_B_IAS_COEFFICIENTS_4313418741







#define TEST_CALC_MULTIPOLES_VIA_B_KOEFS
#ifdef TEST_CALC_MULTIPOLES_VIA_B_KOEFS

    do l = 0, lmax
      U_l(:,:,l) = SVDb_vec(l)%basis
    end do

    do ll = 0, L_max
      !--- Scale the basis vectors for calculating integrals ---!
      xpow = sqrt(((SVDb_vec(1)%xgrid*maxval(rmt(1:nspecies)))**(ll + 2))*SVDb_vec(1)%weights)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
      do i = 1, n
        do j = 0, lmax
          U_l_scaled(:, i, j) = U_l(:, i, j)*xpow
        end do
      end do
      !$OMP END PARALLEL DO

      !--- Calculate the integrals per l,l' pair ---!
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l1,l2)
      do l1 = 0, lmax
        do l2 = 0, lmax
          call dgemm('T', 'N', n, n, Nx, 1._dp, U_l_scaled(:, :, l1), Nx, U_l_scaled(:, :, l2), Nx, 0._dp, I_Llsls(:, :, l1, l2, ll), n)
        end do
      end do
      !$OMP END PARALLEL DO
    end do
    write(*,*) "Maximum abs of I_Llsls: ", maxval(abs(I_Llsls))
    !write(*,*) "Pritnting integrals: ", I_Llsls(:,:,0,0,2)

    write(*,*) "Maximum abs of B_is: ", maxval(abs(B_is))

    write(*,*) "Started calculation of multipoles"
    allocate(q_nnLM_is(nstfv,nstfv,lmmaxvr, natmtot)) ! q_nn_LM_is(n,n',LM, ias)
    q_nnLM_is = cmplx(0._dp, kind=dp)
    !#define old_loop_1412321
#ifdef old_loop_1412321
    !!$OMP PARALLEL DO COLLAPSE(12) PRIVATE(is, ia, ias,ist1,ist2,l,m,lm,l1,m1,lm1,l2,m2,lm2,ibs1,ibs2)
    do is =1, nspecies
      do ia = 1,natoms(is)
        ias = idxas(ia,is)
        do l=0,L_max
          do m=-l,l
            lm=idxlm(l,m)
            do l1=0,lmax
              do m1=-l1,l1
                lm1 = idxlm(l1,m1)
                do l2=0,lmax
                  do m2=-l2,l2
                    lm2 = idxlm(l2,m2)
                    !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(ist1,ist2,ibs1,ibs2)
                    do ist1=1,nstfv
                      do ist2=1,nstfv
                        do ibs1 = 1,n
                          do ibs2 = 1,n
                            q_nnLM_is(ist1,ist2,lm,ias) = q_nnLM_is(ist1,ist2,lm,ias) + gntyyy(lm,lm2,lm1) * B_is(ibs1,lm1,ist1,ias) * conjg(B_is(ibs2,lm2,ist2,ias))*I_Llsls(ibs1,ibs2,l1,l2,l)
                          end do !ibs2
                        end do !ibs1
                      end do !m2
                    end do !l2
                    !$OMP END PARALLEL DO
                  end do  !m1
                end do !l1
              end do !m
            end do !l
          end do !ist2
        end do !ist1
      end do !ia
    end do ! is
    !!$OMP END PARALLEL DOI_Llsls(:
#endif old_loop_1412321

#define new_loop_1312312
#ifdef new_loop_1312312
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ia, ias, l, m, lm, l1, m1, lm1, l2, m2, lm2, ist1, ist2, ibs1, ibs2)
    do is = 1, nspecies
      do ia = 1, natoms(is)
        do l = 0, L_max
          do m = -l, l
            ias = idxas(ia, is)
            lm = idxlm(l, m)
            do l1 = 0, lmax
              do l2 = abs(l-l1),min(lmax,l+l1), 2
                do m1 = max(-l1,m-l2), min(l1,m+l2)
                  lm1 = idxlm(l1, m1)
                  lm2 = idxlm(l2, m1-m)
                  !if (abs(gntyyy(lm1, lm, lm2)) < 1e-14_dp) cycle  ! skip zeros
                  do ist1 = 1, nstfv
                    do ist2 = 1, nstfv
                      do ibs1 = 1, n
                        do ibs2 = 1, n
                          q_nnLM_is(ist1, ist2, lm, ias) = q_nnLM_is(ist1, ist2, lm, ias) &
                            + cmplx(gntyyy(lm, lm2, lm1)) &
                            * B_is(ibs1, lm1, ist1, ias) &
                            * conjg(B_is(ibs2, lm2, ist2, ias)) &
                            * cmplx(I_Llsls(ibs1, ibs2, l1, l2, l))
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO
#endif new_loop_1312312

!#define loop_optimization_1_43141983471283
#ifdef loop_optimization_1_43141983471283

    do is= 1,nspecies
      do ia = 1,natoms(is)
        ias = idxas(ia,is)
        do m = 0,lmax ! we handle positive and negative case seperately
          ! case where m >= 0
          do m2 = -lmax,lmax-m
            m1 = m+m2
            ! now we need to loop over all triples of (l,l1,l2) that have non-zero gaunt at (lm,lm1,lm2)

          end do !m2
        end do !m loop over all magnetic moments
      end do !ia
    end do !is
#endif loop_optimization_1_43141983471283

    write(*,*) "Maximum abs of multipoles: ", maxval(abs(q_nnLM_is))

    do ist1 = 1,nstfv
      do ist2=1,nstfv
        !print '("NEW: For pair ", I0, ", ", I0, " obtained value ", F0.3)', ist1, ist2, maxval(abs(q_nnLM_is(ist1,ist2,:,:)))
        if (ist1 == 1 .and. ist2 == 1) then
          do is = 1, nspecies
            do ia = 1, natoms(is)

              do l = 0,L_max
                do m = -l,l

                  ias = idxas(ia, is)
                  lm = idxlm(l,m)
                  !For l ", I0, ", and m ", I0, " on atom ", I0 ," obtained values ", F0.10,F0.10)', l, m,ias, qlmir(lm,ias)%re,(qlmir(lm,ias)%Im)
                  print '("NEW: For l ", I0, ", and m ", I0, " on atom ", I0 ," obtained value ", F0.10, F0.10)', l, m,ias, (q_nnLM_is(ist1,ist2,lm,ias)%re), (q_nnLM_is(ist1,ist2,lm,ias)%im)
                  !print '(F0.10, " ", F0.10)',dble(q_nnLM_is(ist1,ist2,lm,ias)), dimag(q_nnLM_is(ist1,ist2,lm,ias))
                end do
              end do
            end do !m
          end do !l
        endif
      end do
    end do
    write(*,*) "Finished calculation of multipoles"
#endif TEST_CALC_MULTIPOLES_VIA_B_KOEFS

    write(*,*) "Testing calc finished"

  end subroutine testing_calc










  subroutine testing(wf)
    type(WFType), intent(in) :: wf !! Wave function to be processsed.
    real(dp), parameter :: minx = 0._dp
    real(dp), parameter :: maxx = 1._dp
    real(dp), parameter :: minG = 0.00001_dp
    real(dp), parameter :: maxG = 4._dp
    integer(i32), parameter :: Nx = 4000
    integer(i32), parameter :: NG = 1000
    integer(i32), parameter :: n = 12 !! How many eigenvals & eigenfuncs to keep for representation
    integer(i32), parameter :: lmax = 8 !! Maximum order for spherical Bessel functions
    integer(i32), parameter :: L_max = 5 !! Maximum angular momentum to calculate integrals for
    integer(i32), parameter :: file_unit = 1
    integer :: lm
    real(dp) :: gridx(Nx)
    real(dp) :: weightsx(Nx)
    real(dp) :: gridG(NG)
    real(dp) :: jlxG(0:lmax, Nx, NG)
    real(dp), allocatable :: jlxG_temp(:,:,:), jxG_temp(:,:)
    real(dp) :: jxG_l(Nx, NG) !! Helper object for SVD of sBessel f-ons of fixed order.
    real(dp) :: U(Nx, Nx), S(min(Nx, NG)), VT(NG, NG) !! Result matrices for SVD of sBessels of fixed order.
    real(dp), allocatable :: WORK(:) !! Helper array for algorithms
    integer(i32) :: LWORK !! Size of helper array WORK
    integer(i32) :: INFO !! Diagnostic information about success of subroutine
    real(dp) :: jxG_l_approx(Nx, NG) !! approximately obtained result from SVD
    real(dp) :: jxG_l_approx_tmp(Nx, NG) !! helper for approximately obtained result from SVD
    real(dp) :: I_Llsls(n, n, 0:lmax, 0:lmax, 0:L_max) !! Integral values between all basis elements
    real(dp) :: U_l(Nx, n, 0:lmax) !! sBessel SVD approximation bases for all orders
    real(dp) :: U_l_scaled(Nx, n, 0:lmax) !! Helper in integral calculation
    real(dp) :: xpow(Nx) !! Stores values of form x**p for some p
    real(dp) :: A(Nx, NG) !! Helper matrix for SVD
    real(dp) :: jl1(Nx), jl2(Nx), jl1c(n),  jl_tmp(0:lmax)!,jl2c(n), jl0(Nx), jl0c(n) ! Store example sBessel functions and their expansion coefficients
    real(dp), allocatable :: svdfac(:,:,:) ! Expansion coefficients of spherical Bessel functions in  SVD basis for angular moments, for used lengths of G vectors.
    integer(i32) :: i, j, ll ! loop counter, ll-large l
    integer :: igk, ist,l,m,is ! Loop counters
    complex(dp), allocatable :: B(:,:,:)
    complex(dp), allocatable :: ylmgk(:,:)
    real(dp), allocatable :: check_tmp(:,:)
    complex(dp), allocatable :: Psi_lmn_approx(:,:,:)
    !complex(dp), allocatable ::Psi_lmn(:,:,:)
    integer :: lm1, lm2
    complex(dp), allocatable :: zU_l(:,:,:)
    real(dp) :: omega_inv ! inverse of unit cell volume
    Type(Basis1D) :: SVDb

    omega_inv = sqrt(1/omega)

    gridx = linspace(minx, maxx, Nx)
    call gauss_legendre(gridx, weightsx, interval=[minx, maxx])
    gridG = linspace(minG, maxG, NG)

    !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(i)
    do i = 1, NG
      call sbessel_vec(lmax, gridX*gridG(i), jlxG(:, :, i))
    end do
    !$OMP END PARALLEL DO
    write(*,*) input%groundstate%lmaxvr
    !--- For checking computation results ---!
    ! print *, "G is equal to ", gridG(10)
    ! open(unit=file_unit, file='feld_output.txt', status='replace')
    ! do i=1,Nx
    !   write(file_unit, '(2ES24.16)') gridx(i), jlxG(0,i,10)
    ! end do
    ! close(file_unit)

    !--- Find the singular value decomposition for the spherical Bessel functions
    jxG_l = jlxG(lmax, :, :) ! Gather the sBessel values for fixed order

    ! Determine space needed for SVDreal(dp), allocatable :: svdfac(:,:,:) ! Expansion coefficients of spherical Bessel functions in  SVD basis for angular moments, for used lengths of G vectors.
    allocate (WORK(1))
    LWORK = -1
    call dgesvd('A', 'A', Nx, NG, jxG_l, Nx, S, U, Nx, VT, NG, WORK, LWORK, INFO)
    print *, "INFO from WORK size query returned ", INFO
    LWORK = int(WORK(1))
    deallocate (WORK)
    ! Perform SVD
    allocate (WORK(LWORK))
    call dgesvd('A', 'A', Nx, NG, jxG_l, Nx, S, U, Nx, VT, NG, WORK, LWORK, INFO)

    print *, "INFO from SVD routine returned ", INFO
    jxG_l = jlxG(lmax, :, :) ! Restoring jxG_l, since it is destroyed during SVD

    ! Testing the SVDBasis routine
    lwork = -1
    A = jxG_l
    call SVDBasis(A, work, lwork)

    print *, "SVDBasis performed computation of U with maximum error of: ", maxval(abs(A - U(:, 1:NG)))

    ! Testing Basis1D
    jxG_l = jlxG(lmax, :, :)
    lwork = -1
    A = jxG_l
    call genSVDBasis(Nx,A,weightsx,n,gridx,work, lwork,SVDb)
    write(*,*) "Smalles integration weight is ", minval(weightsx)
    print *, "genSVDBasis performed computation of U with maximum error of: ", maxval(abs(A - U(:, 1:NG)))
    ! Check SVD results
    jxG_l_approx = 0
    do i = 1, size(S)
      jxG_l_approx(i, i) = S(i)
    end do
    call dgemm('N', 'N', Nx, NG, Nx, 1._dp, U, Nx, jxG_l_approx, Nx, 0._dp, jxG_l_approx_tmp, Nx)
    call dgemm('N', 'N', Nx, NG, NG, 1._dp, jxG_l_approx_tmp, Nx, VT, NG, 0._dp, jxG_l_approx, Nx)
    print *, "Size of approximation: ", size(jxG_l_approx)
    print *, "Size of original: ", size(jxG_l)
    print *, "Maximum error in SVD approximation: ", maxval(abs(jxG_l - jxG_l_approx))

    ! Check truncated SVD results
    ! Should assert that n>=min(Nx,NG)
    jxG_l_approx = 0
    do i = 1, n
      jxG_l_approx(i, i) = S(i)
    end do
    call dgemm('N', 'N', Nx, n, n, 1._dp, U, Nx, jxG_l_approx, Nx, 0._dp, jxG_l_approx_tmp, Nx)
    call dgemm('N', 'N', Nx, NG, n, 1._dp, jxG_l_approx_tmp, Nx, VT, NG, 0._dp, jxG_l_approx, Nx)
    print *, "Size of approximation: ", size(jxG_l_approx)
    print *, "Size of original: ", size(jxG_l)
    print *, "Maximum error in truncated SVD approximation: ", maxval(abs(jxG_l - jxG_l_approx))

    !---- Calculating integral values I_L,ls,l's' ----!
    !--- Calculate basis functions for all l
    do l = 0, lmax
      A = jlxG(l, :, :)
      call SVDBasis(A, work, lwork)
      U_l(1:Nx, 1:n, l) = A(1:Nx, 1:n)
    end do

#define alternative_basis
#ifdef alternative_basis
    !----- Calculate weighted basis functions for all l --- !
    do l = 0, lmax
      A = jlxG(l, :, :)
      call genSVDBasis(Nx,A,weightsx,n,gridx,work,lwork,SVDb)
      U_l(1:Nx, 1:n, l) = SVDb%basis(1:Nx,1:n)
    end do
#endif alternative_basis

    do ll = 0, L_max
      !--- Scale the basis vectors for calculating integrals ---!
      xpow = sqrt(gridx**(ll + 2)*weightsx)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
      do i = 1, n
        do j = 0, lmax
          U_l_scaled(:, i, j) = U_l(:, i, j)*xpow
        end do
      end do
      !$OMP END PARALLEL DO

      !--- Calculate the integrals per l,l' pair ---!
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
      do i = 0, lmax
        do j = 0, lmax
          call dgemm('T', 'N', n, n, Nx, 1._dp, U_l_scaled(:, :, i), Nx, U_l_scaled(:, :, j), Nx, 0._dp, I_Llsls(:, :, i, j, ll), n)
        end do
      end do
      !$OMP END PARALLEL DO
    end do

    !--- Check calculated integrals with an example ---!

    !$OMP PARALLEL DO PRIVATE(i, jl_tmp)
    do i = 1, Nx
      call sbessel(lmax, gridx(i)*3.023_dp, jl_tmp(:))
      jl1(i) = jl_tmp(lmax)
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i, jl_tmp)
    do i = 1, Nx
      call sbessel(lmax, gridx(i)*2.423_dp, jl_tmp(:))
      jl2(i) = jl_tmp(lmax - 1)
    end do
    !$OMP END PARALLEL DO

    call dgemv('T', Nx, n, 1._dp, U_l(:, :, lmax), Nx, jl1, 1, 0._dp, jl1c, 1)

    !--- Check projection result
    !call dgemv('N', Nx, n, 1._dp, U_l(:, :, lmax), Nx, jl1c, 1, 0._dp, jl0(:), 1)

    !call dgemv('T', Nx, n, 1._dp, U_l(:, :, lmax - 1), Nx, jl2, 1, 0._dp, jl2c, 1)

    ! Vektori jl1c un jl2c satur izvirzijuma bazes koeficientus.
    ! Varam sarekinat, piemeram, int[x**(L+2)*jl1(x)*jl2(x) dx]
    !call dgemv('N', n, n, 1._dp, I_Llsls(:, :, lmax, lmax - 1, 5), n, jl2c, 1, 0._dp, jl0c, 1)

    !print *, "Integresanas tuvinatais rezultats caur keofiem: ", sum(jl0c*jl1c)

    !--- End of projection result

    !--- Now we need to get the coefficients needed for the next steps.
    !-- For the next step we need to get the radial function structure coefficients


    !---Get the projection coefficients of Bessel functions
    ! Pagaidām sarēķinām tikai pirmajam k vektoram
    ! svdfac(n,ngk)

    allocate(svdfac(n,Gkqset%ngk(1,1),0:lmax))
    allocate(jlxg_temp(0:lmax,Nx,Gkqset%ngk(1,1)))
    allocate(jxG_temp(Nx,Gkqset%ngk(1,1)))
    ! sarēķina sbessel visiem G garumiem, kas tiek lietoti.
    do igk = 1,Gkqset%ngk(1,1)
      call sbessel_vec(lmax, gridx * Gkqset%gkc(igk,1,1),jlxG_temp(:,:,igk))
    end do
    print *, "Maximum length of G that is used: ", maxval(Gkqset%gkc(:,1,1))
    do l = 0,lmax
      ! Jāsarēķina sfēriskās Beseļa funkcijas fiksētam l visiem vajadzīgajiem G.
      ! Savācam vērtības fiksētam l
      jxG_temp = jlxg_temp(l,:,:)
      ! Sarēķin��m izvirzījumu koeficientus.
      call dgemm('T','N',n,Gkqset%ngk(1,1),Nx,1._dp,U_l(:,:,l),Nx,jxG_temp,Nx,0._dp,svdfac(:,:,l),n)

    end do

#define modFeld
#ifdef modFeld
    ! Pārbaudām, ka projekcijas ir pareizas.
    allocate(check_tmp(Nx, Gkqset%ngk(1,1)))
    do l=0,lmax
      call dgemm('N','N',Nx,Gkqset%ngk(1,1),n,1._dp,U_l(:,:,l),Nx,svdfac(:,:,l),n,0._dp,check_tmp,Nx)
      jxG_temp = jlxG_temp(l,:,:)
      print *, "For l =", l, " the maximum projection error is ", maxval(check_tmp - jxG_temp)
    end do
#endif

!#define modFeld1
#ifdef modFeld1
! Atrodam, kurām sbese��a funkcijam ir lielie koeficienti.
    do l=0,lmax
      do igk = 1, Gkqset%ngk(1,1)
        if (maxval(abs(svdfac(:,igk,l)))>1._dp) then
          print *, "Lieli koeficienti pie l=",l," g garuma = ", Gkqset%gkc(igk,1,1),"."
          print *, "Koeficienti: ", svdfac(:,igk,l)
        end if
      end do
    end do
#endif
! Maksimālā kļ��da ir pie l=0, kur tā ir e-11

    write(*,*) Gkqset%ngk(1,1)
    write(*,*) allocated(Gkqset%sfacgk)
    write(*,*) allocated(Gkqset%tpgkc)

    allocate(B(n,lmmaxvr,nstfv))
    allocate(ylmgk(lmmaxvr,Gkqset%ngk(1,1)))
    B = 0._dp

    Do igk = 1, Gkqset%ngk(1,1)
      Call genylm(input%groundstate%lmaxvr, Gkqset%tpgkc(:, igk,1,1), ylmgk(:, igk))
    End Do


    do igk = 1, ngk(1,1)
      do ist = 1, nstfv ! index state
        do l = 0, input%groundstate%lmaxvr
          do m= -l,l
            lm = idxlm(l,m)
            do is=1,n
              B(is,lm,ist)= B(is,lm,ist) + omega_inv*wf%gk(igk,ist)*fourpi* cmplx(0,1, kind=dp) ** l * Gkqset%sfacgk(igk,1,1,1) * conjg(ylmgk(lm,igk))*svdfac(is,igk,l)
            end do !s
          end do ! m
        end do ! l
      end do ! ist
    end do ! igk
    !write(*,*) B
    ! do ist=1,nstfv
    !   print *, B(1,:,ist)
    !   read(*,*)
    ! end do



    ! Psi sfērisko harmoniku koeficienti ir izvirzīti SVD bāzē, kur izvirzījuma koef ir B.
    ! Psi_nlm(r) = Chi_mat * B_vec



    ! call dgemm('N','N',Nx,Gkqset%ngk(1,1),n,1._dp,U_l(:,:,l),Nx,svdfac(:,:,l),n,0._dp,check_tmp,Nx)

    ! Sarēķinām sfērisko harmoniku koeficientus
    allocate(Psi_lmn_approx(Nx,lmmaxvr,nstfv))

    allocate(zU_l(Nx,n,0:lmax))
    zU_l = U_l
    do ist = 1, nstfv
      do l = 0,input%groundstate%lmaxvr
        lm1 = idxlm(l,-l)
        lm2 = idxlm(l,l)
        call zgemm('N','N',Nx,lm2-lm1+1,n,1._dp, zU_l(:,:,l),Nx, B(:,lm1:lm2,ist),n,0._dp,Psi_lmn_approx(:,lm1:lm2,ist),Nx)
      end do ! l
    end do ! ist

    !-----------------------------------------------------------------------
    ! Saglabājam Psi_lmn_approx failā.
    open(unit=10, file='datapsilm.bin', form='unformatted', access='stream', status='replace')

    ! Write dimensions first so the reader is self-contained
    write(10) Nx, lmmaxvr, nstfv

    ! Write the x-grid
    write(10) gridX

    ! Write the full A array (Fortran column-major order)
    write(10) Psi_lmn_approx

    close(10)

    write(*,*) 'Written: Nx=', Nx, ' lmmaxvr=', lmmaxvr, ' nstfv=', nstfv
    !-----------------------------------------------------------------------


    ! ! Paskatāmies, vai ir jau kaut kas sarēķināts dotās viļņa f-as izvirzījumam.
    ! write(*,*) "Is allocated wf%mtrlm:", allocated(wf%mtrlm)

    ! ! sarēķinām to pašu vi��ņu funckijas izvirzījumu ar formulu:
    ! ! psi_nlm(r) = SUM_G psi_nG * 4 * Pi * i**l * j_l(G*r) * Y_lm(G_phth)
    ! write(*,*) "ngk = ", ngk(1,1)
    ! allocate(Psi_lmn(Nx,lmmaxvr,nstfv))
    ! Psi_lmn = cmplx(0,0,kind=dp)
    ! do ist = 1,nstfv
    !   do l = 0, input%groundstate%lmaxvr
    !     jxG_temp = jlxg_temp(l,:,:)
    !     do m = -l, l
    !       lm = idxlm(l,m)
    !       do igk = 1, ngk(1,1)
    !         Psi_lmn(:,lm,ist) = Psi_lmn(:,lm,ist) + wf%gk(igk, ist)* Gkqset%sfacgk(igk,1,1,1)* conjg(ylmgk(lm,igk)) * fourpi * cmplx(0,1, kind=dp) ** l * jxG_temp(:,igk)
    !       end do ! igk
    !     end do ! m
    !   end do ! l
    ! end do ! ist

    ! write(*,*) "Maksimālā kļūda aproksimācijā= ", maxval((abs(Psi_lmn - Psi_lmn_approx)))
  end subroutine testing

  subroutine test_wavefunction(wf)
    !! Checking the state of the wavefunction wf
    type(WFType), intent(in) :: wf
    write(*,*) "wf%gk", allocated(wf%gk)
    write(*,*) "wf%ir",allocated(wf%ir)
    write(*,*) "wf%losize",allocated(wf%losize)
    write(*,*) "wf%mt",allocated(wf%mt)
    write(*,*) "wf%mtmesh",allocated(wf%mtmesh)
    write(*,*) "wf%mtrlm",allocated(wf%mtrlm)

    ! Checking the wavefunction with index 5
    !wf%gk()



  end subroutine test_wavefunction





end module mod_Feldmans
