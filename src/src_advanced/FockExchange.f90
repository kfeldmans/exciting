!
!
!
! Copyright (C) 2021 D. Zavickis, A. Gulans
!
!
Subroutine FockExchange(ikp, q0corr, vnlvv, vxpsiirgk, vxpsimt)
  use modrspace, only: rgrid_max_nmtpoints
  use modbess, only: nfit, zbessi, zbessk, erfc_fit, zilmt, lambda
  use modinteg
  Use modmain
  Use modinput
  Use modgw, only: kqset, Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
  Use potentials, only: coulomb_potential2, coulomb_potential2vFeld
  use weinert, only: poisson_mt_yukawa, pseudocharge_rspace_matrix
  use mod_hybrids, only: gmax_pw_method
  USE OMP_LIB
  use precision, only:dp

  use poterf

  use mod_Feldmans, only: calc_multipoles, hello, testing, &
    testing_calc, set_wf_planewaves, &
    testing_calc_2, psi_sphr_SVDb

  Implicit None
! arguments
  Integer, Intent(In) :: ikp
  Real(8), Intent(In) :: q0corr
  Complex(8), Intent(Out) :: vnlvv(nstsv, nstsv)
  Complex(8), Intent(Out) :: vxpsiirgk(ngkmax, nstsv)
  Complex(8), Intent(Out) :: vxpsimt(lmmaxvr, nrcmtmax, natmtot, nstsv)

! local variables
  Integer :: ngknr, ik, jk, ist1, ist2, ist3
  Integer :: is, ia, ias, ic, m1, m2, lmax, lm
  Integer :: nrc, iq, ig, iv(3), igq0, igk
  Integer :: ilo, loindex
  Integer :: info
  Integer :: ifg, ngvec1, ifit, ifit2
  Logical :: solver, cutoff, handleG0, rpseudo

  Real(8) :: v(3), cfq, ta, tb, t1, norm, uir, x
  Complex(8) zrho01, zrho02, ztmt, zt1, zt2, zt3, zt4, ztir
  Integer :: nr, l, m, io1, lm2, ir, if3, j, lmaxvr, ipt

  Complex(8) :: rhoG0, potG0
  Real(8) :: time_coul, tc, td, time_fft, time_prod, time_rs, time_misc, time_critical, time_pw
  real(dp) :: tf1, tf2
! automatic arrays
  Real(8) :: zn(nspecies)
  Complex(8) sfacgq0(natmtot)
! allocatable arrays
  Real(8), Allocatable :: vgqc(:, :)

  Real(8), Allocatable :: tpgqc(:, :)
  Real(8), Allocatable :: gqc(:)
  Real(8), Allocatable :: jlgqr(:, :, :)
  Real(8), Allocatable :: jlgq0r(:, :, :)
  Real(8), Allocatable :: evalfv(:, :)
  Complex(8), Allocatable :: prodir(:)
  Complex(8), Allocatable :: potir(:)
  Complex(8), Allocatable :: ylmgq(:, :)
  Complex(8), Allocatable :: sfacgq(:, :)
  Complex(8), Allocatable :: vxpsiirtmp(:), vxpsigktmp(:)
  Complex(8), Allocatable :: wfcr1(:, :)
  Complex(8), Allocatable :: wf1ir(:)
  Complex(8), Allocatable :: wf2ir(:)
  Complex(8), Allocatable :: potir0(:)
  Complex(8), Allocatable :: zvclmt(:, :, :, :)
  Complex(8), Allocatable :: zvcltp(:, :)
  Complex(8), Allocatable :: zfmt(:, :), zfmt0(:, :), zrhomt(:, :), potmt0(:, :, :)
  complex(8), allocatable :: q_nnLM_is(:,:,:,:) !! q_nn_LM_is(n,n',LM, ias) in q_nnLM_is(nstfv,nstfv,lmmaxvr, natmtot)
  real(8), allocatable :: jlgqsmallr(:, :, :, :), jlgrtmp(:), rfmt(:)
  complex(dp), allocatable :: qlmir(:,:)
  complex(dp), allocatable :: qlmir_in(:,:)
  complex(dp), allocatable :: qlmir_new(:,:)
  complex(dp), allocatable :: B_ik(:,:,:,:)
  complex(dp), allocatable :: B_jk(:,:,:,:)
  complex(dp), allocatable :: q_nnLM_is_oldnew(:,:,:,:)

  complex(8), Allocatable :: rpseudomat(:, :, :, :) ! (ifit,lm,igr,ias)

  type(WFType) :: wf1, wf2, prod, pot
  integer :: ik_test,ist3_test,jk_test,ist2_test
! external functions
  Complex(8) zfinp, zfmtinp, zfinpir, zfinpmt
  External zfinp, zfmtinp, zfinpir, zfinpmt
  logical :: print_times
  print_times = .true.

  ik_test = 4
  ist3_test = 2
  jk_test = 4
  ist2_test = 2

! allocate local arrays
  Allocate (vgqc(3, ngvec))
  Allocate (tpgqc(2, ngvec))
  Allocate (gqc(ngvec))
  Allocate (jlgqr(0:input%groundstate%lmaxvr + input%groundstate%npsden + 1, ngvec, nspecies))
  Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
  Allocate (ylmgq(lmmaxvr, ngvec))
  Allocate (sfacgq(ngvec, natmtot))
!      Allocate (vxpsiirtmp(ngrtot))   !dealoc
!      Allocate (vxpsigktmp(ngkmax))   !dealoc
  Allocate (wfcr1(ntpll, nrcmtmax))
  ! Allocate (zrhoir(ngrtot))
  Allocate (zvcltp(ntpll, nrcmtmax))
  Allocate (zfmt(lmmaxvr, nrcmtmax), zfmt0(lmmaxvr, nrcmtmax))
  Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot, nstsv)) !dealoc
  !Allocate (zwfir(ngkmax))


  !!write(*,*) "erfcapprox=",input%groundstate%hybrid%erfcapprox
  !write(*,*) "rpseudo=",input%groundstate%hybrid%rpseudo
  rpseudo = input%groundstate%hybrid%rpseudo
  lmaxvr = input%groundstate%lmaxvr
  !write(*,*)"lambda",lambda

  if (allocated(evalfv)) deallocate (evalfv)
  allocate (evalfv(nstfv, kset%nkpt))

  ! test comment for push

  evalfv(:, :) = 0.d0
  Do ik = 1, nkpt
    Call getevalfv(kset%vkl(:, ik), evalfv(:, ik))
  End Do
  call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)

  call WFInit(wf1)

  ik = kset%ikp2ik(ikp) ! 1d reduced index -> 1d non-reduced k-point index
  call genWF(ik, wf1)
  call genWFinMT(wf1)
  call genWFonMesh(wf1)

  call timesec(ta)
  call psi_sphr_SVDb(ik,wf1,B_ik)
  call timesec(tb)
  if (print_times) write(*,*) 'gen expansion for mod_Feld :', tb-ta
  !call hello()
  !call testing(wf1)

  ! ZEROED PARAMETERS HERE FELDMANS EDIT EDIT
  !call set_wf_planewaves(wf1)
  ! ZEROED PARAMETERS HERE FELDMANS EDIT EDIT
  !call testing_calc(ikp,wf1)
  !write(*,*) "Have exited testing_calc_2"
  !stop
  !call testing_calc(ik ,wf1, q_nnLM_is_oldnew,B_ik)

  call WFInit(wf2)

  t1 = 1/sqrt(omega)
! factor for long-range term
  cfq = 0.5d0*(omega/pi)**2
! set the nuclear charges to zero
  zn(:) = 0.d0
  vxpsiirgk(:, :) = 0.d0
  vxpsimt(:, :, :, :) = 0.d0
  vnlvv(:, :) = 0.d0
! calculate the wavefunctions for all states for the input k-point

! if (.true.) then
! start loop over non-reduced k-point set
  write(*,*) "Current non-reduced k-vector value in FockExchange: ", ik
  write(*,*) "Total number of different k-vals to calculate for, based in # if q vectors: ", kqset%nkpt
  Do iq = 1, kqset%nkpt

    write(*,*) "Current iq value in FockExchange: ", iq
    call timesec(ta)
    jk = kqset%kqid(ik, iq) ! ID(k') = ID(k-q)-> ID(k) set

! determine q-vector
    v(:) = kqset%vkc(:, ik) - kqset%vkc(:, jk)

! generate matrix for pseudocharge construction in real-space
    if (rpseudo) then
      if (.not. allocated(rpseudomat)) then
        if ((input%groundstate%hybrid%erfcapprox .eq. "truncatedYukawa") .or. (input%groundstate%hybrid%erfcapprox .eq. "Yukawa")) then
          allocate (rpseudomat(nfit*2 - 1, (lmaxvr + 1)**2, rgrid_max_nmtpoints, natmtot))
        else
          allocate (rpseudomat(1, (lmaxvr + 1)**2, rgrid_max_nmtpoints, natmtot))
        end if
      end if
      if ((input%groundstate%hybrid%erfcapprox .eq. "truncatedYukawa") .or. (input%groundstate%hybrid%erfcapprox .eq. "Yukawa")) then

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ifit,ifit2)
!!$OMP DO
!                do ifit=1,nfit
!                   ifit2=ifit
!                   write(*,*)OMP_GET_THREAD_NUM(),ifit2
!                   call pseudocharge_rspace_matrix (lmaxvr,input%groundstate%npsden,v,rpseudomat(ifit2,:,:,:),&
!                            & yukawa_in=.true.,zlambda=erfc_fit(ifit,2),zilmt=zilmt(ifit,:,:),zbessi=zbessi(:,ifit,:,:))
!                enddo
!!$OMP END DO NOWAIT
!!$OMP DO
!                do ifit=2,nfit
!                   ifit2=nfit+ifit-1
!                   write(*,*)OMP_GET_THREAD_NUM(),ifit2
!                   call pseudocharge_rspace_matrix (lmaxvr,input%groundstate%npsden,v,rpseudomat(ifit2,:,:,:),&
!                            & yukawa_in=.true.,zlambda=conjg(erfc_fit(ifit,2)),zilmt=conjg(zilmt(ifit,:,:)),zbessi=conjg(zbessi(:,ifit,:,:)))
!                enddo
!!$OMP END DO
!!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ifit,ifit2)
!$OMP DO
        do ifit2 = 1, nfit*2 - 1
          if (ifit2 .le. nfit) then
            ifit = ifit2
            call pseudocharge_rspace_matrix(lmaxvr, input%groundstate%npsden, v, rpseudomat(ifit2, :, :, :),&
              & yukawa_in=.true., zlambda=erfc_fit(ifit, 2), zilmt=zilmt(ifit, :, :), zbessi=zbessi(:, ifit, :, :))
          else
            ifit = ifit2 - nfit + 1
            call pseudocharge_rspace_matrix(lmaxvr, input%groundstate%npsden, v, rpseudomat(ifit2, :, :, :),&
              & yukawa_in=.true., zlambda=conjg(erfc_fit(ifit, 2)), zilmt=conjg(zilmt(ifit, :, :)), zbessi=conjg(zbessi(:, ifit, :, :)))
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL

      else
        call pseudocharge_rspace_matrix(lmaxvr, input%groundstate%npsden, v, rpseudomat(1, :, :, :))
      end if !Yukawa
    end if ! rpseudo

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig) SHARED(vgqc,vgc,v,ngvec,gqc,tpgqc,ylmgq,input)
!$OMP DO
    Do ig = 1, ngvec
! determine G+q vectors
      vgqc(:, ig) = vgc(:, ig) + v(:) ! Checked: vgc == Gset%vgc

! G+q-vector length and (theta, phi) coordinates
      Call sphcrd(vgqc(:, ig), gqc(ig), tpgqc(:, ig))
! spherical harmonics for G+q-vector
      Call genylm(input%groundstate%lmaxvr, tpgqc(:, ig), ylmgq(:, ig))
    End Do
!$OMP END DO
!$OMP END PARALLEL

! structure factors for G+q
    Call gensfacgp(ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
    Call findigp0(ngvec, gqc, igq0)
    sfacgq0(:) = sfacgq(igq0, :)
! compute the required spherical Bessel functions
    lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
    Call genjlgpr(lmax, gqc, jlgqr)
    Call genjlgq0r(gqc(igq0), jlgq0r)

!!!variables for the erfcapprox="PW"
    if (input%groundstate%hybrid%erfcapprox .eq. "PW") then
      !write(*,*)"shortest g+q vec",igq0
      gmax_pw_method = 6d0*dsqrt(2d0)*lambda ! 6*sigmas of a Gaussian
!         write(*,*)"min(g+q): G(",igq0,")=", gqc(igq0),"Largest G vector for PW method:",gmax_pw_method
      do ig = 1, ngvec
        if (gqc(ig) .lt. gmax_pw_method) then ! 6*sigmas of a Gaussian
          ngvec1 = ig + 1
!                  write(*,*)"G(",ig,")=",gqc(ig)
          !exit
        end if
      end do
      !write(*,*)"ngvec1=",ngvec1

      if (allocated(jlgqsmallr)) then
        deallocate (jlgqsmallr)
      end if
      if (allocated(jlgrtmp)) then
        deallocate (jlgrtmp)
      end if
      allocate (jlgqsmallr(nrcmtmax, 0:lmaxvr, ngvec1, nspecies))
      allocate (jlgrtmp(0:lmaxvr))

      do is = 1, nspecies
        do ig = 1, ngvec1
          do ir = 1, nrcmt(is)
            x = gqc(ig)*rcmt(ir, is)
            call sbessel(lmaxvr, x, jlgrtmp)
            jlgqsmallr(ir, :, ig, is) = jlgrtmp(:)
          end do
        end do
      end do

    end if

    call timesec(tb)
    if (print_times) write (*, *) 'qpt_init :', tb - ta
! calculate the wavefunctions for occupied states

    call timesec(ta)

    call genWF(jk, wf2)
    call genWFinMT(wf2)
    call genWFonMesh(wf2)


    !--------------------
    ! ZEROED PARAMETERS HERE FELDMANS EDIT EDIT
    !call set_wf_planewaves(wf2)
    !write(*,*) "The volume of unit element: ", omega
    ! ZEROED PARAMETERS HERE FELDMANS EDIT EDIT
    !----------------
    call timesec(tb)
    if (print_times) write (*, *) 'genWFs :', tb - ta

    call timesec(ta)
    call psi_sphr_SVDb(jk,wf2,B_jk)
    call timesec(tb)
    if (print_times) write(*,*) 'gen expansion for mod_Feld :', tb-ta

    solver = (input%groundstate%hybrid%singularity .ne. "exc")

    zvclmt(:, :, :, :) = 0.d0

    if ((input%groundstate%hybrid%erfcapprox .eq. "truncatedYukawa") .or.&
      & ((input%groundstate%hybrid%erfcapprox .eq. "none") .and. (input%groundstate%hybrid%singularity .eq. "exc0d"))) then
      cutoff = .True.
    else
      cutoff = .False.
    end if
    if ((input%groundstate%hybrid%erfcapprox .eq. "truncatedYukawa") .or.&
      & ((input%groundstate%hybrid%erfcapprox .eq. "none") .and. (input%groundstate%hybrid%singularity .eq. "exc0d")) .or.&
      & (input%groundstate%hybrid%erfcapprox .eq. "Yukawa")) then
      handleG0 = .false.
    else
      if (ik .eq. jk) then
        handleG0 = .true.
      else
        handleG0 = .false.
      end if
    end if

    !write(*,*)"cutoff", cutoff,"ik, jk",ik,jk,"handleG0",handleG0

    time_coul = 0d0
    time_fft = 0d0
    time_prod = 0d0
    time_rs = 0d0
    time_misc = 0d0
    time_critical = 0d0
    time_pw = 0d0
    !write(*,*)"B_ik is allocated:", allocated(B_ik)
    !write(*,*)"B_jk is allocated:", allocated(B_jk)
    Allocate (q_nnLM_is(nstfv,nomax,lmmaxvr,natmtot))
    q_nnLM_is = 0._dp
    write(*,*)"q_nnLM_is is allocated:", allocated(q_nnLM_is)
    call timesec(tf1)
    call calc_multipoles(B_ik,B_jk, q_nnLM_is)
    call timesec(tf2)
    write(*,*) "Multipoles calculated in: ", tf2-tf1
!write(*,*)"pirms", OMP_GET_THREAD_NUM()
!write(*,*)"nomax",nomax,"nstfv",nstfv
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qlmir_new,qlmir, ist3, wf1ir,wf2ir,prodir,potir,vxpsiirtmp,vxpsigktmp,potmt0,potir0, igk,ifg,prod,zrho01,pot,j,ifit2,rhoG0,tc,td,ist2,tf1,tf2)REDUCTION(max: time_coul) REDUCTION(max: time_fft) REDUCTION(max: time_prod) REDUCTION(max: time_rs) REDUCTION(max: time_misc) REDUCTION(max: time_critical) REDUCTION(max: time_pw)
    !ist2

    !write(*,*)"pēc", OMP_GET_THREAD_NUM()
    time_coul = 0d0
    time_fft = 0d0
    time_prod = 0d0
    time_rs = 0d0
    time_misc = 0d0
    time_critical = 0d0
    time_pw = 0d0
    call WFInit(prod)
    call WFInit(pot)

    Allocate (pot%mtrlm(lmmaxvr, nrmtmax, natmtot, 1))
    Allocate (potmt0(lmmaxvr, nrmtmax, natmtot))

    Allocate (prod%mtrlm(lmmaxvr, nrmtmax, natmtot, 1))
    Allocate (wf1ir(ngrtot))
    Allocate (wf2ir(ngrtot))
    Allocate (prodir(ngrtot))
    Allocate (potir(ngrtot))
    Allocate (potir0(ngrtot))
    Allocate (vxpsiirtmp(ngrtot))
    Allocate (vxpsigktmp(ngkmax))
    Allocate(qlmir_new(lmmaxvr, natmtot))
    allocate(qlmir(lmmaxvr,natmtot))
    Do ist2 = 1, nomax!1,1!!

      call timesec(tc)
      wf2ir(:) = 0.d0
      Do igk = 1, Gkqset%ngk(1, jk)
        ifg = igfft(Gkqset%igkig(igk, 1, jk))
        wf2ir(ifg) = t1*wf2%gk(igk, ist2)
      End Do
      Call zfftifc(3, ngrid, 1, wf2ir(:))

      call timesec(td)
      time_fft = time_fft + td - tc
!$OMP DO SCHEDULE(DYNAMIC)
      Do ist3 = 1, nstfv!1,1!!
        qlmir_new(:,:) = q_nnLM_is(ist3,ist2,:,:)
        vxpsiirtmp(:) = 0.d0
        call timesec(tc)
        wf1ir(:) = 0.d0
        Do igk = 1, Gkqset%ngk(1, ik)
          ifg = igfft(Gkqset%igkig(igk, 1, ik))
          wf1ir(ifg) = t1*wf1%gk(igk, ist3)
        End Do
        Call zfftifc(3, ngrid, 1, wf1ir(:))
        call timesec(td)
        time_fft = time_fft + td - tc

        ! calculate the complex overlap density
        !-----------------------------------------------------------------------------------
        call timesec(tc)
        call WFprodrs(ist2, wf2, ist3, wf1, prod)

        call timesec(td)
        time_rs = time_rs + td - tc
        call timesec(tc)
        prodir(:) = conjg(wf2ir(:))*wf1ir(:)
        call timesec(td)
        time_prod = time_prod + td - tc

        !if ((ik.eq.jk).and.(.not.solver)) then
        if (handleG0) then
          Call zrhogp(gqc(igq0), jlgq0r, ylmgq(:, &
            & igq0), sfacgq0, prod%mtrlm(:, :, :, 1), prodir(:), rhoG0)
          prodir(:) = prodir(:) - rhoG0
          prod%mtrlm(1, :, :, 1) = prod%mtrlm(1, :, :, 1) - rhoG0/y00
        end if

        call timesec(tc)
        pot%mtrlm(:, :, :, 1) = zzero
        potir = zzero
        if ((input%groundstate%hybrid%erfcapprox .eq. "truncatedYukawa") .or. (input%groundstate%hybrid%erfcapprox .eq. "Yukawa")) then !Yukawa case
          potmt0 = zzero
          potir0 = zzero
          ifit2 = 0
          do j = 1, nfit
            ifit2 = ifit2 + 1
            Call coulomb_potential2(nrcmt, rcmt, ngvec, gqc, igq0, &
              & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:, :, :, 1), &
              & prodir(:), potmt0, potir0, zrho02, v, &
              & cutoff=cutoff, hybrid_in=.true., yukawa_in=.true., &
              & zlambda_in=erfc_fit(j, 2), zbessi=zbessi(:, j, :, :), zbessk=zbessk(:, j, :, :), zilmt=zilmt(j, :, :),&
              & rpseudo_in=rpseudo, rpseudomat=rpseudomat(ifit2, :, :, :))
            pot%mtrlm(:, :, :, 1) = pot%mtrlm(:, :, :, 1) + potmt0*erfc_fit(j, 1)
            potir = potir + potir0*erfc_fit(j, 1)
          end do
          do j = 2, nfit
            ifit2 = ifit2 + 1
            Call coulomb_potential2(nrcmt, rcmt, ngvec, gqc, igq0, &
              & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:, :, :, 1), &
              & prodir(:), potmt0, potir0, zrho02, v, &
              & cutoff=cutoff, hybrid_in=.true., yukawa_in=.true., &
              & zlambda_in=conjg(erfc_fit(j,2)),zbessi=conjg(zbessi(:,j,:,:)),zbessk=conjg(zbessk(:,j,:,:)),zilmt=conjg(zilmt(j,:,:)),&
              & rpseudo_in=rpseudo, rpseudomat=rpseudomat(ifit2, :, :, :))
            pot%mtrlm(:, :, :, 1) = pot%mtrlm(:, :, :, 1) + potmt0*conjg(erfc_fit(j, 1))
            potir = potir + potir0*conjg(erfc_fit(j, 1))
          end do
        end if

        if ((input%groundstate%hybrid%erfcapprox .eq. "PW") .or. (input%groundstate%hybrid%erfcapprox .eq. "none")) then
          potmt0 = zzero
          potir0 = zzero

          Call coulomb_potential2vFeld(qlmir_new, qlmir, nrcmt, rcmt, ngvec, gqc, igq0, &
            & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:, :, :, 1), &
            & prodir(:), pot%mtrlm(:, :, :, 1), potir, zrho02, v, &
            & cutoff=cutoff, hybrid_in=.true.,&
            & rpseudo_in=rpseudo, rpseudomat=rpseudomat(1, :, :, :))

!                  if (input%groundstate%hybrid%erfcapprox.ne."PW")then
!                     pot%mtrlm(:,:,:,1)=potmt0
!                     potir=potir0
!                  endif
          !print '("OLD: For pair ", I0, ", ", I0, " obtained value ", F0.3)', ist2, ist3, maxval(abs(qlmir))
!#define TESTING_PRINTMULTIPOLECOMAPRISON_312312
#ifdef TESTING_PRINTMULTIPOLECOMAPRISON_312312
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              do l=0,input%groundstate%lmaxvr
                do m = -l,l
                  lm = idxlm(l,m)

                  !write(*,*) "Multipole moements for 5,5: ", qlmir
                  !print '("OLD: For l ", I0, ", and m ", I0, " on atom ", I0 ," obtained values ", F0.10,F0.10)', l, m,ias, qlmir(lm,ias)%re,(qlmir(lm,ias)%Im)
                  print '("OLD: For (ik, n) (", I0, ",", I0, ") (jk,m) (", I0,",",I0, ") a ", I0 ," lm ", I0, " obtained values ", F0.10,F0.10)',&
                    ik,ist3,jk, ist2,ias,lm, qlmir(lm,ias)%re,(qlmir(lm,ias)%Im)
                  print '("NEW: For (ik, n) (", I0, ",", I0, ") (jk,m) (", I0,",",I0, ") a ", I0 ," lm ", I0, " obtained values ", F0.10,F0.10)',&
                    ik,ist3,jk, ist2,ias,lm, qlmir_new(lm,ias)%re,(qlmir_new(lm,ias)%Im)
                  !print '(F0.10, " ", F0.10)',dble(qlmir(lm,ias)), dimag(qlmir(lm,ias))
                end do
              end do
            end do
          end do
          !write(*,*) "FockExchange: species radii: ", rmt
#endif TESTING_PRINTMULTIPOLECOMAPRISON_312312

!#define TESTING_7575808138123814824
#ifdef TESTING_7575808138123814824
          do ias = 1, natmtot
            do lm = 1, lmmaxvr
              print '("DIFF: For (ik, n) (", I0, ",", I0, ") (jk,m) (", I0,",",I0, ") a ", I0 ," lm ", I0, " the diff is ", F0.10,F0.10)',&
                ik,ist3,jk, ist2,ias,lm, abs(qlmir_new(lm,ias)-qlmir(lm,ias))
            end do
          end do
#endif TESTING_7575808138123814824

!#define TESTING_41731874817412
#ifdef TESTING_41731874817412
          !print '("Iteration: ", I0," ", I0, " ", I0," ", I0)', ik, ist3, jk, ist2
          if (ik==ik_test .and. ist3 == ist3_test .and. jk == jk_test .and. ist2==ist2_test) then
            do ias = 1,natmtot
              do lm = 1,lmmaxvr
                print '("For (ik, n) (", I0, ",", I0, ") (jk,m) (", I0,",",I0, ") a ", I0 ," lm ", I0, " OLD:", F0.10, " NEW:",F0.10)',&
                  ik,ist3,jk, ist2,ias,lm, abs(qlmir(lm,ias)),abs(qlmir_new(lm,ias))
                !print '("OLD: For (ik, n) (", I0, ",", I0, ") (jk,m) (", I0,",",I0, ") a ", I0 ," lm ", I0, " the diff is ", F0.10,F0.10)',&
                !ik,ist3,jk, ist2,ias,lm
              end do
            end do
          end if
#endif TESTING_41731874817412

!#define TESTING_MULTIPOLE_CALCULATION_8113211203141
#ifdef TESTING_MULTIPOLE_CALCULATION_8113211203141
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              write(*,*) "Scale of multipole: ", maxval(abs(qlmir(:,ias)))
              write(*,*) "Difference output: ", maxval(abs(qlmir(:,ias)-q_nnLM_is(ist3,ist2,:,ias)))
            end do
          end do



#endif TESTING_MULTIPOLE_CALCULATION_8113211203141

        end if
        call timesec(td)
        time_coul = time_coul + td - tc

!if (.false.) then
        if (input%groundstate%hybrid%erfcapprox .eq. "PW") then
          call timesec(tc)
          call poterfpw(ngvec1, prodir, prod%mtrlm(:, :, :, 1), igfft, sfacgq, ylmgq, gqc, jlgqsmallr, potir0, potmt0)
          potir = potir - potir0 !Coulomb - erf
          pot%mtrlm(:, :, :, 1) = -potmt0 + pot%mtrlm(:, :, :, 1)
          call timesec(td)
          time_pw = time_pw + td - tc
        end if

        ! call WFprodrs(ist2,wf2,ist3,wf1,prod)
        ! prodir(:)=conjg(wf2ir(:))*wf1ir(:)

        !if ((ik.eq.jk).and.(.not.solver)) then
        if (handleG0) then
          Call zrhogp(gqc(igq0), jlgq0r, ylmgq(:, &
            & igq0), sfacgq0, pot%mtrlm(:, :, :, 1), potir(:), zrho01)
          potir(:) = potir(:) - zrho01
          pot%mtrlm(1, :, :, 1) = pot%mtrlm(1, :, :, 1) - zrho01/y00
        end if

        if ((input%groundstate%hybrid%erfcapprox .eq. "PW") .and. (handleG0)) then
          potG0 = rhoG0*pi/lambda**2
          potir(:) = potir(:) + potG0
          pot%mtrlm(1, :, :, 1) = pot%mtrlm(1, :, :, 1) + potG0/y00
        end if

        !-----------------------------------------------------------------------------------1

        call timesec(tc)
        call prodshrs(pot%mtrlm(:, :, :, 1), wf2%mtmesh(:, :, :, ist2), prod%mtrlm(:, :, :, 1))
!               call genWFonMeshOne(pot)
!               pot%mtmesh=conjg(pot%mtmesh)
!   call timesec(td)
!   time_misc=time_misc+td-tc

!   call timesec(tc)
!               call WFprodrs(1,pot,ist2,wf2,prod)
        call timesec(td)
        time_misc = time_misc + td - tc
!   time_rs=time_rs+td-tc

        call timesec(tc)
        vxpsiirtmp(:) = potir(:)*wf2ir(:)*cfunir(:)
!               vxpsiirtmp(:) = vxpsiirtmp(:)*cfunir(:)
        call timesec(td)
        time_prod = time_prod + td - tc
        ! ----------------------------------------------------------------------------------
        ! Calculate Fourier transform of vxpsiirtmp(r)
        call timesec(tc)
        Call zfftifc(3, ngrid, -1, vxpsiirtmp)

        Do igk = 1, Gkqset%ngk(1, ik)
          vxpsigktmp(igk) = vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega)*wkptnr(jk)
        End Do

        call timesec(td)
        time_fft = time_fft + td - tc

        call timesec(tc)
!$OMP CRITICAL
        Do igk = 1, Gkqset%ngk(1, ik)
          vxpsiirgk(igk, ist3) = vxpsiirgk(igk, ist3) + vxpsigktmp(igk)
        End Do
        zvclmt(:, :, :, ist3) = zvclmt(:, :, :, ist3) + prod%mtrlm(:, :, :, 1)*wkptnr(jk)
!$OMP END CRITICAL
        call timesec(td)
        time_critical = time_critical + td - tc
!write(*,*) ist2,ist3
      End Do ! ist3
!$OMP END DO NOWAIT
    End Do ! ist2

    call WFRelease(prod)
    call WFRelease(pot)
    Deallocate (wf2ir)
    Deallocate (prodir)
    Deallocate (potir)
    Deallocate (potir0)
    Deallocate (wf1ir)
    Deallocate (potmt0)
    Deallocate (vxpsiirtmp)
    Deallocate (vxpsigktmp)
    deallocate(qlmir)
    deallocate(qlmir_new)
!$OMP END PARALLEL
    deallocate (q_nnLM_is)

    call timesec(ta)
    if (print_times) then
      write (*, *) 'omp_loop :', ta - tb
      write (*, *) '  time_coul :', time_coul
      write (*, *) '  time_fft :', time_fft
      write (*, *) '  time_prod :', time_prod
      write (*, *) '  time_misc :', time_misc
      write (*, *) '  time_rs :', time_rs
      write (*, *) '  time_critical :', time_critical
      write (*, *) '  time_pw :', time_pw
    end if
    vxpsimt = vxpsimt + zvclmt
  End Do ! non-reduced k-point set
!stop
!----------------------------------------------!
!     valence-core-valence contribution        !
!----------------------------------------------!
  call timesec(ta)
  zvclmt(:, :, :, :) = 0.d0

  Allocate (zrhomt(lmmaxvr, nrcmtmax))
  allocate (rfmt(nrcmtmax))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nrc,is,ia,ias,ist2,ist3,l,norm,lm,m,ir,uir,wfcr1,zrhomt,zfmt,ifit,zfmt0,zvcltp,rfmt,ipt)
  Do is = 1, nspecies
    if (spnst(is) .ne. 0) then
      nrc = nrcmt(is)
      Do ir = 1, nrmt(is)
        rfmt(ir) = 1d0/spr(ir, is)
      End Do

      Do ia = 1, natoms(is)
        ias = idxas(ia, is)
!$OMP DO SCHEDULE(DYNAMIC)
        ! Begin loop over occupied   and empty states
        Do ist3 = 1, nstsv
          Do ist2 = 1, spnst(is) !This is essentially ncore(is)
            If (spcore(ist2, is)) Then
              l = spl(ist2, is)
              norm = sqrt(0.5d0*spocc(ist2, is)/dble(2*l + 1))
              Do m = -l, l
                lm = idxlm(l, m)

! calculate the complex overlap density

                Do ir = 1, nrmt(is)
                  uir = norm*rwfcr(ir, 1, ist2, ias)*rfmt(ir)
                  do ipt = 1, ntpll
                    zvcltp(ipt, ir) = uir*conjg(zbshthf(ipt, lm))*wf1%mtmesh(ipt, ir, ias, ist3)
                  end do ! ipt
                End Do ! ir

                Call zgemm('N', 'N', lmmaxvr, nrcmt(is), ntpll, zone, &
                  & zfshthf, lmmaxvr, zvcltp, ntpll, zzero, zrhomt, lmmaxvr)

! calculate the Coulomb potential
                if (input%groundstate%hybrid%erfcapprox .eq. "none") then
                  call poisson_mt_yukawa(lmaxvr, nrc, rcmt(:, is), zrhomt, zfmt, is)
                else
                  zfmt = zzero
                  do ifit = 1, nfit
                    call poisson_mt_yukawa(lmaxvr, nrc, rcmt(:, is), zrhomt, zfmt0, is, &
                      & yukawa_in=.true., zlambda=erfc_fit(ifit, 2), il=zbessi(:, ifit, :, is), kl=zbessk(:, ifit, :, is))
                    zfmt = zfmt + zfmt0*erfc_fit(ifit, 1)
                  end do
                  do ifit = 2, nfit
                    call poisson_mt_yukawa(lmaxvr, nrc, rcmt(:, is), zrhomt, zfmt0, is, &
                      & yukawa_in=.true., zlambda=conjg(erfc_fit(ifit, 2)), il=conjg(zbessi(:, ifit, :, is)), kl=conjg(zbessk(:, ifit, :, is)))
                    zfmt = zfmt + zfmt0*conjg(erfc_fit(ifit, 1))
                  end do
                end if

                Call zgemm('N', 'N', ntpll, nrc, lmmaxvr, &
                  & zone, zbshthf, ntpll, zfmt, lmmaxvr, &
                  & zzero, zvcltp, ntpll) ! Returns zvcltp in SC

! calculate the complex overlap density
                Do ir = 1, nrmt(is)
                  uir = norm*rwfcr(ir, 1, ist2, ias)*rfmt(ir)
                  do ipt = 1, ntpll
                    zvcltp(ipt, ir) = uir*zbshthf(ipt, lm)*zvcltp(ipt, ir)
                  end do ! ipt
                End Do ! ir

                Call zgemm('N', 'N', lmmaxvr, nrcmt(is), ntpll, zone, &
                  & zfshthf, lmmaxvr, zvcltp, ntpll, zone, zvclmt(1, 1, ias, ist3), lmmaxvr)

              End Do ! m
            End If ! spcore(ist2, is)
          End Do ! ist2
        End Do ! ist3
!$OMP END DO
      End Do ! ia
    End If !spnst
  End Do ! is
!$OMP END PARALLEL
  Deallocate (zrhomt, rfmt)
  call timesec(tb)
  if (print_times) write (*, *) 'vcv :', tb - ta

  vxpsimt = vxpsimt + zvclmt
  Allocate (wf1ir(ngrtot))
  Allocate (vxpsiirtmp(ngrtot))
  call timesec(ta)
! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,ist3,ztir,igk,ztmt)
! !$OMP DO collapse(2)
  Do ist1 = 1, nstsv

!!!!!!!!!!!!!!!!!!!!!!!!!! KOREKCIJA ko vajadzētu atslēgt Aux funct method
    if (input%groundstate%hybrid%singularity .eq. "exc") then
      write (*, *) "q0corr", q0corr
      If ((ist1 .le. nomax) .and. (q0corr .ne. 0.d0)) Then
        ! Evaluate wavefunction in real space
        wf1ir(:) = 0.d0
        Do igk = 1, Gkqset%ngk(1, ik)
          ifg = igfft(Gkqset%igkig(igk, 1, ik))
          wf1ir(ifg) = t1*wf1%gk(igk, ist1)
        End Do
        Call zfftifc(3, ngrid, 1, wf1ir(:))

        vxpsimt(:, :, :, ist1) = vxpsimt(:, :, :, ist1) + q0corr*wf1%mtrlm(:, :, :, ist1)

        ! Apply correction to IR part and roll back correction to momentum space
        vxpsiirtmp(:) = q0corr*wf1ir(:)*cfunir(:)
        Call zfftifc(3, ngrid, -1, vxpsiirtmp)
        Do igk = 1, Gkqset%ngk(1, ik)
          vxpsiirgk(igk, ist1) = vxpsiirgk(igk, ist1) + vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega) ! pace IR
        End Do
      end if ! (ist1.le.nomax).and.(q0corr.ne.0.d0)
    End If ! singularity.eq.exc
! !!!!!!!!!!!!!!!!!!!!!!!!!! KOREKCIJA ko vajadzētu atslēgt Aux funct method//
    ! Write(*,*) 'ist1=',ist1
    Do ist3 = 1, nstsv

      ztir = 0.d0

      Do igk = 1, Gkqset%ngk(1, ik)
        ztir = ztir + conjg(wf1%gk(igk, ist1))*vxpsiirgk(igk, ist3)
      End Do
      ztmt = zfinpmt(.True., wf1%mtrlm(:, :, :, ist1), vxpsimt(:, :, :, ist3))

      vnlvv(ist1, ist3) = vnlvv(ist1, ist3) - ztmt - ztir

    End Do ! ist3
  End Do ! ist1
! !$OMP END PARALLEL
  call timesec(tb)

  if (print_times) write (*, *) 'Matrix :', tb - ta

  if (.false.) then
    write (*, *) "element(1,1)", dble(vnlvv(1, 1)), ",", imag(vnlvv(1, 1))
    Write (*, *) "ikp, ik, memopt:", ikp, ik
    write (*, *) 'vnlvv real (1:nstfv,1:nstfv)'
    do ist1 = 1, min(nstfv, 14)
      write (*, '(14F13.9)') dble(vnlvv(ist1, 1:min(nstfv, 14)))
    end do

    write (*, *) 'vnlvv imag (1:nstfv,1:nstfv)--'
    do ist1 = 1, min(nstfv, 14)
      write (*, '(14F13.9)') dimag(vnlvv(ist1, 1:min(nstfv, 14)))
    end do
  end if

  Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
  Deallocate (ylmgq, sfacgq)
  Deallocate (wfcr1)
  Deallocate (wf1ir)
  Deallocate (zvcltp, zfmt, zfmt0)
  Deallocate (vxpsiirtmp)
!      Deallocate (vxpsiirtmp)
!      Deallocate (vxpsigktmp)
  Deallocate (zvclmt)
  if (allocated(jlgqsmallr)) then
    deallocate (jlgqsmallr)
  end if
  if (allocated(jlgrtmp)) then
    deallocate (jlgrtmp)
  end if
  if (allocated(rpseudomat)) Then
    Deallocate (rpseudomat)
  end if

  call WFRelease(wf1)
  call WFRelease(wf2)
  call WFRelease(prod)
  if (allocated(q_nnLM_is)) then
    deallocate(q_nnLM_is)
  end if

  !write(*,*)"FockExchange.f90 stop"
  !stop
  Return
end subroutine FockExchange
!EOC
