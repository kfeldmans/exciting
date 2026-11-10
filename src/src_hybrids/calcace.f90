!
!
!
! Copyright (C) 2021 J. Cimurs, A. Gulans
!
!
Subroutine calcACE (ikp, vnlvv, vxpsiir,wmt)
      Use modmain
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Use modmpi, only : rank
      use mod_fock_mt, only: nfmax, napwf, apw2f, nlof, lo2f
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (In) :: vnlvv (nstsv, nstsv)
      Complex (8), Intent (In) :: vxpsiir (ngkmax, nstsv)
! WP3: muffin-tin part of V_x psi_n, projected onto the sphere functions (FockExchange)
      Complex (8), Intent (In) :: wmt (max(1,nfmax), natmtot, nstsv)

! local variables
      Integer :: ngknr, ik, jk, ist1, ist2, ist3, is4
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm, ilm, irc ! ilm and lm potentially redundant (one of them)
      Integer :: nrc, iq, ig, iv (3), igq0, igk
      Integer :: ilo, loindex
      Integer :: info

      Real (8) :: v (3), cfq, ta,tb
      Complex (8) zrho01, zrho02, zt1, zt2
      Integer :: nr, l, m, io1, lm2, ir, if3
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Integer, Allocatable :: igkignr (:)
      Real (8), Allocatable :: vgklnr (:, :)
      Real (8), Allocatable :: vgkcnr (:, :)
      Real (8), Allocatable :: gkcnr (:)
      Real (8), Allocatable :: tpgkcnr (:, :)
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalsvp (:)
      Real (8), Allocatable :: evalsvnr (:)
      Real (8), Allocatable :: evalfv (:,:)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: sfacgknr (:, :)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :, :, :)
      Complex (8), Allocatable :: wfir1 (:, :, :)
      Complex (8), Allocatable :: wfir2 (:, :, :)
      Complex (8), Allocatable :: wfcr1 (:, :, :)
      Complex (8), Allocatable :: wfcr2 (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvcltp (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: matrixl(:, :)
      Complex (8), Allocatable :: hfxiir(:, :)
      Complex (8), Allocatable :: zfft(:)
      Complex (8), Allocatable :: xiintegral (:, :)
      Complex (8), Allocatable :: apwi(:, :)
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp
! for decomposition
      integer :: lwork, k, len
      complex(8), allocatable :: work(:)
      complex(8), allocatable ::  diag(:,:), matrixl1(:,:)
      real(8), allocatable :: rwork(:), eigvals(:)
! WP3: hfxw = wmt . matrixl^H, the MT part of the ACE vectors in coefficient space
      complex(8), allocatable :: hfxw(:,:,:)
      integer :: ilom




! -- Adaptively Compressed Exchange Operator starts --
      Allocate (matrixl(nstsv,nstsv))
      Allocate (hfxiir(ngkmax,nstsv))

      if (.not.(allocated(pace))) then
        Allocate (pace(nstsv, nmatmax, nkpt))
        pace=zzero
      endif
      Allocate (xiintegral(nstsv,haaijSize))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (apwi(haaijSize,ngkmax))

matrixl = -vnlvv ! create copy

! Remove upper triangular part
Do ist1 = 2, nstsv
    matrixl(ist1-1,ist1:)=0
End Do


        Allocate (diag(nstsv,nstsv))
        Allocate (work(nstsv))
        Allocate (rwork(3*nstsv-2))
        Allocate (eigvals(nstsv))


        call zheev('V', 'L', nstsv, matrixl,nstsv, eigvals, work, -1, rwork, info)
        k = int(work(1))

        deallocate(work)
        allocate(work(k))

        call zheev('V', 'L', nstsv, matrixl, nstsv, eigvals, work, k, rwork, info)
        !write(*,*)info, "eigval info"


        !--sqrt matrix
        diag = 0.d0
        do k =1, nstsv
                if (eigvals(k).gt.(1e-8)) then
                        diag(k,k) = cmplx(1/sqrt(eigvals(k)),0.d0)
                else 
!        write(*,*)"small eigval", eigvals(k)
                        diag(k,k) = cmplx(0.d0, 0.d0)
                end if
        end do


        matrixl = conjg(transpose(matrixl))
        allocate(matrixl1(nstsv, nstsv))
        matrixl1 = 0.d0
        call zgemm('N', 'N', nstsv, nstsv, nstsv,    dcmplx(1.0D0,0.0), diag, nstsv, matrixl, nstsv, dcmplx(0.0D0,0.0), matrixl1, nstsv)   !construct B = sqrt(D^-1)*L^-1
        matrixl = matrixl1
        
        deallocate(work, rwork, eigvals, diag, matrixl1)



        Call zgemm('N','C',ngkmax,nstsv,nstsv,(1.0D0,0.0),vxpsiir,ngkmax,matrixl,nstsv,(0.0D0,0.0),hfxiir,ngkmax) ! compute IR part of ACEO   hfxiir=zvclir*matrixl+


! WP3: hfxw(f,ias,i) = sum_n wmt(f,ias,n) conj(matrixl(i,n)) = int u_f(r) [hfxi_i]_lm(f)(r) r^2 dr,
! the muffin-tin part of the ACE vectors projected onto the sphere functions
Allocate (hfxw(size(wmt,1), natmtot, nstsv))
Call zgemm('N','C', size(wmt,1)*natmtot, nstsv, nstsv, zone, wmt, size(wmt,1)*natmtot, &
           matrixl, nstsv, zzero, hfxw, size(wmt,1)*natmtot)

! -- calculating IR part of <xi|phi>
pace(:,:,ikp)=0.d0
Do ist1 = 1, nstsv
  pace(ist1,1:ngk(1,ikp), ikp) = conjg(hfxiir(1:ngk(1,ikp),ist1))
End Do

! -- calculating MT and LO part of <xi|phi>
! The integrals int u conj(hfxi) r^2 dr are conj(hfxw) (apwfr, lofr real); apw2f and lo2f map
! the APW loop order (l, m, io) and the LO loop order (ilo, m) to the sphere function f.
Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), sfacgk(:, :, 1, ikp), apwalm)
loindex=0
Do is = 1, nspecies
    Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        if3=0
        Do l = 0, input%groundstate%lmaxmat
            Do m = -l, l
                Do io1 = 1, apword (l, is)
                    lm2 = idxlm (l, m)
                    if3=if3+1
                    apwi(if3,1:ngk(1,ikp))=apwalm(1:ngk(1,ikp), io1, lm2, ias)
                End Do
            End Do
        End Do
        if (if3 /= napwf(is) .or. sum(2*lorbl(1:nlorb(is),is)+1) /= nlof(is)) then
            write(*,'("Error(calcACE): APW/LO count differs from mod_fock_mt maps")')
            stop
        end if
        Do if3 = 1, napwf(is)
            xiintegral(:,if3) = conjg(hfxw(apw2f(if3,is),ias,:))
        End Do
        ! LO part; loindex runs on over atoms, lo2f counts per atom
        Do ilom = 1, nlof(is)
            pace(:,loindex+ilom+ngk(1,ikp),ikp) = conjg(hfxw(lo2f(ilom,is),ias,:))
        End Do
        loindex = loindex + nlof(is)
        Call zgemm('N', 'N', nstsv, ngk(1,ikp), napwf(is), zone, xiintegral, nstsv, apwi, haaijSize, &
                   zone, pace(:,:,ikp), nstsv) ! pace = paceMT + pace(IR+LO) = xiintegral*apwi + pace
    End Do
End Do

  Deallocate (matrixl, hfxw)
  Deallocate (hfxiir)
  Deallocate(xiintegral)
  Deallocate (apwi,apwalm)


!-- Adaptively Compressed Exchange Operator ends --

      Return
End Subroutine
!EOC
