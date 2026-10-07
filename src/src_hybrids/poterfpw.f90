module poterf
        
      implicit none
      private
    
      public :: poterfpw, poterfpw_sp
    
      contains

subroutine poterfpw_sp (ngvec1, rhoir,rhomt,igfft,sfacgq,ylmgq,gqc, jlgqsmallr,potir, potmt)
      use modbess, only: lambda
      use mod_atoms, only: nspecies,natoms, natmtot, idxas
      use mod_lattice, only: omega
      use modinput
      use mod_muffin_tin, only: rcmt,nrmtmax,nrmt,idxlm,lmmaxvr
      use modinteg
      use mod_hybrids, only : gmax_pw_method
      use constants, only: fourpi,zil,pi,zzero
      use mod_Gvector, only: cfunir,ngrid,ngrtot,ngvec
      Implicit None
! arguments
      integer, Intent (In)      :: ngvec1
      complex (8), Intent (In)  :: rhoir(:)     !ngrtot
      complex (8), Intent (In)  :: rhomt(:,:,:) !lmmaxvr,nrmtmax,natmtot
      integer, intent(in)       :: igfft(:)     !ngrtot
      complex (8), Intent (In)  :: sfacgq(:,:)  !ngvec, natmtot
      real (8), Intent (In)     :: gqc(:)
      complex (8), Intent (In)  :: ylmgq(:, :)  !lmmaxvr, ngvec
      real (8), Intent (In)     :: jlgqsmallr(:,0:,:,:)  !nrmtmax,0:input%groundstate%lmaxvr,ngvec1,nspecies
      complex (8), Intent (Out) :: potir(:)      ! ngrtot
      complex (8), Intent (Out) :: potmt(:,:,:)  !lmmaxvr,nrmtmax,natmtot
      

      complex (8) :: zt1,zt2,zt3,zt4
      complex (8) :: vig(ngvec1)
      integer :: lmaxvr,is,ia,ias,ig,ifg,l,m,lm,ir,nr,sfld
      real (8) :: fr(nrmtmax)
      real (8) :: refr(nrmtmax,lmmaxvr),imfr(nrmtmax,lmmaxvr)
      real (4) :: sfr(nrmtmax), srefr(nrmtmax,lmmaxvr),simfr(nrmtmax,lmmaxvr), sjl(nrmtmax)
      real (8) :: reint,imint
      real (4) :: sreint,simint
      real (8), external :: ddot
      real (4), external :: sdot


      lmaxvr=input%groundstate%lmaxvr

potir(:) = cfunir(:)*rhoir(:)

Call zfftifc (3, ngrid, -1, potir(:))  !to G space

!!! Obtain Fourier coeficients of the density in the MT and add it in porir
do is=1,nspecies
   nr=nrmt(is)
   do ir=1,nr
     fr(ir)=mt_integw%intw(ir,is) * rcmt(ir,is)**2
   enddo
   do ia=1,natoms(is)
      ias=idxas(ia,is)
      do lm=1,lmmaxvr
        do ir=1,nr
!          refr(ir,lm)=fr(ir) * dble(rhomt(lm,ir,ias))
!          imfr(ir,lm)=fr(ir) * dimag(rhomt(lm,ir,ias))
          srefr(ir,lm)=real(fr(ir) * dble(rhomt(lm,ir,ias)),kind=4)
          simfr(ir,lm)=real(fr(ir) * dimag(rhomt(lm,ir,ias)),kind=4)

        enddo
      enddo
      do ig=1,ngvec1
! ngvec1 is the limit beyong which there are no longer G+q vectors longer than gmax_pw_method.
! But there is no guarantee that the first ngvec1 vectors do not contain too long G+q for non-zero q.
        if (gqc(ig) .lt. gmax_pw_method) then
!        if(.true.) then
          ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq)) 
          do l=0,lmaxvr
            zt3=0d0
            sjl(1:nr)=real(jlgqsmallr(1:nr,l,ig,is),kind=4)
            do m=-l,l
               lm=idxlm(l,m)
               reint=real(sdot(nr,srefr(1,lm),1,sjl,1),kind=8)
               imint=real(sdot(nr,simfr(1,lm),1,sjl,1),kind=8)
!               reint=ddot(nr,refr(1,lm),1,jlgqsmallr(1,l,ig,is),1)
!               imint=ddot(nr,imfr(1,lm),1,jlgqsmallr(1,l,ig,is),1)
               zt3=zt3+dcmplx(reint,imint)*ylmgq(lm,ig)
            enddo ! m
            potir(ifg)=potir(ifg)+zt3*4d0*pi*conjg(sfacgq(ig, ias))/(omega*zil(l)) !!!Fāzes reizinātājs sfacgq(ig, ias) ??
          enddo ! l
        endif
      enddo ! ig
   enddo ! ia 
enddo ! is


   do ig=1, ngvec1
     ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq))   
     if (gqc(ig) .lt. gmax_pw_method) then
!    if(.true.) then
       If (gqc(ig) .Gt. input%structure%epslat) Then  
         potir(ifg)=potir(ifg)*4d0*pi*exp(-gqc(ig)**2/(4d0*lambda**2))/gqc(ig)**2
       else !!!! erfc kernels G=0 *(-1)  
         potir(ifg)=-potir(ifg)*pi/lambda**2
       endif
     else
       potir(ifg)=zzero
     endif
     vig(ig)=potir(ifg)
   enddo

  
  potir=zzero
  do ig=1, ngvec1
    ifg=igfft(ig)
    potir(ifg)=vig(ig)
  enddo
   

!!!obtain radial MT functions from potir and store in pot%mtrlm(:,:,:) (lm,ir,ias)
   do is=1,nspecies
      nr=nrmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
!        refr=0d0
!        imfr=0d0

        srefr=0e0
        simfr=0e0
        do ig=1, ngvec1
          if (gqc(ig) .lt. gmax_pw_method) then
!          if(.true.) then
            zt1=vig(ig)*sfacgq(ig, ias)
            do l=0,lmaxvr
              sjl(1:nr)=real(jlgqsmallr(1:nr,l,ig,is),kind=4)
              do lm=idxlm(l,-l),idxlm(l,l)
                zt2=zt1*conjg(ylmgq(lm,ig))
                sreint=real(zt2,kind=4)
                simint=aimag(zt2)
                call saxpy(nr,sreint,sjl,1,srefr(1,lm),1)
                call saxpy(nr,simint,sjl,1,simfr(1,lm),1)
              enddo
            enddo
          endif
        enddo
        do l=0,lmaxvr
          zt1=4d0*pi*zil(l)
          do m=-l,l
            lm=idxlm(l,m)
            do ir=1,nr
              potmt(lm,ir,ias)=zt1*dcmplx(srefr(ir,lm),simfr(ir,lm))
            enddo
          enddo
        enddo
      enddo
   enddo
   Call zfftifc (3, ngrid, 1, potir(:)) !to realspace 



  end subroutine poterfpw_sp

subroutine poterfpw (ngvec1, rhoir,rhomt,igfft,sfacgq,ylmgq,gqc, jlgqsmallr,potir, potmt)
      use modbess, only: lambda
      use mod_atoms, only: nspecies,natoms, natmtot, idxas
      use mod_lattice, only: omega
      use modinput
      use mod_muffin_tin, only: rcmt,nrmtmax,nrmt,idxlm,lmmaxvr
      use modinteg
      use mod_hybrids, only : gmax_pw_method
      use constants, only: fourpi,zil,pi,zzero
      use mod_Gvector, only: cfunir,ngrid,ngrtot,ngvec
      Implicit None
! arguments
      integer, Intent (In)      :: ngvec1
      complex (8), Intent (In)  :: rhoir(:)     !ngrtot
      complex (8), Intent (In)  :: rhomt(:,:,:) !lmmaxvr,nrmtmax,natmtot
      integer, intent(in)       :: igfft(:)     !ngrtot
      complex (8), Intent (In)  :: sfacgq(:,:)  !ngvec, natmtot
      real (8), Intent (In)     :: gqc(:)
      complex (8), Intent (In)  :: ylmgq(:, :)  !lmmaxvr, ngvec
      real (8), Intent (In)     :: jlgqsmallr(:,0:,:,:)  !nrmtmax,0:input%groundstate%lmaxvr,ngvec1,nspecies
      complex (8), Intent (Out) :: potir(:)      ! ngrtot
      complex (8), Intent (Out) :: potmt(:,:,:)  !lmmaxvr,nrmtmax,natmtot
      

      complex (8) :: zt1,zt2,zt3,zt4
      complex (8) :: vig(ngvec1)
      integer :: lmaxvr,is,ia,ias,ig,ifg,l,m,lm,ir,nr,sfld
      real (8) :: fr(nrmtmax)
      real (8) :: refr(nrmtmax,lmmaxvr),imfr(nrmtmax,lmmaxvr)
      real (8) :: reint,imint
      real (8), external :: ddot


      lmaxvr=input%groundstate%lmaxvr

potir(:) = cfunir(:)*rhoir(:)

Call zfftifc (3, ngrid, -1, potir(:))  !to G space

!!! Obtain Fourier coeficients of the density in the MT and add it in porir
do is=1,nspecies
   nr=nrmt(is)
   do ir=1,nr
     fr(ir)=mt_integw%intw(ir,is) * rcmt(ir,is)**2
   enddo
   do ia=1,natoms(is)
      ias=idxas(ia,is)
      do lm=1,lmmaxvr
        do ir=1,nr
          refr(ir,lm)=fr(ir) * dble(rhomt(lm,ir,ias))
          imfr(ir,lm)=fr(ir) * dimag(rhomt(lm,ir,ias))
        enddo
      enddo

      do ig=1,ngvec1
! ngvec1 is the limit beyong which there are no longer G+q vectors longer than gmax_pw_method.
! But there is no guarantee that the first ngvec1 vectors do not contain too long G+q for non-zero q.
        if (gqc(ig) .lt. gmax_pw_method) then
!        if(.true.) then
          ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq)) 
          do l=0,lmaxvr
            zt3=0d0
            do m=-l,l
               lm=idxlm(l,m)
               reint=ddot(nr,refr(1,lm),1,jlgqsmallr(1,l,ig,is),1)
               imint=ddot(nr,imfr(1,lm),1,jlgqsmallr(1,l,ig,is),1)
               zt3=zt3+dcmplx(reint,imint)*ylmgq(lm,ig)
            enddo ! m
            potir(ifg)=potir(ifg)+zt3*4d0*pi*conjg(sfacgq(ig, ias))/(omega*zil(l)) !!!Fāzes reizinātājs sfacgq(ig, ias) ??
          enddo ! l
        endif
      enddo ! ig
   enddo ! ia 
enddo ! is


   do ig=1, ngvec1
     ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq))   
     if (gqc(ig) .lt. gmax_pw_method) then
!    if(.true.) then
       If (gqc(ig) .Gt. input%structure%epslat) Then  
         potir(ifg)=potir(ifg)*4d0*pi*exp(-gqc(ig)**2/(4d0*lambda**2))/gqc(ig)**2
       else !!!! erfc kernels G=0 *(-1)  
         potir(ifg)=-potir(ifg)*pi/lambda**2
       endif
     else
       potir(ifg)=zzero
     endif
     vig(ig)=potir(ifg)
   enddo

  
  potir=zzero
  do ig=1, ngvec1
    ifg=igfft(ig)
    potir(ifg)=vig(ig)
  enddo
   

!!!obtain radial MT functions from potir and store in pot%mtrlm(:,:,:) (lm,ir,ias)
   do is=1,nspecies
      nr=nrmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        refr=0d0
        imfr=0d0

        do ig=1, ngvec1
          if (gqc(ig) .lt. gmax_pw_method) then
!          if(.true.) then
            zt1=vig(ig)*sfacgq(ig, ias)
            do l=0,lmaxvr
              do lm=idxlm(l,-l),idxlm(l,l)
                zt2=zt1*conjg(ylmgq(lm,ig))
                reint=dble(zt2)
                imint=dimag(zt2)
                call daxpy(nr,reint,jlgqsmallr(1,l,ig,is),1,refr(1,lm),1)
                call daxpy(nr,imint,jlgqsmallr(1,l,ig,is),1,imfr(1,lm),1)
              enddo
            enddo
          endif
        enddo
        do l=0,lmaxvr
          zt1=4d0*pi*zil(l)
          do m=-l,l
            lm=idxlm(l,m)
            do ir=1,nr
              potmt(lm,ir,ias)=zt1*dcmplx(refr(ir,lm),imfr(ir,lm))
            enddo
          enddo
        enddo
      enddo
   enddo
   Call zfftifc (3, ngrid, 1, potir(:)) !to realspace 



  end subroutine poterfpw

end module poterf
