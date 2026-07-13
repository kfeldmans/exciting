!
!BOP
! !ROUTINE: init1
! !INTERFACE:
!
!
Subroutine gengntyyy !(gntyyy)
! !USES:
  Use modinput
  Use modmain
!      Use wigner3j_symbol, only : gaunt_yyy
! !DESCRIPTION:
!  Gaunt coefficients <Y3 | Y2 | Y1 >.
!
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
  Implicit None
! local variables
  Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
  real(8), external :: oldgaunt

! external functions
!      real(8), external :: gaunt_yyy
!
! allocate and generate complex Gaunt coefficient array
!      If (allocated(gntyry)) deallocate (gntyry)
  Allocate (gntyyy(lmmaxvr, lmmaxapw, lmmaxapw))
  Do l1 = 0, input%groundstate%lmaxapw
    Do m1 = - l1, l1
      lm1 = idxlm (l1, m1)
      Do l2 = 0, input%groundstate%lmaxapw
        Do m2 = - l2, l2
          lm2 = idxlm (l2, m2)
          Do l3 = 0, input%groundstate%lmaxapw
            Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)
              gntyyy (lm2, lm3, lm1) = oldgaunt (l1, l2, l3, m1, m2, m3)
            End Do
          End Do
        End Do
      End Do
    End Do
  End Do

  Return
end subroutine gengntyyy
!EOC
