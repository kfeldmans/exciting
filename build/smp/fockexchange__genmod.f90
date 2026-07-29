        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 28 14:17:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FOCKEXCHANGE__genmod
          INTERFACE 
            SUBROUTINE FOCKEXCHANGE(IKP,Q0CORR,VNLVV,VXPSIIRGK,VXPSIMT)
              USE MODMAIN
              INTEGER(KIND=4), INTENT(IN) :: IKP
              REAL(KIND=8), INTENT(IN) :: Q0CORR
              COMPLEX(KIND=8), INTENT(OUT) :: VNLVV(NSTSV,NSTSV)
              COMPLEX(KIND=8), INTENT(OUT) :: VXPSIIRGK(NGKMAX,NSTSV)
              COMPLEX(KIND=8), INTENT(OUT) :: VXPSIMT(LMMAXVR,NRCMTMAX, &
     &NATMTOT,NSTSV)
            END SUBROUTINE FOCKEXCHANGE
          END INTERFACE 
        END MODULE FOCKEXCHANGE__genmod
