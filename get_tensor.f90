MODULE get_tensor
IMPLICIT NONE

CONTAINS
    SUBROUTINE gettensor( Temp, InitD, Dbond, T, FEnergy, TrunError ) 
    !
    !   Purpose:
    !       Get the initial tensor T, FEnergy, and TrunError.
    !   By Yuzhi Liu on 04/28/2013. V1.
    
    USE set_precisions
    !USE bessel
    USE amos
    IMPLICIT NONE

    
    !   Data dictionary: declare calling parameter types & definitions
    REAL(kind=DBL), INTENT(IN) :: Temp                   ! Input temp and  
    INTEGER, INTENT(IN) :: InitD                        ! Input InitD
    INTEGER, INTENT(IN) :: Dbond                        ! Input Dbond
    REAL(kind=DBL), INTENT(OUT), DIMENSION(:,:,:,:) :: T! Output Free Energy
    REAL(kind=DBL), INTENT(OUT) :: FEnergy              ! Output Free Energy
    REAL(kind=DBL), INTENT(OUT) :: TrunError            ! Output errors

    REAL(kind=DBL), DIMENSION(InitD) :: I
    INTEGER ::  j,k,dimen
    INTEGER :: Lid, Kid, Jid, Iid

    !dimen = ( InitD -1 ) / 2

    !WRITE(*,*) 'dimen: ', dimen
    !k = 0
    ! InitD has to be odd number
    IF ( MOD( InitD, 2 ) == 0 ) THEN
        WRITE(*,*) 'InitD has better to be odd number.'
    END IF

    DO j = 1, InitD
        !k = k +1
       I(j) = SQRT( BESSI(j- ( InitD -1 )/2 - 1, 1.0_DBL/Temp) ) 
        !WRITE(*,*) 'k: ', k, 'I" ', I(k)
    END DO

    !WRITE(*,*) I

    DO Lid = 1, InitD
        DO Kid = 1, InitD
            DO Jid = 1, InitD
                DO Iid = 1, InitD
                    T ( Iid, Jid, Kid, Lid ) = I(Iid) * I (Jid )&
                                                * I( Kid ) * I( Lid ) &
                                                * Delta( Iid + Jid, Kid + Lid )
                END DO
            END DO
        END DO
    END DO

    !IF( InitD > Dbond ) THEN
    FEnergy = 0.0_DBL
    TrunError = 0.0_DBL
    
    END SUBROUTINE gettensor
        


    FUNCTION Delta( X1, X2 )
    USE set_precisions
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: X1, X2
    REAL(kind=DBL) :: Delta

    IF ( X1 == X2 ) THEN
        Delta = 1.0_DBL
    ELSE
        Delta = 0.0_DBL
    END IF

    END FUNCTION Delta


    
    FUNCTION BESSI(NN, X)
    !   Purpose:
    !       To calculate the modified Bessel function of the first kind I_NN(X)
    !  for integer nonnegative order NN and real X.
    !       It is based on the AMOS subroutines, Algorithm 644. 
    !       This is realized by setting imaginary part to be zero all the time:
    !   ZI = 0.0 and set the sequence number to be 1 all the time: N = 1. 
    !       Reference: http://www.netlib.org/toms/644
    !   NOTE:  
    !       Will need to check the constants defined in D1MACH, R1MACH, and I1MACH.
    !   By Yuzhi Liu 05/06/2013. V1 
    USE set_precisions
    USE amos
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN
    REAL(kind=DBL), INTENT(IN) :: X
    REAL(kind=DBL) :: BESSI

    REAL(kind=DBL) :: ZR, ZI
    REAL(kind=DBL) :: FNU
    INTEGER :: KODE  
    INTEGER :: N
    REAL(kind=DBL), DIMENSION(1) :: CYR, CYI
    INTEGER :: NZ
    INTEGER :: IERR
    
    ZR = X
    ZI = 0.0_DBL
    FNU = 1.0_DBL*ABS(NN)
    KODE = 1
    N = 1
    CALL ZBESI(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
    IF (IERR /= 0) THEN
        WRITE(*,*) 'Check the BESSI function for errors.'
        WRITE(*,*) 'IERR: ', IERR
    END IF
    
    BESSI = CYR(1)

END FUNCTION BESSI



END MODULE get_tensor
