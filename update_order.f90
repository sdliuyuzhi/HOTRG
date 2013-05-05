MODULE update_order
IMPLICIT NONE

INTERFACE updateorder
    MODULE PROCEDURE updateorder4
END INTERFACE

CONTAINS

    SUBROUTINE updateorder4( A, B, R, L, T, Coef, Comcoef )
    !
    !   Purpose:
    !       Contract tensors A, B, R, and L to form new tensor T.
    !       'A' and 'B' are initial 'T' tensor associated with initial sites.
    !       'R' and 'L' are Gauge tensors calculated in GetGaugeOrder.
    !       'ComCoef' is the input normalization factor.
    !       'T' is the output contracted tensor associated with new sites.
    !       'Coef' is the output normalization factor in order to keep the 
    !   array elements finite.
    !   A up and B down, act R and L, to merge as a single tensor T.
    !   A(1234), B(5647), R(26y), L(15x)
    !   By Yuzhi Liu on 04/28/2013; v1.

    USE set_precisions
    USE generic_contract
    IMPLICIT NONE

    !   Data dictionary: declare calling parameter types & definitions
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:,:) :: A     ! Input array A
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:,:) :: B     ! Input array B
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:)   :: R     ! Input array R
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:)   :: L     ! Input array L
    REAL(kind=DBL), OPTIONAL, INTENT(IN) :: ComCoef         ! Input normalization
    REAL(kind=DBL), INTENT(OUT) :: Coef                     ! Output normalization
    REAL(kind=DBL), INTENT(OUT), DIMENSION(:,:,:,:)  :: T  ! Output array outT

    !   Data dictionary: declare local variable types & definitions
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: AR ! Array A*R
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: ARB! Array A*R*B
    INTEGER :: istat, err 

    !WRITE(*,*) 'The array A is: ',A
    !WRITE(*,*) 'The array R is: ',R
    ! AR(1346y)
    ALLOCATE( AR(SIZE(A,1), SIZE(A,3), SIZE(A,4), SIZE(R,2), SIZE(R,3)),    &
                STAT = istat )
    ifAllocatedAR: IF ( ALLOCATED(AR) ) THEN
        CALL TensorContract( A, R, 4, (/2/), 3, (/1/), error = err, outT = AR )
    ELSE ifAllocatedAR
        WRITE(*,*) 'Warning-Array AR not allocated in update_order.'
    END IF ifAllocatedAR
    !WRITE(*,*) 'The array AR is: ',AR

    ! ARB(13y57)
    ALLOCATE( ARB(SIZE(A,1), SIZE(A,3), SIZE(R,3), SIZE(B,1), SIZE(B,4) ),    &
                STAT = istat )
    ifAllocatedARB: IF ( ALLOCATED(ARB) ) THEN
        CALL TensorContract(AR,B,5,(/3,4/),4,(/3,2/),error=err,outT=ARB)
    ELSE ifAllocatedARB
        WRITE(*,*) 'Warning-Array ARB not allocated in ordate_order.'
    END IF ifAllocatedARB
    !WRITE(*,*) 'The array ARB is: ',ARB

    DEALLOCATE( AR, STAT = istat )

    ! T(xy37)
    CALL TensorContract( ARB, L, 5, (/1,4/), 3, (/1,2/), ABOrder=(/4,2,1,3/), &
                        error = err, outT = T )
    !WRITE(*,*) 'The array T before normalization is: ',T

    ! find the normalization factor Coef
    IF( PRESENT(ComCoef) ) THEN
        Coef = ComCoef
    ELSE
        Coef = MAXVAL( ABS( T ) )
    END IF

    IF( ABS( Coef ) > 1.0E-12_DBL ) THEN
        T = T / Coef
    ELSE
        Coef = 1
    END IF
    !WRITE(*,*) 'The array T after normalization is: ',T

    END SUBROUTINE updateorder4 


END MODULE update_order
