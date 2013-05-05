MODULE renorm_tensor
IMPLICIT NONE

CONTAINS
    SUBROUTINE renormtensor( T, Dbond, outT, Coef, TrunError )
    !
    !   Purpose:
    !       Get new tensor and the normalization factor from the old tensor by
    !   using 'GetGaugeOrder' and 'UpdateOrder' subroutine.
    !       'T' is the input old tensor.
    !       'Dbond' is the bond dimension of the old tensor T
    !       'outT' is the output new tensor
    !       'Coef' is the normalization factor
    !       'TrunError' is the output truncation error.
    !   By Yuzhi Liu on 04/28/2013. V1.
    
    USE set_precisions
    USE get_gauge_order 
    USE update_order
    IMPLICIT NONE

    !   Data dictionary: declare calling parameter types & definitions
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:,:) :: T     ! Input array T
    INTEGER, INTENT(IN) :: Dbond                            ! Input bond dimension
    REAL(kind=DBL), INTENT(OUT), DIMENSION(:,:,:,:) :: outT ! Output array T
    REAL(kind=DBL) :: Coef, TrunError                       ! Output 

    !   Data dictionary: declare local variable types & definitions
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:,:) :: Gauge  ! Gauge tensor
    INTEGER :: err, istat 

    ! get 'Gauge' tensor    
    ALLOCATE( Gauge( SIZE(T,1), SIZE(T,1), Dbond ), STAT = istat )
    ifAllocatedGauge: IF( ALLOCATED(Gauge) ) THEN
        CALL GetGaugeOrder( T, T, Dbond, Gauge, TrunError )
    ELSE ifAllocatedGauge
        WRITE(*,*) 'Warning-Array Gauge not allocated in renorm_tensor.'
    END IF ifAllocatedGauge 

    ! update old tensor T to get the new tensor outT
    CALL UpdateOrder( T, T, Gauge, Gauge, outT, Coef )  

    DEALLOCATE( GAUGE, STAT = istat )

    END SUBROUTINE renormtensor

END MODULE renorm_tensor
