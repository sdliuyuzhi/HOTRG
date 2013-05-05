MODULE coarse_grain
IMPLICIT NONE

CONTAINS
    SUBROUTINE coarsegrain( Temp, InitD, Dbond, SysSize, FEnergy,   &
                             TrunError, DomTr, id )
    !
    !   Purpose:
    !       Perform coarse grain process by doing contraction iteratively.
    !       'Temp' is the input temperature parameter.
    !       'InitD' is the input initial tensor bond dimension.
    !       'Dbond' is bond dimension after the initial cut, which will be used
    !   to perform further cut.
    !       'SysSize' is the input iteration numbers; usually the volume of the
    !   system V= 2^SysSize. 
    !       'FEnergy' is the output free energy density.
    !       'TrunError' is the output truncation error.
    !       'DomTr' is the output 
    !   By Yuzhi Liu on 04/28/2013. V1.

    USE set_precisions
    USE renorm_tensor
    USE get_tensor
    IMPLICIT NONE

    !   Data dictionary: declare calling parameter types & definitions
    REAL(kind=DBL), INTENT(IN) :: Temp                   ! Input temp and  
    INTEGER, INTENT(IN) :: InitD                        ! Input InitD
    INTEGER, INTENT(IN) :: Dbond                        ! Input Dbond
    INTEGER, INTENT(IN) :: SysSize                      ! Input system size
    REAL(kind=DBL), INTENT(OUT) :: FEnergy              ! Output Free Energy
    REAL(kind=DBL), INTENT(OUT) :: TrunError            ! Output errors
    REAL(kind=DBL), INTENT(OUT), DIMENSION(SysSize) :: DomTr  ! Output errors
    INTEGER, INTENT(OUT), DIMENSION(1) :: id                          ! Output dominant id

    !   Data dictionary: declare local variable types & definitions
    REAL(kind=DBL), DIMENSION(Dbond,Dbond,Dbond,Dbond) :: T    ! Input array Ta and Tb=Ta
    REAL(kind=DBL), DIMENSION(Dbond,Dbond,Dbond,Dbond) :: outT ! Temporary input array Ta and Tb=Ta
    REAL(kind=DBL) :: Coef                      ! Normalization factor
    INTEGER :: iterNo
    REAL(kind=DBL) :: FEContrib                ! Free energy density contribution
    REAL(kind=DBL) :: x

!    T = RESHAPE( (/ 2.7621956910836323_DBL, 0.0_DBL,   &
!                   0.0_DBL,    3.6268604078470186_DBL, &
!                   0.0_DBL,    3.6268604078470186_DBL, &
!                   3.6268604078470186_DBL, 0.0_DBL,    &
!                   0.0_DBL,    3.6268604078470186_DBL, &
!                   3.6268604078470186_DBL, 0.0_DBL,    &
!                   3.6268604078470186_DBL, 0.0_DBL,    &
!                   0.0_DBL, 4.7621956910836305_DBL /) &
!       , (/2,2,2,2/) )

    ! Get initial Tensor
    CALL GetTensor( Temp, InitD, Dbond, T, FEnergy, TrunError )  


    !FEnergy = 0.0_DBL
    
    ! initialize DomTr
    DomTr = 0.0_DBL

    ! perform coarse grain
    DO iterNo = 1, SysSize
        CALL RenormTensor( T, Dbond, outT, Coef, TrunError )
        outT = RESHAPE( outT, SHAPE = SHAPE(outT), ORDER = (/3,4,2,1/) )
        T = outT
        DomTr(iterNo) = TrunError
        FEContrib = - Temp * LOG(Coef) * 2.0_DBL**(-iterNo)
        x = ABS( FEContrib ) / MAX(  ABS(FEnergy), EPSILON(1.0_DBL) )
        FEnergy = FEnergy + FEContrib
        WRITE(*,*) 'iterNo: ',iterNo 
        WRITE(*,*) 'Coef: ',Coef
        WRITE(*,*) 'FEnergy:',FEnergy
        WRITE(*,*) 'FEContrib: ',FEContrib
        WRITE(*,*) 'TrunError: ',TrunError
        !WRITE(*,*) 'T: ', T
        ! stop if converges
        IF( x < 1.0E-10 .AND. MOD( iterNo, 4) == 0 ) THEN
            WRITE(*,*) 'Free energy converges'
        END IF
    END DO

    ! find the dominant error and it's positions
    DomTr = MAXVAL( DomTr )
    id = MAXLOC( DomTr )

    END SUBROUTINE coarsegrain

END MODULE coarse_grain
