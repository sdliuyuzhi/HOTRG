PROGRAM test_coarse_grain
!
!   Purpose:
!       To test the subroutine coarsegrain(( Temp, InitD, Dbond, SysSize,i
!    FEnergy,TrunError, DomTr, id )
!   By Yuzhi Liu 04/28/2013. V1
USE set_precisions
USE coarse_grain
IMPLICIT NONE

!   Data dictionary: declare variable types & definitions
REAL(kind=DBL) :: Temp                   ! Input temp and  
INTEGER :: InitD                        ! Input InitD
INTEGER :: Dbond                        ! Input Dbond
INTEGER :: SysSize                      ! Input system size
REAL(kind=DBL) :: FEnergy              ! Output Free Energy
REAL(kind=DBL) :: TrunError            ! Output errors
REAL(kind=DBL), DIMENSION(1) :: DomTr  ! Output errors
INTEGER,DIMENSION(1) :: id                          ! Output dominant id
CHARACTER(len=20) :: filename
INTEGER :: status

!Temp = 0.5_DBL
!InitD = 31
!Dbond = 31
!SysSize = 500

filename = 'param.txt'
OPEN(UNIT=12, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)
openif: IF ( status == 0 ) THEN
    
    ! OPEN was ok. READ values.
    READ(12,*, IOSTAT=status)           ! skip header
    READ(12,*, IOSTAT=status) Temp, InitD, Dbond, SysSize
    readif: IF( status > 0) THEN
        WRITE(*,*) 'An error occured reading parameter file.'
    ELSE readif
        WRITE(*,*) 'Parameters read to the program.'
    END IF readif
ELSE openif
    WRITE(*,*) 'Error opening parameter file: IOSTAT = ', status
END IF openif

CALL CoarseGrain( Temp, InitD, Dbond, SysSize, FEnergy, TrunError, DomTr, id)


END PROGRAM test_coarse_grain
