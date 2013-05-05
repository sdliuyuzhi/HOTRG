MODULE get_gauge_order
IMPLICIT NONE
PRIVATE
PUBLIC :: getgaugeorder

INTERFACE getgaugeorder
    MODULE PROCEDURE getgaugeorder4
END INTERFACE

CONTAINS
    SUBROUTINE getgaugeorder4( Ta, Tb, Dbond, Gauge, TrunError )
    !
    !   Purpose:
    !       Contract tensor Ta and Tb and perform svd decomposition to get
    !   the unitary matrix U1 and U2.
    !       Ta and Tb are input tensors associated with sites. 
    !       Dbond is the (expected) bond dimension of the tensor.
    !       Gauge is the output truncated unitary matrix. 
    !       TrunError is the truncation error.
    !       Details definitions can be found in PRB 86 045139(2012)
    !       http://arxiv.org/abs/1201.1144
    !   By Yuzhi Liu on 04/26/2013; v1.

    USE set_precisions
    USE generic_contract
    IMPLICIT NONE

    !   Data dictionary: declare calling parameter types & definitions
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:,:) :: Ta    ! Input array A
    REAL(kind=DBL), INTENT(IN), DIMENSION(:,:,:,:) :: Tb    ! Input array B
    INTEGER, INTENT(IN) :: Dbond                            ! Bond dimension
    REAL(kind=DBL), INTENT(OUT), DIMENSION(:,:,:) :: Gauge    ! Output U matrix
    REAL(kind=DBL), INTENT(OUT) :: TrunError                ! Truncation error

    !   Data dictionary: declare local variable types & definitions
    INTEGER :: dim1, dim2                                   ! Dim of Ta and Tb
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:,:,:) :: T1   ! Temp arrays
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:,:,:) :: T2   ! Temp arrays
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: M1       ! Temp arrays
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: U1       ! UL arrays
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: U2       ! UR arrays
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:) :: Spectra1   ! A array
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:) :: Spectra2   ! A array
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: V        ! V arrays
    REAL(kind=DBL) :: Tr1, Tr2                              ! Traces of A
    INTEGER :: err, istat                                   ! error flag

    dim1 = SIZE( Ta, 1 )    
    dim2 = SIZE( Tb, 1 )    

    !                                                 
    ! Calculate U2 ( right side)           
    !                                                 

    ! T1(2,4,2',4')
    ALLOCATE( T1(SIZE(Ta,2), SIZE(Ta,4), SIZE(Ta,2), SIZE(Ta,4)), STAT=istat)
    ifAllocatedT1: IF( ALLOCATED(T1) ) THEN
        CALL TensorContract( Ta, Ta, 4, (/1,3/), 4, (/1,3/), error=err,outT=T1)
    ELSE ifAllocatedT1
        WRITE(*,*) 'Warning-Array T1 not allocated in getgaugeorder4.'
    END IF ifAllocatedT1 
    
    ! T2(6,4,6',4')
    ALLOCATE( T2(SIZE(Tb,2), SIZE(Tb,3), SIZE(Tb,2), SIZE(Tb,3)), STAT=istat)
    ifAllocatedT2: IF( ALLOCATED(T2) ) THEN
        CALL TensorContract( Tb, Tb, 4, (/1,4/), 4, (/1,4/), error=err,outT=T2 )
    ELSE ifAllocatedT2
        WRITE(*,*) 'Warning-Array T2 not allocated in getgaugeorder4.'
    END IF ifAllocatedT2 

    ! M1(26,2'6')
    ALLOCATE( M1(SIZE(T1,1)*SIZE(T2,1), SIZE(T1,3)*SIZE(T2,3)), STAT=istat )
    ifAllocatedM1: IF( ALLOCATED(M1) ) THEN
        CALL TensorContract( T1, T2, 4, (/2,4/), 4, (/2,4/), (/1,3,2,4/), &
                            (/2,2/), error=err, outT=M1 )
    ELSE ifAllocatedM1
        WRITE(*,*) 'Warning-Array M1 not allocated in getgaugeorder4.'
    END IF ifAllocatedM1 
    
    DEALLOCATE(T1, STAT=istat)
    DEALLOCATE(T2, STAT=istat)

    ! Reduce the potential round off error
    M1 = ( M1 + TRANSPOSE(M1) ) / 2.0_DBL
    !WRITE(*,*) 'The matrix M1 is: ', M1

    ! SVD decompose M1 to get U2
    ALLOCATE( U2( SIZE(M1,1), SIZE(M1,1) ), STAT = istat )
    ALLOCATE(  V( SIZE(M1,1), SIZE(M1,1) ), STAT = istat )
    ALLOCATE( Spectra2( SIZE(M1,1) ), STAT = istat )
    CALL SVD( M1, U2, Spectra2, V, SIZE(M1,1), SIZE(M1,1) )
    
    ! Notice that the contents of M1 are destroyed.
    !WRITE(*,*) 'The matrix M1 after SVD is: ', M1
    
    DEALLOCATE( V, STAT = istat )
    DEALLOCATE( M1, STAT = istat )

    !WRITE(*,*) 'The matrix U2 is: ', U2
    !WRITE(*,*) 'The singular values Spectra2: ', Spectra2

    !                                                 
    ! Calculate U1 ( left side)           
    !                                                 

    
    ! T1(1,4,1',4')
    ALLOCATE( T1(SIZE(Ta,2), SIZE(Ta,4), SIZE(Ta,2), SIZE(Ta,4)), STAT=istat)
    IF( ALLOCATED(T1) ) THEN
        CALL TensorContract( Ta, Ta, 4, (/2,3/), 4, (/2,3/), error=err,outT=T1)
    ELSE
        WRITE(*,*) 'Warning-Array T1 not allocated in getgaugeorder4.'
    END IF 
    
    ! T2(5,4,5',4')
    ALLOCATE( T2(SIZE(Tb,2), SIZE(Tb,3), SIZE(Tb,2), SIZE(Tb,3)), STAT=istat)
    IF( ALLOCATED(T2) ) THEN
        CALL TensorContract( Tb, Tb, 4, (/2,4/), 4, (/2,4/), error=err,outT=T2 )
    ELSE
        WRITE(*,*) 'Warning-Array T2 not allocated in getgaugeorder4.'
    END IF 

    ! M1(15,1'5')
    ALLOCATE( M1(SIZE(T1,1)*SIZE(T2,1), SIZE(T1,3)*SIZE(T2,3)), STAT=istat )
    IF( ALLOCATED(M1) ) THEN
        CALL TensorContract( T1, T2, 4, (/2,4/), 4, (/2,4/), (/1,3,2,4/), &
                            (/2,2/), error=err, outT=M1 )
    ELSE 
        WRITE(*,*) 'Warning-Array M1 not allocated in getgaugeorder4.'
    END IF
    
    DEALLOCATE(T1, STAT=istat)
    DEALLOCATE(T2, STAT=istat)

    ! Reduce the potential round off error
    M1 = ( M1 + TRANSPOSE(M1) ) / 2.0_DBL
    !WRITE(*,*) 'The matrix M1 is: ', M1

    ! SVD decompose M1 to get U2
    ALLOCATE( U1( SIZE(M1,1), SIZE(M1,1) ), STAT = istat )
    ALLOCATE(  V( SIZE(M1,1), SIZE(M1,1) ), STAT = istat )
    ALLOCATE( Spectra1( SIZE(M1,1) ), STAT = istat )
    CALL SVD( M1, U1, Spectra1, V, SIZE(M1,1), SIZE(M1,1) )
    
    ! Notice that the contents of M1 are destroyed.
    !WRITE(*,*) 'The matrix M1 after SVD is: ', M1
    
    DEALLOCATE( V, STAT = istat )
    DEALLOCATE( M1, STAT = istat )

    !WRITE(*,*) 'The matrix U1 is: ', U1
    !WRITE(*,*) 'The singular values Spectra1: ', Spectra1

    ! Calculate truncation error
    Tr1 = SUM( Spectra1( Dbond+1: ) ) / SUM( Spectra1 ) 
    Tr2 = SUM( Spectra2( Dbond+1: ) ) / SUM( Spectra2 ) 
    TrunError = MIN( Tr1, Tr2 )

    DEALLOCATE( Spectra1, STAT = istat )
    DEALLOCATE( Spectra2, STAT = istat )

    ! Determine the Gauge
    IF ( Tr1 <= Tr2 ) THEN
        GAUGE = RESHAPE( U1(:, 1:Dbond), (/dim1, dim2, Dbond/) )
    ELSE
        GAUGE = RESHAPE( U2(:, 1:Dbond), (/dim1, dim2, Dbond/) )
    END IF

    DEALLOCATE( U1, STAT = istat )
    DEALLOCATE( U2, STAT = istat )

    END SUBROUTINE getgaugeorder4


    SUBROUTINE SVD(A,U,S,V,M,N)
    !
    ! Program computes the matrix singular value decomposition. 
    ! Using Lapack library.
    !
    ! Programmed by sukhbinder Singh
    ! 14th January 2011
    ! Modified by Yuzhi Liu on 04/26/2013.
    !      
    USE set_precisions    
    IMPLICIT NONE


    !   Data dictionary: declare calling parameter types & definitions
    REAL(kind=DBL), INTENT(INOUT), DIMENSION(M,N) :: A
    REAL(kind=DBL), INTENT(OUT), DIMENSION(M,M) :: U
    REAL(kind=DBL), INTENT(OUT), DIMENSION(N,N) :: V
    REAL(kind=DBL), INTENT(OUT), DIMENSION(MIN(M,N))   :: S
 
    !   Data dictionary: declare local variable types & definitions
    REAL(kind=DBL), DIMENSION(N,N) :: VT
    REAL(kind=DBL), ALLOCATABLE, DIMENSION(:) :: WORK
    INTEGER :: LDA, LDU, M, N, LWORK, LDVT, INFO
    CHARACTER :: JOBU, JOBVT
    INTEGER :: istat
    
    JOBU = 'A'
    JOBVT = 'A'
    LDA = M
    LDU = M
    LDVT = N
    
    LWORK = MAX( 1, 3*MIN(M,N) + MAX(M,N), 5*MIN(M,N) )

    ALLOCATE(WORK(LWORK), STAT=istat)

    ifAllocatedWORK: IF ( ALLOCATED(WORK) ) THEN
        CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, &
                    WORK, LWORK, INFO )
        ifInfo: IF ( INFO == 0 ) THEN
            !WRITE(*,*) 'SVD performed successfully!'
        ELSE ifInfo
            WRITE(*,*)  'SVD failed. Error INFO is : ', INFO
        END IF ifInfo
    ELSE ifAllocatedWORK
        WRITE(*,*) 'Warning-Array WORK not allocated in SVD.'
    END IF ifAllocatedWORK
    
    DEALLOCATE(WORK, STAT = istat) 
    V = TRANSPOSE(VT)
    

    END SUBROUTINE SVD
    


END MODULE get_gauge_order
