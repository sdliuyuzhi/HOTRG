MODULE useful_functions
IMPLICIT NONE
PRIVATE
PUBLIC :: setdiff

CONTAINS
    PURE FUNCTION setdiff(A,B)
    !
    !   Purpose:
    !       c=setdiff(A,B) returns the values in A that are not in B. In set theory
    !   terms, c = A - B. This is a simplified version of the setdiff function used
    !   in MATLAB. A and B can be interger arrays only for now.
    !   B is always a subset of A.
    !       If c is empty, then we set c = 0. This may not be optimal but we use it
    !   anyway.
    !
    !   Record of revision:
    !   04/10/2013. Yuzhi Liu v1.1
    !
    IMPLICIT NONE
    
    !   Data dictionary: declare calling parameter types & definitions
    INTEGER, INTENT(IN), DIMENSION(:) :: A              ! Input array A
    INTEGER, INTENT(IN), DIMENSION(:) :: B              ! Input array B
    INTEGER, DIMENSION(abs( size(A) - size(B) )) :: setdiff    ! Output arrary setdiff

    !   Data dictionary: declare local variables
    INTEGER, DIMENSION(size(A)) :: C                    ! Intermediate array C
    INTEGER :: iB                                       ! Loop index: loop over B array
    
    C = A                                               ! Make a copy of array A
    DO iB = 1, size(B)
        C = pack( C, C /= B(iB) ) 
    END DO
    
    setdiff = C
    IF ( size(setdiff) == 0 ) THEN 
        setdiff = 0
    END IF
    
    END FUNCTION setdiff


    
END MODULE useful_functions
