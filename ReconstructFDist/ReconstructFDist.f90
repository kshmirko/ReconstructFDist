!  ReconstructFDist.f90 
!
!  FUNCTIONS:
!  ReconstructFDist - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: ReconstructFDist
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program ReconstructFDist
    use LoadDB, ONLY : load_ds
    implicit none

    INTEGER RANK
    INTEGER(8), ALLOCATABLE :: DIMS(:)
    REAL(4), ALLOCATABLE :: X(:)
    CHARACTER(len=80)    ::  FILENAME = "mieDB_new.h5"
    CHARACTER(Len=10)    ::  dataset = 'mR'
    
    ! read data from hdf file
    CALL LOAD_DS(filename, dataset, x, rank, dims)
    
    write(*,100) dataset, rank
    write(*,101) dims

    if (.not. allocated(X)) STOP "X not allocated"
    print *,X
    deallocate(dims, X)

100 format('Rank(',A,') = ', I2)
101 format('Dimension sizes are: ', <rank>I7)
    end program ReconstructFDist

