MODULE IOMODULE
IMPLICIT NONE

CONTAINS    

    SUBROUTINE PREPAREH5FILE(FNAME, RRI, NRRI, IRI, NIRI, R1, NR1, R2, NR2)
        INTEGER, INTENT(IN) ::  NRRI, NIRI, NR1, NR2
        REAL*8, INTENT(IN)  ::  RRI(NRRI), IRI(NRRI), R1(NR1), R2(NR2)
        CHARACTER*80, INTENT(IN)    :: FNAME
    END SUBROUTINE PREPAREH5FILE

    subroutine createTEST()
        USE hdf5
        CHARACTER(LEN=15), PARAMETER :: filename = "h5ex_d_hyper.h5"
        CHARACTER(LEN=3) , PARAMETER :: dataset  = "DS1"
        INTEGER          , PARAMETER :: dim0     = 6
        INTEGER          , PARAMETER :: dim1     = 8
        INTEGER          , PARAMETER :: RANK = 1

        INTEGER(HID_T)  :: file, space, dset ! Handles
        INTEGER         :: hdferr
        INTEGER(HSIZE_T), DIMENSION(1:RANK) :: dims = (/dim0/)
        INTEGER(HSIZE_T), DIMENSION(1:RANK)   :: start, stride, count, block

        INTEGER, DIMENSION(1:dim0) :: wdata, & ! Write buffer 
                                        rdata    ! Read buffer
        INTEGER :: i, j
        wdata = 1
        !
        ! Initialize FORTRAN interface.
        !
        CALL h5open_f(hdferr)
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)

        CALL h5screate_simple_f(RANK, dims, space, hdferr)
        CALL h5dcreate_f(file, dataset, H5T_NATIVE_INTEGER, space, dset, hdferr)
        
        CALL h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr)
        CALL h5dclose_f(dset , hdferr)
        CALL h5fclose_f(file, hdferr)
        CALL h5close_f(hdferr)
    end subroutine createTEST

END MODULE IOMODULE