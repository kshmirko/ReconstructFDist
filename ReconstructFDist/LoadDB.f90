module LoadDB
implicit none

private
public load_ds
interface load_ds
    MODULE PROCEDURE load_1d_ds,load_2d_ds,load_3d_ds,load_4d_ds
end interface load_ds

contains

subroutine load_1d_ds(filename, dataset, x, rank, dims)
  use h5lt  
  IMPLICIT NONE

  CHARACTER(LEN=*), intent(IN) :: filename
  CHARACTER(LEN=*) , intent(in) :: dataset
  
  INTEGER(HID_T)  :: file, filetype, dset ! Handles
  INTEGER :: hdferr, rank, totsize,i
  INTEGER(HSIZE_T),allocatable   :: dims(:)
  
  real*4, DIMENSION(:), allocatable, intent(inout) :: X
  integer         :: type_class   ! type class
  integer(SIZE_T) :: type_size    ! type size
  
  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  ! 

  ! Now we begin the read section of this example. 
  !
  ! Open file, dataset, and attribute.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  if(hdferr/=0) STOP "Error in h5fopen_f"
  
  call h5ltget_dataset_ndims_f(file, dataset, rank, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_ndims_f"
  
  if(rank/=1) stop "Rank should be equal to 1"
  allocate(dims(rank))
  if(.not. allocated(dims)) STOP "Unable to allocate DIMS"
  
  call h5ltget_dataset_info_f(file, dataset, dims, type_class, type_size, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_info_f"
  
  ALLOCATE( X( dims(1) ) )
  if(.not. allocated(X)) STOP "Unable to allocate X"

  CALL h5ltread_dataset_f(file, dataset, H5T_IEEE_F32LE, X, dims, hdferr)
  if(hdferr/=0) STOP "Error in h5ltread_dataset_f"

  CALL H5Fclose_f (file, hdferr)
  if(hdferr/=0) STOP "Error in h5fclose_f"

  

end subroutine

subroutine load_2d_ds(filename, dataset, x, rank, dims)
  use h5lt  
  IMPLICIT NONE

  CHARACTER(LEN=*), intent(IN) :: filename
  CHARACTER(LEN=*) , intent(in) :: dataset
  
  INTEGER(HID_T)  :: file, filetype, dset ! Handles
  INTEGER :: hdferr, rank, totsize,i
  INTEGER(HSIZE_T),allocatable   :: dims(:)
  
  real*4, DIMENSION(:,:), allocatable, intent(inout) :: X
  integer         :: type_class   ! type class
  integer(SIZE_T) :: type_size    ! type size
  
  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  !

  ! Now we begin the read section of this example. 
  !
  ! Open file, dataset, and attribute.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  if(hdferr/=0) STOP "Error in h5fopen_f"
  
  call h5ltget_dataset_ndims_f(file, dataset, rank, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_ndims_f"
  
  if(rank/=2) stop "Rank should be equal to 2"
  allocate(dims(rank))
  if(.not. allocated(dims)) STOP "Unable to allocate DIMS"
  
  call h5ltget_dataset_info_f(file, dataset, dims, type_class, type_size, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_info_f"
  
  ALLOCATE( X( dims(1), dims(2) ) )
  if(.not. allocated(X)) STOP "Unable to allocate X"

  CALL h5ltread_dataset_f(file, dataset, H5T_IEEE_F32LE, X, dims, hdferr)
  if(hdferr/=0) STOP "Error in h5ltread_dataset_f"

  CALL H5Fclose_f (file, hdferr)
  if(hdferr/=0) STOP "Error in h5fclose_f"

end subroutine

subroutine load_3d_ds(filename, dataset, x, rank, dims)
  use h5lt  
  IMPLICIT NONE

  CHARACTER(LEN=*), intent(IN) :: filename
  CHARACTER(LEN=*) , intent(in) :: dataset
  
  INTEGER(HID_T)  :: file, filetype, dset ! Handles
  INTEGER :: hdferr, rank, totsize,i
  INTEGER(HSIZE_T),allocatable   :: dims(:)
  
  real*4, DIMENSION(:,:,:), allocatable, intent(inout) :: X
  integer         :: type_class   ! type class
  integer(SIZE_T) :: type_size    ! type size
  
  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  !

  ! Now we begin the read section of this example. 
  !
  ! Open file, dataset, and attribute.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  if(hdferr/=0) STOP "Error in h5fopen_f"
  
  call h5ltget_dataset_ndims_f(file, dataset, rank, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_ndims_f"
  
  if(rank/=3) stop "Rank should be equal to 3"
  allocate(dims(rank))
  if(.not. allocated(dims)) STOP "Unable to allocate DIMS"
  
  call h5ltget_dataset_info_f(file, dataset, dims, type_class, type_size, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_info_f"
  
  ALLOCATE( X( dims(1), dims(2), dims(3) ) )
  if(.not. allocated(X)) STOP "Unable to allocate X"

  CALL h5ltread_dataset_f(file, dataset, H5T_IEEE_F32LE, X, dims, hdferr)
  if(hdferr/=0) STOP "Error in h5ltread_dataset_f"

  CALL H5Fclose_f (file, hdferr)
  if(hdferr/=0) STOP "Error in h5fclose_f"

end subroutine

subroutine load_4d_ds(filename, dataset, x, rank, dims)
  use h5lt  
  IMPLICIT NONE

  CHARACTER(LEN=*), intent(IN) :: filename
  CHARACTER(LEN=*) , intent(in) :: dataset
  
  INTEGER(HID_T)  :: file, filetype, dset ! Handles
  INTEGER :: hdferr, rank, totsize,i
  INTEGER(HSIZE_T),allocatable   :: dims(:)
  
  real*4, DIMENSION(:,:,:,:), allocatable, intent(inout) :: X
  integer         :: type_class   ! type class
  integer(SIZE_T) :: type_size    ! type size
  
  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  !

  ! Now we begin the read section of this example. 
  !
  ! Open file, dataset, and attribute.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  if(hdferr/=0) STOP "Error in h5fopen_f"
  
  call h5ltget_dataset_ndims_f(file, dataset, rank, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_ndims_f"
  
  if(rank/=4) stop "Rank should be equal to 4"
  allocate(dims(rank))
  if(.not. allocated(dims)) STOP "Unable to allocate DIMS"
  
  call h5ltget_dataset_info_f(file, dataset, dims, type_class, type_size, hdferr)
  if(hdferr/=0) STOP "Error in h5ltget_dataset_info_f"
  
  ALLOCATE( X( dims(1), dims(2), dims(3), dims(4) ) )
  if(.not. allocated(X)) STOP "Unable to allocate X"

  CALL h5ltread_dataset_f(file, dataset, H5T_IEEE_F32LE, X, dims, hdferr)
  if(hdferr/=0) STOP "Error in h5ltread_dataset_f"

  CALL H5Fclose_f (file, hdferr)
  if(hdferr/=0) STOP "Error in h5fclose_f"

end subroutine


end module LoadDB