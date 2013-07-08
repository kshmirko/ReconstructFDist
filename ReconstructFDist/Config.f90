module config
implicit none

type RefractiveIndex_t
    real*4  ::  mR(3)
    real*4  ::  mI(3)
end type RefractiveIndex_t

type Radius_t
    real*4  ::  Rmin(3)
    real*4  ::  Rmax(3)
end type Radius_t

type SmoothingParams_t
    real*4 K, b
end type SmoothingParams_t

type Discrepancy_t
    real*4 min, max
end type Discrepancy_t

contains

end module config