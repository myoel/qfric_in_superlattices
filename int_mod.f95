! This module contains all dependencies concerning
! the integration. Here, the QUADPACK routine DQAG
! is use with an interface.
module intdef 
 contains
  function integ1(f,lim1,lim2)
    implicit none
    interface
    function f (x) result (f_result)
      real(8), intent (in) :: x
      real(8) :: f_result
    end function f
   end interface
    !----- Integration parameters ----!
   real(8) :: epsrel = 1D-4
   real(8) :: epsabs = -0D0 ! only relative error
   real(8), intent(in) :: lim1,lim2
   real(8) :: integ1
   integer :: key = 6
   integer limit
    parameter ( limit = 10000000)
   integer lenw
    parameter ( lenw = limit * 4 )
   !----- Integration variables -----!
   integer iwork(limit)
   real(8) :: work(lenw)
   real(8) :: res
   real(8) :: abserr
   integer ier
   integer last
   integer neval
   !---------------------------------!
  call dqag( f, lim1, lim2, epsabs, &
            epsrel, key, res,abserr,&
            neval, ier, limit, lenw,&
            last, iwork, work )
  integ1 = res
    if (ier /= 0) then
   print*,'ERROR IN I1'
   print*,'res=',res
!     stop
  endif
  end function

  function integ2(f,lim1,lim2)
    implicit none
    interface
    function f (x) result (f_result)
      real(8), intent (in) :: x
      real(8) :: f_result
    end function f
   end interface
    !----- Integration parameters ----!
   real(8) :: epsrel = 1D-3
   real(8) :: epsabs = -0D0 ! only relative error
   real(8), intent(in) :: lim1,lim2
   real(8) :: integ2
   integer :: key = 6
   integer limit
    parameter ( limit = 10000000)
   integer lenw
    parameter ( lenw = limit * 4 )
   !----- Integration variables -----!
   integer iwork(limit)
   real(8) :: work(lenw)
   real(8) :: res
   real(8) :: abserr
   integer ier
   integer last
   integer neval
   !---------------------------------!
  call dqag( f, lim1, lim2, epsabs, &
            epsrel, key, res,abserr,&
            neval, ier, limit, lenw,&
            last, iwork, work )
  integ2 = res
  if (ier /= 0) then
   print*,'ERROR IN I2'
  endif
  end function

    function integ3(f,lim1,lim2)
    implicit none
    interface
    function f (x) result (f_result)
      real(8), intent (in) :: x
      real(8) :: f_result
    end function f
   end interface
    !----- Integration parameters ----!
   real(8) :: epsrel = 1D-2
   real(8) :: epsabs = -0D0 ! only relative error
   real(8), intent(in) :: lim1,lim2
   real(8) :: integ3
   integer :: key = 6
   integer limit
    parameter ( limit = 10000000)
   integer lenw
    parameter ( lenw = limit * 4 )
   !----- Integration variables -----!
   integer iwork(limit)
   real(8) :: work(lenw)
   real(8) :: res
   real(8) :: abserr
   integer ier
   integer last
   integer neval
   !---------------------------------!
  call dqag( f, lim1, lim2, epsabs, &
            epsrel, key, res,abserr,&
            neval, ier, limit, lenw,&
            last, iwork, work )
  integ3 = res
    if (ier /= 0) then
   print*,'ERROR IN I3'
  endif
  
  end function

end module
