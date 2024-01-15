module test_fast_math
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fast_math
    implicit none

    logical :: verbose = .true. ! change me to .true. if you want to see the results
    interface scramble
        module procedure scramble_sp
        module procedure scramble_dp
    end interface
    interface scramble_l
        module procedure scramble_spl
        module procedure scramble_dpl
    end interface

    character(len=*), parameter :: fmt1 = "('| ',a12,' | <time> [ns/eval] | Speed-Up | relative error' ,'  |')"
    character(len=*), parameter :: fmt2 = "('|--------------|------------------|----------|-----------------| ')"
    character(len=*), parameter :: fmt3 = "('| ',a12,' |        ', f9.4,' | ',f8.2,' |',es16.4,' | ')"

contains

subroutine scramble_sp( x )
    integer, parameter :: wp = sp
    real(wp), intent(inout) :: x(:)
    real(wp) :: u, temp
    integer :: i, j, m

    m = size(x)
    do i = 1, m
        call random_number(u)
        j = 1 + FLOOR(m*u)
        temp = x(j)
        x(j) = x(i)
        x(i) = temp
    end do
end subroutine

subroutine scramble_dp( x )
    integer, parameter :: wp = dp
    real(wp), intent(inout) :: x(:)
    real(wp) :: u, temp
    integer :: i, j, m

    m = size(x)
    do i = 1, m
        call random_number(u)
        j = 1 + FLOOR(m*u)
        temp = x(j)
        x(j) = x(i)
        x(i) = temp
    end do
end subroutine

subroutine scramble_spl( x , l )
    integer, parameter :: wp = sp
    real(wp), intent(inout) :: x(:)
    logical , intent(inout) :: l(:)
    real(wp) :: u, temp
    logical :: ltemp
    integer :: i, j, m

    m = size(x)
    do i = 1, m
        call random_number(u)
        j = 1 + FLOOR(m*u)
        temp = x(j); ltemp = l(j)
        x(j) = x(i); l(j) = l(i)
        x(i) = temp; l(i) = ltemp
    end do
end subroutine

subroutine scramble_dpl( x , l )
    integer, parameter :: wp = dp
    real(wp), intent(inout) :: x(:)
    logical , intent(inout) :: l(:)
    real(wp) :: u, temp
    logical :: ltemp
    integer :: i, j, k, m

    m = size(x)
    do i = 1, m
        call random_number(u)
        j = 1 + FLOOR(m*u)
        temp = x(j); ltemp = l(j)
        x(j) = x(i); l(j) = l(i)
        x(i) = temp; l(i) = ltemp
    end do
end subroutine

real(dp) function timer() result(t)
    integer :: values(8)
    call date_and_time(VALUES=values)
    t = 60*values(6) + values(7) + values(8)*1.d-3
end function

!> Collect all exported unit tests
subroutine collect_suite(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest('fast_sum', test_fast_sum) , &
        new_unittest('fast_dotp', test_fast_dotproduct) , &
        new_unittest('fast_trig', test_fast_trigonometry) , &
        new_unittest('fast_hyper', test_fast_hyperbolic ) , & 
        new_unittest('fast_rsqrt', test_fast_rsqrt ) & !! The Quake III rsqrt implementation here is not realy much faster that what compilers can do.
    ]
end subroutine

subroutine test_fast_sum(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e5, ncalc = 3, niter = 1000
    integer :: iter, i
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    !====================================================================================
    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc), tolerance = epsilon(1._wp)*100

        allocate(x(n) , source = [(8*atan(1._wp)*(real(i,kind=wp)-0.5_wp)/real(n,kind=wp)**2,i=1,n)] )
        
        times_tot(:) = 0
        meanval(:) = 0._wp
        err(:) = 0._wp
        do iter=1,niter
            call scramble(x)
            times(0) = timer()
            xsum(1) = sum(x)       ; times(1) = timer()
            xsum(2) = fsum_kahan(x); times(2) = timer()
            xsum(3) = fsum(x)      ; times(3) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(1:ncalc) = err(1:ncalc) + abs(1._wp-xsum(1:ncalc)/(4*atan(1._wp)))
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        write(*,fmt1) "sum r32"
        write(*,fmt2)
        write(*,fmt3) "intrinsic" , 1e9*times_tot(1)/(niter*n) , times_tot(1)/times_tot(1), err(1)
        write(*,fmt3) "    kahan" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err(2)
        write(*,fmt3) "    chunk" , 1e9*times_tot(3)/(niter*n) , times_tot(1)/times_tot(3), err(3)
        end if

        call check(error, all(err(:)<tolerance) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc), tolerance = epsilon(1._wp)*100

        allocate(x(n) , source = [(8*atan(1._wp)*(real(i,kind=wp)-0.5_wp)/real(n,kind=wp)**2,i=1,n)] )
        
        times_tot(:) = 0
        meanval(:) = 0._wp
        err(:) = 0._wp
        do iter=1,niter
            call scramble(x)
            times(0) = timer()
            xsum(1) = sum(x)       ; times(1) = timer()
            xsum(2) = fsum_kahan(x); times(2) = timer()
            xsum(3) = fsum(x)      ; times(3) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(1:ncalc) = err(1:ncalc) + abs(1._wp-xsum(1:ncalc)/(4*atan(1._wp)))
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        write(*,fmt1) "sum r64"
        write(*,fmt2)
        write(*,fmt3) "intrinsic" , 1e9*times_tot(1)/(niter*n) , times_tot(1)/times_tot(1), err(1)
        write(*,fmt3) "    kahan" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err(2)
        write(*,fmt3) "    chunk" , 1e9*times_tot(3)/(niter*n) , times_tot(1)/times_tot(3), err(3)
        end if

        call check(error, all(err(:)<tolerance) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:)
        logical, allocatable :: mask(:), nmask(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc), tolerance = epsilon(1._wp)*100

        allocate(x(n) , source = [(8*atan(1._wp)*(real(i,kind=wp)-0.5_wp)/real(n,kind=wp)**2,i=1,n)] )
        allocate(mask(n),source=.false.); mask(1:n:2) = .true.
        allocate(nmask(n))

        times_tot(:) = 0
        meanval(:) = 0._wp
        err(:) = 0._wp
        do iter=1,niter
            call scramble_l(x,mask); nmask(:) = .not.mask(:)
            times(0) = timer()
            xsum(1) = sum(x,mask)        + sum(x,nmask)       ; times(1) = timer()
            xsum(2) = fsum_kahan(x,mask) + fsum_kahan(x,nmask); times(2) = timer()
            xsum(3) = fsum(x,mask)       + fsum(x,nmask)      ; times(3) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)

            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(1:ncalc) = err(1:ncalc) + abs(1._wp-(xsum(1:ncalc)/(4*atan(1._wp))))
        end do
        times_tot(:) = times_tot(:) / 2
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        write(*,fmt1) "sum r32 mask"
        write(*,fmt2)
        write(*,fmt3) "intrinsic" , 1e9*times_tot(1)/(niter*n) , times_tot(1)/times_tot(1), err(1)
        write(*,fmt3) "    kahan" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err(2)
        write(*,fmt3) "    chunk" , 1e9*times_tot(3)/(niter*n) , times_tot(1)/times_tot(3), err(3)
        end if

        call check(error, all(err(:)<tolerance) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:)
        logical, allocatable :: mask(:), nmask(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc), tolerance = epsilon(1._wp)*100

        allocate(x(n) , source = [(8*atan(1._wp)*(real(i,kind=wp)-0.5_wp)/real(n,kind=wp)**2,i=1,n)] )
        allocate(mask(n),source=.false.); mask(1:n:2) = .true.
        allocate(nmask(n))

        times_tot(:) = 0
        meanval(:) = 0._wp
        err(:) = 0._wp
        do iter=1,niter
            call scramble_l(x,mask); nmask(:) = .not.mask(:)
            times(0) = timer()
            xsum(1) = sum(x,mask)        + sum(x,nmask)       ; times(1) = timer()
            xsum(2) = fsum_kahan(x,mask) + fsum_kahan(x,nmask); times(2) = timer()
            xsum(3) = fsum(x,mask)       + fsum(x,nmask)      ; times(3) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)

            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(1:ncalc) = err(1:ncalc) + abs(1._wp-(xsum(1:ncalc)/(4*atan(1._wp))))
        end do
        times_tot(:) = times_tot(:) / 2
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        write(*,fmt1) "sum r64 mask"
        write(*,fmt2)
        write(*,fmt3) "intrinsic" , 1e9*times_tot(1)/(niter*n) , times_tot(1)/times_tot(1), err(1)
        write(*,fmt3) "    kahan" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err(2)
        write(*,fmt3) "    chunk" , 1e9*times_tot(3)/(niter*n) , times_tot(1)/times_tot(3), err(3)
        end if

        call check(error, all(err(:)<tolerance) )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_dotproduct(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e5, ncalc = 3, niter = 1000
    integer :: iter, i
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    !====================================================================================
    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc), tolerance = epsilon(1._wp)*100

        allocate(x(n) , source = [(8*atan(1._wp)*(real(i,kind=wp)-0.5_wp)/real(n,kind=wp)**2,i=1,n)] )
        x(:) = sqrt( x(:) )

        times_tot(:) = 0
        meanval(:) = 0._wp
        err(:) = 0._wp
        do iter=1,niter
            call scramble(x)
            times(0) = timer()
            xsum(1) = dot_product(x,x) ; times(1) = timer()
            xsum(2) = fprod_kahan(x,x) ; times(2) = timer()
            xsum(3) = fprod(x,x)       ; times(3) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(1:ncalc) = err(1:ncalc) + abs(1._wp-xsum(1:ncalc)/(4*atan(1._wp)))
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        write(*,fmt1) "dot r32"
        write(*,fmt2)
        write(*,fmt3) "intrinsic" , 1e9*times_tot(1)/(niter*n) , times_tot(1)/times_tot(1), err(1)
        write(*,fmt3) "    kahan", 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err(2)
        write(*,fmt3) "    chunk", 1e9*times_tot(3)/(niter*n) , times_tot(1)/times_tot(3), err(3)
        end if

        call check(error, all(err(:)<tolerance) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc), tolerance = epsilon(1._wp)*100

        allocate(x(n) , source = [(8*atan(1._wp)*(real(i,kind=wp)-0.5_wp)/real(n,kind=wp)**2,i=1,n)] )
        x(:) = sqrt( x(:) )

        times_tot(:) = 0
        meanval(:) = 0._wp
        err(:) = 0._wp
        do iter=1,niter
            call scramble(x)
            times(0) = timer()
            xsum(1) = dot_product(x,x) ; times(1) = timer()
            xsum(2) = fprod_kahan(x,x) ; times(2) = timer()
            xsum(3) = fprod(x,x)       ; times(3) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(1:ncalc) = err(1:ncalc) + abs(1._wp-xsum(1:ncalc)/(4*atan(1._wp)))
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        write(*,fmt1) "dot r64"
        write(*,fmt2)
        write(*,fmt3) "intrinsic" , 1e9*times_tot(1)/(niter*n) , times_tot(1)/times_tot(1), err(1)
        write(*,fmt3) "    kahan", 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err(2)
        write(*,fmt3) "    chunk", 1e9*times_tot(3)/(niter*n) , times_tot(1)/times_tot(3), err(3)
        end if

        call check(error, all(err(:)<tolerance) )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_trigonometry(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 5e5, ncalc = 2, niter = 500
    integer :: i, iter
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    !====================================================================================
    if(verbose)then
        print *,""
        write(*,fmt1) "trigo"
        write(*,fmt2)
    end if
    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = epsilon(1._wp)*500
        !> define a linspace between [-pi,pi]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*acos(-1.0_wp) , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = sin(x); times(1) = timer()
            y = fsin(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "fsin r32" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = epsilon(1._wp)*500
        !> define a linspace between [-pi,pi]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*acos(-1.0_wp) , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = sin(x); times(1) = timer()
            y = fsin(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "fsin r64" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = epsilon(1._wp)*500
        !> define a linspace between [-1,1]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ ((real(i,kind=wp) / n - 0.5_wp)*2 , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = acos(x); times(1) = timer()
            y = facos(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "facos r32" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = epsilon(1._wp)*500
        !> define a linspace between [-1,1]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ ((real(i,kind=wp) / n - 0.5_wp)*2 , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = acos(x); times(1) = timer()
            y = facos(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "facos r64" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 2e-5_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = atan(x); times(1) = timer()
            y = fatan(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "fatan r32" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 2e-5_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = atan(x); times(1) = timer()
            y = fatan(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "fatan r64" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_hyperbolic(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 5e5, ncalc = 2, niter = 500
    integer :: i, iter
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    !====================================================================================
    if(verbose)then
        print *,""
        write(*,fmt1) "hyperb"
        write(*,fmt2)
    end if
    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 1e-5_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = tanh(x); times(1) = timer()
            y = ftanh(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "ftanh r32" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 1e-5_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = tanh(x); times(1) = timer()
            y = ftanh(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "ftanh r64" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 1e-2_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = erf(x); times(1) = timer()
            y = ferf(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "ferf r32" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 1e-2_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = erf(x); times(1) = timer()
            y = ferf(x)  ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / n )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "ferf r64" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_rsqrt(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 5e5, ncalc = 2, niter = 500
    integer :: iter, i
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    !====================================================================================
    if(verbose)then
        print *,""
        write(*,fmt1) "rsqrt"
        write(*,fmt2)
    end if
    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 1e-2_wp
        !> define a log space
        
        allocate( x(n) , y(n), yref(n) )
        call random_number(x)
        x = 10._wp**(-10*(1._wp-x) + 10*x)

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = 1._wp/sqrt(x); times(1) = timer()
            y = frsqrt(x)       ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / sum( yref )**2 )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "frsqrt r32" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: err, tolerance = 1e-2_wp
        !> define a log space
        allocate( x(n) , y(n), yref(n) )
        call random_number(x)
        x = 10._wp**(-200*(1._wp-x) + 200*x)

        times_tot(:) = 0
        err = 0._wp
        do iter=1,niter
            times(0) = timer()
            yref = 1._wp/sqrt(x); times(1) = timer()
            y = frsqrt(x)       ; times(2) = timer()
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            err = err + sqrt( sum( y - yref )**2 / sum( yref )**2 )
        end do
        err = err / niter

        if(verbose) write(*,fmt3) "frsqrt r64" , 1e9*times_tot(2)/(niter*n) , times_tot(1)/times_tot(2), err

        call check(error, err < tolerance )
        if (allocated(error)) return
    end block

end subroutine
    
end module test_fast_math

program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_fast_math, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
  
    stat = 0
  
    testsuites = [ &
      new_testsuite("fast_math", collect_suite) &
      ]
  
    do is = 1, size(testsuites)
      write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if
  
end program tester