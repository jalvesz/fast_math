module test_fast_math
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fast_math
    implicit none

    logical :: verbose = .true. ! change me to .true. if you want to see the results
contains

!> Collect all exported unit tests
subroutine collect_suite(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest('fast_sum', test_fast_sum) , &
        new_unittest('fast_dotp', test_fast_dotproduct) , &
        new_unittest('fast_trig', test_fast_trigonometry) , &
        new_unittest('fast_hyper', test_fast_hyperbolic ) &
    ]
end subroutine

subroutine test_fast_sum(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e6, ncalc = 4, niter = 20
    integer :: iter, i
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    1 format(a10,': <time> = ',f9.4,' ns/eval, speed-up=',f5.2,'X, rel. error=',es16.4)
    !====================================================================================
    call random_seed()
    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = sum(real(x, kind=wp*2)); call cpu_time(times(1))
            xsum(2) = sum(x)       ; call cpu_time(times(2))
            xsum(3) = fsum_pair(x) ; call cpu_time(times(3))
            xsum(4) = fsum(x)      ; call cpu_time(times(4))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," SUM on single precision values "
        print 1, "sum_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fsum_pair" , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        print 1, "fsum_chunk", 1e9*times_tot(4)/(niter*n) , times_tot(2)/times_tot(4), abs(err(4))/abs(meanval(4))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:),y(:)
        logical, allocatable :: xmask(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n),y(n),xmask(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        call random_number(y)
        y = (y - 0.5_wp)*2.0_wp
        xmask = abs(y) > 0.5_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = sum(real(x, kind=wp*2), mask=xmask); call cpu_time(times(1))
            xsum(2) = sum(x, mask=xmask)       ; call cpu_time(times(2))
            xsum(3) = fsum_pair(x, mask=xmask) ; call cpu_time(times(3))
            xsum(4) = fsum(x, mask=xmask)      ; call cpu_time(times(4))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," SUM on single precision values with a mask"
        print 1, "sum_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fsum_pair" , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        print 1, "fsum_chunk", 1e9*times_tot(4)/(niter*n) , times_tot(2)/times_tot(4), abs(err(4))/abs(meanval(4))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = sum(real(x, kind=wp*2)); call cpu_time(times(1))
            xsum(2) = sum(x)       ; call cpu_time(times(2))
            xsum(3) = fsum_pair(x) ; call cpu_time(times(3))
            xsum(4) = fsum(x)      ; call cpu_time(times(4))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," SUM on double precision values "
        print 1, "sum_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fsum_pair" , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        print 1, "fsum_chunk", 1e9*times_tot(4)/(niter*n) , times_tot(2)/times_tot(4), abs(err(4))/abs(meanval(4))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:),y(:)
        logical, allocatable :: xmask(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n),y(n),xmask(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        call random_number(y)
        y = (y - 0.5_wp)*2.0_wp
        xmask = abs(y) > 0.5_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = sum(real(x, kind=wp*2), mask=xmask); call cpu_time(times(1))
            xsum(2) = sum(x, mask=xmask)       ; call cpu_time(times(2))
            xsum(3) = fsum_pair(x, mask=xmask) ; call cpu_time(times(3))
            xsum(4) = fsum(x, mask=xmask)      ; call cpu_time(times(4))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," SUM on double precision values with a mask"
        print 1, "sum_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fsum_pair" , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        print 1, "fsum_chunk", 1e9*times_tot(4)/(niter*n) , times_tot(2)/times_tot(4), abs(err(4))/abs(meanval(4))
        end if 

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_dotproduct(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e6, ncalc = 3, niter = 50
    integer :: iter, i
    real(dp) :: times(0:ncalc), times_tot(ncalc)
    1 format(a10,': <time> = ',f9.4,' ns/eval, speed-up=',f5.2,'X, rel. error=',es16.4)
    !====================================================================================
    call random_seed()
    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)) ; call cpu_time(times(1))
            xsum(2) = dot_product(x,x) ; call cpu_time(times(2))
            xsum(3) = fprod(x,x) ; call cpu_time(times(3))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," dot product on single precision values "
        print 1, "dot_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fprod"     , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)) ; call cpu_time(times(1))
            xsum(2) = dot_product(x,x) ; call cpu_time(times(2))
            xsum(3) = fprod(x,x) ; call cpu_time(times(3))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," dot product on double precision values "
        print 1, "dot_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fprod"     , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block
    
    block
        integer, parameter :: wp=sp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)*real(x, kind=wp*2)) ; call cpu_time(times(1))
            xsum(2) = dot_product(x,x*x) ; call cpu_time(times(2))
            xsum(3) = fprod(x,x,x) ; call cpu_time(times(3))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," weigthed dot product on single precision values "
        print 1, "dot_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fprod"     , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=dp
        real(kind=wp), allocatable :: x(:)
        real(kind=wp) :: xsum(ncalc), meanval(ncalc), err(ncalc)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        allocate(x(n))
        call random_number(x)
        x = (x - 0.5_wp)*2.0_wp
        times_tot(:) = 0
        meanval(:) = 0
        err(:) = 0
        do iter=1,niter
            call cpu_time(times(0))
            xsum(1) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)*real(x, kind=wp*2)) ; call cpu_time(times(1))
            xsum(2) = dot_product(x,x*x) ; call cpu_time(times(2))
            xsum(3) = fprod(x,x,x) ; call cpu_time(times(3))
            times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            meanval(1:ncalc) = meanval(1:ncalc) + xsum(1:ncalc)
            err(2:ncalc) = err(2:ncalc) + abs( xsum(1) - xsum(2:ncalc) )
        end do
        meanval(1:ncalc) = meanval(1:ncalc) / niter
        err(1:ncalc) = err(1:ncalc) / niter 

        if(verbose)then
        print *,""
        print *,"================================================================"
        print *," weigthed dot product on double precision values "
        print 1, "dot_quad"  , 1e9*times_tot(1)/(niter*n) , times_tot(2)/times_tot(1), abs(err(1))/abs(meanval(1))
        print 1, "intrinsic" , 1e9*times_tot(2)/(niter*n) , times_tot(2)/times_tot(2), abs(err(2))/abs(meanval(2))
        print 1, "fprod"     , 1e9*times_tot(3)/(niter*n) , times_tot(2)/times_tot(3), abs(err(3))/abs(meanval(3))
        end if

        call check(error, abs(err(ncalc)) / abs(meanval(1)) < tolerance &
                        & .and. times_tot(ncalc) < times_tot(2) )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_trigonometry(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e6, ncalc = 2
    integer :: i
    real(dp) :: time(0:ncalc), err
    1 format(a10,': <time> = ',f9.4,' ns/eval, speed-up=',f5.2,'X, rel. error=',es16.4)
    !====================================================================================
    if(verbose)then
        print *,""
        print *,"================================================================"
        print *," Fast trigonometric"
    end if
    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        !> define a linspace between [-pi,pi]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*acos(-1.0_wp) , i = 1, n) ]

        call cpu_time(time(0))
        yref = sin(x); call cpu_time(time(1))
        y = fsin(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( sum( y - yref )**2 / n )

        if(verbose) print 1, "fsin r32" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        !> define a linspace between [-pi,pi]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*acos(-1.0_wp) , i = 1, n) ]

        call cpu_time(time(0))
        yref = sin(x); call cpu_time(time(1))
        y = fsin(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( sum( y - yref )**2 / n )

        if(verbose) print 1, "fsin r64" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        !> define a linspace between [-1,1]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ ((real(i,kind=wp) / n - 0.5_wp)*2 , i = 1, n) ]

        call cpu_time(time(0))
        yref = acos(x); call cpu_time(time(1))
        y = facos(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( sum( y - yref )**2 / n )

        if(verbose) print 1, "facos r32" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = epsilon(1._wp)*500
        !> define a linspace between [-1,1]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ ((real(i,kind=wp) / n - 0.5_wp)*2 , i = 1, n) ]

        call cpu_time(time(0))
        yref = acos(x); call cpu_time(time(1))
        y = facos(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( sum( y - yref )**2 / n )

        if(verbose) print 1, "facos r64" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_fast_hyperbolic(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e6, ncalc = 2
    integer :: i
    real(dp) :: time(0:ncalc), err
    1 format(a10,': <time> = ',f9.4,' ns/eval, speed-up=',f5.2,'X, rel. error=',es16.4)
    !====================================================================================
    if(verbose)then
        print *,""
        print *,"================================================================"
        print *," Fast hyperbolic"
    end if
    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = 1e-5_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        call cpu_time(time(0))
        yref = tanh(x); call cpu_time(time(1))
        y = ftanh(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( fsum(( y - yref )**2) ) / sqrt( fsum(( yref )**2) )

        if(verbose) print 1, "ftanh r64" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = 1e-5_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        call cpu_time(time(0))
        yref = tanh(x); call cpu_time(time(1))
        y = ftanh(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( fsum(( y - yref )**2) ) / sqrt( fsum(( yref )**2) )

        if(verbose) print 1, "ftanh r64" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block

    block
        integer, parameter :: wp=sp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = 1e-2_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        call cpu_time(time(0))
        yref = erf(x); call cpu_time(time(1))
        y = ferf(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( fsum(( y - yref )**2) ) / sqrt( fsum(( yref )**2) )

        if(verbose)print 1, "ferf r32" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
        if (allocated(error)) return
    end block
    block
        integer, parameter :: wp=dp
        real(wp), allocatable :: x(:) , y(:), yref(:)
        real(kind=wp) :: tolerance = 1e-2_wp
        !> define a linspace between [-3,3]
        allocate( x(n) , y(n), yref(n) )
        x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*3._wp , i = 1, n) ]

        call cpu_time(time(0))
        yref = erf(x); call cpu_time(time(1))
        y = ferf(x)  ; call cpu_time(time(2))

        time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
        err = sqrt( fsum(( y - yref )**2) ) / sqrt( fsum(( yref )**2) )

        if(verbose) print 1, "ferf r64" , 1e9*time(2)/n, time(1)/time(2), err

        call check(error, err < tolerance .and. time(2) < time(1) )
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
      new_testsuite("fsparse", collect_suite) &
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