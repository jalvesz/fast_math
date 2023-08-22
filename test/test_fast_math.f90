module test_fast_math
    use iso_fortran_env
    use fast_math
    implicit none
    contains

    subroutine test_fast_sum()
        integer, parameter :: n = 1e6, ncalc = 4, niter = 20
        integer :: iter, i
        real(real64) :: times(0:ncalc), times_tot(ncalc)
        character (len=*), parameter :: fmt_cr = "(a10,*(f22.12))", fmt_er = "(a10,*(es22.4))"
        !====================================================================================
        call random_seed()
        print *,""
        print *,"================================================================"
        print *," SUM on single precision values "
        block
            integer, parameter :: wp=real32
            real(kind=wp), allocatable :: x(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = sum(real(x, kind=wp*2)); call cpu_time(times(1))
                xsum(2,iter) = sum(x)       ; call cpu_time(times(2))
                xsum(3,iter) = fsum_pair(x) ; call cpu_time(times(3))
                xsum(4,iter) = fsum(x)      ; call cpu_time(times(4))
                times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "sum_quad", "intrinsic", "fsum_pair", "fsum_chunk"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block

        print *,""
        print *,"================================================================"
        print *," SUM on single precision values with a mask"
        block
            integer, parameter :: wp=real32
            real(kind=wp), allocatable :: x(:),y(:)
            logical, allocatable :: xmask(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n),y(n),xmask(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            call random_number(y)
            y = (y - 0.5_wp)*2.0_wp
            xmask = abs(y) > 0.5_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = sum(real(x, kind=wp*2), mask=xmask); call cpu_time(times(1))
                xsum(2,iter) = sum(x, mask=xmask)       ; call cpu_time(times(2))
                xsum(3,iter) = fsum_pair(x, mask=xmask) ; call cpu_time(times(3))
                xsum(4,iter) = fsum(x, mask=xmask)      ; call cpu_time(times(4))
                times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "sum_quad", "intrinsic", "fsum_pair", "fsum_chunk"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block

        print *,""
        print *,"================================================================"
        print *," SUM on double precision values "
        block
            integer, parameter :: wp=real64
            real(kind=wp), allocatable :: x(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = sum(real(x, kind=wp*2)); call cpu_time(times(1))
                xsum(2,iter) = sum(x)       ; call cpu_time(times(2))
                xsum(3,iter) = fsum_pair(x) ; call cpu_time(times(3))
                xsum(4,iter) = fsum(x)      ; call cpu_time(times(4))
                times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "sum_quad", "intrinsic", "fsum_pair", "fsum_chunk"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block

        print *,""
        print *,"================================================================"
        print *," SUM on double precision values with a mask"
        block
            integer, parameter :: wp=real64
            real(kind=wp), allocatable :: x(:),y(:)
            logical, allocatable :: xmask(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n),y(n),xmask(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            call random_number(y)
            y = (y - 0.5_wp)*2.0_wp
            xmask = abs(y) > 0.5_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = sum(real(x, kind=wp*2), mask=xmask); call cpu_time(times(1))
                xsum(2,iter) = sum(x, mask=xmask)       ; call cpu_time(times(2))
                xsum(3,iter) = fsum_pair(x, mask=xmask) ; call cpu_time(times(3))
                xsum(4,iter) = fsum(x, mask=xmask)      ; call cpu_time(times(4))
                times_tot(:) = times_tot(:) + times(1:ncalc) - times(0:ncalc-1)
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "sum_quad", "intrinsic", "fsum_pair", "fsum_chunk"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block
    end subroutine

    subroutine test_fast_dotproduct()
        integer, parameter :: n = 1e6, ncalc = 3, niter = 100
        integer :: iter, i
        real(real64) :: times(0:ncalc), times_tot(ncalc)
        character (len=*), parameter :: fmt_cr = "(a10,*(f22.12))", fmt_er = "(a10,*(es22.4))"
        !====================================================================================
        call random_seed()
        print *,""
        print *,"================================================================"
        print *," dot product on single precision values "
        block
            integer, parameter :: wp=real32
            real(kind=wp), allocatable :: x(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)) ; call cpu_time(times(1))
                xsum(2,iter) = dot_product(x,x) ; call cpu_time(times(2))
                xsum(3,iter) = fprod(x,x) ; call cpu_time(times(3))
                times_tot(:) = times_tot(:) + ( times(1:ncalc) - times(0:ncalc-1) ) 
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "quad" , "intrinsic", "fprod"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block

        print *,""
        print *,"================================================================"
        print *," dot product on double precision values "
        block
            integer, parameter :: wp=real64
            real(kind=wp), allocatable :: x(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)) ; call cpu_time(times(1))
                xsum(2,iter) = dot_product(x,x) ; call cpu_time(times(2))
                xsum(3,iter) = fprod(x,x) ; call cpu_time(times(3))
                times_tot(:) = times_tot(:) + ( times(1:ncalc) - times(0:ncalc-1) ) 
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "quad" , "intrinsic", "fprod"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block

        print *,""
        print *,"================================================================"
        print *," weigthed dot product on single precision values "
        block
            integer, parameter :: wp=real32
            real(kind=wp), allocatable :: x(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)*real(x, kind=wp*2)) ; call cpu_time(times(1))
                xsum(2,iter) = dot_product(x,x*x) ; call cpu_time(times(2))
                xsum(3,iter) = fprod(x,x,x) ; call cpu_time(times(3))
                times_tot(:) = times_tot(:) + ( times(1:ncalc) - times(0:ncalc-1) ) 
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "quad" , "intrinsic", "fprod"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block

        print *,""
        print *,"================================================================"
        print *," weigthed dot product on double precision values "
        block
            integer, parameter :: wp=real64
            real(kind=wp), allocatable :: x(:)
            real(kind=wp) :: xsum(ncalc,niter), meanval(ncalc), err(ncalc)
            allocate(x(n))
            call random_number(x)
            x = (x - 0.5_wp)*2.0_wp
            times_tot(:) = 0
            do iter=1,niter
                call cpu_time(times(0))
                xsum(1,iter) = dot_product(real(x, kind=wp*2),real(x, kind=wp*2)*real(x, kind=wp*2)) ; call cpu_time(times(1))
                xsum(2,iter) = dot_product(x,x*x) ; call cpu_time(times(2))
                xsum(3,iter) = fprod(x,x,x) ; call cpu_time(times(3))
                times_tot(:) = times_tot(:) + ( times(1:ncalc) - times(0:ncalc-1) ) 
            end do
            meanval(:) = sum( xsum , dim = 2 ) / niter
            do i = 1, ncalc
                err(i) = sum( abs(xsum(i,:) - xsum(1,:)) ) / niter
            end do
            print "(/,a10,*(a22))", "", "quad" , "intrinsic", "fprod"
            print fmt_cr,"Value: ", meanval(:)
            print fmt_er,"Error: ", err(:)
            print fmt_cr,"time : ", times_tot(1:ncalc) / niter
        end block
    end subroutine

    subroutine test_fast_trigonometry()
        integer, parameter :: n = 1e7
        real(real64) :: time(0:2), err
        integer :: i
        character (len=*), parameter :: fmt_cr = "(a10,*(f22.12))", fmt_er = "(a10,*(es22.4))"
        !====================================================================================
        print *,""
        print *,"================================================================"
        print *," Fast trigonometric"
        print "(/,a10,*(a22))", "", "Error", "time intrinsic", "time fast", "Speed up"
        block
            integer, parameter :: wp = real32
            real(wp), allocatable :: x(:) , y(:), yref(:)
            ! -- define a linspace between [-pi,pi]
            allocate( x(n) , y(n), yref(n) )
            x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*acos(-1.0_wp) , i = 1, n) ]

            call cpu_time(time(0))
            yref = sin(x); call cpu_time(time(1))
            y = fsin(x)  ; call cpu_time(time(2))

            time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
            err = sqrt( sum( y - yref )**2 / n )
            
            print fmt_er, "real32 sin", err, time(1), time(2), time(1)/time(2)
        end block
        block
            integer, parameter :: wp = real64
            real(wp), allocatable :: x(:) , y(:), yref(:)
            ! -- define a linspace between [-pi,pi]
            allocate( x(n) , y(n), yref(n) )
            x(:) = [ (2*(real(i,kind=wp) / n - 0.5_wp)*acos(-1.0_wp) , i = 1, n) ]

            call cpu_time(time(0))
            yref = sin(x); call cpu_time(time(1))
            y = fsin(x)  ; call cpu_time(time(2))

            time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
            err = sqrt( sum( y - yref )**2 / n )

            print fmt_er, "real64 sin", err, time(1), time(2), time(1)/time(2)
        end block

        block
            integer, parameter :: wp = real32
            real(wp), allocatable :: x(:) , y(:), yref(:)
            ! -- define a linspace between [-1,1]
            allocate( x(n) , y(n), yref(n) )
            x(:) = [ ((real(i,kind=wp) / n - 0.5_wp)*2 , i = 1, n) ]

            call cpu_time(time(0))
            yref = acos(x); call cpu_time(time(1))
            y = facos(x)  ; call cpu_time(time(2))

            time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
            err = sqrt( sum( y - yref )**2 / n )
            
            print fmt_er, "real32 acos", err, time(1), time(2), time(1)/time(2)
        end block
        block
            integer, parameter :: wp = real64
            real(wp), allocatable :: x(:) , y(:), yref(:)
            ! -- define a linspace between [-1,1]
            allocate( x(n) , y(n), yref(n) )
            x(:) = [ ((real(i,kind=wp) / n - 0.5_wp)*2 , i = 1, n) ]

            call cpu_time(time(0))
            yref = acos(x); call cpu_time(time(1))
            y = facos(x)  ; call cpu_time(time(2))

            time(2:1:-1) = time(2:1:-1) - time(1:0:-1)
            err = sqrt( sum( y - yref )**2 / n )

            print fmt_er, "real64 acos", err, time(1), time(2), time(1)/time(2)
        end block
    end subroutine

end module