module fast_dotp
    !! A faster and more accurate implementation of the dot_product intrinsic. 
    !! It uses the same principle as fsum_chunk but considering local multiplications that can be vectorized for faster summation.
    use iso_fortran_env
    implicit none
    private
    
    public :: fprod
    integer, parameter :: chunk64 = 32
    integer, parameter :: chunk32 = 64
    
    interface fprod
        module procedure fprod_r32
        module procedure fprod_r32_weighted
        module procedure fprod_r64
        module procedure fprod_r64_weighted
    end interface
    
    contains

    pure function fprod_r32(a,b) result(p)
        integer, parameter :: wp = real32
        integer, parameter :: chunk = chunk32
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
  
    pure function fprod_r64(a,b) result(p)
        integer, parameter :: wp = real64
        integer, parameter :: chunk = chunk64
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
    
    pure function fprod_r32_weighted(a,b,w) result(p)
        integer, parameter :: wp = real32
        integer, parameter :: chunk = chunk32
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(in) :: w(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)*w(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)*w(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
  
    pure function fprod_r64_weighted(a,b,w) result(p)
        integer, parameter :: wp = real64
        integer, parameter :: chunk = chunk64
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(in) :: w(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)*w(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)*w(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
  end module fast_dotp