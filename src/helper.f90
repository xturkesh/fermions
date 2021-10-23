module helper
use, intrinsic :: iso_fortran_env, dp=> real64
implicit none 

complex(dp), parameter :: ii=(0._dp,1._dp)
complex(dp), parameter :: zero=(0._dp,0._dp)
complex(dp), parameter :: one=(1._dp,0._dp)

contains

subroutine eigvalsh(n,mat,eigs)
integer :: n, lda, lwork,info
complex(dp) :: mat(1:n,1:n), w(1:2*n*n)
real(dp) :: rw(1:3*n-2)
real(dp) :: eigs(1:n)
lwork =2*n*n
call zheev('N','U',n,mat,n,eigs,w,lwork,rw,info)
end subroutine eigvalsh

subroutine  qr(m,n,mat)
integer :: m,n, info, lwork
complex(dp) :: mat(1:m,1:n), tau(1:n), work(n*(3+n/2))
real(dp) :: rwork(3*n-2)

lwork = n*(3+n/2)
call zgeqrf(m,n,mat,m,tau,work,lwork,info)
call zungqr(m,n,n,mat,m,tau,work,lwork,info)
end subroutine qr

subroutine expm(t,A,expA)
implicit none
complex(dp), dimension(:,:) :: A, expA
real(dp) :: t 

integer :: ideg,m,ldh,lwsp,iexph,ns,iflag 
integer, allocatable, dimension(:) :: ipiv
complex(dp), allocatable, dimension(:) :: wsp

ideg = 6
m = size(A(1,:))
ldh = m
lwsp = 4*m*m+ideg+1
allocate(ipiv(m),wsp(lwsp))
call zgpadm(ideg,m,t,A,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
expA = reshape(wsp(iexph:iexph+m*m),(/m,m/))
deallocate(ipiv,wsp)
end subroutine expm

subroutine init_random_seed()
use iso_fortran_env, only: int64
implicit none
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid,getpid
integer(int64) :: t
            
call random_seed(size = n)
allocate(seed(n))
! First try if the OS provides a random number generator
open(newunit=un, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
    read(un) seed
    close(un)
else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
    call system_clock(t)
    if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                        + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                        + dt(3) * 24_int64 * 60 * 60 * 1000 &
                        + dt(5) * 60 * 60 * 1000 &
                        + dt(6) * 60 * 1000 + dt(7) * 1000 &
                        + dt(8)
    end if
    pid = getpid()
    t = ieor(t, int(pid, kind(t)))
    do i = 1, n
        seed(i) = lcg(t)
    end do
end if
call random_seed(put=seed)
contains
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
    function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
        s = 104729
    else
        s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
end subroutine init_random_seed

end module helper