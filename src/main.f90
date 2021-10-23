program main
use, intrinsic :: iso_fortran_env, dp=> real64
use free_fermions
use helper
implicit none 

integer :: lsys, nt, ti, ndis, ri, ell 
real(dp) :: eps, delta, tmax, tmin, ent, tt  
real(dp), allocatable, dimension(:) :: trange, jhop, hmag  
character(80) :: arg1, arg2, arg3, fname
type(fermions) :: state

call get_command_argument(1,arg1)
call get_command_argument(2,arg2)
call get_command_argument(3,arg3)
call get_command_argument(4,fname)

read(arg1,*) lsys
read(arg2,*) ell 
read(arg3,*) delta 

call init_random_seed

nt = 300
ndis = 1000
tmax = 1e12_dp
tmin = 1e-3_dp
eps = exp((log(tmax)-log(tmin))/real(nt,dp))

allocate(trange(nt),jhop(lsys),hmag(lsys))
do ti=1,nt 
    trange(ti) = tmin*(eps**real(ti,dp))
end do 

open(unit=88,file=fname,action='write',position='append')
call state%init(lsys,lsys/2,.true.)
do ri = 1,ndis 
    call random_number(jhop)
    jhop = jhop**delta 
    hmag = 0._dp 

    call state%build_ham(jhop,hmag)
    call state%init_state()
    do ti = 1,nt 
        tt = trange(ti)
        call state%exp_ham(tt)
        call state%entropy(1,ell,ent)
        write(88,'(I6,F32.15,F20.8)') ri, tt, ent
    end do 
end do
call state%dest
close(88)

deallocate(trange,jhop,hmag)
end program 