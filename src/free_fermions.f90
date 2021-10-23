module free_fermions
use, intrinsic :: iso_fortran_env, dp=> real64
use helper

implicit none 

type fermions
    integer :: sites, num 
    complex(dp), allocatable, dimension(:,:) :: umat, cmat, umat0
    complex(dp), allocatable, dimension(:,:) :: hmat, exph
    complex(dp), allocatable, dimension(:) :: hj, jj 
    logical :: bc 
contains 
    procedure :: init 
    procedure :: init_state
    procedure :: build_ham 
    procedure :: exp_ham 
    procedure :: normalize 
    procedure :: entropy 
    procedure :: negativity 
    procedure :: correlation 
    procedure :: dest  
end type fermions


contains 

subroutine init(state,lsys,num,bc)
class(fermions) :: state 
integer :: lsys, num 
logical, optional :: bc 
if (present(bc)) then 
    state%bc = bc 
else 
    state%bc = .false.
end if 
state%sites = lsys 
state%num = num 
if (allocated(state%cmat)) deallocate(state%cmat)
allocate(state%cmat(1:lsys,1:lsys))
if (allocated(state%umat)) deallocate(state%umat)
allocate(state%umat(1:lsys,1:num))
if (allocated(state%umat0)) deallocate(state%umat0)
allocate(state%umat0(1:lsys,1:num))
if (allocated(state%hmat)) deallocate(state%hmat)
allocate(state%hmat(1:lsys,1:lsys))
if (allocated(state%exph)) deallocate(state%exph)
allocate(state%exph(1:lsys,1:lsys))
if (allocated(state%jj)) deallocate(state%jj)
allocate(state%jj(1:lsys))
if (allocated(state%hj)) deallocate(state%hj)
allocate(state%hj(1:lsys))
end subroutine 

subroutine dest(state)
class(fermions) :: state
if (allocated(state%cmat)) deallocate(state%cmat)
if (allocated(state%umat)) deallocate(state%umat)
if (allocated(state%umat0)) deallocate(state%umat0)
if (allocated(state%hmat)) deallocate(state%hmat)
if (allocated(state%exph)) deallocate(state%exph)
if (allocated(state%jj)) deallocate(state%jj)
if (allocated(state%hj)) deallocate(state%hj)
end subroutine 

subroutine build_ham(state,jcoup,hmag)
class(fermions) :: state
real(dp), dimension(:) :: jcoup, hmag
integer :: i
associate(ls=> state%sites, hmat => state%hmat, jj=>state%jj, hj=> state%hj)
hmat = zero
do i =1,ls
    jj(i) = jcoup(i)
    hj(i) = hmag(i)
    hmat(i,i+1) = -ii*jcoup(i)
    hmat(i+1,i) = -ii*jcoup(i)
    hmat(i,i) = -ii*hmag(i)
end do
hmat(ls,ls) = -ii*hmag(ls)
if (state%bc) then
    hmat(ls,1) = -ii*jcoup(ls)
    hmat(1,ls) = -ii*jcoup(ls)
end if
end associate
end subroutine 

subroutine init_state(state)
class(fermions) :: state
integer :: i 
associate(ls => state%sites, ns=> state%num, u => state%umat0)
u = zero
do i=1,ns
    if (mod(i,2)==0) then
        u(2*i,i) = one
    end if
end do
end associate
end subroutine 

subroutine exp_ham(state,  t)
class(fermions) :: state 
real(dp) :: t
associate(hmat=>state%hmat,exph=>state%exph,uin=>state%umat0,uout=>state%umat)
call expm(t,hmat,exph)
uout = matmul(exph,uin)
end associate
end subroutine 

subroutine normalize(state)
class(fermions) :: state
associate(ls=>state%sites,num=>state%num,u=>state%umat)
call qr(ls,num,u)
end associate
end subroutine 

subroutine correlation(state) 
class(fermions) :: state 
associate(u => state%umat, corr => state%cmat)
corr = matmul(u,conjg(transpose(u)))
end associate
end subroutine

subroutine entropy(state,a1,a2,ent)
class(fermions) :: state 
real(dp) :: ent, ee  
integer :: a1, a2, la, i 
real(dp), allocatable, dimension(:) :: eigs
complex(dp), allocatable, dimension(:,:) :: wmat 
la = a2-a1+1 
allocate(eigs(la),wmat(la,la))
call state%correlation 
wmat = state%cmat(a1:a2,a1:a2)
call eigvalsh(la,wmat,eigs)
ent = 0._dp 
do i=1,la 
    ee = eigs(i)
    if (ee>0._dp .and. ee<1._dp) then 
        ent = ent - (1._dp - ee)*log(1._dp-ee) - ee*log(ee)
    end if
end do
deallocate(eigs,wmat)
end subroutine 

subroutine negativity(state)
class(fermions) :: state

end subroutine negativity


end module free_fermions