#!/usr/bin/env python3
"""Compile the SciFort optional FITPACK adapter against a synthetic API stub.

This checks wrapper syntax only, NOT the upstream numerical code. For genuine
FITPACK tests, use check_optional_fitpack_local.py with the source ZIP.
"""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
MOCK = r'''
module fitpack
    use iso_fortran_env, only: real64
    implicit none
    integer, parameter :: FP_REAL=real64, FP_FLAG=kind(0)
    type :: fitpack_curve
        integer :: knots=0, order=3
        real(real64), allocatable :: t(:),c(:)
    contains
        procedure :: new_fit => curve_new_fit
        procedure :: mse => curve_mse
        procedure :: destroy => curve_destroy
    end type
    type, extends(fitpack_curve) :: fitpack_periodic_curve
    end type
    type :: fitpack_parametric_curve
        integer :: knots=0, order=3, idim=0
        real(real64), allocatable :: t(:),c(:),u(:)
    contains
        procedure :: new_fit => param_new_fit
        procedure :: mse => param_mse
        procedure :: destroy => param_destroy
    end type
    type, extends(fitpack_parametric_curve) :: fitpack_closed_curve
    end type
    type :: fitpack_surface
        integer :: m=0,order(2)=3,nest(2)=0,nmax=0,knots(2)=0
        integer :: iopt=0,lwrk=0,liwrk=0,lwrk2=0
        integer, allocatable :: iwrk(:)
        real(real64), allocatable :: t(:,:),c(:),wrk(:),wrk2(:)
    contains
        procedure :: new_fit => surface_new_fit
        procedure :: new_points => surface_new_points
        procedure :: fit => surface_fit
        procedure :: mse => surface_mse
        procedure :: destroy => surface_destroy
    end type
contains
    pure logical function FITPACK_SUCCESS(code)
        integer(FP_FLAG), intent(in) :: code
        FITPACK_SUCCESS = code <= 0
    end function
    integer(FP_FLAG) function curve_new_fit(self,x,y,w,smoothing,order) result(status)
        class(fitpack_curve), intent(inout) :: self
        real(real64), intent(in) :: x(:),y(:)
        real(real64), intent(in), optional :: w(:),smoothing
        integer, intent(in), optional :: order
        integer :: degree
        degree=3
        if(present(order)) degree=order
        self%order=degree
        self%knots=2*(degree+1)
        self%t=[spread(x(1),1,degree+1),spread(x(size(x)),1,degree+1)]
        self%c=spread(0.0_real64,1,self%knots)
        status=0
    end function
    real(real64) function curve_mse(self) result(value)
        class(fitpack_curve), intent(in) :: self
        value=0.0_real64
    end function
    subroutine curve_destroy(self)
        class(fitpack_curve), intent(inout) :: self
        if(allocated(self%t)) deallocate(self%t)
        if(allocated(self%c)) deallocate(self%c)
    end subroutine
    integer(FP_FLAG) function param_new_fit(self,x,u,w,smoothing,order) result(status)
        class(fitpack_parametric_curve), intent(inout) :: self
        real(real64), intent(in) :: x(:,:)
        real(real64), intent(in), optional :: u(:),w(:),smoothing
        integer, intent(in), optional :: order
        integer :: degree, j
        degree=3
        if(present(order)) degree=order
        self%order=degree
        self%idim=size(x,1)
        self%knots=2*(degree+1)
        self%t=[spread(0.0_real64,1,degree+1),spread(1.0_real64,1,degree+1)]
        self%c=spread(0.0_real64,1,self%knots*self%idim)
        self%u=[(real(j-1,real64)/real(max(1,size(x,2)-1),real64),j=1,size(x,2))]
        if(present(u)) self%u=u
        status=0
    end function
    real(real64) function param_mse(self) result(v)
        class(fitpack_parametric_curve), intent(in) :: self
        v=0.0_real64
    end function
    subroutine param_destroy(self)
        class(fitpack_parametric_curve), intent(inout) :: self
        if(allocated(self%t)) deallocate(self%t)
        if(allocated(self%c)) deallocate(self%c)
        if(allocated(self%u)) deallocate(self%u)
    end subroutine
    subroutine surface_new_points(self,x,y,z,w)
        class(fitpack_surface), intent(inout) :: self
        real(real64),intent(in) :: x(:),y(:),z(:)
        real(real64),intent(in),optional :: w(:)
        self%m=size(x)
        self%nmax=2*(self%m+4)
        self%nest=self%nmax
        allocate(self%t(self%nmax,2),source=0.0_real64)
        allocate(self%c(1),source=0.0_real64)
        allocate(self%wrk(1),source=0.0_real64)
        allocate(self%wrk2(1),source=0.0_real64)
        allocate(self%iwrk(1),source=0)
    end subroutine
    integer(FP_FLAG) function surface_new_fit(self,x,y,z,w,smoothing,order) result(status)
        class(fitpack_surface), intent(inout) :: self
        real(real64),intent(in) :: x(:),y(:),z(:)
        real(real64),intent(in),optional :: w(:),smoothing
        integer,intent(in),optional :: order
        call self%new_points(x,y,z,w)
        status=self%fit(smoothing=smoothing,order=order)
    end function
    integer(FP_FLAG) function surface_fit(self,smoothing,order,keep_knots) result(status)
        class(fitpack_surface), intent(inout) :: self
        real(real64),intent(in),optional :: smoothing
        integer,intent(in),optional :: order
        logical,intent(in),optional :: keep_knots
        integer :: degree
        degree=3
        if(present(order)) degree=order
        self%order=degree
        self%knots=2*(degree+1)
        self%t(:,1)=0.0_real64
        self%t(:,2)=0.0_real64
        self%t(degree+2:self%knots(1),1)=1.0_real64
        self%t(degree+2:self%knots(2),2)=1.0_real64
        status=0
    end function
    real(real64) function surface_mse(self) result(v)
        class(fitpack_surface),intent(in) :: self
        v=0.0_real64
    end function
    subroutine surface_destroy(self)
        class(fitpack_surface),intent(inout) :: self
        if(allocated(self%t))deallocate(self%t)
        if(allocated(self%c))deallocate(self%c)
        if(allocated(self%wrk))deallocate(self%wrk)
        if(allocated(self%wrk2))deallocate(self%wrk2)
        if(allocated(self%iwrk))deallocate(self%iwrk)
    end subroutine
end module
'''
TEST = r'''
program verify_optional_adapter_mock
    use scifort_kinds, only: dp
    use scifort_optional_fitpack, only: univariate_spline,fitpack_smooth_success, &
        fitpack_smooth_invalid_input
    use scifort_interpolate_fitpack_compat, only: fitpack_tck_type
    implicit none
    type(fitpack_tck_type) :: tck
    real(dp) :: x(10),y(10),resid
    integer :: status,i
    do i=1,10
        x(i)=0.1_dp*real(i-1,dp)
        y(i)=sin(x(i))
    end do
    tck=univariate_spline(x,y,s=0.1_dp,status=status,residual=resid)
    if(status /= fitpack_smooth_success .or. .not.tck%valid()) error stop 'adapter mock'
    if(size(tck%t) /= 8 .or. size(tck%c) /= 8) error stop 'tuple sizing mock'
    tck=univariate_spline(x,y,s=-1.0_dp,status=status)
    if(status /= fitpack_smooth_invalid_input) error stop 'invalid target mock'
    print *, 'Synthetic FITPACK adapter interface PASS (NOT upstream verification)'
end program
'''

def cmd(args, cwd):
    p = subprocess.run(args,cwd=cwd,capture_output=True,text=True)
    if p.returncode:
        raise RuntimeError(' '.join(args)+'\n'+p.stdout[-1000:]+'\n'+p.stderr[-2500:])
    return p.stdout.strip()

def main():
    with tempfile.TemporaryDirectory() as tmp:
        work=Path(tmp)
        (work/'fitpack_mock.f90').write_text(MOCK)
        (work/'test_mock.f90').write_text(TEST)
        compiler=['gfortran','-std=f2018','-ffree-line-length-none','-O0','-g','-fcheck=all']
        sources=(ROOT/'src/kinds.f90',ROOT/'src/interpolate_bspline.f90',
                 ROOT/'src/interpolate_fitpack_bivariate.f90',
                 ROOT/'src/interpolate_fitpack_compat.f90',work/'fitpack_mock.f90',
                 ROOT/'optional/fitpack_backend/src/scifort_optional_fitpack.f90')
        for i,source in enumerate(sources):
            cmd([*compiler,'-J'+str(work),'-I'+str(work),'-c',str(source),
                 '-o',str(work/f'{i}.o')],work)
        exe=work/'test_fitpack_mock'
        cmd([*compiler,'-I'+str(work),str(work/'test_mock.f90'),
             *(str(work/f'{i}.o') for i in range(len(sources))),'-o',str(exe)],work)
        print(cmd([str(exe)],work))
if __name__=='__main__':main()
