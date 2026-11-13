! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Machine-readable comparison driver; built/run by check_bvp.py, outside fpm test.
program bvp_reference_driver
    use scifort_integrate_bvp, only : solve_bvp,solve_bvp_result,solve_bvp_complex_result
    use scifort_kinds, only : dp
    implicit none
    type(solve_bvp_result) :: r
    type(solve_bvp_complex_result) :: z
    integer :: family,branch,j,nu
    real(dp) :: coefficient,tol,x(5),y(2,5),t(81),out(2,81),s(2,2)
    complex(dp) :: cy(2,5),cout(2,81),cs(2,2)
    character(80) :: arg
    call get_command_argument(1,arg)
    read(arg,*) family
    call get_command_argument(2,arg)
    read(arg,*) coefficient
    call get_command_argument(3,arg)
    read(arg,*) tol
    call get_command_argument(4,arg)
    read(arg,*) branch
    x=[(real(j-1,dp)/4.0_dp,j=1,5)]
    t=[(real(j-1,dp)/80.0_dp,j=1,81)]
    y=0.0_dp
    s=0.0_dp
    s(2,2)=-2.0_dp
    select case(family)
    case(0,2)
        call solve_bvp(rhs,bc,x,y,r,tol=tol,max_nodes=3000)
    case(1)
        if (branch==1) y(1,:)=3.0_dp
        call solve_bvp(rhs,bc,x,y,r,tol=tol,max_nodes=3000)
    case(3)
        y(1,:)=1.0_dp
        call solve_bvp(rhs,bc,x,y,r,s=s,tol=tol,max_nodes=3000)
    case(4)
        y(1,:)=exp(0.8_dp*coefficient*x)
        y(2,:)=2.0_dp*exp(-1.5_dp*coefficient*x)
        call solve_bvp(rhs,bc,x,y,r,p=[0.8_dp*coefficient,-1.5_dp*coefficient], &
            tol=tol,max_nodes=3000)
    case(5,6)
        cy=cmplx(1.0_dp,0.0_dp,dp)
        if (family==6) then
            cy(2,:)=cmplx(0.0_dp,0.0_dp,dp)
            cs=cmplx(s,0.0_dp,dp)
            call solve_bvp(crhs,cbc,x,cy,z,s=cs,tol=tol,max_nodes=3000)
        else
            call solve_bvp(crhs,cbc,x,cy,z,tol=tol,max_nodes=3000)
        end if
    end select
    if (family<5) then
        print *, r%status,size(r%x),r%niter,maxval(r%rms_residuals),maxval(abs(r%bc_residuals))
        print *, size(r%p),r%p
        do nu=0,2
            call r%sol%evaluate(t,out,nu=nu)
            do j=1,81
                print '(5es26.16e3)',t(j),out(:,j),0.0_dp,0.0_dp
            end do
        end do
    else
        print *, z%status,size(z%x),z%niter,maxval(z%rms_residuals),maxval(abs(z%bc_residuals))
        print *, 0
        do nu=0,2
            call z%sol%evaluate(t,cout,nu=nu)
            do j=1,81
                print '(5es26.16e3)',t(j),real(cout(:,j),dp),aimag(cout(:,j))
            end do
        end do
    end if
contains
    subroutine rhs(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: f(:,:) !! differential equations
        select case(family)
        case(0)
            f(1,:)=coefficient*v(1,:)
            f(2,:)=-coefficient*v(2,:)
        case(1)
            f(1,:)=v(2,:)
            f(2,:)=-coefficient*exp(v(1,:))
        case(2)
            f(1,:)=v(2,:)
            f(2,:)=-coefficient*v(2,:)
        case(3)
            f(1,:)=v(2,:)
            f(2,:)=-coefficient**2*v(1,:)
        case(4)
            f(1,:)=p(1)*v(1,:)
            f(2,:)=p(2)*v(2,:)
        end select
    end subroutine rhs
    subroutine bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! boundary constraints
        select case(family)
        case(0)
            r=[a(1)+b(2)-1.0_dp-2.0_dp*exp(-coefficient), &
                b(1)+2.0_dp*a(2)-exp(coefficient)-4.0_dp]
        case(1)
            r=[a(1),b(1)]
        case(2)
            r=[a(1),b(1)-1.0_dp+exp(-coefficient)]
        case(3)
            r=[a(2),b(1)-sin(coefficient)/coefficient]
        case(4)
            r=[a(1)-1.0_dp,a(2)-2.0_dp,b(1)-exp(coefficient), &
                b(2)-2.0_dp*exp(-2.0_dp*coefficient)]
        end select
    end subroutine bc
    subroutine crhs(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:) !! coordinates
        complex(dp), intent(in) :: v(:,:),p(:) !! complex states and parameters
        complex(dp), intent(out) :: f(:,:) !! differential equations
        complex(dp) :: q
        if (family==5) then
            f(1,:)=cmplx(0.0_dp,coefficient,dp)*v(1,:)
            f(2,:)=-coefficient*v(2,:)
        else
            q=coefficient*cmplx(1.0_dp,0.3_dp,dp)
            f(1,:)=v(2,:)
            f(2,:)=-q*q*v(1,:)
        end if
    end subroutine crhs
    subroutine cbc(a,b,p,r)
        implicit none
        complex(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        complex(dp), intent(out) :: r(:) !! complex endpoint constraints
        complex(dp) :: q
        if (family==5) then
            r=[a(1)+b(2)-1.0_dp-2.0_dp*exp(-coefficient), &
                b(1)+2.0_dp*a(2)-exp(cmplx(0.0_dp,coefficient,dp))-4.0_dp]
        else
            q=coefficient*cmplx(1.0_dp,0.3_dp,dp)
            r=[a(2),b(1)-sin(q)/q]
        end if
    end subroutine cbc
end program bvp_reference_driver

