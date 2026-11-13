! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program aaa_reference_driver
    use scifort_interpolate_aaa, only : aaa,complex_aaa_type
    use scifort_kinds, only : dp
    implicit none
    type(complex_aaa_type) :: r
    complex(dp) :: x(201),f(201),q(121),a,b,value
    complex(dp), allocatable :: support(:),poles(:),roots(:),residues(:)
    real(dp) :: coefficient,tol,t,ctol
    integer :: family,budget,clean,i,status,info,ps,rs,ds
    character(80) :: arg
    call get_command_argument(1,arg)
    read(arg,*) family
    call get_command_argument(2,arg)
    read(arg,*) coefficient
    call get_command_argument(3,arg)
    read(arg,*) tol
    call get_command_argument(4,arg)
    read(arg,*) budget
    call get_command_argument(5,arg)
    read(arg,*) clean
    do i=1,201
        t=-2.0_dp+4.0_dp*real(i-1,dp)/200.0_dp
        x(i)=cmplx(t,0.0_dp,dp)
    end do
    do i=1,121
        q(i)=cmplx(-1.99_dp+3.98_dp*real(i-1,dp)/120.0_dp,0.0_dp,dp)
    end do
    select case(family)
    case(0)
        f=(x+0.2_dp)/(x*x+0.4_dp*x+1+coefficient**2)
    case(1)
        f=exp(coefficient*x)
    case(2)
        f=tanh(coefficient*x)
    case(3)
        f=sin(coefficient*x)
    case(4)
        x=exp(-coefficient*(2.0_dp-real(x,dp))/4.0_dp)
        q=exp(-coefficient*(2.0_dp-real(q,dp))/4.0_dp)
        f=log(x)
    case(5)
        x=2.0_dp*exp(cmplx(0.0_dp,acos(-1.0_dp),dp)*(x+2.0_dp)/2.01_dp)
        q=1.4_dp*exp(cmplx(0.0_dp,acos(-1.0_dp),dp)*(q+2.0_dp)/2.01_dp)
        a=coefficient*cmplx(0.3_dp,0.4_dp,dp)
        b=coefficient*cmplx(-0.2_dp,0.6_dp,dp)
        f=(x+0.2_dp)/((x-a)*(x-b))
    case(6)
        x=0.75_dp*x+0.003_dp
        q=0.75_dp*q+0.003_dp
        do i=1,201
            f(i)=cmplx(gamma(real(x(i),dp)),0.0_dp,dp)
        end do
    case(7)
        f=1.0_dp/(x+3.0_dp)+coefficient/(x-0.313_dp)
    case(8)
        f=1.0_dp/(x+3.0_dp)+1e-8_dp*sin(coefficient*x)
    end select
    ctol=1e-13_dp
    if (clean==2) ctol=1e-6_dp
    r=aaa(x,f,rtol=tol,max_terms=budget,clean_up=clean/=0,clean_up_tol=ctol,status=status)
    support=r%support_points()
    print *, status,r%fit_status(),r%cleanup_status(),size(support),r%max_error()
    if (.not.r%valid()) stop
    do i=1,121
        value=r%value(q(i),info)
        print '(4es26.17e3,1x,i0)',real(q(i),dp),aimag(q(i)),real(value,dp),aimag(value),info
    end do
    poles=r%poles(ps)
    roots=r%roots(rs)
    residues=r%residues(ds)
    print *, ps,rs,ds,size(poles),size(roots),size(residues)
    do i=1,size(poles)
        print '(2es26.17e3)',real(poles(i),dp),aimag(poles(i))
    end do
    do i=1,size(roots)
        print '(2es26.17e3)',real(roots(i),dp),aimag(roots(i))
    end do
    do i=1,size(residues)
        print '(2es26.17e3)',real(residues(i),dp),aimag(residues(i))
    end do
end program aaa_reference_driver
