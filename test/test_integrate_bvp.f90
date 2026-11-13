! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_integrate_bvp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_is_nan
    use scifort_integrate, only : dp, solve_bvp, solve_bvp_result, solve_bvp_complex_result, &
        bvp_status_invalid_argument, bvp_status_max_nodes, bvp_status_singular, &
        bvp_status_numerical_failure, bvp_status_bc_tolerance
    implicit none
    type(solve_bvp_result) :: r,r2
    type(solve_bvp_complex_result) :: z
    real(dp) :: x(5),y(2,5),xx(101),yy(2,101),val(2),bad(1),s(2,2),nan,inf,pi
    complex(dp) :: cy(1,5),cv(1,101),cp(1)
    integer :: j,info,checks
    checks=0
    pi=acos(-1.0_dp)
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    x=[(real(j-1,dp)/4.0_dp,j=1,5)]
    xx=[(real(j-1,dp)/100.0_dp,j=1,101)]
    y=0.0_dp

    call solve_bvp(linear,linear_bc,x,y,r,tol=1e-8_dp)
    call check(r%success,'linear convergence')
    call r%sol%evaluate(xx,yy,info)
    call check(info==0,'vector dense evaluation')
    call check(maxval(abs(yy(1,:)-exp(xx)))<2e-8_dp,'exponential dense solution')
    call check(maxval(abs(yy(2,:)-2.0_dp*exp(-xx)))<2e-8_dp,'nonseparated BC solution')
    call check(maxval(r%rms_residuals)<=r%tol,'residual tolerance')
    call check(maxval(abs(r%bc_residuals))<1e-10_dp,'boundary tolerance')
    call check(size(r%x)>5 .and. r%niter>1,'adaptive refinement')
    call solve_bvp(linear,linear_bc,x,y,r2,fun_jac=linear_jac,bc_jac=linear_bc_jac,tol=1e-8_dp)
    call check(r2%success,'analytic Jacobians')
    call r2%sol%evaluate(xx,yy)
    call check(maxval(abs(yy(1,:)-exp(xx)))<2e-8_dp,'analytic dense solution')
    call r%sol%evaluate(0.3_dp,val,info,nu=1)
    call check(info==0 .and. abs(val(1)-exp(0.3_dp))<2e-7_dp,'dense first derivative')
    call r%sol%evaluate(0.3_dp,val,info,nu=2)
    call check(info==0 .and. abs(val(1)-exp(0.3_dp))<2e-4_dp,'dense second derivative')
    call r%sol%evaluate(0.3_dp,val,info,nu=3)
    call check(info==0 .and. abs(val(1)-exp(0.3_dp))<0.02_dp,'dense third derivative')
    call r%sol%evaluate(0.3_dp,val,info,nu=4)
    call check(info==0 .and. all(val==0.0_dp),'higher derivatives vanish')
    call r%sol%evaluate([0.8_dp,0.1_dp],yy(:,:2),info)
    call check(info==0 .and. abs(yy(1,1)-exp(0.8_dp))<2e-8_dp,'unsorted queries')
    r%y=0.0_dp
    call r%sol%evaluate(0.3_dp,val)
    call check(abs(val(1)-exp(0.3_dp))<2e-8_dp,'dense object owns independent state')
    call r%sol%evaluate(-0.01_dp,val,info)
    call check(info==0 .and. abs(val(1)-exp(-0.01_dp))<1e-6_dp,'extrapolation')
    call r%sol%evaluate(-0.01_dp,val,info,extrapolate=.false.)
    call check(info==-1 .and. all(ieee_is_nan(val)),'disabled extrapolation')
    call r%sol%evaluate(nan,val,info)
    call check(info==-1 .and. all(ieee_is_nan(val)),'nonfinite query')
    call r%sol%evaluate(0.3_dp,bad,info)
    call check(info==-1 .and. all(ieee_is_nan(bad)),'invalid output shape')
    call r%sol%evaluate(0.3_dp,val,info,nu=-1)
    call check(info==-1,'invalid derivative order')

    ! Bratu y''=-exp(y), zero endpoint values: two independent analytic branches.
    do j=1,2
        y=0.0_dp
        if (j==2) y(1,:)=3.0_dp
        call solve_bvp(bratu,zero_bc,x,y,r,tol=1e-7_dp)
        call check(r%success,'Bratu convergence')
        call r%sol%evaluate(xx,yy)
        if (j==1) then
            call check(maxval(abs(yy(1,:)-bratu_exact(xx,0.3792911497626886_dp))) &
                <2e-7_dp,'Bratu lower branch')
        else
            call check(maxval(abs(yy(1,:)-bratu_exact(xx,2.7346756930305267_dp))) &
                <2e-7_dp,'Bratu upper branch')
        end if
    end do

    y(1,:)=sin(pi*x)
    y(2,:)=pi*cos(pi*x)
    call solve_bvp(eigen_rhs,eigen_bc,x,y,r,p=[3.0_dp],tol=1e-8_dp)
    call check(r%success,'unknown parameter convergence')
    call check(abs(r%p(1)-pi)<2e-8_dp,'eigenvalue parameter')
    call r%sol%evaluate(xx,yy)
    call check(maxval(abs(yy(1,:)-sin(pi*xx)))<2e-8_dp,'normalized eigenfunction')
    call solve_bvp(eigen_rhs,eigen_bc,x,y,r,p=[3.0_dp],tol=1e-8_dp, &
        fun_jac=eigen_jac,bc_jac=eigen_bc_jac)
    call check(r%success .and. abs(r%p(1)-pi)<2e-8_dp,'analytic parameter Jacobians')

    y=0.0_dp
    call solve_bvp(layer,layer_bc,x,y,r,tol=1e-6_dp)
    call check(r%success,'boundary layer convergence')
    call check(maxval(r%x(2:)-r%x(:size(r%x)-1))> &
        2.0_dp*minval(r%x(2:)-r%x(:size(r%x)-1)),'nonuniform adaptive mesh')
    call r%sol%evaluate(xx,yy)
    call check(maxval(abs(yy(1,:)-(1.0_dp-exp(-100.0_dp*xx))))<3e-7_dp,'boundary layer')

    ! Lane-Emden n=1: y''+2*y'/x=-y, y(0)=1, y(1)=sin(1).
    y(1,:)=1.0_dp
    y(2,:)=0.0_dp
    s=0.0_dp
    s(2,2)=-2.0_dp
    call solve_bvp(emden,emden_bc,x,y,r,s=s,tol=1e-7_dp,fun_jac=emden_jac)
    call check(r%success,'singular endpoint convergence')
    call check(abs(r%y(2,1))<1e-12_dp,'singular endpoint constraint')
    call r%sol%evaluate(xx,yy)
    call check(abs(yy(1,1)-1.0_dp)<2e-7_dp,'Emden origin')
    call check(maxval(abs(yy(1,2:)-sin(xx(2:))/xx(2:)))<2e-7_dp,'Emden analytic solution')
    call solve_bvp(emden,emden_bc,x,y,r,s=s,tol=1e-7_dp)
    call check(r%success,'singular endpoint finite-difference Jacobian')

    cy=cmplx(1.0_dp,0.0_dp,dp)
    cp=cmplx(0.2_dp,0.9_dp,dp)
    call solve_bvp(complex_rhs,complex_bc,x,cy,z,p=cp,tol=1e-7_dp)
    call check(z%success,'complex parameter convergence')
    call check(abs(z%p(1)-cmplx(0.0_dp,1.0_dp,dp))<3e-7_dp,'complex parameter')
    call z%sol%evaluate(xx,cv,info)
    call check(info==0 .and. maxval(abs(cv(1,:)-exp(cmplx(0.0_dp,xx,dp)))) &
        <3e-7_dp,'complex dense solution')
    call solve_bvp(complex_rhs,complex_bc,x,cy,z,p=cp,tol=1e-7_dp, &
        fun_jac=complex_jac,bc_jac=complex_bc_jac)
    call check(z%success,'complex analytic Jacobians')

    y=0.0_dp
    call solve_bvp(linear,linear_bc,x,y,r,max_nodes=5,tol=1e-10_dp)
    call check(r%status==bvp_status_max_nodes .and. .not.r%success,'mesh budget')
    call check(allocated(r%rms_residuals),'failure residual diagnostics')
    call solve_bvp(linear,singular_bc,x,y,r)
    call check(r%status==bvp_status_singular,'singular boundary Jacobian')
    call solve_bvp(linear,impossible_bc,x,y,r,bc_jac=impossible_jac)
    call check(r%status==bvp_status_bc_tolerance,'unattainable boundary tolerance')
    call solve_bvp(bad_rhs,linear_bc,x,y,r)
    call check(r%status==bvp_status_numerical_failure .and. .not.r%success,'NaN callback')
    call solve_bvp(linear,linear_bc,x,y,r,fun_jac=bad_jac)
    call check(r%status==bvp_status_numerical_failure,'NaN Jacobian')
    call solve_bvp(linear,linear_bc,x,y,r,tol=0.0_dp)
    call check(r%status==bvp_status_invalid_argument,'invalid tolerance')
    call solve_bvp(linear,linear_bc,x,y,r,bc_tol=inf)
    call check(r%status==bvp_status_invalid_argument,'infinite boundary tolerance')
    call solve_bvp(linear,linear_bc,x,y,r,p=[nan])
    call check(r%status==bvp_status_invalid_argument,'NaN parameter')
    call solve_bvp(linear,linear_bc,x,y,r,max_nodes=4)
    call check(r%status==bvp_status_invalid_argument,'invalid budget')
    call solve_bvp(linear,linear_bc,[0.0_dp,0.5_dp,0.5_dp,0.7_dp,1.0_dp],y,r)
    call check(r%status==bvp_status_invalid_argument,'repeated mesh coordinate')
    y(1,1)=inf
    call solve_bvp(linear,linear_bc,x,y,r)
    call check(r%status==bvp_status_invalid_argument,'infinite state')
    call r%sol%evaluate(0.1_dp,val,info)
    call check(info==-1,'uninitialized dense solution')
    print '(a,i0,a)', 'BVP: ',checks,' checks passed.'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! assertion result
        character(*), intent(in) :: label !! assertion context
        if (.not.ok) then
            print *, 'FAIL: ',label
            error stop 1
        end if
        checks=checks+1
    end subroutine check
    function bratu_exact(t,a) result(out)
        implicit none
        real(dp), intent(in) :: t(:),a !! coordinates and positive Bratu root
        real(dp) :: out(size(t))
        out=-2.0_dp*log(cosh(a*(2.0_dp*t-1.0_dp))/cosh(a))
    end function bratu_exact
    subroutine linear(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: f(:,:) !! differential equations
        f(1,:)=v(1,:)
        f(2,:)=-v(2,:)
    end subroutine linear
    subroutine linear_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! coupled endpoint constraints
        r=[a(1)+b(2)-1.0_dp-2.0_dp/exp(1.0_dp),b(1)+2.0_dp*a(2)-exp(1.0_dp)-4.0_dp]
    end subroutine linear_bc
    subroutine linear_jac(t,v,p,a,b)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: a(:,:,:),b(:,:,:) !! state and parameter Jacobians
        a=0.0_dp
        a(1,1,:)=1.0_dp
        a(2,2,:)=-1.0_dp
        b=0.0_dp
    end subroutine linear_jac
    subroutine linear_bc_jac(a,b,p,da,db,dpb)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: da(:,:),db(:,:),dpb(:,:) !! boundary Jacobians
        da=0.0_dp
        db=0.0_dp
        dpb=0.0_dp
        da(1,1)=1.0_dp
        da(2,2)=2.0_dp
        db(1,2)=1.0_dp
        db(2,1)=1.0_dp
    end subroutine linear_bc_jac
    subroutine bratu(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: f(:,:) !! first-order Bratu equations
        f(1,:)=v(2,:)
        f(2,:)=-exp(v(1,:))
    end subroutine bratu
    subroutine zero_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! homogeneous endpoint values
        r=[a(1),b(1)]
    end subroutine zero_bc
    subroutine eigen_rhs(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, eigenvalue
        real(dp), intent(out) :: f(:,:) !! oscillation equations
        f(1,:)=v(2,:)
        f(2,:)=-p(1)**2*v(1,:)
    end subroutine eigen_rhs
    subroutine eigen_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and eigenvalue
        real(dp), intent(out) :: r(:) !! endpoints and normalization
        r=[a(1),b(1),a(2)-p(1)]
    end subroutine eigen_bc
    subroutine eigen_jac(t,v,p,a,b)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, eigenvalue
        real(dp), intent(out) :: a(:,:,:),b(:,:,:) !! state and parameter derivatives
        a=0.0_dp
        b=0.0_dp
        a(1,2,:)=1.0_dp
        a(2,1,:)=-p(1)**2
        b(2,1,:)=-2.0_dp*p(1)*v(1,:)
    end subroutine eigen_jac
    subroutine eigen_bc_jac(a,b,p,da,db,dpb)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and eigenvalue
        real(dp), intent(out) :: da(:,:),db(:,:),dpb(:,:) !! boundary derivatives
        da=0.0_dp
        db=0.0_dp
        dpb=0.0_dp
        da(1,1)=1.0_dp
        db(2,1)=1.0_dp
        da(3,2)=1.0_dp
        dpb(3,1)=-1.0_dp
    end subroutine eigen_bc_jac
    subroutine layer(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: f(:,:) !! thin boundary layer equations
        f(1,:)=v(2,:)
        f(2,:)=-100.0_dp*v(2,:)
    end subroutine layer
    subroutine layer_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! prescribed layer endpoints
        r=[a(1),b(1)-1.0_dp+exp(-100.0_dp)]
    end subroutine layer_bc
    subroutine emden(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: f(:,:) !! regular part of Lane-Emden equations
        f(1,:)=v(2,:)
        f(2,:)=-v(1,:)
    end subroutine emden
    subroutine emden_jac(t,v,p,a,b)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: a(:,:,:),b(:,:,:) !! regular-part Jacobians
        a=0.0_dp
        a(1,2,:)=1.0_dp
        a(2,1,:)=-1.0_dp
        b=0.0_dp
    end subroutine emden_jac
    subroutine emden_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! regularity and outer value
        r=[a(2),b(1)-sin(1.0_dp)]
    end subroutine emden_bc
    subroutine complex_rhs(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:) !! coordinates
        complex(dp), intent(in) :: v(:,:),p(:) !! complex state and parameter
        complex(dp), intent(out) :: f(:,:) !! holomorphic equations
        f=p(1)*v
    end subroutine complex_rhs
    subroutine complex_bc(a,b,p,r)
        implicit none
        complex(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        complex(dp), intent(out) :: r(:) !! complex endpoint constraints
        r=[a(1)-1.0_dp,b(1)-exp(cmplx(0.0_dp,1.0_dp,dp))]
    end subroutine complex_bc
    subroutine complex_jac(t,v,p,a,b)
        implicit none
        real(dp), intent(in) :: t(:) !! coordinates
        complex(dp), intent(in) :: v(:,:),p(:) !! state and parameters
        complex(dp), intent(out) :: a(:,:,:),b(:,:,:) !! holomorphic Jacobians
        a(1,1,:)=p(1)
        b(1,1,:)=v(1,:)
    end subroutine complex_jac
    subroutine complex_bc_jac(a,b,p,da,db,dpb)
        implicit none
        complex(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        complex(dp), intent(out) :: da(:,:),db(:,:),dpb(:,:) !! boundary Jacobians
        da=0.0_dp
        db=0.0_dp
        dpb=0.0_dp
        da(1,1)=1.0_dp
        db(2,1)=1.0_dp
    end subroutine complex_bc_jac
    subroutine singular_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! redundant constraints
        r=[a(1),a(1)-1.0_dp]
    end subroutine singular_bc
    subroutine impossible_bc(a,b,p,r)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: r(:) !! intentionally unattainable constraints
        r=1.0_dp
    end subroutine impossible_bc
    subroutine impossible_jac(a,b,p,da,db,dpb)
        implicit none
        real(dp), intent(in) :: a(:),b(:),p(:) !! endpoint states and parameters
        real(dp), intent(out) :: da(:,:),db(:,:),dpb(:,:) !! deliberately inconsistent derivative
        da=0.0_dp
        db=0.0_dp
        dpb=0.0_dp
        da(1,1)=1.0_dp
        da(2,2)=1.0_dp
    end subroutine impossible_jac
    subroutine bad_rhs(t,v,p,f)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: f(:,:) !! deliberate nonfinite output
        f=ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine bad_rhs
    subroutine bad_jac(t,v,p,a,b)
        implicit none
        real(dp), intent(in) :: t(:),v(:,:),p(:) !! coordinates, states, parameters
        real(dp), intent(out) :: a(:,:,:),b(:,:,:) !! deliberate nonfinite derivatives
        a=ieee_value(0.0_dp,ieee_quiet_nan)
        b=0.0_dp
    end subroutine bad_jac
end program test_integrate_bvp

