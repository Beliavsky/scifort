! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
!
! Adaptive Lobatto IIIA collocation adapted from SciPy v1.15.3 integrate/_bvp.py.
! Full upstream terms: LICENSES/scipy-bvp-BSD-3-Clause.txt.
! Native sparse LU, callback/result interfaces and failure checks are SciFort code.
module scifort_integrate_bvp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_integrate_sparse_jacobian, only : ivp_sparse_matrix, ivp_complex_sparse_matrix
    use scifort_integrate_sparse_lu, only : ivp_sparse_lu_real, ivp_sparse_lu_complex, &
        sparse_lu_factor
    use scifort_linalg_svd, only : pinv
    implicit none
    private
    integer, parameter, public :: bvp_status_success = 0
    integer, parameter, public :: bvp_status_max_nodes = 1
    integer, parameter, public :: bvp_status_singular = 2
    integer, parameter, public :: bvp_status_bc_tolerance = 3
    integer, parameter, public :: bvp_status_invalid_argument = -1
    integer, parameter, public :: bvp_status_numerical_failure = -2
    integer, parameter, public :: bvp_status_setup_failure = -3
    public :: solve_bvp

    type, public :: bvp_solution
        private
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:,:), f(:,:)
    contains
        procedure, private :: evaluate_scalar => evaluate_scalar_real
        procedure, private :: evaluate_vector => evaluate_vector_real
        generic :: evaluate => evaluate_scalar, evaluate_vector
    end type bvp_solution

    type, public :: solve_bvp_result
        real(dp), allocatable :: x(:) !! final mesh
        real(dp), allocatable :: y(:,:), yp(:,:), p(:) !! states, derivatives, fitted parameters
        real(dp), allocatable :: rms_residuals(:) !! normalized RMS defect per interval
        real(dp), allocatable :: bc_residuals(:) !! final boundary residuals
        type(bvp_solution) :: sol !! independent C1 piecewise-cubic dense solution
        integer :: status = bvp_status_invalid_argument
        integer :: niter = 0
        integer :: nfev = 0
        integer :: njev = 0
        logical :: success = .false.
        real(dp) :: tol = 0.0_dp
        character(:), allocatable :: message
    end type solve_bvp_result
    public :: bvp_fun_real, bvp_bc_real, bvp_fun_jac_real, bvp_bc_jac_real
    abstract interface
        subroutine bvp_fun_real(x,y,p,f)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent mesh coordinates
            real(dp), intent(in) :: y(:,:) !! state columns at x, shape (n,m)
            real(dp), intent(in) :: p(:) !! unknown parameters, possibly zero length
            real(dp), intent(out) :: f(:,:) !! dy/dx columns, same shape as y
        end subroutine bvp_fun_real
        subroutine bvp_bc_real(ya,yb,p,residual)
            import :: dp
            implicit none
            real(dp), intent(in) :: ya(:) !! left endpoint state
            real(dp), intent(in) :: yb(:) !! right endpoint state
            real(dp), intent(in) :: p(:) !! unknown parameters, possibly zero length
            real(dp), intent(out) :: residual(:) !! n+k two-point boundary residuals
        end subroutine bvp_bc_real
        subroutine bvp_fun_jac_real(x,y,p,df_dy,df_dp)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent coordinates
            real(dp), intent(in) :: y(:,:) !! state columns, shape (n,m)
            real(dp), intent(in) :: p(:) !! unknown parameters
            real(dp), intent(out) :: df_dy(:,:,:) !! state Jacobian, shape (n,n,m)
            real(dp), intent(out) :: df_dp(:,:,:) !! parameter Jacobian, shape (n,k,m)
        end subroutine bvp_fun_jac_real
        subroutine bvp_bc_jac_real(ya,yb,p,dya,dyb,dp_bc)
            import :: dp
            implicit none
            real(dp), intent(in) :: ya(:) !! left endpoint state
            real(dp), intent(in) :: yb(:) !! right endpoint state
            real(dp), intent(in) :: p(:) !! unknown parameters
            real(dp), intent(out) :: dya(:,:) !! left-state boundary Jacobian, shape (n+k,n)
            real(dp), intent(out) :: dyb(:,:) !! right-state boundary Jacobian, shape (n+k,n)
            real(dp), intent(out) :: dp_bc(:,:) !! parameter boundary Jacobian, shape (n+k,k)
        end subroutine bvp_bc_jac_real
    end interface

    type, public :: bvp_complex_solution
        private
        real(dp), allocatable :: x(:)
        complex(dp), allocatable :: y(:,:), f(:,:)
    contains
        procedure, private :: evaluate_scalar => evaluate_scalar_complex
        procedure, private :: evaluate_vector => evaluate_vector_complex
        generic :: evaluate => evaluate_scalar, evaluate_vector
    end type bvp_complex_solution

    type, public :: solve_bvp_complex_result
        real(dp), allocatable :: x(:) !! final mesh
        complex(dp), allocatable :: y(:,:), yp(:,:), p(:) !! states, derivatives, parameters
        real(dp), allocatable :: rms_residuals(:) !! normalized RMS defect per interval
        complex(dp), allocatable :: bc_residuals(:) !! final boundary residuals
        type(bvp_complex_solution) :: sol !! independent C1 piecewise-cubic dense solution
        integer :: status = bvp_status_invalid_argument
        integer :: niter = 0
        integer :: nfev = 0
        integer :: njev = 0
        logical :: success = .false.
        real(dp) :: tol = 0.0_dp
        character(:), allocatable :: message
    end type solve_bvp_complex_result
    public :: bvp_fun_complex, bvp_bc_complex, bvp_fun_jac_complex, bvp_bc_jac_complex
    abstract interface
        subroutine bvp_fun_complex(x,y,p,f)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent mesh coordinates
            complex(dp), intent(in) :: y(:,:) !! state columns at x, shape (n,m)
            complex(dp), intent(in) :: p(:) !! unknown parameters, possibly zero length
            complex(dp), intent(out) :: f(:,:) !! dy/dx columns, same shape as y
        end subroutine bvp_fun_complex
        subroutine bvp_bc_complex(ya,yb,p,residual)
            import :: dp
            implicit none
            complex(dp), intent(in) :: ya(:) !! left endpoint state
            complex(dp), intent(in) :: yb(:) !! right endpoint state
            complex(dp), intent(in) :: p(:) !! unknown parameters, possibly zero length
            complex(dp), intent(out) :: residual(:) !! n+k two-point boundary residuals
        end subroutine bvp_bc_complex
        subroutine bvp_fun_jac_complex(x,y,p,df_dy,df_dp)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent coordinates
            complex(dp), intent(in) :: y(:,:) !! state columns, shape (n,m)
            complex(dp), intent(in) :: p(:) !! unknown parameters
            complex(dp), intent(out) :: df_dy(:,:,:) !! state Jacobian, shape (n,n,m)
            complex(dp), intent(out) :: df_dp(:,:,:) !! parameter Jacobian, shape (n,k,m)
        end subroutine bvp_fun_jac_complex
        subroutine bvp_bc_jac_complex(ya,yb,p,dya,dyb,dp_bc)
            import :: dp
            implicit none
            complex(dp), intent(in) :: ya(:) !! left endpoint state
            complex(dp), intent(in) :: yb(:) !! right endpoint state
            complex(dp), intent(in) :: p(:) !! unknown parameters
            complex(dp), intent(out) :: dya(:,:) !! left-state boundary Jacobian, shape (n+k,n)
            complex(dp), intent(out) :: dyb(:,:) !! right-state boundary Jacobian, shape (n+k,n)
            complex(dp), intent(out) :: dp_bc(:,:) !! parameter boundary Jacobian, shape (n+k,k)
        end subroutine bvp_bc_jac_complex
    end interface

    interface solve_bvp
        module procedure solve_bvp_real, solve_bvp_complex
    end interface solve_bvp
    interface finite
        module procedure finite_r1, finite_r2, finite_r3, finite_c1, finite_c2, finite_c3
    end interface finite
contains

    subroutine solve_bvp_real(fun,bc,x,y,res,p,s,fun_jac,bc_jac,tol,max_nodes,bc_tol)
        implicit none
        procedure(bvp_fun_real) :: fun !! vectorized differential equations
        procedure(bvp_bc_real) :: bc !! n+k boundary equations
        real(dp), intent(in) :: x(:) !! strictly increasing initial mesh, at least two nodes
        real(dp), intent(in) :: y(:,:) !! finite initial guess, shape (n,size(x))
        type(solve_bvp_result), intent(out) :: res !! solution and convergence diagnostics
        real(dp), intent(in), optional :: p(:) !! initial guess for unknown parameters
        real(dp), intent(in), optional :: s(:,:) !! singular term S*y/(x-x(1)), shape (n,n)
        procedure(bvp_fun_jac_real), optional :: fun_jac !! analytic ODE derivatives
        procedure(bvp_bc_jac_real), optional :: bc_jac !! analytic boundary derivatives
        real(dp), intent(in), optional :: tol !! positive RMS defect tolerance, default 1e-3
        integer, intent(in), optional :: max_nodes !! mesh budget, default 1000
        real(dp), intent(in), optional :: bc_tol !! positive absolute boundary tolerance
        procedure(bvp_fun_jac_real), pointer :: user_fun_jac
        procedure(bvp_bc_jac_real), pointer :: user_bc_jac
        real(dp), allocatable :: mesh(:),h(:),xm(:),rowmax(:),rms(:),newx(:)
        real(dp), allocatable :: v(:,:),par(:),f(:,:),mid(:,:),fm(:,:),col(:,:),br(:),r(:)
        real(dp), allocatable :: a(:,:,:),ap(:,:,:),am(:,:,:),amp(:,:,:)
        real(dp), allocatable :: left(:,:,:),right(:,:,:),param(:,:,:),ba(:,:),bb(:,:),bp(:,:)
        real(dp), allocatable :: step(:),trialstep(:),vnew(:,:),pnew(:),rnew(:),bnew(:)
        real(dp), allocatable :: ft(:,:),mt(:,:),fmt(:,:),ct(:,:),vmesh(:,:)
        real(dp), allocatable :: project(:,:),d(:,:),eye(:,:),sinv(:,:)
        type(ivp_sparse_matrix) :: jac
        type(ivp_sparse_lu_real) :: lu
        real(dp) :: tolerance,btol,alpha,cost,newcost,scale,tn,newtn
        integer :: n,m,k,budget,info,it,j,l,q,idx,pos,additions,outer,istat
        integer(int64) :: total64,nnz64,add64
        logical :: accepted,converged,valid,has_s,have_solution

        have_solution=.false.
        user_fun_jac=>null()
        user_bc_jac=>null()
        if (present(fun_jac)) user_fun_jac=>fun_jac
        if (present(bc_jac)) user_bc_jac=>bc_jac
        res%message='Invalid boundary-value problem arguments.'
        n=size(y,1)
        m=size(x)
        k=0
        if (present(p)) k=size(p)
        if (n<1 .or. m<2 .or. size(y,2)/=m) return
        if (.not.finite(x) .or. .not.finite(y)) return
        if (any(x(2:)<=x(:m-1))) return
        if (.not.finite(x(2:)-x(:m-1))) return
        if (present(p)) then
            if (.not.finite(p)) return
        end if
        tolerance=1e-3_dp
        if (present(tol)) tolerance=tol
        if (.not.ieee_is_finite(tolerance)) return
        if (tolerance<=0.0_dp) return
        tolerance=max(tolerance,100.0_dp*epsilon(1.0_dp))
        btol=tolerance
        if (present(bc_tol)) btol=bc_tol
        if (.not.ieee_is_finite(btol)) return
        if (btol<=0.0_dp) return
        budget=1000
        if (present(max_nodes)) budget=max_nodes
        if (budget<m) return
        total64=int(n,int64)*int(budget,int64)+int(k,int64)
        if (total64>=int(huge(1),int64)) return
        has_s=present(s)
        if (has_s) then
            if (any(shape(s)/=[n,n])) return
            if (.not.finite(s)) return
        end if
        allocate(mesh(m),v(n,m),par(k),eye(n,n),project(n,n),d(n,n),sinv(n,n),stat=info)
        if (info/=0) then
            res%status=bvp_status_setup_failure
            res%message='Could not allocate boundary-value setup arrays.'
            return
        end if
        mesh=x
        v=y
        if (present(p)) par=p
        eye=0.0_dp
        do j=1,n
            eye(j,j)=1.0_dp
        end do
        project=eye
        d=eye
        if (has_s) then
            call pinv(s,sinv,info)
            if (info/=0) then
                res%status=bvp_status_setup_failure
                res%message='Singular-term pseudoinverse did not converge.'
                return
            end if
            project=eye-matmul(sinv,s)
            call pinv(eye-s,d,info)
            if (info/=0) then
                res%status=bvp_status_setup_failure
                res%message='Singular-term endpoint pseudoinverse did not converge.'
                return
            end if
            v(:,1)=matmul(project,v(:,1))
        end if
        res%tol=tolerance
        res%status=bvp_status_numerical_failure
        outer=0
        do
            outer=outer+1
            res%niter=outer
            m=size(mesh)
            allocate(h(m-1),xm(m-1),f(n,m),mid(n,m-1),fm(n,m-1),col(n,m-1), &
                br(n+k),r(n*m+k),step(n*m+k),trialstep(n*m+k),rowmax(n*m+k), &
                vnew(n,m),pnew(k),rnew(n*m+k),bnew(n+k),ft(n,m),mt(n,m-1), &
                fmt(n,m-1),ct(n,m-1),rms(m-1),a(n,n,m),ap(n,k,m),am(n,n,m-1), &
                amp(n,k,m-1),left(n,n,m-1),right(n,n,m-1),param(n,k,m-1), &
                ba(n+k,n),bb(n+k,n),bp(n+k,k),stat=info)
            if (info/=0) then
                res%status=bvp_status_setup_failure
                exit
            end if
            f=ieee_value(0.0_dp,ieee_quiet_nan)
            br=ieee_value(0.0_dp,ieee_quiet_nan)
            rms=ieee_value(0.0_dp,ieee_quiet_nan)
            accepted=.false.
            h=mesh(2:)-mesh(:m-1)
            xm=mesh(:m-1)+0.5_dp*h
            have_solution=.false.
            call residual(v,par,f,mid,fm,col,br,r,valid)
            if (.not.valid) exit
            have_solution=.true.
            converged=.false.
            do it=1,8
                accepted=.false.
                call build_jacobian(valid)
                if (.not.valid) exit
                call sparse_lu_factor(jac,0.0_dp,1.0_dp,lu,info)
                if (info/=0) then
                    res%status=bvp_status_singular
                    exit
                end if
                call lu%solve(r/rowmax,step,info)
                if (info/=0 .or. .not.finite(step)) exit
                alpha=1.0_dp
                accepted=.false.
                do l=1,5
                    vnew=v-alpha*reshape(step(:n*m),[n,m])
                    pnew=par-alpha*step(n*m+1:)
                    if (has_s) vnew(:,1)=matmul(project,vnew(:,1))
                    call residual(vnew,pnew,ft,mt,fmt,ct,bnew,rnew,valid)
                    if (valid) then
                        call lu%solve(rnew/rowmax,trialstep,info)
                        if (info==0 .and. finite(trialstep)) then
                            scale=max(maxval(abs(step)),maxval(abs(trialstep)),tiny(1.0_dp))
                            cost=sum(abs(step/scale)**2)
                            newcost=sum(abs(trialstep/scale)**2)
                            if (newcost<(1.0_dp-0.4_dp*alpha)*cost .or. cost==0.0_dp) then
                                accepted=.true.
                                exit
                            end if
                            ! The last finite trial is retained, as in SciPy's Newton iteration.
                            if (l==5) accepted=.true.
                        end if
                    end if
                    if (l<5) alpha=alpha*0.5_dp
                end do
                if (.not.accepted) exit
                v=vnew
                par=pnew
                f=ft
                mid=mt
                fm=fmt
                col=ct
                br=bnew
                r=rnew
                converged=all(abs(col)<spread(h*(0.05_dp*tolerance*2.0_dp/3.0_dp),1,n)* &
                    (1.0_dp+abs(fm))) .and. maxval(abs(br))<=btol
                if (converged) exit
            end do
            if (res%status==bvp_status_singular) exit
            if (.not.accepted) exit
            call estimate_rms(valid)
            if (.not.valid) exit
            add64=int(count(rms>tolerance .and. rms/100.0_dp<tolerance),int64)+ &
                2_int64*count(rms/100.0_dp>=tolerance)
            if (add64>int(budget-m,int64)) then
                res%status=bvp_status_max_nodes
                exit
            end if
            additions=int(add64)
            if (additions==0) then
                if (maxval(abs(br))<=btol) then
                    res%status=bvp_status_success
                    exit
                else if (outer>=10) then
                    res%status=bvp_status_bc_tolerance
                    exit
                end if
            else
                if (additions>budget-m) then
                    res%status=bvp_status_max_nodes
                    exit
                end if
                allocate(newx(m+additions),vmesh(n,m+additions),stat=info)
                if (info/=0) then
                    res%status=bvp_status_setup_failure
                    exit
                end if
                idx=1
                newx(1)=mesh(1)
                vmesh(:,1)=v(:,1)
                do j=1,m-1
                    q=0
                    if (rms(j)>tolerance) q=1
                    if (rms(j)/100.0_dp>=tolerance) q=2
                    do l=1,q
                        tn=real(l,dp)/real(q+1,dp)
                        newtn=mesh(j)+tn*h(j)
                        if (newtn<=newx(idx) .or. newtn>=mesh(j+1)) exit
                        idx=idx+1
                        newx(idx)=newtn
                        call hermite_real(v(:,j),v(:,j+1),f(:,j),f(:,j+1),h(j),tn,0,vmesh(:,idx))
                    end do
                    if (l<=q) exit
                    idx=idx+1
                    newx(idx)=mesh(j+1)
                    vmesh(:,idx)=v(:,j+1)
                end do
                if (idx/=m+additions) exit
                call move_alloc(newx,mesh)
                call move_alloc(vmesh,v)
            end if
            deallocate(h,xm,f,mid,fm,col,br,r,step,trialstep,rowmax,vnew,pnew,rnew,bnew, &
                ft,mt,fmt,ct,rms,a,ap,am,amp,left,right,param,ba,bb,bp)
        end do
        res%success=res%status==bvp_status_success
        select case(res%status)
        case(bvp_status_success)
            res%message='Boundary-value solution converged.'
        case(bvp_status_max_nodes)
            res%message='Mesh refinement would exceed max_nodes.'
        case(bvp_status_singular)
            res%message='Collocation Jacobian is numerically singular.'
        case(bvp_status_bc_tolerance)
            res%message='Boundary tolerance not reached in ten iterations.'
        case(bvp_status_setup_failure)
            res%message='Could not allocate solver work arrays.'
        case default
            res%message='Nonfinite evaluation, failed Newton step, or mesh stagnation.'
        end select
        res%x=mesh
        res%y=v
        res%p=par
        if (have_solution .and. allocated(f)) then
            if (finite(f) .and. size(f,2)==size(mesh)) then
                res%yp=f
                res%sol%x=mesh
                res%sol%y=v
                res%sol%f=f
            end if
        end if
        if (allocated(br)) res%bc_residuals=br
        if (allocated(rms)) then
            if (res%status==bvp_status_success .or. res%status==bvp_status_max_nodes .or. &
                res%status==bvp_status_bc_tolerance) res%rms_residuals=rms
        end if
    contains
        subroutine rhs(xx,yy,pp,ff)
            implicit none
            real(dp), intent(in) :: xx(:) !! evaluation coordinates
            real(dp), intent(in) :: yy(:,:),pp(:) !! state and parameters
            real(dp), intent(out) :: ff(:,:) !! modified differential equations
            integer :: jj
            ff=ieee_value(0.0_dp,ieee_quiet_nan)
            call fun(xx,yy,pp,ff)
            res%nfev=res%nfev+1
            if (.not.finite(ff)) return
            if (has_s) then
                do jj=1,size(xx)
                    if (xx(jj)==x(1)) then
                        ff(:,jj)=matmul(d,ff(:,jj))
                    else
                        ff(:,jj)=ff(:,jj)+matmul(s,yy(:,jj))/(xx(jj)-x(1))
                    end if
                end do
            end if
        end subroutine rhs

        subroutine residual(yy,pp,ff,ym,ffm,cc,bbres,rr,ok)
            implicit none
            real(dp), intent(in) :: yy(:,:),pp(:) !! current iterate
            real(dp), intent(out) :: ff(:,:),ym(:,:),ffm(:,:) !! mesh and midpoint RHS work
            real(dp), intent(out) :: cc(:,:),bbres(:),rr(:) !! collocation and boundary residuals
            logical, intent(out) :: ok !! all callback outputs finite
            ok=.false.
            ff=ieee_value(0.0_dp,ieee_quiet_nan)
            bbres=ieee_value(0.0_dp,ieee_quiet_nan)
            rr=ieee_value(0.0_dp,ieee_quiet_nan)
            if (.not.finite(yy) .or. .not.finite(pp)) return
            call rhs(mesh,yy,pp,ff)
            if (.not.finite(ff)) return
            ym=0.5_dp*(yy(:,2:)+yy(:,:m-1))-0.125_dp*spread(h,1,n)* &
                (ff(:,2:)-ff(:,:m-1))
            if (.not.finite(ym)) return
            call rhs(xm,ym,pp,ffm)
            if (.not.finite(ffm)) return
            cc=yy(:,2:)-yy(:,:m-1)-spread(h/6.0_dp,1,n)* &
                (ff(:,:m-1)+ff(:,2:)+4.0_dp*ffm)
            call bc(yy(:,1),yy(:,m),pp,bbres)
            if (.not.finite(cc) .or. .not.finite(bbres)) return
            rr(:n*(m-1))=reshape(cc,[n*(m-1)])
            rr(n*(m-1)+1:)=bbres
            ok=.true.
        end subroutine residual

        subroutine derivatives(xx,yy,pp,ff,aa,aa_p)
            implicit none
            real(dp), intent(in) :: xx(:) !! coordinates
            real(dp), intent(in) :: yy(:,:),pp(:),ff(:,:) !! iterate and known RHS
            real(dp), intent(out) :: aa(:,:,:),aa_p(:,:,:) !! differential Jacobians
            real(dp) :: yt(size(yy,1),size(yy,2)),pt(size(pp)), &
                fout(size(yy,1),size(yy,2)),dh(size(yy,2)),ph
            integer :: jj,ll
            aa=ieee_value(0.0_dp,ieee_quiet_nan)
            aa_p=ieee_value(0.0_dp,ieee_quiet_nan)
            res%njev=res%njev+1
            if (associated(user_fun_jac)) then
                call user_fun_jac(xx,yy,pp,aa,aa_p)
                if (.not.finite(aa) .or. .not.finite(aa_p)) return
                if (has_s) then
                    do jj=1,size(xx)
                        if (xx(jj)==x(1)) then
                            aa(:,:,jj)=matmul(d,aa(:,:,jj))
                            aa_p(:,:,jj)=matmul(d,aa_p(:,:,jj))
                        else
                            aa(:,:,jj)=aa(:,:,jj)+s/(xx(jj)-x(1))
                        end if
                    end do
                end if
            else
                do ll=1,n
                    yt=yy
                    yt(ll,:)=yy(ll,:)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(yy(ll,:)))
                    dh=yt(ll,:)-yy(ll,:)
                    call rhs(xx,yt,pp,fout)
                    if (.not.finite(fout)) return
                    aa(:,ll,:)=(fout-ff)/spread(dh,1,n)
                end do
                do ll=1,k
                    pt=pp
                    pt(ll)=pp(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(pp(ll)))
                    ph=pt(ll)-pp(ll)
                    call rhs(xx,yy,pt,fout)
                    if (.not.finite(fout)) return
                    aa_p(:,ll,:)=(fout-ff)/ph
                end do
            end if
        end subroutine derivatives

        subroutine build_jacobian(ok)
            implicit none
            logical, intent(out) :: ok !! finite nonsingular assembly succeeded
            real(dp) :: ya(n),yb(n),pp(k),bout(n+k),dh
            integer :: jj,ll,ii,rr,cc,nnz
            ok=.false.
            call derivatives(mesh,v,par,f,a,ap)
            call derivatives(xm,mid,par,fm,am,amp)
            if (.not.finite(a) .or. .not.finite(ap)) return
            if (.not.finite(am) .or. .not.finite(amp)) return
            ba=ieee_value(0.0_dp,ieee_quiet_nan)
            bb=ieee_value(0.0_dp,ieee_quiet_nan)
            bp=ieee_value(0.0_dp,ieee_quiet_nan)
            if (associated(user_bc_jac)) then
                call user_bc_jac(v(:,1),v(:,m),par,ba,bb,bp)
            else
                do ll=1,n
                    ya=v(:,1)
                    ya(ll)=ya(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(ya(ll)))
                    dh=ya(ll)-v(ll,1)
                    bout=ieee_value(0.0_dp,ieee_quiet_nan)
                    call bc(ya,v(:,m),par,bout)
                    if (.not.finite(bout)) return
                    ba(:,ll)=(bout-br)/dh
                    yb=v(:,m)
                    yb(ll)=yb(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(yb(ll)))
                    dh=yb(ll)-v(ll,m)
                    bout=ieee_value(0.0_dp,ieee_quiet_nan)
                    call bc(v(:,1),yb,par,bout)
                    if (.not.finite(bout)) return
                    bb(:,ll)=(bout-br)/dh
                end do
                do ll=1,k
                    pp=par
                    pp(ll)=pp(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(pp(ll)))
                    dh=pp(ll)-par(ll)
                    bout=ieee_value(0.0_dp,ieee_quiet_nan)
                    call bc(v(:,1),v(:,m),pp,bout)
                    if (.not.finite(bout)) return
                    bp(:,ll)=(bout-br)/dh
                end do
            end if
            if (.not.finite(ba) .or. .not.finite(bb) .or. .not.finite(bp)) return
            do jj=1,m-1
                left(:,:,jj)=-eye-h(jj)/6.0_dp*(a(:,:,jj)+2.0_dp*am(:,:,jj))- &
                    matmul(h(jj)*am(:,:,jj),h(jj)*a(:,:,jj))/12.0_dp
                right(:,:,jj)=eye-h(jj)/6.0_dp*(a(:,:,jj+1)+2.0_dp*am(:,:,jj))+ &
                    matmul(h(jj)*am(:,:,jj),h(jj)*a(:,:,jj+1))/12.0_dp
                param(:,:,jj)=-h(jj)/6.0_dp*(ap(:,:,jj)+ap(:,:,jj+1)+ &
                    4.0_dp*(amp(:,:,jj)+matmul(h(jj)*am(:,:,jj), &
                    ap(:,:,jj)-ap(:,:,jj+1))/8.0_dp))
            end do
            nnz64=2_int64*(m-1)*int(n,int64)*n+2_int64*(n+k)*n+ &
                (int(n,int64)*m+k)*k
            if (nnz64>=huge(1)) then
                res%status=bvp_status_setup_failure
                return
            end if
            nnz=int(nnz64)
            if (allocated(jac%data)) deallocate(jac%data,jac%indices,jac%indptr)
            allocate(jac%data(nnz),jac%indices(nnz),jac%indptr(n*m+k+1),stat=istat)
            if (istat/=0) then
                res%status=bvp_status_setup_failure
                return
            end if
            jac%nrows=n*m+k
            jac%ncols=jac%nrows
            pos=1
            cc=0
            do jj=1,m
                do ll=1,n
                    cc=cc+1
                    jac%indptr(cc)=pos
                    if (jj>1) call insert_block((jj-2)*n,right(:,ll,jj-1))
                    if (jj<m) call insert_block((jj-1)*n,left(:,ll,jj))
                    if (jj==1) call insert_block((m-1)*n,ba(:,ll))
                    if (jj==m) call insert_block((m-1)*n,bb(:,ll))
                end do
            end do
            do ll=1,k
                cc=cc+1
                jac%indptr(cc)=pos
                do jj=1,m-1
                    call insert_block((jj-1)*n,param(:,ll,jj))
                end do
                call insert_block((m-1)*n,bp(:,ll))
            end do
            jac%indptr(cc+1)=pos
            if (.not.finite(jac%data)) return
            rowmax=0.0_dp
            do ii=1,nnz
                rr=jac%indices(ii)
                rowmax(rr)=max(rowmax(rr),abs(jac%data(ii)))
            end do
            if (any(rowmax==0.0_dp)) then
                res%status=bvp_status_singular
                return
            end if
            do ii=1,nnz
                jac%data(ii)=jac%data(ii)/rowmax(jac%indices(ii))
            end do
            ok=.true.
        end subroutine build_jacobian

        subroutine insert_block(offset,values)
            implicit none
            integer, intent(in) :: offset !! preceding residual rows
            real(dp), intent(in) :: values(:) !! one column block
            integer :: ii
            do ii=1,size(values)
                jac%indices(pos)=offset+ii
                jac%data(pos)=values(ii)
                pos=pos+1
            end do
        end subroutine insert_block

        subroutine estimate_rms(ok)
            implicit none
            logical, intent(out) :: ok !! quadrature and callbacks finite
            real(dp) :: yy(n,m-1),dy(n,m-1),ff(n,m-1),defect(n,m-1)
            real(dp) :: xx(m-1),t,accum(m-1)
            integer :: jj,ll
            ok=.false.
            accum=32.0_dp/45.0_dp*sum(abs(1.5_dp*col/ &
                spread(h,1,n)/(1.0_dp+abs(fm)))**2,dim=1)
            do ll=1,2
                t=0.5_dp+(2*ll-3)*0.5_dp*sqrt(3.0_dp/7.0_dp)
                xx=mesh(:m-1)+t*h
                do jj=1,m-1
                    call hermite_real(v(:,jj),v(:,jj+1),f(:,jj),f(:,jj+1),h(jj),t,0,yy(:,jj))
                    call hermite_real(v(:,jj),v(:,jj+1),f(:,jj),f(:,jj+1),h(jj),t,1,dy(:,jj))
                end do
                if (.not.finite(yy) .or. .not.finite(dy)) return
                call rhs(xx,yy,par,ff)
                if (.not.finite(ff)) return
                defect=(dy-ff)/(1.0_dp+abs(ff))
                accum=accum+49.0_dp/90.0_dp*sum(abs(defect)**2,dim=1)
            end do
            rms=sqrt(0.5_dp*accum)
            ok=finite(rms)
        end subroutine estimate_rms
    end subroutine solve_bvp_real

    subroutine hermite_real(y0,y1,f0,f1,h,t,nu,out)
        implicit none
        real(dp), intent(in) :: y0(:),y1(:),f0(:),f1(:) !! endpoint states and derivatives
        real(dp), intent(in) :: h,t !! interval width and normalized coordinate
        integer, intent(in) :: nu !! derivative order, nonnegative
        real(dp), intent(out) :: out(:) !! interpolated derivative
        real(dp) :: b(size(y0)),c(size(y0))
        b=3.0_dp*(y1-y0)-h*(2.0_dp*f0+f1)
        c=2.0_dp*(y0-y1)+h*(f0+f1)
        select case(nu)
        case(0)
            out=y0+t*(h*f0+t*(b+t*c))
        case(1)
            out=f0+t*(2.0_dp*b+3.0_dp*t*c)/h
        case(2)
            out=(2.0_dp*b+6.0_dp*t*c)/h/h
        case(3)
            out=6.0_dp*c/h/h/h
        case default
            out=0.0_dp
        end select
    end subroutine hermite_real

    subroutine evaluate_scalar_real(self,x,y,status,nu,extrapolate)
        implicit none
        class(bvp_solution), intent(in) :: self !! stored dense solution
        real(dp), intent(in) :: x !! query coordinate
        real(dp), intent(out) :: y(:) !! state or derivative, size n
        integer, intent(out), optional :: status !! zero success, -1 invalid query
        integer, intent(in), optional :: nu !! derivative order, default zero
        logical, intent(in), optional :: extrapolate !! extend endpoint cubics, default true
        real(dp) :: yy(size(y),1)
        integer :: info
        call evaluate_vector_real(self,[x],yy,info,nu,extrapolate)
        y=yy(:,1)
        if (present(status)) status=info
    end subroutine evaluate_scalar_real

    subroutine evaluate_vector_real(self,x,y,status,nu,extrapolate)
        implicit none
        class(bvp_solution), intent(in) :: self !! stored dense solution
        real(dp), intent(in) :: x(:) !! query coordinates in any order
        real(dp), intent(out) :: y(:,:) !! state or derivatives, shape (n,size(x))
        integer, intent(out), optional :: status !! zero success, -1 invalid query
        integer, intent(in), optional :: nu !! nonnegative derivative order, default zero
        logical, intent(in), optional :: extrapolate !! extend endpoint cubics, default true
        integer :: j,lo,hi,mid,order,m
        logical :: extend
        real(dp) :: h,t
        y=ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status=bvp_status_invalid_argument
        if (.not.allocated(self%x)) return
        if (size(y,1)/=size(self%y,1) .or. size(y,2)/=size(x)) return
        if (.not.finite(x)) return
        order=0
        if (present(nu)) order=nu
        if (order<0) return
        extend=.true.
        if (present(extrapolate)) extend=extrapolate
        m=size(self%x)
        if (.not.extend) then
            if (any(x<self%x(1)) .or. any(x>self%x(m))) return
        end if
        do j=1,size(x)
            lo=1
            hi=m
            do while(hi-lo>1)
                mid=lo+(hi-lo)/2
                if (x(j)<self%x(mid)) then
                    hi=mid
                else
                    lo=mid
                end if
            end do
            h=self%x(lo+1)-self%x(lo)
            t=(x(j)-self%x(lo))/h
            call hermite_real(self%y(:,lo),self%y(:,lo+1),self%f(:,lo), &
                self%f(:,lo+1),h,t,order,y(:,j))
        end do
        if (present(status)) status=bvp_status_success
    end subroutine evaluate_vector_real

    subroutine solve_bvp_complex(fun,bc,x,y,res,p,s,fun_jac,bc_jac,tol,max_nodes,bc_tol)
        implicit none
        procedure(bvp_fun_complex) :: fun !! vectorized differential equations
        procedure(bvp_bc_complex) :: bc !! n+k boundary equations
        real(dp), intent(in) :: x(:) !! strictly increasing initial mesh, at least two nodes
        complex(dp), intent(in) :: y(:,:) !! finite initial guess, shape (n,size(x))
        type(solve_bvp_complex_result), intent(out) :: res !! solution and convergence diagnostics
        complex(dp), intent(in), optional :: p(:) !! initial guess for unknown parameters
        complex(dp), intent(in), optional :: s(:,:) !! singular term S*y/(x-x(1)), shape (n,n)
        procedure(bvp_fun_jac_complex), optional :: fun_jac !! analytic ODE derivatives
        procedure(bvp_bc_jac_complex), optional :: bc_jac !! analytic boundary derivatives
        real(dp), intent(in), optional :: tol !! positive RMS defect tolerance, default 1e-3
        integer, intent(in), optional :: max_nodes !! mesh budget, default 1000
        real(dp), intent(in), optional :: bc_tol !! positive absolute boundary tolerance
        procedure(bvp_fun_jac_complex), pointer :: user_fun_jac
        procedure(bvp_bc_jac_complex), pointer :: user_bc_jac
        real(dp), allocatable :: mesh(:),h(:),xm(:),rowmax(:),rms(:),newx(:)
        complex(dp), allocatable :: v(:,:),par(:),f(:,:),mid(:,:),fm(:,:),col(:,:),br(:),r(:)
        complex(dp), allocatable :: a(:,:,:),ap(:,:,:),am(:,:,:),amp(:,:,:)
        complex(dp), allocatable :: left(:,:,:),right(:,:,:),param(:,:,:),ba(:,:),bb(:,:),bp(:,:)
        complex(dp), allocatable :: step(:),trialstep(:),vnew(:,:),pnew(:),rnew(:),bnew(:)
        complex(dp), allocatable :: ft(:,:),mt(:,:),fmt(:,:),ct(:,:),vmesh(:,:)
        complex(dp), allocatable :: project(:,:),d(:,:),eye(:,:),sinv(:,:)
        type(ivp_complex_sparse_matrix) :: jac
        type(ivp_sparse_lu_complex) :: lu
        real(dp) :: tolerance,btol,alpha,cost,newcost,scale,tn,newtn
        integer :: n,m,k,budget,info,it,j,l,q,idx,pos,additions,outer,istat
        integer(int64) :: total64,nnz64,add64
        logical :: accepted,converged,valid,has_s,have_solution

        have_solution=.false.
        user_fun_jac=>null()
        user_bc_jac=>null()
        if (present(fun_jac)) user_fun_jac=>fun_jac
        if (present(bc_jac)) user_bc_jac=>bc_jac
        res%message='Invalid boundary-value problem arguments.'
        n=size(y,1)
        m=size(x)
        k=0
        if (present(p)) k=size(p)
        if (n<1 .or. m<2 .or. size(y,2)/=m) return
        if (.not.finite(x) .or. .not.finite(y)) return
        if (any(x(2:)<=x(:m-1))) return
        if (.not.finite(x(2:)-x(:m-1))) return
        if (present(p)) then
            if (.not.finite(p)) return
        end if
        tolerance=1e-3_dp
        if (present(tol)) tolerance=tol
        if (.not.ieee_is_finite(tolerance)) return
        if (tolerance<=0.0_dp) return
        tolerance=max(tolerance,100.0_dp*epsilon(1.0_dp))
        btol=tolerance
        if (present(bc_tol)) btol=bc_tol
        if (.not.ieee_is_finite(btol)) return
        if (btol<=0.0_dp) return
        budget=1000
        if (present(max_nodes)) budget=max_nodes
        if (budget<m) return
        total64=int(n,int64)*int(budget,int64)+int(k,int64)
        if (total64>=int(huge(1),int64)) return
        has_s=present(s)
        if (has_s) then
            if (any(shape(s)/=[n,n])) return
            if (.not.finite(s)) return
        end if
        allocate(mesh(m),v(n,m),par(k),eye(n,n),project(n,n),d(n,n),sinv(n,n),stat=info)
        if (info/=0) then
            res%status=bvp_status_setup_failure
            res%message='Could not allocate boundary-value setup arrays.'
            return
        end if
        mesh=x
        v=y
        if (present(p)) par=p
        eye=cmplx(0.0_dp,0.0_dp,dp)
        do j=1,n
            eye(j,j)=cmplx(1.0_dp,0.0_dp,dp)
        end do
        project=eye
        d=eye
        if (has_s) then
            call pinv(s,sinv,info)
            if (info/=0) then
                res%status=bvp_status_setup_failure
                res%message='Singular-term pseudoinverse did not converge.'
                return
            end if
            project=eye-matmul(sinv,s)
            call pinv(eye-s,d,info)
            if (info/=0) then
                res%status=bvp_status_setup_failure
                res%message='Singular-term endpoint pseudoinverse did not converge.'
                return
            end if
            v(:,1)=matmul(project,v(:,1))
        end if
        res%tol=tolerance
        res%status=bvp_status_numerical_failure
        outer=0
        do
            outer=outer+1
            res%niter=outer
            m=size(mesh)
            allocate(h(m-1),xm(m-1),f(n,m),mid(n,m-1),fm(n,m-1),col(n,m-1), &
                br(n+k),r(n*m+k),step(n*m+k),trialstep(n*m+k),rowmax(n*m+k), &
                vnew(n,m),pnew(k),rnew(n*m+k),bnew(n+k),ft(n,m),mt(n,m-1), &
                fmt(n,m-1),ct(n,m-1),rms(m-1),a(n,n,m),ap(n,k,m),am(n,n,m-1), &
                amp(n,k,m-1),left(n,n,m-1),right(n,n,m-1),param(n,k,m-1), &
                ba(n+k,n),bb(n+k,n),bp(n+k,k),stat=info)
            if (info/=0) then
                res%status=bvp_status_setup_failure
                exit
            end if
            f=ieee_value(0.0_dp,ieee_quiet_nan)
            br=ieee_value(0.0_dp,ieee_quiet_nan)
            rms=ieee_value(0.0_dp,ieee_quiet_nan)
            accepted=.false.
            h=mesh(2:)-mesh(:m-1)
            xm=mesh(:m-1)+0.5_dp*h
            have_solution=.false.
            call residual(v,par,f,mid,fm,col,br,r,valid)
            if (.not.valid) exit
            have_solution=.true.
            converged=.false.
            do it=1,8
                accepted=.false.
                call build_jacobian(valid)
                if (.not.valid) exit
                call sparse_lu_factor(jac,cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp),lu,info)
                if (info/=0) then
                    res%status=bvp_status_singular
                    exit
                end if
                call lu%solve(r/rowmax,step,info)
                if (info/=0 .or. .not.finite(step)) exit
                alpha=1.0_dp
                accepted=.false.
                do l=1,5
                    vnew=v-alpha*reshape(step(:n*m),[n,m])
                    pnew=par-alpha*step(n*m+1:)
                    if (has_s) vnew(:,1)=matmul(project,vnew(:,1))
                    call residual(vnew,pnew,ft,mt,fmt,ct,bnew,rnew,valid)
                    if (valid) then
                        call lu%solve(rnew/rowmax,trialstep,info)
                        if (info==0 .and. finite(trialstep)) then
                            scale=max(maxval(abs(step)),maxval(abs(trialstep)),tiny(1.0_dp))
                            cost=sum(abs(step/scale)**2)
                            newcost=sum(abs(trialstep/scale)**2)
                            if (newcost<(1.0_dp-0.4_dp*alpha)*cost .or. cost==0.0_dp) then
                                accepted=.true.
                                exit
                            end if
                            ! The last finite trial is retained, as in SciPy's Newton iteration.
                            if (l==5) accepted=.true.
                        end if
                    end if
                    if (l<5) alpha=alpha*0.5_dp
                end do
                if (.not.accepted) exit
                v=vnew
                par=pnew
                f=ft
                mid=mt
                fm=fmt
                col=ct
                br=bnew
                r=rnew
                converged=all(abs(col)<spread(h*(0.05_dp*tolerance*2.0_dp/3.0_dp),1,n)* &
                    (1.0_dp+abs(fm))) .and. maxval(abs(br))<=btol
                if (converged) exit
            end do
            if (res%status==bvp_status_singular) exit
            if (.not.accepted) exit
            call estimate_rms(valid)
            if (.not.valid) exit
            add64=int(count(rms>tolerance .and. rms/100.0_dp<tolerance),int64)+ &
                2_int64*count(rms/100.0_dp>=tolerance)
            if (add64>int(budget-m,int64)) then
                res%status=bvp_status_max_nodes
                exit
            end if
            additions=int(add64)
            if (additions==0) then
                if (maxval(abs(br))<=btol) then
                    res%status=bvp_status_success
                    exit
                else if (outer>=10) then
                    res%status=bvp_status_bc_tolerance
                    exit
                end if
            else
                if (additions>budget-m) then
                    res%status=bvp_status_max_nodes
                    exit
                end if
                allocate(newx(m+additions),vmesh(n,m+additions),stat=info)
                if (info/=0) then
                    res%status=bvp_status_setup_failure
                    exit
                end if
                idx=1
                newx(1)=mesh(1)
                vmesh(:,1)=v(:,1)
                do j=1,m-1
                    q=0
                    if (rms(j)>tolerance) q=1
                    if (rms(j)/100.0_dp>=tolerance) q=2
                    do l=1,q
                        tn=real(l,dp)/real(q+1,dp)
                        newtn=mesh(j)+tn*h(j)
                        if (newtn<=newx(idx) .or. newtn>=mesh(j+1)) exit
                        idx=idx+1
                        newx(idx)=newtn
                        call hermite_complex(v(:,j),v(:,j+1),f(:,j),f(:,j+1),h(j),tn,0,vmesh(:,idx))
                    end do
                    if (l<=q) exit
                    idx=idx+1
                    newx(idx)=mesh(j+1)
                    vmesh(:,idx)=v(:,j+1)
                end do
                if (idx/=m+additions) exit
                call move_alloc(newx,mesh)
                call move_alloc(vmesh,v)
            end if
            deallocate(h,xm,f,mid,fm,col,br,r,step,trialstep,rowmax,vnew,pnew,rnew,bnew, &
                ft,mt,fmt,ct,rms,a,ap,am,amp,left,right,param,ba,bb,bp)
        end do
        res%success=res%status==bvp_status_success
        select case(res%status)
        case(bvp_status_success)
            res%message='Boundary-value solution converged.'
        case(bvp_status_max_nodes)
            res%message='Mesh refinement would exceed max_nodes.'
        case(bvp_status_singular)
            res%message='Collocation Jacobian is numerically singular.'
        case(bvp_status_bc_tolerance)
            res%message='Boundary tolerance not reached in ten iterations.'
        case(bvp_status_setup_failure)
            res%message='Could not allocate solver work arrays.'
        case default
            res%message='Nonfinite evaluation, failed Newton step, or mesh stagnation.'
        end select
        res%x=mesh
        res%y=v
        res%p=par
        if (have_solution .and. allocated(f)) then
            if (finite(f) .and. size(f,2)==size(mesh)) then
                res%yp=f
                res%sol%x=mesh
                res%sol%y=v
                res%sol%f=f
            end if
        end if
        if (allocated(br)) res%bc_residuals=br
        if (allocated(rms)) then
            if (res%status==bvp_status_success .or. res%status==bvp_status_max_nodes .or. &
                res%status==bvp_status_bc_tolerance) res%rms_residuals=rms
        end if
    contains
        subroutine rhs(xx,yy,pp,ff)
            implicit none
            real(dp), intent(in) :: xx(:) !! evaluation coordinates
            complex(dp), intent(in) :: yy(:,:),pp(:) !! state and parameters
            complex(dp), intent(out) :: ff(:,:) !! modified differential equations
            integer :: jj
            ff=ieee_value(0.0_dp,ieee_quiet_nan)
            call fun(xx,yy,pp,ff)
            res%nfev=res%nfev+1
            if (.not.finite(ff)) return
            if (has_s) then
                do jj=1,size(xx)
                    if (xx(jj)==x(1)) then
                        ff(:,jj)=matmul(d,ff(:,jj))
                    else
                        ff(:,jj)=ff(:,jj)+matmul(s,yy(:,jj))/(xx(jj)-x(1))
                    end if
                end do
            end if
        end subroutine rhs

        subroutine residual(yy,pp,ff,ym,ffm,cc,bbres,rr,ok)
            implicit none
            complex(dp), intent(in) :: yy(:,:),pp(:) !! current iterate
            complex(dp), intent(out) :: ff(:,:),ym(:,:),ffm(:,:) !! mesh and midpoint RHS work
            complex(dp), intent(out) :: cc(:,:),bbres(:),rr(:) !! collocation and boundary residuals
            logical, intent(out) :: ok !! all callback outputs finite
            ok=.false.
            ff=ieee_value(0.0_dp,ieee_quiet_nan)
            bbres=ieee_value(0.0_dp,ieee_quiet_nan)
            rr=ieee_value(0.0_dp,ieee_quiet_nan)
            if (.not.finite(yy) .or. .not.finite(pp)) return
            call rhs(mesh,yy,pp,ff)
            if (.not.finite(ff)) return
            ym=0.5_dp*(yy(:,2:)+yy(:,:m-1))-0.125_dp*spread(h,1,n)* &
                (ff(:,2:)-ff(:,:m-1))
            if (.not.finite(ym)) return
            call rhs(xm,ym,pp,ffm)
            if (.not.finite(ffm)) return
            cc=yy(:,2:)-yy(:,:m-1)-spread(h/6.0_dp,1,n)* &
                (ff(:,:m-1)+ff(:,2:)+4.0_dp*ffm)
            call bc(yy(:,1),yy(:,m),pp,bbres)
            if (.not.finite(cc) .or. .not.finite(bbres)) return
            rr(:n*(m-1))=reshape(cc,[n*(m-1)])
            rr(n*(m-1)+1:)=bbres
            ok=.true.
        end subroutine residual

        subroutine derivatives(xx,yy,pp,ff,aa,aa_p)
            implicit none
            real(dp), intent(in) :: xx(:) !! coordinates
            complex(dp), intent(in) :: yy(:,:),pp(:),ff(:,:) !! iterate and known RHS
            complex(dp), intent(out) :: aa(:,:,:),aa_p(:,:,:) !! differential Jacobians
            complex(dp) :: yt(size(yy,1),size(yy,2)),pt(size(pp)), &
                fout(size(yy,1),size(yy,2)),dh(size(yy,2)),ph
            integer :: jj,ll
            aa=ieee_value(0.0_dp,ieee_quiet_nan)
            aa_p=ieee_value(0.0_dp,ieee_quiet_nan)
            res%njev=res%njev+1
            if (associated(user_fun_jac)) then
                call user_fun_jac(xx,yy,pp,aa,aa_p)
                if (.not.finite(aa) .or. .not.finite(aa_p)) return
                if (has_s) then
                    do jj=1,size(xx)
                        if (xx(jj)==x(1)) then
                            aa(:,:,jj)=matmul(d,aa(:,:,jj))
                            aa_p(:,:,jj)=matmul(d,aa_p(:,:,jj))
                        else
                            aa(:,:,jj)=aa(:,:,jj)+s/(xx(jj)-x(1))
                        end if
                    end do
                end if
            else
                do ll=1,n
                    yt=yy
                    yt(ll,:)=yy(ll,:)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(yy(ll,:)))
                    dh=yt(ll,:)-yy(ll,:)
                    call rhs(xx,yt,pp,fout)
                    if (.not.finite(fout)) return
                    aa(:,ll,:)=(fout-ff)/spread(dh,1,n)
                end do
                do ll=1,k
                    pt=pp
                    pt(ll)=pp(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(pp(ll)))
                    ph=pt(ll)-pp(ll)
                    call rhs(xx,yy,pt,fout)
                    if (.not.finite(fout)) return
                    aa_p(:,ll,:)=(fout-ff)/ph
                end do
            end if
        end subroutine derivatives

        subroutine build_jacobian(ok)
            implicit none
            logical, intent(out) :: ok !! finite nonsingular assembly succeeded
            complex(dp) :: ya(n),yb(n),pp(k),bout(n+k),dh
            integer :: jj,ll,ii,rr,cc,nnz
            ok=.false.
            call derivatives(mesh,v,par,f,a,ap)
            call derivatives(xm,mid,par,fm,am,amp)
            if (.not.finite(a) .or. .not.finite(ap)) return
            if (.not.finite(am) .or. .not.finite(amp)) return
            ba=ieee_value(0.0_dp,ieee_quiet_nan)
            bb=ieee_value(0.0_dp,ieee_quiet_nan)
            bp=ieee_value(0.0_dp,ieee_quiet_nan)
            if (associated(user_bc_jac)) then
                call user_bc_jac(v(:,1),v(:,m),par,ba,bb,bp)
            else
                do ll=1,n
                    ya=v(:,1)
                    ya(ll)=ya(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(ya(ll)))
                    dh=ya(ll)-v(ll,1)
                    bout=ieee_value(0.0_dp,ieee_quiet_nan)
                    call bc(ya,v(:,m),par,bout)
                    if (.not.finite(bout)) return
                    ba(:,ll)=(bout-br)/dh
                    yb=v(:,m)
                    yb(ll)=yb(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(yb(ll)))
                    dh=yb(ll)-v(ll,m)
                    bout=ieee_value(0.0_dp,ieee_quiet_nan)
                    call bc(v(:,1),yb,par,bout)
                    if (.not.finite(bout)) return
                    bb(:,ll)=(bout-br)/dh
                end do
                do ll=1,k
                    pp=par
                    pp(ll)=pp(ll)+sqrt(epsilon(1.0_dp))*(1.0_dp+abs(pp(ll)))
                    dh=pp(ll)-par(ll)
                    bout=ieee_value(0.0_dp,ieee_quiet_nan)
                    call bc(v(:,1),v(:,m),pp,bout)
                    if (.not.finite(bout)) return
                    bp(:,ll)=(bout-br)/dh
                end do
            end if
            if (.not.finite(ba) .or. .not.finite(bb) .or. .not.finite(bp)) return
            do jj=1,m-1
                left(:,:,jj)=-eye-h(jj)/6.0_dp*(a(:,:,jj)+2.0_dp*am(:,:,jj))- &
                    matmul(h(jj)*am(:,:,jj),h(jj)*a(:,:,jj))/12.0_dp
                right(:,:,jj)=eye-h(jj)/6.0_dp*(a(:,:,jj+1)+2.0_dp*am(:,:,jj))+ &
                    matmul(h(jj)*am(:,:,jj),h(jj)*a(:,:,jj+1))/12.0_dp
                param(:,:,jj)=-h(jj)/6.0_dp*(ap(:,:,jj)+ap(:,:,jj+1)+ &
                    4.0_dp*(amp(:,:,jj)+matmul(h(jj)*am(:,:,jj), &
                    ap(:,:,jj)-ap(:,:,jj+1))/8.0_dp))
            end do
            nnz64=2_int64*(m-1)*int(n,int64)*n+2_int64*(n+k)*n+ &
                (int(n,int64)*m+k)*k
            if (nnz64>=huge(1)) then
                res%status=bvp_status_setup_failure
                return
            end if
            nnz=int(nnz64)
            if (allocated(jac%data)) deallocate(jac%data,jac%indices,jac%indptr)
            allocate(jac%data(nnz),jac%indices(nnz),jac%indptr(n*m+k+1),stat=istat)
            if (istat/=0) then
                res%status=bvp_status_setup_failure
                return
            end if
            jac%nrows=n*m+k
            jac%ncols=jac%nrows
            pos=1
            cc=0
            do jj=1,m
                do ll=1,n
                    cc=cc+1
                    jac%indptr(cc)=pos
                    if (jj>1) call insert_block((jj-2)*n,right(:,ll,jj-1))
                    if (jj<m) call insert_block((jj-1)*n,left(:,ll,jj))
                    if (jj==1) call insert_block((m-1)*n,ba(:,ll))
                    if (jj==m) call insert_block((m-1)*n,bb(:,ll))
                end do
            end do
            do ll=1,k
                cc=cc+1
                jac%indptr(cc)=pos
                do jj=1,m-1
                    call insert_block((jj-1)*n,param(:,ll,jj))
                end do
                call insert_block((m-1)*n,bp(:,ll))
            end do
            jac%indptr(cc+1)=pos
            if (.not.finite(jac%data)) return
            rowmax=0.0_dp
            do ii=1,nnz
                rr=jac%indices(ii)
                rowmax(rr)=max(rowmax(rr),abs(jac%data(ii)))
            end do
            if (any(rowmax==0.0_dp)) then
                res%status=bvp_status_singular
                return
            end if
            do ii=1,nnz
                jac%data(ii)=jac%data(ii)/rowmax(jac%indices(ii))
            end do
            ok=.true.
        end subroutine build_jacobian

        subroutine insert_block(offset,values)
            implicit none
            integer, intent(in) :: offset !! preceding residual rows
            complex(dp), intent(in) :: values(:) !! one column block
            integer :: ii
            do ii=1,size(values)
                jac%indices(pos)=offset+ii
                jac%data(pos)=values(ii)
                pos=pos+1
            end do
        end subroutine insert_block

        subroutine estimate_rms(ok)
            implicit none
            logical, intent(out) :: ok !! quadrature and callbacks finite
            complex(dp) :: yy(n,m-1),dy(n,m-1),ff(n,m-1),defect(n,m-1)
            real(dp) :: xx(m-1),t,accum(m-1)
            integer :: jj,ll
            ok=.false.
            accum=32.0_dp/45.0_dp*sum(abs(1.5_dp*col/ &
                spread(h,1,n)/(1.0_dp+abs(fm)))**2,dim=1)
            do ll=1,2
                t=0.5_dp+(2*ll-3)*0.5_dp*sqrt(3.0_dp/7.0_dp)
                xx=mesh(:m-1)+t*h
                do jj=1,m-1
                    call hermite_complex(v(:,jj),v(:,jj+1),f(:,jj),f(:,jj+1),h(jj),t,0,yy(:,jj))
                    call hermite_complex(v(:,jj),v(:,jj+1),f(:,jj),f(:,jj+1),h(jj),t,1,dy(:,jj))
                end do
                if (.not.finite(yy) .or. .not.finite(dy)) return
                call rhs(xx,yy,par,ff)
                if (.not.finite(ff)) return
                defect=(dy-ff)/(1.0_dp+abs(ff))
                accum=accum+49.0_dp/90.0_dp*sum(abs(defect)**2,dim=1)
            end do
            rms=sqrt(0.5_dp*accum)
            ok=finite(rms)
        end subroutine estimate_rms
    end subroutine solve_bvp_complex

    subroutine hermite_complex(y0,y1,f0,f1,h,t,nu,out)
        implicit none
        complex(dp), intent(in) :: y0(:),y1(:),f0(:),f1(:) !! endpoint states and derivatives
        real(dp), intent(in) :: h,t !! interval width and normalized coordinate
        integer, intent(in) :: nu !! derivative order, nonnegative
        complex(dp), intent(out) :: out(:) !! interpolated derivative
        complex(dp) :: b(size(y0)),c(size(y0))
        b=3.0_dp*(y1-y0)-h*(2.0_dp*f0+f1)
        c=2.0_dp*(y0-y1)+h*(f0+f1)
        select case(nu)
        case(0)
            out=y0+t*(h*f0+t*(b+t*c))
        case(1)
            out=f0+t*(2.0_dp*b+3.0_dp*t*c)/h
        case(2)
            out=(2.0_dp*b+6.0_dp*t*c)/h/h
        case(3)
            out=6.0_dp*c/h/h/h
        case default
            out=cmplx(0.0_dp,0.0_dp,dp)
        end select
    end subroutine hermite_complex

    subroutine evaluate_scalar_complex(self,x,y,status,nu,extrapolate)
        implicit none
        class(bvp_complex_solution), intent(in) :: self !! stored dense solution
        real(dp), intent(in) :: x !! query coordinate
        complex(dp), intent(out) :: y(:) !! state or derivative, size n
        integer, intent(out), optional :: status !! zero success, -1 invalid query
        integer, intent(in), optional :: nu !! derivative order, default zero
        logical, intent(in), optional :: extrapolate !! extend endpoint cubics, default true
        complex(dp) :: yy(size(y),1)
        integer :: info
        call evaluate_vector_complex(self,[x],yy,info,nu,extrapolate)
        y=yy(:,1)
        if (present(status)) status=info
    end subroutine evaluate_scalar_complex

    subroutine evaluate_vector_complex(self,x,y,status,nu,extrapolate)
        implicit none
        class(bvp_complex_solution), intent(in) :: self !! stored dense solution
        real(dp), intent(in) :: x(:) !! query coordinates in any order
        complex(dp), intent(out) :: y(:,:) !! state or derivatives, shape (n,size(x))
        integer, intent(out), optional :: status !! zero success, -1 invalid query
        integer, intent(in), optional :: nu !! nonnegative derivative order, default zero
        logical, intent(in), optional :: extrapolate !! extend endpoint cubics, default true
        integer :: j,lo,hi,mid,order,m
        logical :: extend
        real(dp) :: h,t
        y=ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status=bvp_status_invalid_argument
        if (.not.allocated(self%x)) return
        if (size(y,1)/=size(self%y,1) .or. size(y,2)/=size(x)) return
        if (.not.finite(x)) return
        order=0
        if (present(nu)) order=nu
        if (order<0) return
        extend=.true.
        if (present(extrapolate)) extend=extrapolate
        m=size(self%x)
        if (.not.extend) then
            if (any(x<self%x(1)) .or. any(x>self%x(m))) return
        end if
        do j=1,size(x)
            lo=1
            hi=m
            do while(hi-lo>1)
                mid=lo+(hi-lo)/2
                if (x(j)<self%x(mid)) then
                    hi=mid
                else
                    lo=mid
                end if
            end do
            h=self%x(lo+1)-self%x(lo)
            t=(x(j)-self%x(lo))/h
            call hermite_complex(self%y(:,lo),self%y(:,lo+1),self%f(:,lo), &
                self%f(:,lo+1),h,t,order,y(:,j))
        end do
        if (present(status)) status=bvp_status_success
    end subroutine evaluate_vector_complex

    logical function finite_r1(a) result(ok)
        implicit none
        real(dp), intent(in) :: a(:) !! values to check for NaN or infinity
        ok=all(ieee_is_finite(a))
    end function finite_r1

    logical function finite_r2(a) result(ok)
        implicit none
        real(dp), intent(in) :: a(:,:) !! values to check for NaN or infinity
        ok=all(ieee_is_finite(a))
    end function finite_r2

    logical function finite_r3(a) result(ok)
        implicit none
        real(dp), intent(in) :: a(:,:,:) !! values to check for NaN or infinity
        ok=all(ieee_is_finite(a))
    end function finite_r3

    logical function finite_c1(a) result(ok)
        implicit none
        complex(dp), intent(in) :: a(:) !! values to check for NaN or infinity
        ok=all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_c1

    logical function finite_c2(a) result(ok)
        implicit none
        complex(dp), intent(in) :: a(:,:) !! values to check for NaN or infinity
        ok=all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_c2

    logical function finite_c3(a) result(ok)
        implicit none
        complex(dp), intent(in) :: a(:,:,:) !! values to check for NaN or infinity
        ok=all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_c3
end module scifort_integrate_bvp
