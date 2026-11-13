! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2017, The Chancellor, Masters and Scholars of the University
! of Oxford, and the Chebfun Developers. All rights reserved.
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! AAA selection, Loewner fitting, arrowhead spectra and cleanup adapted from
! SciPy v1.18.1 interpolate/_bary_rational.py. Full terms:
! LICENSES/scipy-aaa-BSD-3-Clause.txt. See CODE_PROVENANCE.md for modifications.
module scifort_interpolate_aaa
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_is_nan, &
        ieee_value,ieee_quiet_nan,ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_linalg_svd, only : svd
    use scifort_linalg_generalized_eigen, only : eigvals
    use scifort_linalg_matrix_functions, only : matrix_function_max_order
    implicit none
    private
    integer, parameter, public :: aaa_success=0,aaa_invalid_input=1,aaa_max_terms=2
    integer, parameter, public :: aaa_numerical_failure=3,aaa_spectral_size_limit=4
    integer, parameter, public :: aaa_setup_failure=5,aaa_undefined=6
    public :: aaa

    type :: aaa_core
        complex(dp), allocatable :: points(:),samples(:),z(:),f(:),s(:),v(:),w(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: errors(:)
        complex(dp) :: center=cmplx(0.0_dp,0.0_dp,dp)
        real(dp) :: xscale=1.0_dp,yscale=1.0_dp,rtol=0.0_dp,error=0.0_dp
        integer :: fit_status=aaa_invalid_input,cleanup_status=aaa_success
        logical :: converged=.false.
    end type aaa_core

    type, public :: aaa_type
        private
        type(aaa_core) :: core
    contains
        procedure :: valid => real_valid
        procedure :: converged => real_converged
        procedure :: fit_status => real_fit_status
        procedure :: cleanup_status => real_cleanup_status
        procedure :: max_error => real_max_error
        procedure :: support_points => real_support_points
        procedure :: support_values => real_support_values
        procedure :: weights => real_weights
        procedure :: errors => real_errors
        procedure, private :: value_real => real_value_real
        procedure, private :: value_complex => real_value_complex
        generic :: value => value_real,value_complex
        procedure, private :: values_real => real_values_real
        procedure, private :: values_complex => real_values_complex
        generic :: values => values_real,values_complex
        procedure :: poles => real_poles
        procedure :: roots => real_roots
        procedure :: residues => real_residues
        procedure :: clean_up => real_clean_up
    end type aaa_type

    type, public :: complex_aaa_type
        private
        type(aaa_core) :: core
    contains
        procedure :: valid => complex_valid
        procedure :: converged => complex_converged
        procedure :: fit_status => complex_fit_status
        procedure :: cleanup_status => complex_cleanup_status
        procedure :: max_error => complex_max_error
        procedure :: support_points => complex_support_points
        procedure :: support_values => complex_support_values
        procedure :: weights => complex_weights
        procedure :: errors => complex_errors
        procedure, private :: value_real => complex_value_real
        procedure, private :: value_complex => complex_value_complex
        generic :: value => value_real,value_complex
        procedure, private :: values_real => complex_values_real
        procedure, private :: values_complex => complex_values_complex
        generic :: values => values_real,values_complex
        procedure :: poles => complex_poles
        procedure :: roots => complex_roots
        procedure :: residues => complex_residues
        procedure :: clean_up => complex_clean_up
    end type complex_aaa_type

    interface aaa
        module procedure construct_rr,construct_rc,construct_cr,construct_cc
    end interface aaa
contains
    pure elemental logical function finite(z) result(ok)
        implicit none
        complex(dp), intent(in) :: z !! complex value to check
        ok=ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite

    pure complex(dp) function nan_value() result(v)
        implicit none
        v=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),ieee_value(0.0_dp,ieee_quiet_nan),dp)
    end function nan_value

    pure logical function less(a,b) result(ok)
        implicit none
        complex(dp), intent(in) :: a,b !! finite coordinates in lexicographic order
        ok=real(a,dp)<real(b,dp)
        if (real(a,dp)==real(b,dp)) ok=aimag(a)<aimag(b)
    end function less

    recursive subroutine sort_indices(x,index,work,lo,hi)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite sample coordinates
        integer, intent(inout) :: index(:),work(:) !! stable permutation and merge workspace
        integer, intent(in) :: lo,hi !! inclusive permutation range
        integer :: mid,i,j,k
        if (hi<=lo) return
        mid=lo+(hi-lo)/2
        call sort_indices(x,index,work,lo,mid)
        call sort_indices(x,index,work,mid+1,hi)
        i=lo
        j=mid+1
        do k=lo,hi
            if (i>mid) then
                work(k)=index(j)
                j=j+1
            else if (j>hi) then
                work(k)=index(i)
                i=i+1
            else if (less(x(index(j)),x(index(i)))) then
                work(k)=index(j)
                j=j+1
            else
                work(k)=index(i)
                i=i+1
            end if
        end do
        index(lo:hi)=work(lo:hi)
    end subroutine sort_indices

    subroutine build(x,y,core,status,rtol,max_terms,clean_up,clean_up_tol)
        implicit none
        complex(dp), intent(in) :: x(:),y(:) !! coordinates and scalar samples, equal length
        type(aaa_core), intent(out) :: core !! fitted normalized representation
        integer, intent(out) :: status !! aaa_* construction outcome
        real(dp), intent(in), optional :: rtol !! finite nonnegative tolerance, default eps**0.75
        integer, intent(in), optional :: max_terms !! positive support budget, default 100
        logical, intent(in), optional :: clean_up !! automatic doublet removal, default true
        real(dp), intent(in), optional :: clean_up_tol !! cleanup threshold >= 0, default 1e-13
        integer, allocatable :: index(:),work(:),selected(:)
        logical, allocatable :: mask(:),nonzero(:)
        real(dp), allocatable :: err(:),history(:)
        complex(dp), allocatable :: weights(:),zs(:),fs(:)
        complex(dp) :: mean,v
        real(dp) :: tolerance,ctol,normf,big
        integer :: n,i,j,k,m,budget,info,removed
        logical :: do_cleanup
        status=aaa_invalid_input
        n=size(x)
        if (n<1 .or. size(y)/=n) return
        if (.not.all(finite(x))) return
        tolerance=epsilon(1.0_dp)**0.75_dp
        if (present(rtol)) tolerance=rtol
        if (.not.ieee_is_finite(tolerance)) return
        if (tolerance<0.0_dp) return
        ctol=1e-13_dp
        if (present(clean_up_tol)) ctol=clean_up_tol
        if (.not.ieee_is_finite(ctol)) return
        if (ctol<0.0_dp) return
        budget=100
        if (present(max_terms)) budget=max_terms
        if (budget<1) return
        do_cleanup=.true.
        if (present(clean_up)) do_cleanup=clean_up
        allocate(index(n),work(n),stat=info)
        if (info/=0) then
            status=aaa_setup_failure
            return
        end if
        k=0
        do i=1,n
            if (.not.finite(y(i))) cycle
            k=k+1
            index(k)=i
        end do
        if (k==0) return
        call sort_indices(x,index,work,1,k)
        j=1
        do i=2,k
            if (x(index(i))==x(index(j))) cycle
            j=j+1
            index(j)=index(i)
        end do
        n=j
        budget=min(budget,n)
        allocate(core%points(n),core%samples(n),core%z(n),core%f(n),selected(budget), &
            mask(n),err(n),history(budget),zs(budget),fs(budget),stat=info)
        if (info/=0) then
            status=aaa_setup_failure
            return
        end if
        core%points=x(index(:n))
        core%samples=y(index(:n))
        core%center=cmplx(0.5_dp*minval(real(core%points,dp))+ &
            0.5_dp*maxval(real(core%points,dp)),0.5_dp*minval(aimag(core%points))+ &
            0.5_dp*maxval(aimag(core%points)),dp)
        core%z=core%points-core%center
        core%xscale=max(maxval(abs(real(core%z,dp))),maxval(abs(aimag(core%z))))
        if (core%xscale==0.0_dp) core%xscale=1.0_dp
        core%z=core%z/core%xscale
        ! Distinct inputs must remain distinguishable after affine normalization.
        if (any(core%z(2:)==core%z(:n-1))) then
            status=aaa_numerical_failure
            return
        end if
        core%yscale=max(maxval(abs(real(y(index(:n)),dp))),maxval(abs(aimag(y(index(:n))))))
        if (core%yscale==0.0_dp) core%yscale=1.0_dp
        core%f=y(index(:n))/core%yscale
        core%rtol=tolerance
        normf=maxval(abs(core%f))
        mean=sum(core%f)/real(n,dp)
        err=abs(core%f-mean)
        mask=.true.
        status=aaa_max_terms
        do m=1,budget
            j=maxloc(err,dim=1,mask=mask)
            selected(m)=j
            zs(m)=core%z(j)
            fs(m)=core%f(j)
            mask(j)=.false.
            call fit_weights(core%z,core%f,zs(:m),fs(:m),weights,info)
            if (info/=aaa_success) then
                status=info
                return
            end if
            do i=1,n
                v=bary_value(zs(:m),fs(:m),weights,core%z(i),info)
                if (info==0) then
                    err(i)=abs(core%f(i)-v)
                else
                    err(i)=ieee_value(0.0_dp,ieee_positive_inf)
                end if
            end do
            big=maxval(err)
            history(m)=absolute_error(big,core%yscale)
            ! Divide the error by the data norm to avoid overflowing rtol*normf.
            core%converged=big==0.0_dp
            if (normf>0.0_dp .and. ieee_is_finite(big)) &
                core%converged=big/normf<=tolerance
            if (core%converged) then
                status=aaa_success
                exit
            end if
        end do
        m=min(m,budget)
        nonzero=abs(weights)>0.0_dp
        core%s=pack(zs(:m),nonzero)
        core%v=pack(fs(:m),nonzero)
        core%w=pack(weights,nonzero)
        core%indices=pack(selected(:m),nonzero)
        core%errors=history(:m)
        core%error=history(m)
        core%fit_status=status
        if (do_cleanup .and. size(core%w)>1) then
            removed=cleanup(core,ctol,info)
            if (info/=0) then
                status=info
            else
                status=core%fit_status
            end if
        end if
    end subroutine build

    subroutine fit_weights(z,f,s,v,w,status)
        implicit none
        complex(dp), intent(in) :: z(:),f(:) !! normalized data coordinates and samples
        complex(dp), intent(in) :: s(:),v(:) !! selected support coordinates and samples
        complex(dp), allocatable, intent(out) :: w(:) !! normalized barycentric weights
        integer, intent(out) :: status !! zero on successful SVD fit
        complex(dp), allocatable :: a(:,:),u(:,:),vh(:,:)
        real(dp), allocatable :: sigma(:),col_norm(:)
        logical, allocatable :: mask(:)
        real(dp) :: scale,norm,tolerance
        integer :: n,m,rows,i,j,k,nu,ns,info,rank,first
        logical :: equilibrate
        n=size(z)
        m=size(s)
        status=aaa_setup_failure
        allocate(mask(n),w(m),stat=info)
        if (info/=0) return
        mask=.true.
        do j=1,m
            where(z==s(j)) mask=.false.
        end do
        rows=count(mask)
        w=cmplx(1.0_dp/sqrt(real(m,dp)),0.0_dp,dp)
        status=aaa_success
        if (rows==0) return
        ns=min(rows,m)
        nu=ns
        allocate(a(rows,m),u(rows,nu),vh(m,m),sigma(ns),col_norm(m),stat=info)
        if (info/=0) then
            status=aaa_setup_failure
            return
        end if
        status=aaa_numerical_failure
        k=0
        do i=1,n
            if (.not.mask(i)) cycle
            k=k+1
            do j=1,m
                a(k,j)=divide_value(f(i)-v(j),z(i)-s(j),info)
                if (info/=0) return
            end do
        end do
        status=aaa_numerical_failure
        if (.not.all(finite(a))) return
        scale=max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
        if (scale>0.0_dp) a=a/scale
        call svd(a,u,sigma,vh,info,full_matrices=rows<m)
        if (info/=0) return
        col_norm=1.0_dp
        equilibrate=.false.
        if (rows>=m) equilibrate=sigma(m)<=3.0_dp*epsilon(1.0_dp)*sigma(1)
        if (equilibrate) then
            do j=1,m
                scale=maxval(abs(a(:,j)))
                if (scale==0.0_dp) cycle
                col_norm(j)=scale*sqrt(sum(abs(a(:,j)/scale)**2))
                a(:,j)=a(:,j)/col_norm(j)
            end do
            call svd(a,u,sigma,vh,info)
            if (info/=0) return
        end if
        if (rows>=m) then
            first=minloc(sigma,dim=1)
            w=cmplx(0.0_dp,0.0_dp,dp)
            k=0
            do j=1,m
                if (sigma(j)/=sigma(first)) cycle
                w=w+conjg(vh(j,:))
                k=k+1
            end do
            ! A common scale cancels from the barycentric quotient. Applying it
            ! before undoing column scaling avoids reciprocal overflow.
            w=(w/sqrt(real(k,dp)))*(minval(col_norm)/col_norm)
        else
            tolerance=real(max(rows,m),dp)*epsilon(1.0_dp)*sigma(1)
            rank=count(sigma>tolerance)
            w=sum(conjg(vh(rank+1:m,:)),dim=1)/sqrt(real(m-rank,dp))
        end if
        scale=maxval(abs(w))
        if (.not.ieee_is_finite(scale) .or. scale==0.0_dp) return
        w=w/scale
        norm=sqrt(sum(abs(w)**2))
        w=w/norm
        status=aaa_success
    end subroutine fit_weights

    function divide_value(num,den,status) result(value)
        implicit none
        complex(dp), intent(in) :: num,den !! finite numerator and nonzero denominator
        integer, intent(out) :: status !! zero if the quotient is representable
        complex(dp) :: value,q
        real(dp) :: sn,sd,factor,parts(2),e,logs
        integer :: j
        value=nan_value()
        status=aaa_numerical_failure
        if (.not.finite(num) .or. .not.finite(den)) return
        sn=max(abs(real(num,dp)),abs(aimag(num)))
        sd=max(abs(real(den,dp)),abs(aimag(den)))
        if (sd==0.0_dp) return
        if (sn==0.0_dp) then
            value=cmplx(0.0_dp,0.0_dp,dp)
            status=aaa_success
            return
        end if
        q=(num/sn)/(den/sd)
        if (sd>=1.0_dp) then
            factor=sn/sd
        else if (sn<=huge(1.0_dp)*sd) then
            factor=sn/sd
        else
            factor=0.0_dp
        end if
        if (factor>0.0_dp) then
            value=scale_value(q,factor,status)
            return
        end if
        ! Extreme exponent ratios use componentwise logarithmic rescaling.
        parts=[real(q,dp),aimag(q)]
        logs=log(sn)-log(sd)
        do j=1,2
            if (parts(j)==0.0_dp) cycle
            e=log(abs(parts(j)))+logs
            if (e>log(huge(1.0_dp))) return
            if (e<log(tiny(1.0_dp))+log(epsilon(1.0_dp))-log(2.0_dp)) then
                parts(j)=0.0_dp
            else
                parts(j)=sign(exp(e),parts(j))
            end if
        end do
        value=cmplx(parts(1),parts(2),dp)
        status=aaa_success
    end function divide_value

    function bary_value(s,v,w,q,status) result(value)
        implicit none
        complex(dp), intent(in) :: s(:),v(:),w(:) !! support points, normalized values, weights
        complex(dp), intent(in) :: q !! scalar query coordinate, possibly infinite or NaN
        integer, intent(out) :: status !! zero on finite evaluation, otherwise aaa_* error
        complex(dp) :: value,num,den
        complex(dp) :: d(size(s)),terms(size(s))
        real(dp) :: scale,near
        integer :: i
        value=nan_value()
        status=aaa_invalid_input
        if (size(w)==0) return
        if (ieee_is_nan(real(q,dp)) .or. ieee_is_nan(aimag(q))) return
        status=aaa_numerical_failure
        if (.not.finite(q)) then
            den=sum(w)
            if (abs(den)==0.0_dp) return
            value=divide_value(sum(w*v),den,status)
        else
            do i=1,size(s)
                if (q==s(i) .and. abs(w(i))>0.0_dp) then
                    value=v(i)
                    status=aaa_success
                    return
                end if
            end do
            scale=max(maxval(abs(real(s,dp))),maxval(abs(aimag(s))), &
                abs(real(q,dp)),abs(aimag(q)))
            if (scale==0.0_dp) scale=1.0_dp
            d=q/scale-s/scale
            near=huge(1.0_dp)
            do i=1,size(s)
                if (abs(w(i))==0.0_dp) cycle
                if (abs(d(i))==0.0_dp) return
                near=min(near,abs(d(i)))
            end do
            terms=cmplx(0.0_dp,0.0_dp,dp)
            do i=1,size(s)
                if (abs(w(i))==0.0_dp) cycle
                terms(i)=w(i)*(near/d(i))
            end do
            num=sum(terms*v)
            den=sum(terms)
            if (abs(den)==0.0_dp) return
            value=divide_value(num,den,status)
        end if
        if (finite(value)) status=aaa_success
    end function bary_value

    pure real(dp) function absolute_error(error,scale) result(value)
        implicit none
        real(dp), intent(in) :: error,scale !! nonnegative normalized error and positive data scale
        value=ieee_value(0.0_dp,ieee_positive_inf)
        if (.not.ieee_is_finite(error)) return
        if (scale>1.0_dp) then
            if (error>huge(1.0_dp)/scale) return
        end if
        value=error*scale
    end function absolute_error

    function scale_value(value,scale,status) result(out)
        implicit none
        complex(dp), intent(in) :: value !! normalized value to convert to physical units
        real(dp), intent(in) :: scale !! finite positive scale
        integer, intent(out) :: status !! zero if the result is representable
        complex(dp) :: out
        out=nan_value()
        status=aaa_numerical_failure
        if (.not.finite(value)) return
        if (scale>1.0_dp) then
            if (max(abs(real(value,dp)),abs(aimag(value)))>huge(1.0_dp)/scale) return
        end if
        out=value*scale
        status=aaa_success
    end function scale_value

    function add_value(a,b,status) result(value)
        implicit none
        complex(dp), intent(in) :: a,b !! finite complex summands
        integer, intent(out) :: status !! zero if both result components are representable
        complex(dp) :: value
        real(dp) :: aa(2),bb(2)
        integer :: j
        value=nan_value()
        status=aaa_numerical_failure
        if (.not.finite(a) .or. .not.finite(b)) return
        aa=[real(a,dp),aimag(a)]
        bb=[real(b,dp),aimag(b)]
        do j=1,2
            if (bb(j)>0.0_dp) then
                if (aa(j)>huge(1.0_dp)-bb(j)) return
            else
                if (aa(j)<-huge(1.0_dp)-bb(j)) return
            end if
        end do
        value=a+b
        status=aaa_success
    end function add_value

    function scale_two(value,first,second,status) result(out)
        implicit none
        complex(dp), intent(in) :: value !! normalized value
        real(dp), intent(in) :: first,second !! positive physical scale factors
        integer, intent(out) :: status !! zero for a representable result
        complex(dp) :: out
        real(dp) :: factor,small,large,parts(2),e
        integer :: j
        small=min(first,second)
        large=max(first,second)
        out=nan_value()
        status=aaa_numerical_failure
        if (.not.finite(value)) return
        factor=0.0_dp
        if (large<=1.0_dp) then
            factor=small*large
        else if (small<=huge(1.0_dp)/large) then
            factor=small*large
        end if
        if (factor>=tiny(1.0_dp)) then
            out=scale_value(value,factor,status)
            return
        end if
        parts=[real(value,dp),aimag(value)]
        do j=1,2
            if (parts(j)==0.0_dp) cycle
            e=log(abs(parts(j)))+log(first)+log(second)
            if (e>log(huge(1.0_dp))) return
            if (e<log(tiny(1.0_dp))+log(epsilon(1.0_dp))-log(2.0_dp)) then
                parts(j)=0.0_dp
            else
                parts(j)=sign(exp(e),parts(j))
            end if
        end do
        out=cmplx(parts(1),parts(2),dp)
        status=aaa_success
    end function scale_two

    function evaluate(core,q,status) result(value)
        implicit none
        type(aaa_core), intent(in) :: core !! retained fit
        complex(dp), intent(in) :: q !! physical query coordinate
        integer, intent(out) :: status !! scalar evaluation outcome
        complex(dp) :: value,temp,query
        integer :: j,info
        value=nan_value()
        status=aaa_invalid_input
        if (.not.allocated(core%w)) return
        if (finite(q)) then
            do j=1,size(core%w)
                if (q/=core%points(core%indices(j))) cycle
                value=core%samples(core%indices(j))
                status=aaa_success
                return
            end do
        end if
        info=aaa_numerical_failure
        if (finite(q)) then
            temp=add_value(q,-core%center,info)
            if (info==0) query=divide_value(temp,cmplx(core%xscale,0.0_dp,dp),info)
        end if
        if (info==0) then
            temp=bary_value(core%s,core%v,core%w,query,status)
        else
            temp=bary_value(core%points(core%indices),core%v,core%w,q,status)
        end if
        if (status/=0) return
        value=scale_value(temp,core%yscale,status)
    end function evaluate

    function spectrum(core,zeros,status,normalized) result(values)
        implicit none
        type(aaa_core), intent(in) :: core !! retained fit
        logical, intent(in) :: zeros !! numerator zeros if true; denominator poles otherwise
        integer, intent(out) :: status !! zero, invalid object, size limit, or numerical failure
        logical, intent(in), optional :: normalized !! keep normalized coordinates, default false
        complex(dp), allocatable :: values(:),a(:,:),b(:,:),eigen(:),u(:),candidate(:)
        complex(dp), allocatable :: alpha(:),beta(:)
        complex(dp) :: p,delta,dv,der,temp
        complex(dp) :: d(size_if_valid(core)),terms(size_if_valid(core))
        real(dp) :: near,bound
        integer :: n,m,i,j,k,info,it
        logical :: norm_output
        norm_output=.false.
        if (present(normalized)) norm_output=normalized
        allocate(values(0))
        status=aaa_invalid_input
        if (.not.allocated(core%w)) return
        m=size(core%w)
        if (zeros .and. all(abs(core%v)==0.0_dp)) then
            status=aaa_undefined
            return
        end if
        status=aaa_success
        if (m==1) return
        if (m+1>matrix_function_max_order) then
            status=aaa_spectral_size_limit
            return
        end if
        n=m+1
        allocate(a(n,n),b(n,n),eigen(n),u(m),candidate(n),alpha(n),beta(n),stat=info)
        if (info/=0) then
            status=aaa_setup_failure
            return
        end if
        u=core%w
        if (zeros) u=u*core%v
        if (maxval(abs(u))>0.0_dp) u=u/maxval(abs(u))
        a=cmplx(0.0_dp,0.0_dp,dp)
        b=a
        a(1,2:)=u
        a(2:,1)=cmplx(1.0_dp,0.0_dp,dp)
        do i=1,m
            a(i+1,i+1)=core%s(i)
            b(i+1,i+1)=cmplx(1.0_dp,0.0_dp,dp)
        end do
        call eigvals(a,b,eigen,info,alpha=alpha,beta=beta)
        status=aaa_numerical_failure
        if (info/=0) return
        ! Resolve near-infinite projective roots at binary64 precision.
        where(abs(beta)<=64.0_dp*epsilon(1.0_dp)*abs(alpha)) eigen=nan_value()
        k=count(finite(eigen))
        if (k>m-1) return
        k=0
        do i=1,n
            if (.not.finite(eigen(i))) cycle
            p=eigen(i)
            ! Exact numerator factors at zero-valued support points need no division.
            j=minloc(abs(p-core%s),dim=1)
            if (abs(u(j))==0.0_dp .and. abs(p-core%s(j))<1e-8_dp) then
                p=core%s(j)
            else
                do it=1,8
                    d=p-core%s
                    if (any(abs(d)==0.0_dp)) return
                    near=minval(abs(d))
                    terms=u*(near/d)
                    dv=sum(terms)
                    der=-sum(terms/d)
                    if (abs(der)==0.0_dp) exit
                    delta=divide_value(dv,der,info)
                    if (info/=0) return
                    if (.not.finite(delta)) return
                    if (abs(delta)>0.5_dp*max(1.0_dp,abs(p))) exit
                    p=p-delta
                    if (abs(delta)<=16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(p))) exit
                end do
                d=p-core%s
                if (any(abs(d)==0.0_dp)) return
                terms=u*(minval(abs(d))/d)
                bound=sum(abs(terms))
                if (abs(sum(terms))>1e-7_dp*bound) return
            end if
            if (norm_output) then
                temp=p
            else
                temp=scale_value(p,core%xscale,info)
                if (info/=0) return
                temp=add_value(temp,core%center,info)
                if (info/=0) return
            end if
            k=k+1
            candidate(k)=temp
        end do
        values=candidate(:k)
        status=aaa_success
    end function spectrum

    pure integer function size_if_valid(core) result(n)
        implicit none
        type(aaa_core), intent(in) :: core !! fitted object, possibly uninitialized
        n=0
        if (allocated(core%w)) n=size(core%w)
    end function size_if_valid

    function residue_values(core,status,normalized) result(values)
        implicit none
        type(aaa_core), intent(in) :: core !! retained fit with simple finite poles
        integer, intent(out) :: status !! zero or spectral/numerical failure
        logical, intent(in), optional :: normalized !! keep normalized residues, default false
        complex(dp), allocatable :: values(:),poles(:),work(:)
        complex(dp) :: p,num,den,temp
        complex(dp) :: d(size_if_valid(core)),terms(size_if_valid(core))
        integer :: i,info
        logical :: norm_output
        norm_output=.false.
        if (present(normalized)) norm_output=normalized
        allocate(values(0))
        poles=spectrum(core,.false.,status,.true.)
        if (status/=0) return
        allocate(work(size(poles)))
        status=aaa_numerical_failure
        do i=1,size(poles)
            p=poles(i)
            d=p-core%s
            if (any(abs(d)==0.0_dp)) return
            terms=core%w*(minval(abs(d))/d)
            num=sum(terms*core%v)
            den=-sum(terms/d)
            if (abs(den)==0.0_dp) return
            temp=divide_value(num,den,info)
            if (info/=0) return
            if (norm_output) then
                work(i)=temp
            else
                work(i)=scale_two(temp,core%xscale,core%yscale,info)
                if (info/=0) return
            end if
        end do
        values=work
        status=aaa_success
    end function residue_values

    integer function cleanup(core,tolerance,status) result(removed)
        implicit none
        type(aaa_core), intent(inout) :: core !! fit changed only after a successful refit
        real(dp), intent(in) :: tolerance !! finite nonnegative relative residue threshold
        integer, intent(out) :: status !! zero or cleanup failure, retained on the object
        complex(dp), allocatable :: poles(:),residues(:),w(:),s(:),v(:)
        logical, allocatable :: keep(:)
        integer, allocatable :: indices(:)
        complex(dp) :: val
        real(dp) :: geom,dist,ratio,big,normf
        integer :: i,j,info,ni
        removed=0
        status=aaa_invalid_input
        core%cleanup_status=status
        if (.not.allocated(core%w)) return
        if (.not.ieee_is_finite(tolerance)) return
        if (tolerance<0.0_dp) return
        status=aaa_success
        core%cleanup_status=status
        if (size(core%w)==1 .or. tolerance==0.0_dp) return
        poles=spectrum(core,.false.,status,.true.)
        if (status/=0) then
            core%cleanup_status=status
            return
        end if
        residues=residue_values(core,status,.true.)
        if (status/=0) then
            core%cleanup_status=status
            return
        end if
        if (any(abs(core%f)==0.0_dp)) return
        geom=exp(sum(log(abs(core%f)))/real(size(core%f),dp))
        allocate(keep(size(core%w)))
        keep=.true.
        do i=1,size(poles)
            ! Use normalized coordinates and residues to avoid scale-dependent decisions.
            val=poles(i)
            dist=minval(abs(val-core%z))
            if (dist==0.0_dp) cycle
            ratio=abs(residues(i))/dist
            if (ratio/geom>=tolerance) cycle
            j=minloc(abs(val-core%s),dim=1)
            keep(j)=.false.
        end do
        ni=count(.not.keep)
        if (ni==0) return
        status=aaa_numerical_failure
        if (.not.any(keep)) then
            core%cleanup_status=status
            return
        end if
        s=pack(core%s,keep)
        v=pack(core%v,keep)
        indices=pack(core%indices,keep)
        call fit_weights(core%z,core%f,s,v,w,status)
        if (status/=0) then
            core%cleanup_status=status
            return
        end if
        big=0.0_dp
        do i=1,size(core%z)
            val=bary_value(s,v,w,core%z(i),info)
            if (info/=0) then
                status=aaa_numerical_failure
                core%cleanup_status=status
                return
            end if
            big=max(big,abs(core%f(i)-val))
        end do
        keep=abs(w)>0.0_dp
        core%s=pack(s,keep)
        core%v=pack(v,keep)
        core%w=pack(w,keep)
        core%indices=pack(indices,keep)
        core%error=absolute_error(big,core%yscale)
        normf=maxval(abs(core%f))
        core%converged=big==0.0_dp
        if (normf>0.0_dp) core%converged=big/normf<=core%rtol
        core%fit_status=aaa_max_terms
        if (core%converged) core%fit_status=aaa_success
        core%cleanup_status=aaa_success
        status=aaa_success
        removed=ni
    end function cleanup

    function construct_rr(x,y,rtol,max_terms,clean_up,clean_up_tol,status) result(res)
        implicit none
        real(dp), intent(in) :: x(:) !! finite sample coordinates, nonempty
        real(dp), intent(in) :: y(:) !! scalar function samples, same length as x
        real(dp), intent(in), optional :: rtol !! nonnegative relative tolerance, default eps**0.75
        integer, intent(in), optional :: max_terms !! positive support budget, default 100
        logical, intent(in), optional :: clean_up !! remove spurious pole-zero pairs, default true
        real(dp), intent(in), optional :: clean_up_tol !! cleanup threshold, default 1e-13
        integer, intent(out), optional :: status !! aaa_* fit or automatic-cleanup outcome
        type(aaa_type) :: res
        integer :: info
        call build(cmplx(x,0.0_dp,dp),cmplx(y,0.0_dp,dp),res%core,info, &
            rtol,max_terms,clean_up,clean_up_tol)
        if (.not.allocated(res%core%w)) res%core%fit_status=info
        if (present(status)) status=info
    end function construct_rr

    function construct_rc(x,y,rtol,max_terms,clean_up,clean_up_tol,status) result(res)
        implicit none
        real(dp), intent(in) :: x(:) !! finite sample coordinates, nonempty
        complex(dp), intent(in) :: y(:) !! scalar function samples, same length as x
        real(dp), intent(in), optional :: rtol !! nonnegative relative tolerance, default eps**0.75
        integer, intent(in), optional :: max_terms !! positive support budget, default 100
        logical, intent(in), optional :: clean_up !! remove spurious pole-zero pairs, default true
        real(dp), intent(in), optional :: clean_up_tol !! cleanup threshold, default 1e-13
        integer, intent(out), optional :: status !! aaa_* fit or automatic-cleanup outcome
        type(complex_aaa_type) :: res
        integer :: info
        call build(cmplx(x,0.0_dp,dp),y,res%core,info,rtol,max_terms,clean_up,clean_up_tol)
        if (.not.allocated(res%core%w)) res%core%fit_status=info
        if (present(status)) status=info
    end function construct_rc

    function construct_cr(x,y,rtol,max_terms,clean_up,clean_up_tol,status) result(res)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite sample coordinates, nonempty
        real(dp), intent(in) :: y(:) !! scalar function samples, same length as x
        real(dp), intent(in), optional :: rtol !! nonnegative relative tolerance, default eps**0.75
        integer, intent(in), optional :: max_terms !! positive support budget, default 100
        logical, intent(in), optional :: clean_up !! remove spurious pole-zero pairs, default true
        real(dp), intent(in), optional :: clean_up_tol !! cleanup threshold, default 1e-13
        integer, intent(out), optional :: status !! aaa_* fit or automatic-cleanup outcome
        type(complex_aaa_type) :: res
        integer :: info
        call build(x,cmplx(y,0.0_dp,dp),res%core,info,rtol,max_terms,clean_up,clean_up_tol)
        if (.not.allocated(res%core%w)) res%core%fit_status=info
        if (present(status)) status=info
    end function construct_cr

    function construct_cc(x,y,rtol,max_terms,clean_up,clean_up_tol,status) result(res)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite sample coordinates, nonempty
        complex(dp), intent(in) :: y(:) !! scalar function samples, same length as x
        real(dp), intent(in), optional :: rtol !! nonnegative relative tolerance, default eps**0.75
        integer, intent(in), optional :: max_terms !! positive support budget, default 100
        logical, intent(in), optional :: clean_up !! remove spurious pole-zero pairs, default true
        real(dp), intent(in), optional :: clean_up_tol !! cleanup threshold, default 1e-13
        integer, intent(out), optional :: status !! aaa_* fit or automatic-cleanup outcome
        type(complex_aaa_type) :: res
        integer :: info
        call build(x,y,res%core,info,rtol,max_terms,clean_up,clean_up_tol)
        if (.not.allocated(res%core%w)) res%core%fit_status=info
        if (present(status)) status=info
    end function construct_cc

    logical function real_valid(self) result(ok)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        ok=.false.
        if (allocated(self%core%w)) ok=size(self%core%w)>0
    end function real_valid

    logical function real_converged(self) result(ok)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        ok=self%core%converged .and. self%valid()
    end function real_converged

    integer function real_fit_status(self) result(code)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        code=self%core%fit_status
    end function real_fit_status

    integer function real_cleanup_status(self) result(code)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        code=self%core%cleanup_status
    end function real_cleanup_status

    real(dp) function real_max_error(self) result(error)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        error=ieee_value(0.0_dp,ieee_quiet_nan)
        if (self%valid()) error=self%core%error
    end function real_max_error

    function real_support_points(self) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        real(dp), allocatable :: values(:)
        if (self%valid()) then
            values=real(self%core%points(self%core%indices),dp)
        else
            allocate(values(0))
        end if
    end function real_support_points

    function real_support_values(self) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        real(dp), allocatable :: values(:)
        if (self%valid()) then
            values=real(self%core%samples(self%core%indices),dp)
        else
            allocate(values(0))
        end if
    end function real_support_values

    function real_weights(self) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        real(dp), allocatable :: values(:)
        if (self%valid()) then
            values=real(self%core%w,dp)
        else
            allocate(values(0))
        end if
    end function real_weights

    function real_errors(self) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        real(dp), allocatable :: values(:)
        allocate(values(0))
        if (allocated(self%core%errors)) values=self%core%errors
    end function real_errors

    function real_value_real(self,q,status) result(value)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        real(dp), intent(in) :: q !! scalar query coordinate
        integer, intent(out), optional :: status !! aaa_* evaluation outcome
        real(dp) :: value
        complex(dp) :: temp
        integer :: info
        temp=evaluate(self%core,cmplx(q,0.0_dp,dp),info)
        value=real(temp,dp)
        if (present(status)) status=info
    end function real_value_real

    function real_values_real(self,q,status) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        real(dp), intent(in) :: q(:) !! query coordinates, in any order
        integer, intent(out), optional :: status !! zero on success; first query error otherwise
        real(dp), allocatable :: values(:)
        integer :: i,info,code
        allocate(values(size(q)))
        code=aaa_invalid_input
        if (self%valid()) code=aaa_success
        do i=1,size(q)
            values(i)=self%value(q(i),info)
            if (code==0 .and. info/=0) code=info
        end do
        if (present(status)) status=code
    end function real_values_real

    function real_value_complex(self,q,status) result(value)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        complex(dp), intent(in) :: q !! scalar query coordinate
        integer, intent(out), optional :: status !! aaa_* evaluation outcome
        complex(dp) :: value
        complex(dp) :: temp
        integer :: info
        temp=evaluate(self%core,q,info)
        value=temp
        if (present(status)) status=info
    end function real_value_complex

    function real_values_complex(self,q,status) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        complex(dp), intent(in) :: q(:) !! query coordinates, in any order
        integer, intent(out), optional :: status !! zero on success; first query error otherwise
        complex(dp), allocatable :: values(:)
        integer :: i,info,code
        allocate(values(size(q)))
        code=aaa_invalid_input
        if (self%valid()) code=aaa_success
        do i=1,size(q)
            values(i)=self%value(q(i),info)
            if (code==0 .and. info/=0) code=info
        end do
        if (present(status)) status=code
    end function real_values_complex

    function real_poles(self,status) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        integer, intent(out), optional :: status !! aaa_* spectral calculation outcome
        complex(dp), allocatable :: values(:)
        integer :: info
        values=spectrum(self%core,.false.,info)
        if (present(status)) status=info
    end function real_poles

    function real_roots(self,status) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        integer, intent(out), optional :: status !! aaa_* spectral calculation outcome
        complex(dp), allocatable :: values(:)
        integer :: info
        values=spectrum(self%core,.true.,info)
        if (present(status)) status=info
    end function real_roots

    function real_residues(self,status) result(values)
        implicit none
        class(aaa_type), intent(in) :: self !! approximation object
        integer, intent(out), optional :: status !! aaa_* spectral calculation outcome
        complex(dp), allocatable :: values(:)
        integer :: info
        values=residue_values(self%core,info)
        if (present(status)) status=info
    end function real_residues

    integer function real_clean_up(self,cleanup_tol,status) result(removed)
        implicit none
        class(aaa_type), intent(inout) :: self !! fit updated only on successful cleanup
        real(dp), intent(in), optional :: cleanup_tol !! nonnegative threshold, default 1e-13
        integer, intent(out), optional :: status !! aaa_* cleanup outcome
        real(dp) :: tolerance
        integer :: info
        tolerance=1e-13_dp
        if (present(cleanup_tol)) tolerance=cleanup_tol
        removed=cleanup(self%core,tolerance,info)
        if (present(status)) status=info
    end function real_clean_up

    logical function complex_valid(self) result(ok)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        ok=.false.
        if (allocated(self%core%w)) ok=size(self%core%w)>0
    end function complex_valid

    logical function complex_converged(self) result(ok)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        ok=self%core%converged .and. self%valid()
    end function complex_converged

    integer function complex_fit_status(self) result(code)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        code=self%core%fit_status
    end function complex_fit_status

    integer function complex_cleanup_status(self) result(code)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        code=self%core%cleanup_status
    end function complex_cleanup_status

    real(dp) function complex_max_error(self) result(error)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        error=ieee_value(0.0_dp,ieee_quiet_nan)
        if (self%valid()) error=self%core%error
    end function complex_max_error

    function complex_support_points(self) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        complex(dp), allocatable :: values(:)
        if (self%valid()) then
            values=self%core%points(self%core%indices)
        else
            allocate(values(0))
        end if
    end function complex_support_points

    function complex_support_values(self) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        complex(dp), allocatable :: values(:)
        if (self%valid()) then
            values=self%core%samples(self%core%indices)
        else
            allocate(values(0))
        end if
    end function complex_support_values

    function complex_weights(self) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        complex(dp), allocatable :: values(:)
        if (self%valid()) then
            values=self%core%w
        else
            allocate(values(0))
        end if
    end function complex_weights

    function complex_errors(self) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        real(dp), allocatable :: values(:)
        allocate(values(0))
        if (allocated(self%core%errors)) values=self%core%errors
    end function complex_errors

    function complex_value_real(self,q,status) result(value)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        real(dp), intent(in) :: q !! scalar query coordinate
        integer, intent(out), optional :: status !! aaa_* evaluation outcome
        complex(dp) :: value
        complex(dp) :: temp
        integer :: info
        temp=evaluate(self%core,cmplx(q,0.0_dp,dp),info)
        value=temp
        if (present(status)) status=info
    end function complex_value_real

    function complex_values_real(self,q,status) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        real(dp), intent(in) :: q(:) !! query coordinates, in any order
        integer, intent(out), optional :: status !! zero on success; first query error otherwise
        complex(dp), allocatable :: values(:)
        integer :: i,info,code
        allocate(values(size(q)))
        code=aaa_invalid_input
        if (self%valid()) code=aaa_success
        do i=1,size(q)
            values(i)=self%value(q(i),info)
            if (code==0 .and. info/=0) code=info
        end do
        if (present(status)) status=code
    end function complex_values_real

    function complex_value_complex(self,q,status) result(value)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        complex(dp), intent(in) :: q !! scalar query coordinate
        integer, intent(out), optional :: status !! aaa_* evaluation outcome
        complex(dp) :: value
        complex(dp) :: temp
        integer :: info
        temp=evaluate(self%core,q,info)
        value=temp
        if (present(status)) status=info
    end function complex_value_complex

    function complex_values_complex(self,q,status) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        complex(dp), intent(in) :: q(:) !! query coordinates, in any order
        integer, intent(out), optional :: status !! zero on success; first query error otherwise
        complex(dp), allocatable :: values(:)
        integer :: i,info,code
        allocate(values(size(q)))
        code=aaa_invalid_input
        if (self%valid()) code=aaa_success
        do i=1,size(q)
            values(i)=self%value(q(i),info)
            if (code==0 .and. info/=0) code=info
        end do
        if (present(status)) status=code
    end function complex_values_complex

    function complex_poles(self,status) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        integer, intent(out), optional :: status !! aaa_* spectral calculation outcome
        complex(dp), allocatable :: values(:)
        integer :: info
        values=spectrum(self%core,.false.,info)
        if (present(status)) status=info
    end function complex_poles

    function complex_roots(self,status) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        integer, intent(out), optional :: status !! aaa_* spectral calculation outcome
        complex(dp), allocatable :: values(:)
        integer :: info
        values=spectrum(self%core,.true.,info)
        if (present(status)) status=info
    end function complex_roots

    function complex_residues(self,status) result(values)
        implicit none
        class(complex_aaa_type), intent(in) :: self !! approximation object
        integer, intent(out), optional :: status !! aaa_* spectral calculation outcome
        complex(dp), allocatable :: values(:)
        integer :: info
        values=residue_values(self%core,info)
        if (present(status)) status=info
    end function complex_residues

    integer function complex_clean_up(self,cleanup_tol,status) result(removed)
        implicit none
        class(complex_aaa_type), intent(inout) :: self !! fit updated only on successful cleanup
        real(dp), intent(in), optional :: cleanup_tol !! nonnegative threshold, default 1e-13
        integer, intent(out), optional :: status !! aaa_* cleanup outcome
        real(dp) :: tolerance
        integer :: info
        tolerance=1e-13_dp
        if (present(cleanup_tol)) tolerance=cleanup_tol
        removed=cleanup(self%core,tolerance,info)
        if (present(status)) status=info
    end function complex_clean_up
end module scifort_interpolate_aaa
