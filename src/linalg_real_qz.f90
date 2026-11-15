! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Real generalized Schur adapters for the attributed private LAPACK kernels.
module scifort_linalg_real_qz
    use scifort_kinds, only: dp
    use scifort_linalg, only: linalg_status_success,linalg_status_invalid_shape, &
        linalg_status_no_convergence,linalg_status_singular
    use scifort_linalg_matrix_functions, only: matrix_function_max_order
    use scifort_linalg_qr, only: qr
    use scifort_linalg_qz_lapack, only: dgghrd,dhgeqz,dlag2
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: real_qz,real_qz_eigenvalues,real_qz_check
contains
    subroutine real_qz(a,b,aa,bb,q,z,status,alpha,beta,method)
        implicit none
        real(dp), intent(in) :: a(:,:) !! Finite first square matrix, order 0..40.
        real(dp), intent(in) :: b(:,:) !! Finite second square matrix, possibly singular.
        real(dp), intent(out) :: aa(:,:) !! Real quasi-upper-triangular first Schur factor.
        real(dp), intent(out) :: bb(:,:) !! Real upper-triangular second factor, nonnegative diagonal.
        real(dp), intent(out) :: q(:,:) !! Left real orthogonal Schur vectors.
        real(dp), intent(out) :: z(:,:) !! Right real orthogonal Schur vectors.
        integer, intent(out) :: status !! Existing linalg status code; zero means verified success.
        complex(dp), intent(out), optional :: alpha(:) !! Homogeneous numerators, size n.
        real(dp), intent(out), optional :: beta(:) !! Nonnegative homogeneous denominators, size n.
        character(*), intent(in), optional :: method !! auto (default) or iterative; both real QZ.
        real(dp), allocatable :: h(:,:),t(:,:),left(:,:),right(:,:),r(:,:),ar(:),ai(:),be(:),work(:)
        real(dp) :: sa,sb
        integer :: n,info,i
        aa=0.0_dp
        bb=0.0_dp
        q=0.0_dp
        z=0.0_dp
        if (present(alpha)) alpha=cmplx(0.0_dp,0.0_dp,dp)
        if (present(beta)) beta=0.0_dp
        status=linalg_status_invalid_shape
        n=size(a,1)
        if (n>matrix_function_max_order .or. size(a,2)/=n) return
        if (any(shape(b)/=[n,n]) .or. any(shape(aa)/=[n,n])) return
        if (any(shape(bb)/=[n,n]) .or. any(shape(q)/=[n,n])) return
        if (any(shape(z)/=[n,n])) return
        if (present(alpha)) then
            if (size(alpha)/=n) return
        end if
        if (present(beta)) then
            if (size(beta)/=n) return
        end if
        if (present(method)) then
            if (trim(adjustl(method))/='auto' .and. trim(adjustl(method))/='iterative') return
        end if
        status=linalg_status_no_convergence
        if (.not.all(ieee_is_finite(a)) .or. .not.all(ieee_is_finite(b))) return
        if (n==0) then
            status=linalg_status_success
            return
        end if
        sa=maxval(abs(a))
        sb=maxval(abs(b))
        if (sa==0.0_dp) sa=1.0_dp
        if (sb==0.0_dp) sb=1.0_dp
        allocate(h(n,n),t(n,n),left(n,n),right(n,n),r(n,n),ar(n),ai(n),be(n),work(n))
        h=a/sa
        t=b/sb
        if (any(a/=0.0_dp .and. h==0.0_dp)) return
        if (any(b/=0.0_dp .and. t==0.0_dp)) return
        call qr(t,left,r,info)
        if (info/=linalg_status_success) return
        t=r
        h=matmul(transpose(left),h)
        call dgghrd('v','i',n,1,n,h,n,t,n,left,n,right,n,info)
        if (info/=0) return
        call dhgeqz('s','v','v',n,1,n,h,n,t,n,ar,ai,be,left,n,right,n,work,n,info)
        if (info/=0) return
        if (any(ar==0.0_dp .and. ai==0.0_dp .and. be==0.0_dp)) then
            status=linalg_status_singular
            return
        end if
        call real_qz_check(a/sa,b/sb,h,t,left,right,info)
        if (info/=linalg_status_success) return
        if (sa>1.0_dp) then
            if (maxval(abs(h))>huge(sa)/sa) return
            if (maxval(hypot(ar,ai))>huge(sa)/sa) return
        end if
        if (sb>1.0_dp) then
            if (maxval(abs(t))>huge(sb)/sb .or. maxval(abs(be))>huge(sb)/sb) return
        end if
        if (any(h/=0.0_dp .and. h*sa==0.0_dp)) return
        if (any(t/=0.0_dp .and. t*sb==0.0_dp)) return
        if (any(ar/=0.0_dp .and. ar*sa==0.0_dp)) return
        if (any(ai/=0.0_dp .and. ai*sa==0.0_dp)) return
        if (any(be/=0.0_dp .and. be*sb==0.0_dp)) return
        h=h*sa
        t=t*sb
        ar=ar*sa
        ai=ai*sa
        be=be*sb
        if (.not.all(ieee_is_finite(h)) .or. .not.all(ieee_is_finite(t))) return
        if (.not.all(ieee_is_finite(ar)) .or. .not.all(ieee_is_finite(ai))) return
        if (.not.all(ieee_is_finite(be))) return
        aa=h
        bb=t
        q=left
        z=right
        if (present(alpha)) alpha=cmplx(ar,ai,dp)
        if (present(beta)) beta=be
        status=linalg_status_success
    end subroutine real_qz

    subroutine real_qz_eigenvalues(aa,bb,alpha,beta,status)
        implicit none
        real(dp), intent(in) :: aa(:,:) !! Finite real quasi-triangular first Schur factor.
        real(dp), intent(in) :: bb(:,:) !! Finite real triangular second factor.
        complex(dp), intent(out) :: alpha(:) !! Homogeneous numerators in block order.
        real(dp), intent(out) :: beta(:) !! Real homogeneous denominators in block order.
        integer, intent(out) :: status !! Zero success; four range/structure failure; five indeterminate.
        real(dp) :: sa,sb,a2(2,2),b2(2,2),s1,s2,wr1,wr2,wi
        integer :: i,n
        status=linalg_status_no_convergence
        n=size(aa,1)
        if (size(aa,2)/=n .or. any(shape(bb)/=[n,n])) return
        if (size(alpha)/=n .or. size(beta)/=n) return
        if (.not.all(ieee_is_finite(aa)) .or. .not.all(ieee_is_finite(bb))) return
        i=1
        do while (i<=n)
            if (i<n) then
                if (aa(i+1,i)/=0.0_dp) then
                    sa=maxval(abs(aa(i:i+1,i:i+1)))
                    sb=maxval(abs(bb(i:i+1,i:i+1)))
                    if (sa==0.0_dp .or. sb==0.0_dp) return
                    a2=aa(i:i+1,i:i+1)/sa
                    b2=bb(i:i+1,i:i+1)/sb
                    call dlag2(a2,2,b2,2,tiny(1.0_dp),s1,s2,wr1,wr2,wi)
                    if (wi<=0.0_dp) return
                    if (max(hypot(wr1,wi),hypot(wr2,wi))>huge(sa)/max(1.0_dp,sa)) return
                    if (max(abs(s1),abs(s2))>huge(sb)/max(1.0_dp,sb)) return
                    alpha(i)=cmplx(wr1*sa,wi*sa,dp)
                    alpha(i+1)=cmplx(wr2*sa,-wi*sa,dp)
                    beta(i)=s1*sb
                    beta(i+1)=s2*sb
                    i=i+2
                    cycle
                end if
            end if
            alpha(i)=cmplx(aa(i,i),0.0_dp,dp)
            beta(i)=bb(i,i)
            if (alpha(i)==cmplx(0.0_dp,0.0_dp,dp) .and. beta(i)==0.0_dp) then
                status=linalg_status_singular
                return
            end if
            i=i+1
        end do
        status=linalg_status_success
    end subroutine real_qz_eigenvalues

    subroutine real_qz_check(a,b,aa,bb,q,z,status)
        implicit none
        real(dp), intent(in) :: a(:,:) !! Normalized first original matrix.
        real(dp), intent(in) :: b(:,:) !! Normalized second original matrix.
        real(dp), intent(in) :: aa(:,:) !! First candidate Schur factor in matching units.
        real(dp), intent(in) :: bb(:,:) !! Second candidate Schur factor in matching units.
        real(dp), intent(in) :: q(:,:) !! Candidate left orthogonal basis.
        real(dp), intent(in) :: z(:,:) !! Candidate right orthogonal basis.
        integer, intent(out) :: status !! Zero verified; four nonfinite or failed reconstruction.
        real(dp), allocatable :: identity(:,:),error(:,:)
        real(dp) :: tol
        integer :: n,i,j
        status=linalg_status_no_convergence
        n=size(a,1)
        if (n==0) then
            status=linalg_status_success
            return
        end if
        if (.not.all(ieee_is_finite(aa)) .or. .not.all(ieee_is_finite(bb))) return
        if (.not.all(ieee_is_finite(q)) .or. .not.all(ieee_is_finite(z))) return
        tol=4096*epsilon(1.0_dp)*real(max(1,n),dp)
        if (maxval(abs(q))>1.0_dp+tol .or. maxval(abs(z))>1.0_dp+tol) return
        identity=0.0_dp*a
        do i=1,n
            identity(i,i)=1.0_dp
        end do
        if (maxval(abs(matmul(transpose(q),q)-identity))>tol) return
        if (maxval(abs(matmul(transpose(z),z)-identity))>tol) return
        error=matmul(q,matmul(aa,transpose(z)))-a
        if (.not.all(ieee_is_finite(error)) .or. maxval(abs(error))>tol) return
        error=matmul(q,matmul(bb,transpose(z)))-b
        if (.not.all(ieee_is_finite(error)) .or. maxval(abs(error))>tol) return
        do j=1,n
            if (bb(j,j)<0.0_dp) return
            do i=j+1,n
                if (bb(i,j)/=0.0_dp) return
                if (i>j+1 .and. aa(i,j)/=0.0_dp) return
            end do
            if (j<n-1) then
                if (aa(j+1,j)/=0.0_dp .and. aa(j+2,j+1)/=0.0_dp) return
            end if
        end do
        status=linalg_status_success
    end subroutine real_qz_check
end module scifort_linalg_real_qz
