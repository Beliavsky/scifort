! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-k QR updates and row/column insertion and deletion.
! Reuse existing QR factors. The low-rank update is factored in an augmented
! orthogonal basis; column deletion refactors only the small R core.
! Row edits and column insertion re-factor the explicitly edited matrix.
module scifort_linalg_qr_updates
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: qr_update, qr_insert, qr_delete

    interface qr_update
        module procedure qr_update_real_vector
        module procedure qr_update_real_matrix
        module procedure qr_update_complex_vector
        module procedure qr_update_complex_matrix
    end interface
    interface qr_insert
        module procedure qr_insert_real_vector
        module procedure qr_insert_real_matrix
        module procedure qr_insert_complex_vector
        module procedure qr_insert_complex_matrix
    end interface
    interface qr_delete
        module procedure qr_delete_real
        module procedure qr_delete_complex
    end interface
contains
    pure logical function finite_real(a) result(ok)
        real(dp), intent(in) :: a(:, :) !! input values
        ok = all(ieee_is_finite(a))
    end function finite_real

    pure logical function valid_factors_real(q, r) result(ok)
        real(dp), intent(in) :: q(:, :) !! existing orthonormal QR columns
        real(dp), intent(in) :: r(:, :) !! existing upper trapezoidal factor
        integer :: m, h, n
        m = size(q, 1)
        h = size(q, 2)
        n = size(r, 2)
        ok = m >= 1 .and. n >= 1 .and. &
            (h == m .or. h == min(m,n)) .and. size(r,1) == h
    end function valid_factors_real

    subroutine qr_update_real_vector(q, r, u, v, qnew, rnew, status)
        real(dp), intent(in) :: q(:, :) !! m-by-h full or economic Q
        real(dp), intent(in) :: r(:, :) !! h-by-n upper trapezoidal R
        real(dp), intent(in) :: u(:) !! m-vector for rank-one update
        real(dp), intent(in) :: v(:) !! n-vector; update with v transpose/conjugate transpose
        real(dp), intent(out) :: qnew(:, :) !! updated m-by-h Q
        real(dp), intent(out) :: rnew(:, :) !! updated h-by-n R
        integer, intent(out) :: status !! linalg_status_* code
        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (size(u) /= size(q,1) .or. size(v) /= size(r,2)) return
        call qr_update_real_matrix(q, r, reshape(u,[size(u),1]), &
            reshape(v,[size(v),1]), qnew, rnew, status)
    end subroutine qr_update_real_vector

    subroutine qr_update_real_matrix(q, r, u, v, qnew, rnew, status)
        real(dp), intent(in) :: q(:, :) !! m-by-h full or economic Q
        real(dp), intent(in) :: r(:, :) !! h-by-n upper trapezoidal R
        real(dp), intent(in) :: u(:, :) !! m-by-k rank-k left update
        real(dp), intent(in) :: v(:, :) !! n-by-k rank-k right update
        real(dp), intent(out) :: qnew(:, :) !! updated m-by-h Q
        real(dp), intent(out) :: rnew(:, :) !! updated h-by-n R
        integer, intent(out) :: status !! linalg_status_* code
        real(dp), allocatable :: w(:, :), residual(:, :), basis(:, :)
        real(dp), allocatable :: augmented(:, :), core(:, :), qcore(:, :), rcore(:, :)
        real(dp), allocatable :: trial(:), proj(:)
        integer :: m, n, h, k, maxnew, used, j, i, pass
        real(dp) :: length, original_length

        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_factors_real(q,r)) return
        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        k = size(u,2)
        if (size(qnew,1) /= m .or. size(qnew,2) /= h) return
        if (size(rnew,1) /= h .or. size(rnew,2) /= n) return
        if (size(u,1) /= m .or. size(v,1) /= n .or. size(v,2) /= k) return
        status = linalg_status_no_convergence
        if (.not. finite_real(q) .or. .not. finite_real(r)) return
        if (.not. finite_real(u) .or. .not. finite_real(v)) return
        if (k == 0) then
            qnew = q
            rnew = r
            status = linalg_status_success
            return
        end if
        ! The update lives in span(Q) plus at most k additional directions.
        ! Remove the Q projection, then construct a reorthogonalized basis.
        w = matmul(transpose(q),u)
        residual = u - matmul(q,w)
        maxnew = min(m-h,k)
        allocate(basis(m,maxnew),trial(m),proj(h))
        basis = 0.0_dp
        used = 0
        do j = 1, k
            if (used == maxnew) exit
            trial = residual(:,j)
            original_length = sqrt(sum(abs(u(:,j))**2))
            do pass = 1, 2
                proj = matmul(transpose(q),trial)
                trial = trial - matmul(q,proj)
                do i = 1, used
                    trial = trial - dot_product(basis(:,i),trial)*basis(:,i)
                end do
            end do
            length = sqrt(sum(abs(trial)**2))
            if (length <= 64.0_dp*epsilon(1.0_dp)*original_length) cycle
            used = used + 1
            basis(:,used) = trial/length
        end do
        allocate(augmented(m,h+used),core(h+used,n),qcore(h+used,h),rcore(h,n))
        augmented(:,1:h) = q
        if (used > 0) augmented(:,h+1:) = basis(:,1:used)
        core = 0.0_dp
        core(1:h,:) = r + matmul(w,transpose(v))
        if (used > 0) core(h+1:,:) = matmul(matmul(transpose(basis(:,1:used)),u),transpose(v))
        call qr(core,qcore,rcore,status)
        if (status /= linalg_status_success) return
        qnew = matmul(augmented,qcore)
        rnew = rcore
        if (.not. finite_real(qnew) .or. .not. finite_real(rnew)) then
            qnew = 0.0_dp
            rnew = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine qr_update_real_matrix

    ! Column insertions can be handled by a reduced QR on R augmented with
    ! projections of u. If the new columns do not expand the thin Q basis
    ! enough, fall back to a complete (but still correct) refactorization.
    subroutine qr_insert_cols_real(q, r, u, k, qnew, rnew, status)
        real(dp), intent(in) :: q(:, :) !! original Q (already validated)
        real(dp), intent(in) :: r(:, :) !! original R (already validated)
        real(dp), intent(in) :: u(:, :) !! m-by-p columns to insert
        integer, intent(in) :: k !! one-based insertion location
        real(dp), intent(out) :: qnew(:, :) !! m-by-hnew updated Q
        real(dp), intent(out) :: rnew(:, :) !! hnew-by-(n+p) updated R
        integer, intent(out) :: status !! linalg_status_* code
        real(dp), allocatable :: w(:, :), residual(:, :), basis(:, :), trial(:), proj(:)
        real(dp), allocatable :: augmented(:, :), core(:, :), qcore(:, :), rcore(:, :)
        real(dp), allocatable :: original(:, :), edited(:, :)
        integer :: m, n, h, p, hnew, maxnew, used, j, i, pass
        real(dp) :: length, original_length

        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        p = size(u,2)
        hnew = size(qnew,2)
        w = matmul(transpose(q),u)
        residual = u-matmul(q,w)
        maxnew = min(m-h,p)
        allocate(basis(m,maxnew),trial(m),proj(h))
        basis = 0.0_dp
        used = 0
        do j = 1,p
            if (used == maxnew) exit
            trial = residual(:,j)
            original_length = sqrt(sum(abs(u(:,j))**2))
            do pass = 1,2
                proj = matmul(transpose(q),trial)
                trial = trial-matmul(q,proj)
                do i = 1,used
                    trial = trial-dot_product(basis(:,i),trial)*basis(:,i)
                end do
            end do
            length = sqrt(sum(abs(trial)**2))
            if (length <= 64.0_dp*epsilon(1.0_dp)*original_length) cycle
            used = used+1
            basis(:,used) = trial/length
        end do
        if (h+used < hnew) then
            ! The inserted columns are dependent or near-dependent on Q.
            original = matmul(q,r)
            allocate(edited(m,n+p))
            edited(:,1:k-1) = original(:,1:k-1)
            edited(:,k:k+p-1) = u
            edited(:,k+p:n+p) = original(:,k:n)
            call qr(edited,qnew,rnew,status)
            return
        end if
        allocate(augmented(m,h+used),core(h+used,n+p), &
            qcore(h+used,hnew),rcore(hnew,n+p))
        augmented(:,1:h) = q
        if (used > 0) augmented(:,h+1:) = basis(:,1:used)
        core = 0.0_dp
        core(1:h,1:k-1) = r(:,1:k-1)
        core(1:h,k:k+p-1) = w
        core(1:h,k+p:n+p) = r(:,k:n)
        if (used > 0) core(h+1:,k:k+p-1) = matmul(transpose(basis(:,1:used)),u)
        call qr(core,qcore,rcore,status)
        if (status /= linalg_status_success) return
        qnew = matmul(augmented,qcore)
        rnew = rcore
        if (.not. finite_real(qnew) .or. .not. finite_real(rnew)) then
            qnew = 0.0_dp
            rnew = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine qr_insert_cols_real

    subroutine qr_insert_real_vector(q, r, u, k, qnew, rnew, status, which)
        real(dp), intent(in) :: q(:, :) !! m-by-h original Q
        real(dp), intent(in) :: r(:, :) !! h-by-n original R
        real(dp), intent(in) :: u(:) !! one row (n entries) or column (m entries)
        integer, intent(in) :: k !! one-based insertion position, 1 through dimension+1
        real(dp), intent(out) :: qnew(:, :) !! updated orthonormal factor
        real(dp), intent(out) :: rnew(:, :) !! updated upper trapezoidal factor
        integer, intent(out) :: status !! linalg_status_* code
        character(len=*), intent(in), optional :: which !! 'row' default or 'col'
        character(len=:), allocatable :: direction
        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        direction = 'row'
        if (present(which)) direction = trim(which)
        select case(direction)
        case('row')
            if (size(u) /= size(r,2)) return
            call qr_insert_real_matrix(q,r,reshape(u,[1,size(u)]),k,qnew,rnew,status,which=direction)
        case('col')
            if (size(u) /= size(q,1)) return
            call qr_insert_real_matrix(q,r,reshape(u,[size(u),1]),k,qnew,rnew,status,which=direction)
        end select
    end subroutine qr_insert_real_vector

    subroutine qr_insert_real_matrix(q, r, u, k, qnew, rnew, status, which)
        real(dp), intent(in) :: q(:, :) !! m-by-h original Q
        real(dp), intent(in) :: r(:, :) !! h-by-n original R
        real(dp), intent(in) :: u(:, :) !! p-by-n row block or m-by-p column block
        integer, intent(in) :: k !! one-based insertion position
        real(dp), intent(out) :: qnew(:, :) !! new full or economic Q
        real(dp), intent(out) :: rnew(:, :) !! new matching R
        integer, intent(out) :: status !! linalg_status_* code
        character(len=*), intent(in), optional :: which !! 'row' default or 'col'
        real(dp), allocatable :: original(:, :), edited(:, :)
        character(len=:), allocatable :: direction
        integer :: m, n, h, p, mnew, nnew, hnew
        logical :: full

        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_factors_real(q,r)) return
        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        full = h == m
        direction = 'row'
        if (present(which)) direction = trim(which)
        select case(direction)
        case('row')
            if (size(u,2) /= n) return
            p = size(u,1)
            if (k < 1 .or. k > m+1) return
            mnew = m+p
            nnew = n
        case('col')
            if (size(u,1) /= m) return
            p = size(u,2)
            if (k < 1 .or. k > n+1) return
            mnew = m
            nnew = n+p
        case default
            return
        end select
        hnew = min(mnew,nnew)
        if (full) hnew = mnew
        if (size(qnew,1) /= mnew .or. size(qnew,2) /= hnew) return
        if (size(rnew,1) /= hnew .or. size(rnew,2) /= nnew) return
        status = linalg_status_no_convergence
        if (.not. finite_real(q) .or. .not. finite_real(r)) return
        if (.not. finite_real(u)) return
        if (p == 0) then
            qnew = q
            rnew = r
            status = linalg_status_success
            return
        end if
        if (direction == 'col') then
            call qr_insert_cols_real(q,r,u,k,qnew,rnew,status)
            return
        end if
        original = matmul(q,r)
        allocate(edited(mnew,nnew))
        select case(direction)
        case('row')
            edited(1:k-1,:) = original(1:k-1,:)
            edited(k:k+p-1,:) = u
            edited(k+p:mnew,:) = original(k:m,:)
        case('col')
            edited(:,1:k-1) = original(:,1:k-1)
            edited(:,k:k+p-1) = u
            edited(:,k+p:nnew) = original(:,k:n)
        end select
        call qr(edited,qnew,rnew,status)
    end subroutine qr_insert_real_matrix

    subroutine qr_delete_real(q, r, k, qnew, rnew, status, p, which)
        real(dp), intent(in) :: q(:, :) !! m-by-h original Q
        real(dp), intent(in) :: r(:, :) !! h-by-n original R
        integer, intent(in) :: k !! one-based first row or column to remove
        real(dp), intent(out) :: qnew(:, :) !! updated full or economic Q
        real(dp), intent(out) :: rnew(:, :) !! updated upper trapezoidal R
        integer, intent(out) :: status !! linalg_status_* code
        integer, intent(in), optional :: p !! number of rows/columns to delete, default 1
        character(len=*), intent(in), optional :: which !! 'row' default or 'col'
        real(dp), allocatable :: edited(:, :), qcore(:, :), rcore(:, :)
        character(len=:), allocatable :: direction
        integer :: m, n, h, mnew, nnew, hnew, count
        logical :: full

        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_factors_real(q,r)) return
        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        full = h == m
        direction = 'row'
        if (present(which)) direction = trim(which)
        count = 1
        if (present(p)) count = p
        if (count < 1 .or. k < 1) return
        select case(direction)
        case('row')
            if (k+count-1 > m) return
            mnew = m-count
            nnew = n
        case('col')
            if (k+count-1 > n) return
            mnew = m
            nnew = n-count
        case default
            return
        end select
        hnew = min(mnew,nnew)
        if (full) hnew = mnew
        if (size(qnew,1) /= mnew .or. size(qnew,2) /= hnew) return
        if (size(rnew,1) /= hnew .or. size(rnew,2) /= nnew) return
        status = linalg_status_no_convergence
        if (.not. finite_real(q) .or. .not. finite_real(r)) return
        if (mnew == 0 .or. nnew == 0) then
            if (direction == 'col' .and. full) qnew = q
            status = linalg_status_success
            return
        end if
        if (direction == 'col') then
            allocate(edited(h,nnew),qcore(h,hnew),rcore(hnew,nnew))
            edited(:,1:k-1) = r(:,1:k-1)
            edited(:,k:nnew) = r(:,k+count:n)
            call qr(edited,qcore,rcore,status)
            if (status /= linalg_status_success) return
            qnew = matmul(q,qcore)
            rnew = rcore
        else
            allocate(edited(mnew,n))
            edited(1:k-1,:) = matmul(q(1:k-1,:),r)
            edited(k:mnew,:) = matmul(q(k+count:m,:),r)
            call qr(edited,qnew,rnew,status)
        end if
    end subroutine qr_delete_real

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! input values
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    pure logical function valid_factors_complex(q, r) result(ok)
        complex(dp), intent(in) :: q(:, :) !! existing orthonormal QR columns
        complex(dp), intent(in) :: r(:, :) !! existing upper trapezoidal factor
        integer :: m, h, n
        m = size(q, 1)
        h = size(q, 2)
        n = size(r, 2)
        ok = m >= 1 .and. n >= 1 .and. &
            (h == m .or. h == min(m,n)) .and. size(r,1) == h
    end function valid_factors_complex

    subroutine qr_update_complex_vector(q, r, u, v, qnew, rnew, status)
        complex(dp), intent(in) :: q(:, :) !! m-by-h full or economic Q
        complex(dp), intent(in) :: r(:, :) !! h-by-n upper trapezoidal R
        complex(dp), intent(in) :: u(:) !! m-vector for rank-one update
        complex(dp), intent(in) :: v(:) !! n-vector; update with v transpose/conjugate transpose
        complex(dp), intent(out) :: qnew(:, :) !! updated m-by-h Q
        complex(dp), intent(out) :: rnew(:, :) !! updated h-by-n R
        integer, intent(out) :: status !! linalg_status_* code
        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (size(u) /= size(q,1) .or. size(v) /= size(r,2)) return
        call qr_update_complex_matrix(q, r, reshape(u,[size(u),1]), &
            reshape(v,[size(v),1]), qnew, rnew, status)
    end subroutine qr_update_complex_vector

    subroutine qr_update_complex_matrix(q, r, u, v, qnew, rnew, status)
        complex(dp), intent(in) :: q(:, :) !! m-by-h full or economic Q
        complex(dp), intent(in) :: r(:, :) !! h-by-n upper trapezoidal R
        complex(dp), intent(in) :: u(:, :) !! m-by-k rank-k left update
        complex(dp), intent(in) :: v(:, :) !! n-by-k rank-k right update
        complex(dp), intent(out) :: qnew(:, :) !! updated m-by-h Q
        complex(dp), intent(out) :: rnew(:, :) !! updated h-by-n R
        integer, intent(out) :: status !! linalg_status_* code
        complex(dp), allocatable :: w(:, :), residual(:, :), basis(:, :)
        complex(dp), allocatable :: augmented(:, :), core(:, :), qcore(:, :), rcore(:, :)
        complex(dp), allocatable :: trial(:), proj(:)
        integer :: m, n, h, k, maxnew, used, j, i, pass
        real(dp) :: length, original_length

        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_factors_complex(q,r)) return
        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        k = size(u,2)
        if (size(qnew,1) /= m .or. size(qnew,2) /= h) return
        if (size(rnew,1) /= h .or. size(rnew,2) /= n) return
        if (size(u,1) /= m .or. size(v,1) /= n .or. size(v,2) /= k) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(q) .or. .not. finite_complex(r)) return
        if (.not. finite_complex(u) .or. .not. finite_complex(v)) return
        if (k == 0) then
            qnew = q
            rnew = r
            status = linalg_status_success
            return
        end if
        ! The update lives in span(Q) plus at most k additional directions.
        ! Remove the Q projection, then construct a reorthogonalized basis.
        w = matmul(transpose(conjg(q)),u)
        residual = u - matmul(q,w)
        maxnew = min(m-h,k)
        allocate(basis(m,maxnew),trial(m),proj(h))
        basis = 0.0_dp
        used = 0
        do j = 1, k
            if (used == maxnew) exit
            trial = residual(:,j)
            original_length = sqrt(sum(abs(u(:,j))**2))
            do pass = 1, 2
                proj = matmul(transpose(conjg(q)),trial)
                trial = trial - matmul(q,proj)
                do i = 1, used
                    trial = trial - dot_product(basis(:,i),trial)*basis(:,i)
                end do
            end do
            length = sqrt(sum(abs(trial)**2))
            if (length <= 64.0_dp*epsilon(1.0_dp)*original_length) cycle
            used = used + 1
            basis(:,used) = trial/length
        end do
        allocate(augmented(m,h+used),core(h+used,n),qcore(h+used,h),rcore(h,n))
        augmented(:,1:h) = q
        if (used > 0) augmented(:,h+1:) = basis(:,1:used)
        core = 0.0_dp
        core(1:h,:) = r + matmul(w,transpose(conjg(v)))
        if (used > 0) core(h+1:,:) = matmul(matmul(transpose(conjg(basis(:,1:used))),u),transpose(conjg(v)))
        call qr(core,qcore,rcore,status)
        if (status /= linalg_status_success) return
        qnew = matmul(augmented,qcore)
        rnew = rcore
        if (.not. finite_complex(qnew) .or. .not. finite_complex(rnew)) then
            qnew = 0.0_dp
            rnew = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine qr_update_complex_matrix

    ! Column insertions can be handled by a reduced QR on R augmented with
    ! projections of u. If the new columns do not expand the thin Q basis
    ! enough, fall back to a complete (but still correct) refactorization.
    subroutine qr_insert_cols_complex(q, r, u, k, qnew, rnew, status)
        complex(dp), intent(in) :: q(:, :) !! original Q (already validated)
        complex(dp), intent(in) :: r(:, :) !! original R (already validated)
        complex(dp), intent(in) :: u(:, :) !! m-by-p columns to insert
        integer, intent(in) :: k !! one-based insertion location
        complex(dp), intent(out) :: qnew(:, :) !! m-by-hnew updated Q
        complex(dp), intent(out) :: rnew(:, :) !! hnew-by-(n+p) updated R
        integer, intent(out) :: status !! linalg_status_* code
        complex(dp), allocatable :: w(:, :), residual(:, :), basis(:, :), trial(:), proj(:)
        complex(dp), allocatable :: augmented(:, :), core(:, :), qcore(:, :), rcore(:, :)
        complex(dp), allocatable :: original(:, :), edited(:, :)
        integer :: m, n, h, p, hnew, maxnew, used, j, i, pass
        real(dp) :: length, original_length

        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        p = size(u,2)
        hnew = size(qnew,2)
        w = matmul(transpose(conjg(q)),u)
        residual = u-matmul(q,w)
        maxnew = min(m-h,p)
        allocate(basis(m,maxnew),trial(m),proj(h))
        basis = 0.0_dp
        used = 0
        do j = 1,p
            if (used == maxnew) exit
            trial = residual(:,j)
            original_length = sqrt(sum(abs(u(:,j))**2))
            do pass = 1,2
                proj = matmul(transpose(conjg(q)),trial)
                trial = trial-matmul(q,proj)
                do i = 1,used
                    trial = trial-dot_product(basis(:,i),trial)*basis(:,i)
                end do
            end do
            length = sqrt(sum(abs(trial)**2))
            if (length <= 64.0_dp*epsilon(1.0_dp)*original_length) cycle
            used = used+1
            basis(:,used) = trial/length
        end do
        if (h+used < hnew) then
            ! The inserted columns are dependent or near-dependent on Q.
            original = matmul(q,r)
            allocate(edited(m,n+p))
            edited(:,1:k-1) = original(:,1:k-1)
            edited(:,k:k+p-1) = u
            edited(:,k+p:n+p) = original(:,k:n)
            call qr(edited,qnew,rnew,status)
            return
        end if
        allocate(augmented(m,h+used),core(h+used,n+p), &
            qcore(h+used,hnew),rcore(hnew,n+p))
        augmented(:,1:h) = q
        if (used > 0) augmented(:,h+1:) = basis(:,1:used)
        core = 0.0_dp
        core(1:h,1:k-1) = r(:,1:k-1)
        core(1:h,k:k+p-1) = w
        core(1:h,k+p:n+p) = r(:,k:n)
        if (used > 0) core(h+1:,k:k+p-1) = matmul(transpose(conjg(basis(:,1:used))),u)
        call qr(core,qcore,rcore,status)
        if (status /= linalg_status_success) return
        qnew = matmul(augmented,qcore)
        rnew = rcore
        if (.not. finite_complex(qnew) .or. .not. finite_complex(rnew)) then
            qnew = 0.0_dp
            rnew = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine qr_insert_cols_complex

    subroutine qr_insert_complex_vector(q, r, u, k, qnew, rnew, status, which)
        complex(dp), intent(in) :: q(:, :) !! m-by-h original Q
        complex(dp), intent(in) :: r(:, :) !! h-by-n original R
        complex(dp), intent(in) :: u(:) !! one row (n entries) or column (m entries)
        integer, intent(in) :: k !! one-based insertion position, 1 through dimension+1
        complex(dp), intent(out) :: qnew(:, :) !! updated orthonormal factor
        complex(dp), intent(out) :: rnew(:, :) !! updated upper trapezoidal factor
        integer, intent(out) :: status !! linalg_status_* code
        character(len=*), intent(in), optional :: which !! 'row' default or 'col'
        character(len=:), allocatable :: direction
        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        direction = 'row'
        if (present(which)) direction = trim(which)
        select case(direction)
        case('row')
            if (size(u) /= size(r,2)) return
            call qr_insert_complex_matrix(q,r,reshape(u,[1,size(u)]),k,qnew,rnew,status,which=direction)
        case('col')
            if (size(u) /= size(q,1)) return
            call qr_insert_complex_matrix(q,r,reshape(u,[size(u),1]),k,qnew,rnew,status,which=direction)
        end select
    end subroutine qr_insert_complex_vector

    subroutine qr_insert_complex_matrix(q, r, u, k, qnew, rnew, status, which)
        complex(dp), intent(in) :: q(:, :) !! m-by-h original Q
        complex(dp), intent(in) :: r(:, :) !! h-by-n original R
        complex(dp), intent(in) :: u(:, :) !! p-by-n row block or m-by-p column block
        integer, intent(in) :: k !! one-based insertion position
        complex(dp), intent(out) :: qnew(:, :) !! new full or economic Q
        complex(dp), intent(out) :: rnew(:, :) !! new matching R
        integer, intent(out) :: status !! linalg_status_* code
        character(len=*), intent(in), optional :: which !! 'row' default or 'col'
        complex(dp), allocatable :: original(:, :), edited(:, :)
        character(len=:), allocatable :: direction
        integer :: m, n, h, p, mnew, nnew, hnew
        logical :: full

        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_factors_complex(q,r)) return
        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        full = h == m
        direction = 'row'
        if (present(which)) direction = trim(which)
        select case(direction)
        case('row')
            if (size(u,2) /= n) return
            p = size(u,1)
            if (k < 1 .or. k > m+1) return
            mnew = m+p
            nnew = n
        case('col')
            if (size(u,1) /= m) return
            p = size(u,2)
            if (k < 1 .or. k > n+1) return
            mnew = m
            nnew = n+p
        case default
            return
        end select
        hnew = min(mnew,nnew)
        if (full) hnew = mnew
        if (size(qnew,1) /= mnew .or. size(qnew,2) /= hnew) return
        if (size(rnew,1) /= hnew .or. size(rnew,2) /= nnew) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(q) .or. .not. finite_complex(r)) return
        if (.not. finite_complex(u)) return
        if (p == 0) then
            qnew = q
            rnew = r
            status = linalg_status_success
            return
        end if
        if (direction == 'col') then
            call qr_insert_cols_complex(q,r,u,k,qnew,rnew,status)
            return
        end if
        original = matmul(q,r)
        allocate(edited(mnew,nnew))
        select case(direction)
        case('row')
            edited(1:k-1,:) = original(1:k-1,:)
            edited(k:k+p-1,:) = u
            edited(k+p:mnew,:) = original(k:m,:)
        case('col')
            edited(:,1:k-1) = original(:,1:k-1)
            edited(:,k:k+p-1) = u
            edited(:,k+p:nnew) = original(:,k:n)
        end select
        call qr(edited,qnew,rnew,status)
    end subroutine qr_insert_complex_matrix

    subroutine qr_delete_complex(q, r, k, qnew, rnew, status, p, which)
        complex(dp), intent(in) :: q(:, :) !! m-by-h original Q
        complex(dp), intent(in) :: r(:, :) !! h-by-n original R
        integer, intent(in) :: k !! one-based first row or column to remove
        complex(dp), intent(out) :: qnew(:, :) !! updated full or economic Q
        complex(dp), intent(out) :: rnew(:, :) !! updated upper trapezoidal R
        integer, intent(out) :: status !! linalg_status_* code
        integer, intent(in), optional :: p !! number of rows/columns to delete, default 1
        character(len=*), intent(in), optional :: which !! 'row' default or 'col'
        complex(dp), allocatable :: edited(:, :), qcore(:, :), rcore(:, :)
        character(len=:), allocatable :: direction
        integer :: m, n, h, mnew, nnew, hnew, count
        logical :: full

        qnew = 0.0_dp
        rnew = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_factors_complex(q,r)) return
        m = size(q,1)
        h = size(q,2)
        n = size(r,2)
        full = h == m
        direction = 'row'
        if (present(which)) direction = trim(which)
        count = 1
        if (present(p)) count = p
        if (count < 1 .or. k < 1) return
        select case(direction)
        case('row')
            if (k+count-1 > m) return
            mnew = m-count
            nnew = n
        case('col')
            if (k+count-1 > n) return
            mnew = m
            nnew = n-count
        case default
            return
        end select
        hnew = min(mnew,nnew)
        if (full) hnew = mnew
        if (size(qnew,1) /= mnew .or. size(qnew,2) /= hnew) return
        if (size(rnew,1) /= hnew .or. size(rnew,2) /= nnew) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(q) .or. .not. finite_complex(r)) return
        if (mnew == 0 .or. nnew == 0) then
            if (direction == 'col' .and. full) qnew = q
            status = linalg_status_success
            return
        end if
        if (direction == 'col') then
            allocate(edited(h,nnew),qcore(h,hnew),rcore(hnew,nnew))
            edited(:,1:k-1) = r(:,1:k-1)
            edited(:,k:nnew) = r(:,k+count:n)
            call qr(edited,qcore,rcore,status)
            if (status /= linalg_status_success) return
            qnew = matmul(q,qcore)
            rnew = rcore
        else
            allocate(edited(mnew,n))
            edited(1:k-1,:) = matmul(q(1:k-1,:),r)
            edited(k:mnew,:) = matmul(q(k+count:m,:),r)
            call qr(edited,qnew,rnew,status)
        end if
    end subroutine qr_delete_complex

end module scifort_linalg_qr_updates
