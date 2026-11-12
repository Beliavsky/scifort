! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real state-feedback assignment for partially actuated MIMO systems.
! This implementation constructs candidate right eigenvectors from the null
! spaces of [A - lambda*I, -B]. Real KNV0 rank-one and Yang-Tits
! real-pair / complex-pair determinant updates improve the transfer matrix.
! The group sweep order and initialization differ from SciPy's reference.
module scifort_signal_place_poles_mimo
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : symmetric_eigen_jacobi, linalg_status_success
    implicit none
    private
    public :: place_poles_mimo_nullspace

    type :: pole_space
        complex(dp) :: lambda = (0.0_dp, 0.0_dp)
        complex(dp), allocatable :: x(:,:) !! orthonormal feasible eigenvector columns
        complex(dp), allocatable :: u(:,:) !! corresponding input vectors
        integer :: nvec = 0 !! usable nullspace dimension
        integer :: first = 0 !! initial column in the real transfer matrix
        integer :: width = 1 !! 1 for real, 2 for conjugate pair
    end type pole_space
contains
    subroutine place_poles_mimo_nullspace(a, b, poles, gain, status, method, &
        target_rtol, max_iterations, iterations, relative_change, det_transfer)
        real(dp), intent(in) :: a(:,:) !! state matrix (n,n)
        real(dp), intent(in) :: b(:,:) !! input matrix (n,m)
        complex(dp), intent(in) :: poles(:) !! requested poles, conjugate paired
        real(dp), allocatable, intent(out) :: gain(:,:) !! feedback K, m by n
        integer, intent(out) :: status !! 0 success, 2 uncontrollable, 4 numerical failure
        character(len=*), intent(in) :: method !! "YT" or "KNV0"
        real(dp), intent(in) :: target_rtol !! target relative determinant change
        integer, intent(in) :: max_iterations !! maximum determinant-improvement sweeps
        integer, intent(out) :: iterations !! number of sweeps performed
        real(dp), intent(out) :: relative_change !! achieved relative determinant change
        real(dp), intent(out) :: det_transfer !! absolute determinant of real transfer matrix
        type(pole_space), allocatable :: group(:)
        real(dp), allocatable :: x(:,:), u(:,:), solved(:,:), next_x(:,:), next_u(:,:)
        logical, allocatable :: used(:)
        logical :: ok
        integer :: i, j, ng, n, m, k, mate, sweep, rank_b
        real(dp) :: tol, residual, scale, before, after
        complex(dp) :: lambda

        status = 4
        iterations = 0
        relative_change = 0.0_dp
        det_transfer = 0.0_dp
        n = size(a, 1)
        m = size(b, 2)
        tol = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(poles)))
        rank_b = numerical_rank(b)
        if (rank_b == 0) then
            status = 2
            return
        end if
        do i = 1, n
            if (count(abs(poles - poles(i)) <= tol) > rank_b) then
                status = 1
                return
            end if
        end do

        allocate(group(n), used(n), x(n,n), u(m,n), next_x(n,2), next_u(m,2))
        used = .false.
        x = 0.0_dp
        u = 0.0_dp
        ng = 0
        k = 1
        do i = 1, n
            if (used(i)) cycle
            used(i) = .true.
            ng = ng + 1
            group(ng)%first = k
            lambda = poles(i)
            if (abs(aimag(lambda)) > tol) then
                mate = 0
                do j = i+1,n
                    if (used(j)) cycle
                    if (abs(poles(j) - conjg(lambda)) <= tol) then
                        mate = j
                        exit
                    end if
                end do
                if (mate == 0) then
                    status = 1
                    return
                end if
                used(mate) = .true.
                group(ng)%width = 2
            end if
            group(ng)%lambda = lambda
            call feasible_eigenspace(a, b, lambda, group(ng)%x, group(ng)%u, &
                group(ng)%nvec, ok)
            if (.not. ok .or. group(ng)%nvec == 0) then
                status = 2
                return
            end if
            call select_vectors(group(ng), x(:,1:k-1), next_x(:,1:group(ng)%width), &
                next_u(:,1:group(ng)%width), ok)
            if (.not. ok) then
                status = 2
                return
            end if
            x(:,k:k+group(ng)%width-1) = next_x(:,1:group(ng)%width)
            u(:,k:k+group(ng)%width-1) = next_u(:,1:group(ng)%width)
            k = k + group(ng)%width
        end do

        ! The rank-one KNV0 step projects the excluded transfer column onto
        ! its feasible eigenspace. YT also updates pairs of real columns with
        ! a dominant singular-vector pair and complex conjugate columns with
        ! a Hermitian eigenvector. This is a determinant-ascent iteration;
        ! unlike SciPy, groups are visited in their original input order.
        do sweep = 1, max_iterations
            before = transfer_log_volume(x)
            if (method == 'KNV0') then
                do i = 1, ng
                    call update_knv0(group(i), x, u)
                end do
            else
                call yt_sweep(group(:ng), x, u)
            end if
            after = transfer_log_volume(x)
            iterations = sweep
            relative_change = abs(1.0_dp-exp(min(0.0_dp, before-after)))
            if (relative_change <= target_rtol .and. target_rtol >= 0.0_dp) exit
        end do
        det_transfer = exp(max(log(tiny(1.0_dp)), transfer_log_volume(x)))

        call solve_real_system(transpose(x), transpose(u), solved, ok)
        if (.not. ok) then
            status = 4
            return
        end if
        gain = transpose(solved)
        if (.not. all(ieee_is_finite(gain))) then
            deallocate(gain)
            return
        end if

        ! Check the eigenvector relation without relying on an eigensolver.
        ! For a complex pair z=x+i*y, the feedback is real and the real
        ! transfer block obeys F*x = a*x-b*y, F*y=b*x+a*y.
        residual = 0.0_dp
        do i = 1, ng
            k = group(i)%first
            lambda = group(i)%lambda
            if (group(i)%width == 1) then
                residual = max(residual, maxval(abs(matmul(a,x(:,k)) - &
                    matmul(b,u(:,k)) - real(lambda,dp)*x(:,k))))
            else
                residual = max(residual, maxval(abs(matmul(a,x(:,k)) - &
                    matmul(b,u(:,k)) - real(lambda,dp)*x(:,k) + aimag(lambda)*x(:,k+1))))
                residual = max(residual, maxval(abs(matmul(a,x(:,k+1)) - &
                    matmul(b,u(:,k+1)) - real(lambda,dp)*x(:,k+1) - aimag(lambda)*x(:,k))))
            end if
        end do
        scale = max(1.0_dp, maxval(abs(a)), maxval(abs(b))*maxval(abs(gain)), &
            maxval(abs(poles)))
        if (.not. ieee_is_finite(residual) .or. residual > 1.0e-8_dp*scale) then
            deallocate(gain)
            return
        end if
        status = 0
    end subroutine place_poles_mimo_nullspace

    real(dp) function transfer_log_volume(x) result(log_volume)
        real(dp), intent(in) :: x(:,:) !! real transfer-matrix columns
        real(dp), allocatable :: q(:,:), v(:)
        real(dp) :: h, nv
        integer :: n, i, j, rep
        n = size(x, 1)
        allocate(q(n,n),v(n))
        q = 0.0_dp
        log_volume = 0.0_dp
        do i = 1,n
            v = x(:,i)
            do rep = 1,2
                do j = 1,i-1
                    h = dot_product(q(:,j),v)
                    v = v-h*q(:,j)
                end do
            end do
            nv = sqrt(sum(v**2))
            if (nv <= tiny(1.0_dp)) then
                log_volume = log(tiny(1.0_dp))
                return
            end if
            log_volume = log_volume+log(nv)
            q(:,i) = v/nv
        end do
    end function transfer_log_volume

    subroutine complement_two(x, first, width, q, nr)
        real(dp), intent(in) :: x(:,:) !! current real transfer matrix
        integer, intent(in) :: first !! first omitted column
        integer, intent(in) :: width !! number of omitted columns
        real(dp), intent(out) :: q(:,:) !! orthogonal complement, n by 2
        integer, intent(out) :: nr !! independent complement vectors
        real(dp), allocatable :: fixed(:,:), basis(:,:), v(:)
        real(dp) :: nv
        integer :: n, rank, i, j
        n = size(x,1)
        fixed = pack_other_columns(x, first, width)
        allocate(basis(n,n),v(n))
        call orthonormal_columns(fixed, basis(:,1:size(fixed,2)), rank)
        nr = 0
        q = 0.0_dp
        do i = 1,n
            v = 0.0_dp
            v(i) = 1.0_dp
            block
                real(dp) :: projected(n)
                call project_orthogonal(v,basis(:,1:rank),projected)
                v = projected
            end block
            do j = 1,nr
                v = v - dot_product(v,q(:,j))*q(:,j)
            end do
            nv = sqrt(sum(v**2))
            if (nv <= 1.0e-10_dp) cycle
            nr = nr+1
            q(:,nr) = v/nv
            if (nr == size(q,2)) exit
        end do
    end subroutine complement_two

    subroutine update_knv0(space, x, u)
        type(pole_space), intent(in) :: space !! real eigenspace
        real(dp), intent(inout) :: x(:,:) !! real transfer matrix
        real(dp), intent(inout) :: u(:,:) !! real control-input matrix
        real(dp), allocatable :: other(:,:), basis(:,:), v(:), q(:), coeff(:)
        real(dp) :: normv, before
        integer :: n, i, nr, k, p
        if (space%width /= 1) return
        n = size(x,1)
        p = space%nvec
        k = space%first
        other = pack_other_columns(x,k,1)
        allocate(basis(n,n-1),v(n),q(n),coeff(p))
        call orthonormal_columns(other,basis,nr)
        if (nr /= n-1) return
        q = 0.0_dp
        do i = 1,n
            v = 0.0_dp
            v(i) = 1.0_dp
            call project_orthogonal(v,basis,q)
            normv = sqrt(sum(q**2))
            if (normv > 1.0e-9_dp) exit
        end do
        if (normv <= 1.0e-9_dp) return
        q = q/normv
        do i = 1,p
            coeff(i) = dot_product(q,real(space%x(:,i),dp))
        end do
        v = matmul(real(space%x(:,1:p),dp),coeff)
        normv = sqrt(sum(v**2))
        if (normv <= 1.0e-12_dp) return
        before = transfer_log_volume(x)
        block
            real(dp), allocatable :: trial(:,:)
            trial = x
            trial(:,k) = v/normv
            if (transfer_log_volume(trial) > before + 1.0e-13_dp) then
                x(:,k) = trial(:,k)
                u(:,k) = real(matmul(space%u(:,1:p),cmplx(coeff,0.0_dp,dp)),dp)/normv
            end if
        end block
    end subroutine update_knv0

    subroutine complement_except_pair(x, i, j, q, nr)
        real(dp), intent(in) :: x(:,:) !! transfer columns
        integer, intent(in) :: i !! first excluded column
        integer, intent(in) :: j !! second excluded column
        real(dp), intent(out) :: q(:,:) !! orthonormal complement n by 2
        integer, intent(out) :: nr !! dimension of returned complement
        real(dp), allocatable :: fixed(:,:), basis(:,:), v(:)
        real(dp) :: nv
        integer :: n, k, c, rank, h
        n = size(x,1)
        allocate(fixed(n,n-2),basis(n,n-2),v(n))
        c = 0
        do k=1,n
            if (k == i .or. k == j) cycle
            c = c+1
            fixed(:,c) = x(:,k)
        end do
        call orthonormal_columns(fixed,basis,rank)
        nr = 0
        q = 0.0_dp
        do k=1,n
            v = 0.0_dp
            v(k) = 1.0_dp
            block
                real(dp) :: projected(n)
                call project_orthogonal(v,basis(:,1:rank),projected)
                v = projected
            end block
            do h=1,nr
                v = v-dot_product(v,q(:,h))*q(:,h)
            end do
            nv = sqrt(sum(v**2))
            if (nv <= 1.0e-10_dp) cycle
            nr = nr+1
            q(:,nr) = v/nv
            if (nr == 2) exit
        end do
    end subroutine complement_except_pair

    subroutine yt_sweep(group, x, u)
        type(pole_space), intent(in) :: group(:) !! feasible real and complex pole groups
        real(dp), intent(inout) :: x(:,:) !! real transfer matrix
        real(dp), intent(inout) :: u(:,:) !! matching input columns
        integer, allocatable :: re(:), cx(:)
        integer :: nr, nc, half, r, i, j, k, tmp
        nr = count(group%width == 1)
        nc = count(group%width == 2)
        allocate(re(nr),cx(nc))
        nr = 0
        nc = 0
        do k=1,size(group)
            if (group(k)%width == 1) then
                nr = nr+1
                re(nr) = k
            else
                nc = nc+1
                cx(nc) = k
            end if
        end do
        ! The real/complex pole grouping and update sequence follows
        ! the real-pair and conjugate-pair schedule of Tits--Yang.
        do i=2,nr
            tmp = re(i)
            j = i-1
            do while (j >= 1)
                if (real(group(re(j))%lambda,dp) <= real(group(tmp)%lambda,dp)) exit
                re(j+1) = re(j)
                j = j-1
            end do
            re(j+1) = tmp
        end do
        half = nr/2
        if (nr > 0) call real_pair(nr,1)
        do r=1, half+mod(nr,2)-1
            call real_pair(2*r,2*r+1)
        end do
        call complex_pairs()
        do r=1,half
            call real_pair(2*r-1,2*r)
        end do
        call single_real()
        call complex_pairs()
        do j=2, half+mod(nr,2)-1
            do i=1,half
                call real_pair(i,i+j)
            end do
        end do
        call single_real()
        call complex_pairs()
        do j=2,half+mod(nr,2)-1
            do i=half+1,nr
                k = i+j
                if (k > nr) k = k-nr
                call real_pair(i,k)
            end do
        end do
        call single_real()
        call complex_pairs()
        do i=1,half
            call real_pair(i,i+half)
        end do
        call single_real()
        call complex_pairs()
    contains
        subroutine real_pair(i,j)
            integer, intent(in) :: i,j !! one-based positions in sorted real poles
            if (i == j) then
                call update_knv0(group(re(i)),x,u)
            else
                call update_yt_real(group(re(i)),group(re(j)),x,u)
            end if
        end subroutine real_pair
        subroutine complex_pairs()
            integer :: h
            do h=1,nc
                call update_yt_complex(group(cx(h)),x,u)
            end do
        end subroutine complex_pairs
        subroutine single_real()
            if (nr == 1) call update_knv0(group(re(1)),x,u)
        end subroutine single_real
    end subroutine yt_sweep

    subroutine update_yt_real(space_i, space_j, x, u)
        type(pole_space), intent(in) :: space_i !! first real-pole space
        type(pole_space), intent(in) :: space_j !! second real-pole space
        real(dp), intent(inout) :: x(:,:) !! real transfer matrix
        real(dp), intent(inout) :: u(:,:) !! real control-input matrix
        real(dp), allocatable :: q(:,:), mat(:,:), gram(:,:), vals(:), vectors(:,:), &
            nu(:), mu(:), xi(:), xj(:), ui(:), uj(:)
        real(dp) :: s, old_volume
        integer :: n, p, r, i, j, nr, status
        n = size(x,1)
        p = space_i%nvec
        r = space_j%nvec
        i = space_i%first
        j = space_j%first
        allocate(q(n,2),mat(p,r),gram(r,r),vals(r),vectors(r,r),nu(r),mu(p), &
            xi(n),xj(n),ui(size(u,1)),uj(size(u,1)))
        call complement_except_pair(x,i,j,q,nr)
        if (nr /= 2) return
        mat = 0.0_dp
        do j = 1,r
            do i = 1,p
                mat(i,j) = dot_product(real(space_i%x(:,i),dp),q(:,1))* &
                    dot_product(real(space_j%x(:,j),dp),q(:,2)) - &
                    dot_product(real(space_i%x(:,i),dp),q(:,2))* &
                    dot_product(real(space_j%x(:,j),dp),q(:,1))
            end do
        end do
        gram = matmul(transpose(mat),mat)
        call symmetric_eigen_jacobi(gram,vals,vectors,status)
        if (status /= linalg_status_success) return
        s = sqrt(max(0.0_dp,vals(r)))
        if (s <= 1.0e-13_dp) return
        nu = vectors(:,r)
        mu = matmul(mat,nu)/s
        xi = matmul(real(space_i%x(:,1:p),dp),mu)
        xj = matmul(real(space_j%x(:,1:r),dp),nu)
        ui = matmul(real(space_i%u(:,1:p),dp),mu)
        uj = matmul(real(space_j%u(:,1:r),dp),nu)
        i = space_i%first
        j = space_j%first
        old_volume = transfer_log_volume(x)
        ! Keep the original columns until the prospective pair passes the
        ! determinant-ascent check.
        block
            real(dp), allocatable :: trial(:,:)
            trial = x
            trial(:,i) = xi
            trial(:,j) = xj
            if (transfer_log_volume(trial) > old_volume + 1.0e-13_dp) then
                x(:,i) = xi
                x(:,j) = xj
                u(:,i) = ui
                u(:,j) = uj
            end if
        end block
    end subroutine update_yt_real

    subroutine update_yt_complex(space, x, u)
        type(pole_space), intent(in) :: space !! conjugate-pole eigenspace
        real(dp), intent(inout) :: x(:,:) !! real transfer matrix
        real(dp), intent(inout) :: u(:,:) !! real input-transfer matrix
        real(dp), allocatable :: q(:,:), gr(:,:), vals(:), eigvec(:,:), augmented(:,:)
        complex(dp), allocatable :: herm(:,:), coeff(:), z(:), w(:)
        real(dp) :: vnorm, old_volume
        integer :: p, n, nr, a, b, j, istat, imax
        if (space%width /= 2) return
        n = size(x,1)
        p = space%nvec
        j = space%first
        allocate(q(n,2),herm(p,p),gr(2*p,2*p),vals(2*p),eigvec(2*p,2*p), &
            coeff(p),z(n),w(size(u,1)))
        call complement_two(x,j,2,q,nr)
        if (nr /= 2) return
        do b=1,p
            do a=1,p
                herm(a,b) = cmplx(0.0_dp,4.0_dp,dp)* &
                    (conjg(dot_product(cmplx(q(:,2),0.0_dp,dp),space%x(:,a))) * &
                     dot_product(cmplx(q(:,1),0.0_dp,dp),space%x(:,b)) - &
                     conjg(dot_product(cmplx(q(:,1),0.0_dp,dp),space%x(:,a))) * &
                     dot_product(cmplx(q(:,2),0.0_dp,dp),space%x(:,b)))
            end do
        end do
        gr(:p,:p) = real(herm,dp)
        gr(:p,p+1:) = -aimag(herm)
        gr(p+1:,:p) = aimag(herm)
        gr(p+1:,p+1:) = real(herm,dp)
        call symmetric_eigen_jacobi(gr,vals,eigvec,istat)
        if (istat /= linalg_status_success) return
        imax = maxloc(abs(vals),dim=1)
        if (abs(vals(imax)) <= 1.0e-13_dp) return
        coeff = cmplx(eigvec(:p,imax),eigvec(p+1:,imax),dp)
        vnorm = sqrt(sum(abs(coeff)**2))
        if (vnorm <= 1.0e-13_dp) return
        coeff = coeff/vnorm
        z = matmul(space%x(:,1:p),coeff)
        w = matmul(space%u(:,1:p),coeff)
        old_volume = transfer_log_volume(x)
        ! Only accept determinant-increasing updates, to protect against
        ! eigenvalue degeneracy in the 2*p real symmetric representation.
        augmented = x
        augmented(:,j) = real(z,dp)
        augmented(:,j+1) = aimag(z)
        if (transfer_log_volume(augmented) > old_volume + 1.0e-13_dp) then
            x(:,j:j+1) = augmented(:,j:j+1)
            u(:,j) = real(w,dp)
            u(:,j+1) = aimag(w)
        end if
    end subroutine update_yt_complex

    function pack_other_columns(x, first, width) result(other)
        real(dp), intent(in) :: x(:,:) !! complete transfer matrix
        integer, intent(in) :: first !! first excluded column
        integer, intent(in) :: width !! number of excluded columns
        real(dp), allocatable :: other(:,:) !! non-excluded columns
        integer :: n
        n = size(x,2)
        allocate(other(size(x,1), n-width))
        if (first > 1) other(:,1:first-1) = x(:,1:first-1)
        if (first+width <= n) other(:,first:) = x(:,first+width:n)
    end function pack_other_columns

    subroutine feasible_eigenspace(a, b, lambda, vx, vu, nr, ok)
        real(dp), intent(in) :: a(:,:) !! n-by-n state matrix
        real(dp), intent(in) :: b(:,:) !! n-by-m control matrix
        complex(dp), intent(in) :: lambda !! eigenvalue
        complex(dp), allocatable, intent(out) :: vx(:,:) !! orthonormal x-vectors
        complex(dp), allocatable, intent(out) :: vu(:,:) !! matching u-vectors
        integer, intent(out) :: nr !! number of usable independent vectors
        logical, intent(out) :: ok !! full row rank pencil found
        complex(dp), allocatable :: mat(:,:), tmp_row(:), full(:), v(:), w(:)
        complex(dp) :: t
        real(dp) :: largest, eps, norms
        integer, allocatable :: perm(:)
        integer :: n, m, nc, i, j, k, ip, jp, p, free_col

        ok = .false.
        nr = 0
        n = size(a,1)
        m = size(b,2)
        nc = n+m
        allocate(mat(n,nc), tmp_row(nc), full(nc), v(n), w(m), perm(nc))
        mat = (0.0_dp,0.0_dp)
        mat(:,1:n) = cmplx(a,0.0_dp,dp)
        do i = 1, n
            mat(i,i) = mat(i,i)-lambda
        end do
        mat(:,n+1:nc) = cmplx(-b,0.0_dp,dp)
        largest = maxval(abs(mat))
        eps = max(256.0_dp*epsilon(1.0_dp)*real(n,dp)*largest, tiny(1.0_dp))
        do i = 1, nc
            perm(i) = i
        end do
        ! Complete pivoting on the n-by-(n+m) constraint pencil.
        do k = 1, n
            largest = 0.0_dp
            ip = k
            jp = k
            do j = k, nc
                do i = k, n
                    if (abs(mat(i,j)) > largest) then
                        largest = abs(mat(i,j))
                        ip = i
                        jp = j
                    end if
                end do
            end do
            if (largest <= eps) return
            if (ip /= k) then
                tmp_row = mat(k,:)
                mat(k,:) = mat(ip,:)
                mat(ip,:) = tmp_row
            end if
            if (jp /= k) then
                tmp_row(1:n) = mat(:,k)
                mat(:,k) = mat(:,jp)
                mat(:,jp) = tmp_row(1:n)
                p = perm(k)
                perm(k) = perm(jp)
                perm(jp) = p
            end if
            t = mat(k,k)
            mat(k,:) = mat(k,:)/t
            do i = 1, n
                if (i == k) cycle
                t = mat(i,k)
                mat(i,:) = mat(i,:) - t*mat(k,:)
            end do
        end do

        allocate(vx(n,m), vu(m,m))
        vx = (0.0_dp,0.0_dp)
        vu = (0.0_dp,0.0_dp)
        do free_col = 1, m
            tmp_row = (0.0_dp,0.0_dp)
            tmp_row(1:n) = -mat(:,n+free_col)
            tmp_row(n+free_col) = (1.0_dp,0.0_dp)
            full = (0.0_dp,0.0_dp)
            do j = 1,nc
                full(perm(j)) = tmp_row(j)
            end do
            v = full(1:n)
            w = full(n+1:nc)
            ! Orthogonalize the x-component while preserving the equation
            ! (A-lambda*I)*x = B*u by applying the same updates to u.
            do p = 1, nr
                t = dot_product(vx(:,p),v)
                v = v-t*vx(:,p)
                w = w-t*vu(:,p)
            end do
            norms = sqrt(sum(abs(v)**2))
            if (norms <= 1.0e-11_dp*max(1.0_dp,sqrt(sum(abs(full)**2)))) cycle
            nr = nr+1
            vx(:,nr) = v/norms
            vu(:,nr) = w/norms
        end do
        ok = .true.
    end subroutine feasible_eigenspace

    subroutine select_vectors(space, other, x, u, ok)
        type(pole_space), intent(in) :: space !! feasible eigenvectors
        real(dp), intent(in) :: other(:,:) !! fixed transfer columns
        real(dp), intent(out) :: x(:,:) !! selected real transfer columns
        real(dp), intent(out) :: u(:,:) !! corresponding control input columns
        logical, intent(out) :: ok !! found independent columns
        real(dp), allocatable :: q(:,:), residual1(:), residual2(:), trial_x(:,:), trial_u(:,:), qbest(:), unit(:)
        complex(dp), allocatable :: coeff(:), z(:), w(:)
        real(dp) :: score, best, normz, det, t, weights(2), qscore
        integer :: n, m, p, i, j, it, nc, k, trials
        logical :: real_pole
        n = size(x,1)
        m = size(u,1)
        p = space%nvec
        nc = space%width
        real_pole = nc == 1
        allocate(q(n,size(other,2)), residual1(n), residual2(n), &
            trial_x(n,nc), trial_u(m,nc), coeff(p), z(n), w(m), qbest(n), unit(n))
        call orthonormal_columns(other,q,k)
        best = -1.0_dp
        x = 0.0_dp
        u = 0.0_dp
        ok = .false.

        trials = max(96,32*p)
        if (real_pole .and. k == n-1) then
            ! KNV0 rank-one projection: use the vector orthogonal to all
            ! other transfer columns to select the best feasible real x.
            qscore = -1.0_dp
            do i = 1,n
                unit = 0.0_dp
                unit(i) = 1.0_dp
                call project_orthogonal(unit,q(:,1:k),residual1)
                t = sqrt(sum(residual1**2))
                if (t > qscore) then
                    qscore = t
                    qbest = residual1
                end if
            end do
            if (qscore > 1.0e-12_dp) qbest = qbest/qscore
            trials = trials+1
        end if
        do it = 1,trials
            coeff = (0.0_dp,0.0_dp)
            if (it <= p) then
                coeff(it) = (1.0_dp,0.0_dp)
            else if (real_pole .and. k == n-1 .and. it == trials) then
                do j=1,p
                    coeff(j) = cmplx(dot_product(qbest,real(space%x(:,j),dp)),0.0_dp,dp)
                end do
            else
                do j=1,p
                    t = real(it-p,dp)
                    if (real_pole) then
                        coeff(j) = cmplx(sin((t+11*j)*0.731_dp) + &
                            cos((t*j+3*j)*0.419_dp),0.0_dp,dp)
                    else
                        coeff(j) = cmplx(sin((t+11*j)*0.731_dp)+cos((t*j+3*j)*0.419_dp), &
                            sin((t*j+7*j)*0.823_dp)+cos((t+13*j)*0.647_dp),dp)
                    end if
                end do
            end if
            z = matmul(space%x(:,1:p),coeff)
            w = matmul(space%u(:,1:p),coeff)
            normz = sqrt(sum(abs(z)**2))
            if (normz <= 1.0e-14_dp) cycle
            z = z/normz
            w = w/normz
            trial_x(:,1) = real(z,dp)
            trial_u(:,1) = real(w,dp)
            if (.not. real_pole) then
                trial_x(:,2) = aimag(z)
                trial_u(:,2) = aimag(w)
            end if
            call project_orthogonal(trial_x(:,1),q(:,1:k),residual1)
            if (real_pole) then
                score = sqrt(sum(residual1**2))
            else
                call project_orthogonal(trial_x(:,2),q(:,1:k),residual2)
                weights = [dot_product(residual1,residual1), dot_product(residual2,residual2)]
                det = weights(1)*weights(2)-dot_product(residual1,residual2)**2
                score = sqrt(max(0.0_dp,det))
            end if
            if (score > best) then
                best = score
                x = trial_x
                u = trial_u
            end if
        end do
        if (best > 1.0e-10_dp) ok = .true.
    end subroutine select_vectors

    subroutine orthonormal_columns(x, q, nr)
        real(dp), intent(in) :: x(:,:) !! input column set
        real(dp), intent(out) :: q(:,:) !! independent orthonormal columns
        integer, intent(out) :: nr !! number of independent columns
        real(dp) :: v(size(x,1)), h, nv
        integer :: i, j, rep
        nr = 0
        q = 0.0_dp
        do i = 1,size(x,2)
            v = x(:,i)
            do rep = 1,2
                do j = 1,nr
                    h = dot_product(q(:,j),v)
                    v = v-h*q(:,j)
                end do
            end do
            nv = sqrt(sum(v**2))
            if (nv <= 1.0e-12_dp*max(1.0_dp,sqrt(sum(x(:,i)**2)))) cycle
            nr = nr+1
            q(:,nr) = v/nv
        end do
    end subroutine orthonormal_columns

    subroutine project_orthogonal(v,q,out)
        real(dp), intent(in) :: v(:) !! column to project
        real(dp), intent(in) :: q(:,:) !! orthonormal basis for excluded subspace
        real(dp), intent(out) :: out(:) !! residual orthogonal vector
        integer :: rep, j
        out = v
        do rep = 1,2
            do j = 1,size(q,2)
                out = out-dot_product(q(:,j),out)*q(:,j)
            end do
        end do
    end subroutine project_orthogonal

    integer function numerical_rank(b) result(rank)
        real(dp), intent(in) :: b(:,:) !! real matrix
        real(dp), allocatable :: q(:,:), v(:)
        real(dp) :: nv, tol
        integer :: i, j, n, m
        n = size(b,1)
        m = size(b,2)
        allocate(q(n,min(n,m)),v(n))
        rank = 0
        q = 0.0_dp
        tol = 256.0_dp*epsilon(1.0_dp)*max(n,m)*maxval(abs(b))
        do j=1,m
            v = b(:,j)
            do i=1,rank
                v = v-dot_product(v,q(:,i))*q(:,i)
            end do
            nv = sqrt(sum(v**2))
            if (nv <= tol) cycle
            rank = rank+1
            if (rank > min(n,m)) exit
            q(:,rank) = v/nv
        end do
    end function numerical_rank

    subroutine solve_real_system(a, rhs, x, ok)
        real(dp), intent(in) :: a(:,:) !! n-by-n coefficient matrix
        real(dp), intent(in) :: rhs(:,:) !! n-by-p RHS matrix
        real(dp), allocatable, intent(out) :: x(:,:) !! solution
        logical, intent(out) :: ok !! nonsingular solve
        real(dp), allocatable :: mat(:,:), work(:,:), temp(:)
        real(dp) :: t, tol
        integer :: n, j, k, pivot, p
        ok = .false.
        n = size(a,1)
        p = size(rhs,2)
        if (size(a,2) /= n .or. size(rhs,1) /= n) return
        allocate(mat(n,n),work(n,p),temp(max(n,p)))
        mat = a
        work = rhs
        tol = 128.0_dp*epsilon(1.0_dp)*real(n,dp)*maxval(abs(a))
        do k=1,n
            pivot = k-1+maxloc(abs(mat(k:n,k)),dim=1)
            if (abs(mat(pivot,k)) <= tol) return
            if (pivot /= k) then
                temp(:n) = mat(k,:)
                mat(k,:) = mat(pivot,:)
                mat(pivot,:) = temp(:n)
                temp(:p) = work(k,:)
                work(k,:) = work(pivot,:)
                work(pivot,:) = temp(:p)
            end if
            t = mat(k,k)
            mat(k,:) = mat(k,:)/t
            work(k,:) = work(k,:)/t
            do j=1,n
                if (j == k) cycle
                t = mat(j,k)
                mat(j,:) = mat(j,:) - t*mat(k,:)
                work(j,:) = work(j,:) - t*work(k,:)
            end do
        end do
        if (.not. all(ieee_is_finite(work))) return
        x = work
        ok = .true.
    end subroutine solve_real_system
end module scifort_signal_place_poles_mimo
