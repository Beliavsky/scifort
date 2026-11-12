! SPDX-License-Identifier: MIT
program test_linalg_qr_updates
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr
    use scifort_linalg_qr_updates, only : qr_update, qr_insert, qr_delete
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp), allocatable :: a(:, :), q(:, :), r(:, :), q1(:, :), r1(:, :), u(:, :), v(:, :), expected(:, :)
    complex(dp), allocatable :: ac(:, :), qc(:, :), rc(:, :), qc1(:, :), rc1(:, :), uc(:, :), vc(:, :), ec(:, :)
    real(dp) :: tol, worst, saved
    integer :: m, n, h, m2, n2, h2, mode, op, p, k, i, j, s, count
    logical :: full
    character(len=3) :: which

    tol = 2.0e-9_dp
    worst = 0.0_dp
    count = 0
    do m = 2, 10
        do n = 2, 10
            allocate(a(m,n),ac(m,n))
            do j = 1, n
                do i = 1, m
                    a(i,j) = sin(real(5*i+7*j+i*j,dp)) + 0.3_dp*cos(real(i*13+j,dp))
                    ac(i,j) = cmplx(a(i,j),0.4_dp*cos(real(i*3+j*8,dp)),dp)
                end do
            end do
            do mode = 1, 2
                full = mode == 2
                h = min(m,n)
                if (full) h = m
                allocate(q(m,h),r(h,n),qc(m,h),rc(h,n))
                call qr(a,q,r,s)
                if (s /= linalg_status_success) error stop 'qr real factorization'
                call qr(ac,qc,rc,s)
                if (s /= linalg_status_success) error stop 'qr complex factorization'
                do op = 1, 6
                    select case(op)
                    case(1,2) ! rank-one and rank-k updates
                        p = 1
                        if (op == 2) p = 3
                        allocate(u(m,p),v(n,p),uc(m,p),vc(n,p))
                        do k = 1, p
                            do i = 1, m
                                u(i,k) = cos(real(i*5+k*7,dp))
                                uc(i,k) = cmplx(u(i,k),sin(real(i+k*3,dp)),dp)
                            end do
                            do j = 1, n
                                v(j,k) = sin(real(j*9+k*3,dp))
                                vc(j,k) = cmplx(v(j,k),cos(real(j*7+k,dp)),dp)
                            end do
                        end do
                        allocate(q1(m,h),r1(h,n),qc1(m,h),rc1(h,n))
                        expected = a + matmul(u,transpose(v))
                        ec = ac + matmul(uc,transpose(conjg(vc)))
                        if (op == 1) then
                            call qr_update(q,r,u(:,1),v(:,1),q1,r1,s)
                        else
                            call qr_update(q,r,u,v,q1,r1,s)
                        end if
                        if (s /= linalg_status_success) then
                            print *, 'qr_update real status',m,n,mode,op,s
                            error stop 'real update'
                        end if
                        call check_real(q1,r1,expected,tol,worst,count)
                        if (op == 1) then
                            call qr_update(qc,rc,uc(:,1),vc(:,1),qc1,rc1,s)
                        else
                            call qr_update(qc,rc,uc,vc,qc1,rc1,s)
                        end if
                        if (s /= linalg_status_success) then
                            print *, 'qr_update complex status',m,n,mode,op,s
                            error stop 'complex update'
                        end if
                        call check_complex(qc1,rc1,ec,tol,worst,count)
                        deallocate(u,v,uc,vc)
                    case(3,4) ! two inserted rows / columns
                        p = 2
                        k = 1 + n/2
                        which = 'col'
                        m2 = m
                        n2 = n+p
                        if (op == 3) then
                            which = 'row'
                            k = 1+m/2
                            m2 = m+p
                            n2 = n
                        end if
                        h2 = min(m2,n2)
                        if (h == m) h2 = m2
                        allocate(q1(m2,h2),r1(h2,n2),qc1(m2,h2),rc1(h2,n2))
                        if (op == 3) then
                            allocate(u(p,n),uc(p,n))
                            do j = 1, n
                                do i = 1, p
                                    u(i,j) = sin(real(3*i+j*13,dp))
                                    uc(i,j) = cmplx(u(i,j),cos(real(11*i+j,dp)),dp)
                                end do
                            end do
                            allocate(expected(m2,n2),ec(m2,n2))
                            expected(1:k-1,:) = a(1:k-1,:)
                            expected(k:k+p-1,:) = u
                            expected(k+p:,:) = a(k:,:)
                            ec(1:k-1,:) = ac(1:k-1,:)
                            ec(k:k+p-1,:) = uc
                            ec(k+p:,:) = ac(k:,:)
                        else
                            allocate(u(m,p),uc(m,p))
                            do j = 1, p
                                do i = 1, m
                                    u(i,j) = sin(real(3*i+j*13,dp))
                                    uc(i,j) = cmplx(u(i,j),cos(real(11*i+j,dp)),dp)
                                end do
                            end do
                            allocate(expected(m2,n2),ec(m2,n2))
                            expected(:,1:k-1) = a(:,1:k-1)
                            expected(:,k:k+p-1) = u
                            expected(:,k+p:) = a(:,k:)
                            ec(:,1:k-1) = ac(:,1:k-1)
                            ec(:,k:k+p-1) = uc
                            ec(:,k+p:) = ac(:,k:)
                        end if
                        call qr_insert(q,r,u,k,q1,r1,s,which)
                        if (s /= linalg_status_success) then
                            print *, 'real insertion status', m,n,mode,which,s,shape(q1),shape(r1)
                            error stop 'real insertion'
                        end if
                        call check_real(q1,r1,expected,tol,worst,count)
                        call qr_insert(qc,rc,uc,k,qc1,rc1,s,which)
                        if (s /= linalg_status_success) error stop 'complex insertion'
                        call check_complex(qc1,rc1,ec,tol,worst,count)
                        deallocate(u,uc)
                    case(5,6) ! delete two rows / columns
                        p = 2
                        k = 1
                        which = 'col'
                        m2 = m
                        n2 = n-p
                        if (op == 5) then
                            which = 'row'
                            m2 = m-p
                            n2 = n
                        end if
                        h2 = min(m2,n2)
                        if (h == m) h2 = m2
                        allocate(q1(m2,h2),r1(h2,n2),qc1(m2,h2),rc1(h2,n2))
                        allocate(expected(m2,n2),ec(m2,n2))
                        if (op == 5) then
                            expected = a(k+p:,:)
                            ec = ac(k+p:,:)
                        else
                            expected = a(:,k+p:)
                            ec = ac(:,k+p:)
                        end if
                        call qr_delete(q,r,k,q1,r1,s,p,which)
                        if (s /= linalg_status_success) error stop 'real deletion'
                        call check_real(q1,r1,expected,tol,worst,count)
                        call qr_delete(qc,rc,k,qc1,rc1,s,p,which)
                        if (s /= linalg_status_success) error stop 'complex deletion'
                        call check_complex(qc1,rc1,ec,tol,worst,count)
                    end select
                    if (allocated(q1)) deallocate(q1,r1,qc1,rc1)
                    if (allocated(expected)) deallocate(expected,ec)
                end do
                deallocate(q,r,qc,rc)
            end do
            deallocate(a,ac)
        end do
    end do
    ! Edge and validation tests
    allocate(a(5,3),q(5,5),r(5,3))
    do j = 1, 3
        do i = 1, 5
            a(i,j) = sin(real(i*j,dp))
        end do
    end do
    call qr(a,q,r,s)
    if (s /= linalg_status_success) error stop 'initial factor'
    allocate(q1(5,5),r1(5,0),expected(5,0))
    call qr_delete(q,r,1,q1,r1,s,p=3,which='col')
    if (s /= linalg_status_success) error stop 'all column deletion'
    if (maxval(abs(matmul(transpose(q1),q1)-eye(5))) > tol) error stop 'all column Q'
    deallocate(q1,r1,expected)
    allocate(q1(0,0),r1(0,3))
    call qr_delete(q,r,1,q1,r1,s,p=5,which='row')
    if (s /= linalg_status_success) error stop 'all row deletion'
    call qr_delete(q,r,1,q1,r1,s,p=6,which='row')
    if (s /= linalg_status_invalid_shape) error stop 'invalid row deletion'
    deallocate(q1,r1)
    allocate(q1(5,5),r1(5,3))
    call qr_delete(q,r,1,q1,r1,s,p=1,which='wrong')
    if (s /= linalg_status_invalid_shape) error stop 'invalid which'
    a(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call qr_insert(q,r,a(1,:),2,q1,r1,s,which='row')
    if (s /= linalg_status_invalid_shape) error stop 'wrong insertion output shape'
    ! A rank-zero update preserves the factors exactly.
    allocate(u(5,0),v(3,0))
    call qr_update(q,r,u,v,q1,r1,s)
    if (s /= linalg_status_success .or. any(q1 /= q) .or. any(r1 /= r)) &
        error stop 'zero-rank QR update'
    deallocate(u,v)
    ! A rank-one update that stays in Q's column space is well-defined.
    call qr_update(q,r,q(:,1),[0.2_dp,0.1_dp,-0.4_dp],q1,r1,s)
    if (s /= linalg_status_success) error stop 'dependent update'
    if (maxval(abs(matmul(q1,r1)-matmul(q,r)- &
        matmul(reshape(q(:,1),[5,1]),reshape([0.2_dp,0.1_dp,-0.4_dp],[1,3])))) > tol) &
        error stop 'dependent update reconstruction'
    ! Null inserts should copy the input factors.
    allocate(u(0,3))
    call qr_insert(q,r,u,2,q1,r1,s)
    if (s /= linalg_status_success .or. any(q1 /= q) .or. any(r1 /= r)) &
        error stop 'empty row insertion'
    deallocate(u)
    ! Unusable sizes / invalid indices / nonfinite inputs return a status.
    call qr_delete(q,r,0,q1,r1,s)
    if (s /= linalg_status_invalid_shape) error stop 'zero index accepted'
    call qr_update(q,r,[1.0_dp],[1.0_dp,1.0_dp,1.0_dp],q1,r1,s)
    if (s /= linalg_status_invalid_shape) error stop 'bad update shape accepted'
    saved = q(1,1)
    q(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call qr_update(q,r,[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
        [1.0_dp,2.0_dp,3.0_dp],q1,r1,s)
    if (s /= linalg_status_no_convergence) error stop 'nonfinite Q update'
    call qr_delete(q,r,2,q1,r1,s,p=0)
    if (s /= linalg_status_invalid_shape) error stop 'nonpositive deletion count'
    q(1,1) = saved
    ! Full-Q output dimensions are required after inserting a row.
    deallocate(q1,r1)
    allocate(q1(6,6),r1(6,3))
    call qr_insert(q,r,[1.0_dp,2.0_dp,3.0_dp],6,q1,r1,s)
    if (s /= linalg_status_success) error stop 'insert final row'
    ! A dependent column in economic mode exercises the dense QR fallback.
    deallocate(q,r,q1,r1,a)
    allocate(a(5,3),q(5,3),r(3,3),q1(5,4),r1(4,4))
    do j = 1,3
        do i = 1,5
            a(i,j) = sin(real(i*7+j*11,dp))
        end do
    end do
    call qr(a,q,r,s)
    if (s /= linalg_status_success) error stop 'thin QR before dependent insertion'
    call qr_insert(q,r,q(:,1),2,q1,r1,s,which='col')
    if (s /= linalg_status_success) error stop 'dependent column insertion'
    if (maxval(abs(matmul(q1,r1(:,1:4))- &
        reshape([a(:,1),q(:,1),a(:,2),a(:,3)],[5,4]))) > tol) &
        error stop 'dependent inserted column reconstruction'
    print '(A,I0,A,ES12.4)', 'qr_update/insert/delete native cases: ',count,', max error ',worst
contains
    function eye(n) result(x)
        integer, intent(in) :: n
        real(dp) :: x(n,n)
        integer :: idx
        x = 0.0_dp
        do idx = 1,n
            x(idx,idx) = 1.0_dp
        end do
    end function eye

    subroutine check_real(f,g,expected,tolerance,worst,count)
        real(dp), intent(in) :: f(:,:),g(:,:),expected(:,:),tolerance
        real(dp), intent(inout) :: worst
        integer, intent(inout) :: count
        real(dp) :: err, scale, ortho
        if (size(expected) == 0) return
        scale = max(1.0_dp,maxval(abs(expected)))
        err = maxval(abs(matmul(f,g)-expected))/scale
        ortho = maxval(abs(matmul(transpose(f),f)-eye(size(f,2))))
        worst = max(worst,err,ortho)
        if (max(err,ortho) > tolerance) then
            print *, 'real mismatch:',err,ortho,shape(f),shape(g),shape(expected)
            error stop 'real QR update results'
        end if
        count = count+1
    end subroutine check_real

    subroutine check_complex(f,g,expected,tolerance,worst,count)
        complex(dp), intent(in) :: f(:,:),g(:,:),expected(:,:)
        real(dp), intent(in) :: tolerance
        real(dp), intent(inout) :: worst
        integer, intent(inout) :: count
        real(dp) :: err, scale, ortho
        if (size(expected) == 0) return
        scale = max(1.0_dp,maxval(abs(expected)))
        err = maxval(abs(matmul(f,g)-expected))/scale
        ortho = maxval(abs(matmul(transpose(conjg(f)),f)-cmplx(eye(size(f,2)),0.0_dp,dp)))
        worst = max(worst,err,ortho)
        if (max(err,ortho) > tolerance) then
            print *, 'complex mismatch:',err,ortho,shape(f),shape(g),shape(expected)
            error stop 'complex QR update results'
        end if
        count = count+1
    end subroutine check_complex
end program test_linalg_qr_updates
