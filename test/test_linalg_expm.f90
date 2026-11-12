! SPDX-License-Identifier: MIT
program test_linalg_expm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_expm, only : expm
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp), allocatable :: a(:,:), x(:,:), y(:,:), neg(:,:), eye(:,:)
    real(dp), allocatable :: batch(:,:,:), result(:,:,:)
    complex(dp), allocatable :: z(:,:), w(:,:), zn(:,:), zeye(:,:)
    complex(dp), allocatable :: cbatch(:,:,:), cresult(:,:,:)
    real(dp) :: angle, error, nan_value
    integer :: n, v, i, j, s, count, b, status

    count = 0
    ! Zero-dimensional matrices have a valid empty exponential.
    allocate(a(0,0), x(0,0), z(0,0), w(0,0))
    call expm(a,x,status)
    if (status /= linalg_status_success) error stop 'empty real'
    call expm(z,w,status)
    if (status /= linalg_status_success) error stop 'empty complex'
    count = count+2
    deallocate(a,x,z,w)

    do n = 1, 12
        allocate(a(n,n),x(n,n),y(n,n),neg(n,n),eye(n,n))
        allocate(z(n,n),w(n,n),zn(n,n),zeye(n,n))
        eye = 0.0_dp
        zeye = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, n
            eye(i,i) = 1.0_dp
            zeye(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        do v = 1, 40
            a = 0.0_dp
            z = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, n
                do i = 1, n
                    a(i,j) = 0.2_dp*sin(real(31*i+7*j+11*v,dp))/sqrt(real(n,dp))
                    z(i,j) = cmplx(a(i,j), &
                        0.17_dp*cos(real(13*i+9*j+v,dp))/sqrt(real(n,dp)),dp)
                end do
            end do
            select case (mod(v,5))
            case (1)
                do i = 1, n
                    a(i,i) = a(i,i) + 0.3_dp
                    z(i,i) = z(i,i) + cmplx(0.3_dp,-0.1_dp,dp)
                end do
            case (2)
                a = a * 3.0_dp
                z = z * 3.0_dp
            case (3)
                a = a * 0.001_dp
                z = z * 0.001_dp
            case (4)
                a = a - transpose(a)
                z = z - transpose(conjg(z))
            end select
            call expm(a, x, status, s)
            if (status /= linalg_status_success) error stop 'real status'
            call expm(-a, neg, status)
            if (status /= linalg_status_success) error stop 'real inverse status'
            error = maxval(abs(matmul(x,neg)-eye))
            if (error > 1.0e-9_dp) then
                print *, 'real', n, v, error
                error stop 'exp(A)*exp(-A) not identity'
            end if
            call expm(z, w, status)
            if (status /= linalg_status_success) error stop 'complex status'
            call expm(-z,zn,status)
            if (status /= linalg_status_success) error stop 'complex inverse status'
            error = maxval(abs(matmul(w,zn)-zeye))
            if (error > 1.0e-9_dp) then
                print *, 'complex', n, v, error
                error stop 'complex exp(A)*exp(-A) not identity'
            end if
            count = count+2
        end do
        deallocate(a,x,y,neg,eye,z,w,zn,zeye)
    end do

    allocate(a(3,3),x(3,3),y(3,3),z(3,3),w(3,3))
    a = 0.0_dp
    a(1,2) = 2.0_dp
    a(2,3) = 3.0_dp
    call expm(a,x,status)
    if (status /= 0) error stop 'nilpotent status'
    y = 0.0_dp
    y(1,1)=1.0_dp
    y(2,2)=1.0_dp
    y(3,3)=1.0_dp
    y(1,2)=2.0_dp
    y(2,3)=3.0_dp
    y(1,3)=3.0_dp
    if (maxval(abs(x-y)) > 1.0e-13_dp) error stop 'nilpotent exponential'
    count = count+1
    deallocate(a,x,y,z,w)

    allocate(a(2,2),x(2,2),y(2,2),z(2,2),w(2,2))
    do i = 1, 8
        angle = 10.0_dp**(real(i-4,dp)/1.3_dp)
        a = 0.0_dp
        a(1,2) = -angle
        a(2,1) = angle
        call expm(a,x,status,s)
        if (status /= linalg_status_success) error stop 'rotation exponential'
        y(1,1) = cos(angle)
        y(1,2) = -sin(angle)
        y(2,1) = sin(angle)
        y(2,2) = cos(angle)
        if (maxval(abs(x-y)) > 2.0e-11_dp) then
            print *, 'rotation',angle,maxval(abs(x-y)),s
            error stop 'rotation inaccurate'
        end if
        z = cmplx(a, 0.5_dp*a,dp)
        call expm(z,w,status)
        if (status /= linalg_status_success) error stop 'complex rotation'
        count = count+2
    end do
    a = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,2) = -3.0_dp
    call expm(a,x,status)
    if (status /= linalg_status_success) error stop 'diagonal status'
    if (abs(x(1,1)-exp(2.0_dp)) > 1.0e-14_dp .or. &
            abs(x(2,2)-exp(-3.0_dp)) > 1.0e-14_dp) error stop 'diagonal exponential'
    z = cmplx(0.0_dp,0.0_dp,dp)
    z(1,1) = cmplx(0.1_dp,1.0_dp,dp)
    z(2,2) = cmplx(-0.2_dp,-2.0_dp,dp)
    call expm(z,w,status)
    if (status /= linalg_status_success) error stop 'complex diagonal status'
    if (abs(w(1,1)-exp(z(1,1))) > 1.0e-14_dp .or. &
            abs(w(2,2)-exp(z(2,2))) > 1.0e-14_dp) error stop 'complex diagonal'
    count = count+2

    ! Confirm one leading batch dimension, including zero batches.
    allocate(batch(7,2,2),result(7,2,2),cbatch(7,2,2),cresult(7,2,2))
    do b=1,7
        batch(b,:,:) = a*real(b,dp)*0.1_dp
        cbatch(b,:,:) = z*real(b,dp)*0.1_dp
    end do
    call expm(batch,result,status)
    if (status /= linalg_status_success) error stop 'batch real'
    call expm(cbatch,cresult,status)
    if (status /= linalg_status_success) error stop 'batch complex'
    do b=1,7
        call expm(batch(b,:,:),x,status)
        if (maxval(abs(x-result(b,:,:))) > 1.0e-14_dp) error stop 'batch real contents'
        call expm(cbatch(b,:,:),w,status)
        if (maxval(abs(w-cresult(b,:,:))) > 1.0e-14_dp) error stop 'batch complex contents'
        count = count+2
    end do
    deallocate(batch,result,cbatch,cresult)

    allocate(batch(0,2,2),result(0,2,2),cbatch(0,2,2),cresult(0,2,2))
    call expm(batch,result,status)
    if (status /= linalg_status_success) error stop 'zero batch real'
    call expm(cbatch,cresult,status)
    if (status /= linalg_status_success) error stop 'zero batch complex'
    deallocate(batch,result,cbatch,cresult)
    count = count+2

    ! Reject shape and nonfinite failures, and initialize output.
    allocate(neg(2,1))
    call expm(a,neg,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid output shape'
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    a(1,2) = nan_value
    call expm(a,x,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite input'
    if (any(abs(x) > 0.0_dp)) error stop 'nonfinite output not initialized'
    count = count+2

    print '(a,i0)', 'linalg_expm native checks passed: ',count
end program test_linalg_expm
