! SPDX-License-Identifier: MIT
! Stream protocol for independent scipy.linalg.convolution_matrix comparisons.
program driver_convolution_matrix
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_convolution_matrix, only : convolution_matrix, convolution_matrix_rows
    use scifort_linalg_matmul_convolution_matrix, only : matmul_convolution_matrix
    implicit none
    integer :: m, n, nrhs, coeff_type, rhs_type, status, ios, k, cases, s
    real(dp) :: max_constructor_error, max_product_error, norm, err
    character(5) :: mode
    complex(dp), allocatable :: ca(:), cx(:,:), cref(:,:), cyref(:,:)
    complex(dp), allocatable :: cm(:,:), cy(:,:)
    real(dp), allocatable :: ra(:), rx(:,:), rm(:,:), ry(:,:)

    cases=0
    max_constructor_error=0.0_dp
    max_product_error=0.0_dp
    do
        read(*,*,iostat=ios) m,n,nrhs,mode,coeff_type,rhs_type
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid reference case header'
        k=convolution_matrix_rows(m,n,trim(mode),status)
        if (status /= linalg_status_success) error stop 'invalid reference geometry'
        allocate(ca(m),cx(n,nrhs),cref(k,n),cyref(k,nrhs))
        read(*,*,iostat=ios) ca
        if (ios /= 0) error stop 'cannot read coefficients'
        read(*,*,iostat=ios) cx
        if (ios /= 0) error stop 'cannot read RHS'
        read(*,*,iostat=ios) cref
        if (ios /= 0) error stop 'cannot read expected convolution matrix'
        read(*,*,iostat=ios) cyref
        if (ios /= 0) error stop 'cannot read expected products'
        allocate(cm(k,n),cy(k,nrhs),ra(m),rx(n,nrhs),rm(k,n),ry(k,nrhs))
        ra = real(ca,dp)
        rx = real(cx,dp)

        if (coeff_type == 1) then
            call convolution_matrix(ra,rm,s,trim(mode))
            if (s /= linalg_status_success) error stop 'real constructor status'
            err = maxval(abs(cmplx(rm,kind=dp)-cref))
        else
            call convolution_matrix(ca,cm,s,trim(mode))
            if (s /= linalg_status_success) error stop 'complex constructor status'
            err = maxval(abs(cm-cref))
        end if
        norm = max(1.0e-200_dp,maxval(abs(cref)))
        max_constructor_error = max(max_constructor_error,err/norm)
        if (err/norm > 2.0e-13_dp) then
            print *, 'constructor mismatch', m,n,nrhs,mode,coeff_type,rhs_type,err/norm
            error stop 'SciPy constructor mismatch'
        end if

        if (coeff_type == 1 .and. rhs_type == 1) then
            call matmul_convolution_matrix(ra,rx,ry,s,trim(mode))
            if (s /= linalg_status_success) error stop 'real product status'
            err = maxval(abs(cmplx(ry,kind=dp)-cyref))
        else if (coeff_type == 1 .and. rhs_type == 2) then
            call matmul_convolution_matrix(ra,cx,cy,s,trim(mode))
            if (s /= linalg_status_success) error stop 'real-complex product status'
            err = maxval(abs(cy-cyref))
        else if (coeff_type == 2 .and. rhs_type == 1) then
            call matmul_convolution_matrix(ca,rx,cy,s,trim(mode))
            if (s /= linalg_status_success) error stop 'complex-real product status'
            err = maxval(abs(cy-cyref))
        else
            call matmul_convolution_matrix(ca,cx,cy,s,trim(mode))
            if (s /= linalg_status_success) error stop 'complex product status'
            err = maxval(abs(cy-cyref))
        end if
        norm = max(1.0e-200_dp,maxval(abs(cyref)))
        max_product_error = max(max_product_error,err/norm)
        if (err/norm > 1.0e-10_dp) then
            print *, 'product mismatch',m,n,nrhs,mode,coeff_type,rhs_type,err/norm
            error stop 'SciPy product mismatch'
        end if
        cases=cases+1
        deallocate(ca,cx,cref,cyref,cm,cy,ra,rx,rm,ry)
    end do
    print '(a,i0,a,es13.5,a,es13.5)', 'PASS cases=',cases, &
        ' max_constructor_relative_error=',max_constructor_error, &
        ' max_product_relative_error=',max_product_error
end program driver_convolution_matrix
