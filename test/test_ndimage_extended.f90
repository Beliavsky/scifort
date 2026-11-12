! SPDX-License-Identifier: MIT
! New ndimage full-kernel, order-statistic, morphology, and distance tests.
program test_ndimage_extended
    use scifort_kinds, only : dp
    use scifort_ndimage, only : correlate_nd, convolve_nd, &
        median_filter_nd, rank_filter_nd, percentile_filter_nd, &
        minimum_filter_nd, maximum_filter_nd, prewitt_nd, sobel_nd, &
        binary_opening_nd, binary_closing_nd, binary_erosion_nd, &
        binary_dilation_nd, binary_fill_holes_nd, binary_propagation_nd, &
        distance_transform_edt_nd, ndimage_status_success, &
        ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    real(dp) :: x(5),result(5),a(5,6),r(5,6),r2(5,6),edt(5), &
        edt2(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: cz(2),cw(1),cr(2)
    logical :: mask(5,6),tmp(5,6),work(5,6),holes(5,6),seed(5),allowed(5), &
        reached(5),mask15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: status,i,j

    x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    call correlate_nd(x,[10.0_dp,1.0_dp],result,status,mode='constant')
    if (status/=0 .or. any(abs(result-[1.0_dp,12.0_dp,23.0_dp,34.0_dp,45.0_dp])>1.e-13_dp)) &
        error stop 'asymmetric even correlation'
    call convolve_nd(x,[10.0_dp,1.0_dp],result,status,mode='constant')
    if (status/=0 .or. any(abs(result-[21.0_dp,32.0_dp,43.0_dp,54.0_dp,5.0_dp])>1.e-13_dp)) &
        error stop 'asymmetric even convolution'
    cz=cmplx([1.0_dp,2.0_dp],[1.0_dp,-1.0_dp],dp)
    cw=[cmplx(2.0_dp,3.0_dp,dp)]
    call correlate_nd(cz,cw,cr,status)
    if (status/=0 .or. any(abs(cr-cz*conjg(cw(1)))>1.e-13_dp)) &
        error stop 'conjugated complex correlation weights'
    call convolve_nd(cz,cw,cr,status)
    if (status/=0 .or. any(abs(cr-cz*cw(1))>1.e-13_dp)) &
        error stop 'complex convolution weights'
    do j=1,6
        do i=1,5
            a(i,j)=real(modulo(i*7+j*11+i*j,17),dp)
        end do
    end do
    call rank_filter_nd(a,0,[3,2],r,status)
    if (status/=ndimage_status_success) error stop 'rank filter status'
    call minimum_filter_nd(a,[3,2],r2,status)
    if (status/=0 .or. any(r/=r2)) error stop 'rank zero is minimum'
    call rank_filter_nd(a,-1,[2,3],r,status)
    if (status/=0) error stop 'rank negative status'
    call maximum_filter_nd(a,[2,3],r2,status)
    if (status/=0 .or. any(r/=r2)) error stop 'rank minus one is maximum'
    call percentile_filter_nd(a,100.0_dp,[2,3],r,status)
    if (status/=0 .or. any(r/=r2)) error stop '100th percentile maximum'
    call median_filter_nd(a,3,r,status)
    if (status/=0) error stop 'median status'
    call percentile_filter_nd(a,50.0_dp,3,r2,status)
    if (status/=0 .or. any(r/=r2)) error stop 'median percentile identity'
    call percentile_filter_nd(a,101.0_dp,3,r,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid percentile'
    call rank_filter_nd(a,9,3,r,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid rank'
    call rank_filter_nd(a,0,[100,100],r,status)
    if (status/=ndimage_status_success) error stop 'wide rank window'
    call rank_filter_nd(a,0,[230,230],r,status)
    if (status/=ndimage_status_too_expensive) error stop 'oversized rank window'
    call prewitt_nd(a,1,r,status)
    if (status/=0) error stop 'prewitt status'
    call sobel_nd(a,1,r2,status)
    if (status/=0) error stop 'sobel status'
    call prewitt_nd(a,0,r,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid prewitt axis'

    mask=.false.
    mask(2:4,2:5)=.true.
    call binary_erosion_nd(mask,tmp,status)
    if (status/=0) error stop 'erosion status'
    call binary_dilation_nd(tmp,work,status)
    if (status/=0) error stop 'dilation status'
    call binary_opening_nd(mask,holes,status)
    if (status/=0 .or. any(holes.neqv.work)) error stop 'opening composition'
    call binary_dilation_nd(mask,tmp,status)
    if (status/=0) error stop 'dilation status'
    call binary_erosion_nd(tmp,work,status)
    if (status/=0) error stop 'erosion status'
    call binary_closing_nd(mask,holes,status)
    if (status/=0 .or. any(holes.neqv.work)) error stop 'closing composition'
    holes=.true.
    holes(3,3)=.false.
    call binary_fill_holes_nd(holes,work,status)
    if (status/=0 .or. .not.all(work)) error stop 'fill one interior hole'
    holes(1,1)=.false.
    call binary_fill_holes_nd(holes,work,status)
    if (status/=0 .or. work(1,1) .or. .not.work(3,3)) error stop 'preserve exterior'
    seed=[.false.,.true.,.false.,.false.,.false.]
    allowed=[.false.,.false.,.true.,.true.,.true.]
    call binary_propagation_nd(seed,allowed,reached,status)
    if (status/=0 .or. any(reached.neqv.[.false.,.true.,.true.,.true.,.true.])) &
        error stop 'mask propagation with seed outside mask'
    call binary_erosion_nd(mask,tmp,status,iterations=0)
    if (status/=0) error stop 'converged opening erosion'
    call binary_dilation_nd(tmp,work,status,iterations=0)
    if (status/=0) error stop 'converged opening dilation'
    call binary_opening_nd(mask,holes,status,iterations=0)
    if (status/=0 .or. any(holes.neqv.work)) error stop 'converged opening composition'
    call binary_propagation_nd(seed,allowed,reached,status,connectivity=0)
    if (status/=ndimage_status_invalid_input) error stop 'invalid propagation connectivity'

    seed=[.true.,.true.,.false.,.true.,.true.]
    call distance_transform_edt_nd(seed,edt,status)
    if (status/=0 .or. any(abs(edt-[2.0_dp,1.0_dp,0.0_dp,1.0_dp,2.0_dp])>1.e-13_dp)) &
        error stop 'distance transform 1d'
    call distance_transform_edt_nd(seed,edt,status,sampling=[1.5_dp])
    if (status/=0 .or. any(abs(edt-[3.0_dp,1.5_dp,0.0_dp,1.5_dp,3.0_dp])>1.e-13_dp)) &
        error stop 'sampled distance transform'
    call distance_transform_edt_nd(seed,edt,status,sampling=[-1.0_dp])
    if (status/=ndimage_status_invalid_input) error stop 'invalid sampling'
    call distance_transform_edt_nd(seed,edt,status,sampling=[1.0_dp,1.0_dp])
    if (status/=ndimage_status_invalid_input) error stop 'sampling dimension mismatch'
    mask15=.true.
    mask15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    call distance_transform_edt_nd(mask15,edt2,status)
    if (status/=0 .or. abs(edt2(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-1.0_dp)>1.e-13_dp) &
        error stop 'rank15 distance transform'
    print '(a)', 'PASS: extended ndimage filters, morphology, EDT and validation'
end program test_ndimage_extended
