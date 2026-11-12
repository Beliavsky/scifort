! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_spatial
    use scifort_kinds, only : dp
    use scifort_signal_spatial, only : order_filter, sepfir2d, &
        spatial_status_success, spatial_status_invalid_input
    implicit none
    integer :: code, i
    integer :: x(5), image(3,3), cube(3,3,3)
    real(dp) :: rx(5), expected(3,3)
    real(dp), allocatable :: y(:), z(:,:), zr(:,:,:)
    integer, allocatable :: iy(:), iz(:,:), iz3(:,:,:)
    complex(dp), allocatable :: cz(:,:)

    x = [1,2,3,4,5]
    iy = order_filter(x, [1,1,1], 0, status=code)
    if (code /= spatial_status_success) error stop 'integer order status'
    if (any(iy /= [0,1,2,3,0])) error stop 'integer order rank 0'
    iy = order_filter(x, [1,1,1], 2)
    if (any(iy /= [2,3,4,5,5])) error stop 'integer order rank 2'
    iy = order_filter(x, [.true., .false., .true.], 1)
    if (any(iy /= [2,3,4,5,4])) error stop 'sparse logical domain'
    iy = order_filter([7,6,1,5,4,3,2], [0,0,0,1,1,1,0], 1)
    if (any(iy /= [6,5,4,4,3,2,0])) error stop 'sparse median footprint'
    rx = real(x, dp)
    y = order_filter(rx, [0.25_dp, 0.0_dp, -1.0_dp], 0)
    if (maxval(abs(y - [0.0_dp,1.0_dp,2.0_dp,3.0_dp,0.0_dp])) > 1.0e-12_dp) &
        error stop 'nonbinary real footprint'

    image = reshape([(i,i=1,9)], [3,3])
    iz = order_filter(image, reshape([1,0,0,0,1,0,0,0,1], [3,3]), 0)
    if (iz(2,2) /= 1 .or. iz(1,1) /= 0 .or. iz(3,3) /= 0) error stop 'diagonal footprint'
    cube = reshape([(i,i=1,27)], [3,3,3])
    iz3 = order_filter(cube, reshape([(.true., i=1,27)], [3,3,3]), 13)
    if (iz3(2,2,2) /= 14) error stop '3D order median'
    iz = order_filter(image(1:3:2,1:3:2), reshape([1,1,1,1,1,1,1,1,1], [3,3]), 4)
    if (any(shape(iz) /= [2,2])) error stop 'noncontiguous input'

    iy = order_filter(x, [1,1], 0, status=code)
    if (code /= spatial_status_invalid_input .or. any(iy /= 0)) error stop 'even footprint'
    iy = order_filter(x, [0,0,0], 0, status=code)
    if (code /= spatial_status_invalid_input) error stop 'empty footprint'
    iy = order_filter(x, [1,1,1], 3, status=code)
    if (code /= spatial_status_invalid_input) error stop 'out-of-range rank'
    iy = order_filter(x, [1,1,1], -1, status=code)
    if (code /= spatial_status_invalid_input) error stop 'negative rank'

    expected = reshape([4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp,9.0_dp,7.0_dp,8.0_dp,9.0_dp], [3,3])
    z = sepfir2d(real(image,dp), [1.0_dp,0.0_dp,0.0_dp], [0.0_dp,1.0_dp,0.0_dp], code)
    if (code /= spatial_status_success) error stop 'sepfir2d status'
    if (maxval(abs(z-expected)) > 1.0e-12_dp) error stop 'sepfir2d row mirror'
    z = sepfir2d(image, [0.0_dp,1.0_dp,0.0_dp], [1.0_dp,0.0_dp,0.0_dp])
    expected = reshape([2.0_dp,3.0_dp,3.0_dp,5.0_dp,6.0_dp,6.0_dp,8.0_dp,9.0_dp,9.0_dp], [3,3])
    if (maxval(abs(z-expected)) > 1.0e-12_dp) error stop 'sepfir2d col mirror'

    cz = sepfir2d(cmplx(real(image,dp),real(image,dp),dp), &
        [cmplx(0.0_dp,0.0_dp,dp), cmplx(0.0_dp,1.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp)], &
        [0.0_dp,2.0_dp,0.0_dp])
    if (maxval(abs(cz - cmplx(-2.0_dp,2.0_dp,dp)*real(image,dp))) > 1.0e-12_dp) &
        error stop 'complex taps and image'
    cz = sepfir2d(image, [cmplx(1.0_dp,1.0_dp,dp)], [cmplx(1.0_dp,0.0_dp,dp)])
    if (maxval(abs(cz - cmplx(1.0_dp,1.0_dp,dp)*real(image,dp))) > 1.0e-12_dp) &
        error stop 'integer input complex taps'
    z = sepfir2d(reshape([2.0_dp], [1,1]), [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
        [1.0_dp,2.0_dp,3.0_dp], code)
    if (code /= spatial_status_success) error stop 'overlong filter status'
    if (abs(z(1,1)-180.0_dp) > 1.0e-12_dp) error stop 'repeated reflection'
    z = sepfir2d(real(image,dp), [1.0_dp,1.0_dp], [1.0_dp], code)
    if (code /= spatial_status_invalid_input .or. any(z /= 0.0_dp)) error stop 'even-length taps'
    z = sepfir2d(real(image,dp), [real(dp) ::], [1.0_dp], code)
    if (code /= spatial_status_invalid_input) error stop 'empty taps'
    z = sepfir2d(reshape([real(dp) ::], [0,2]), [1.0_dp], [1.0_dp], code)
    if (code /= spatial_status_invalid_input) error stop 'empty image'
    print *, 'PASS test_signal_spatial'
end program test_signal_spatial
