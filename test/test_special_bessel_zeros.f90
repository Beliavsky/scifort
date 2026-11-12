! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_zeros
    use scifort_special, only : bessel_zeros_status_invalid_argument, &
        bessel_zeros_status_success, dp, jn_zeros, jnjnp_zeros, jnp_zeros, jnyn_zeros, &
        jv, jvp, y0_zeros, y1_zeros, y1p_zeros, yn_zeros, ynp_zeros, yv, yvp
    implicit none

    real(dp), allocatable :: jz(:)
    real(dp), allocatable :: jpz(:)
    real(dp), allocatable :: yz(:)
    real(dp), allocatable :: ypz(:)
    real(dp), allocatable :: companions(:)
    real(dp), allocatable :: merged(:)
    integer, allocatable :: orders(:)
    integer, allocatable :: serials(:)
    integer, allocatable :: kinds(:)
    integer :: status

    real(dp), parameter :: j0_ref(5) = [ &
        2.4048255576957724_dp, 5.520078110286311_dp, 8.653727912911013_dp, &
        11.791534439014281_dp, 14.930917708487787_dp]
    real(dp), parameter :: jp0_ref(5) = [ &
        3.8317059702075125_dp, 7.015586669815619_dp, 10.173468135062722_dp, &
        13.323691936314223_dp, 16.470630050877630_dp]
    real(dp), parameter :: y0_ref(5) = [ &
        0.8935769662791675_dp, 3.957678419314858_dp, 7.086051060301773_dp, &
        10.222345043496416_dp, 13.361097473872764_dp]
    real(dp), parameter :: yp0_ref(5) = [ &
        2.197141326031017_dp, 5.429681040794135_dp, 8.596005868331169_dp, &
        11.749154830839881_dp, 14.897442128336726_dp]
    real(dp), parameter :: j2_ref(5) = [ &
        5.135622301840683_dp, 8.417244140399866_dp, 11.619841172149060_dp, &
        14.795951782351262_dp, 17.959819494987826_dp]
    real(dp), parameter :: jp2_ref(5) = [ &
        3.0542369282271404_dp, 6.706133194158459_dp, 9.969467823087596_dp, &
        13.170370856016124_dp, 16.347522318321783_dp]
    real(dp), parameter :: y2_ref(5) = [ &
        3.384241767149594_dp, 6.793807513268268_dp, 10.023477979360038_dp, &
        13.209986710206417_dp, 16.378966558947457_dp]
    real(dp), parameter :: yp2_ref(5) = [ &
        5.002582931446065_dp, 8.350724701413078_dp, 11.574195465217645_dp, &
        14.760909306207676_dp, 17.931285939466857_dp]
    real(dp), parameter :: j100_ref(3) = [ &
        108.83616589840977_dp, 115.73935123918876_dp, 121.57533101701064_dp]
    real(dp), parameter :: jp100_ref(3) = [ &
        103.76837768254228_dp, 112.38667144602339_dp, 118.69505734057492_dp]
    real(dp), parameter :: y100_ref(3) = [ &
        104.38020425686610_dp, 112.48596483026810_dp, 118.74455848574802_dp]
    real(dp), parameter :: yp100_ref(3) = [ &
        108.65470181072872_dp, 115.67276549379478_dp, 121.53617203885214_dp]
    real(dp), parameter :: jnjnp_ref(20) = [ &
        0.0_dp, 1.8411837813406593_dp, 2.404825557695773_dp, &
        3.054236928227140_dp, 3.831705970207512_dp, 3.8317059702075125_dp, &
        4.201188941210528_dp, 5.135622301840682_dp, 5.317553126083994_dp, &
        5.3314427735250325_dp, 5.520078110286311_dp, 6.380161895923984_dp, &
        6.415616375700240_dp, 6.706133194158459_dp, 7.015586669815619_dp, &
        7.015586669815619_dp, 7.501266144684147_dp, 7.5883424345038035_dp, &
        8.015236598375951_dp, 8.417244140399864_dp]
    integer, parameter :: jnjnp_n_ref(20) = [ &
        0, 1, 0, 2, 1, 0, 3, 2, 4, 1, 0, 3, 5, 2, 1, 0, 6, 4, 3, 2]
    integer, parameter :: jnjnp_m_ref(20) = [ &
        0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 2, 2]
    integer, parameter :: jnjnp_t_ref(20) = [ &
        1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0]

    real(dp), parameter :: y0_special_ref(4) = [ &
        0.8935769662791675_dp, 3.957678419314858_dp, &
        7.086051060301773_dp, 10.222345043496416_dp]
    real(dp), parameter :: y0_derivative_ref(4) = [ &
        0.8794208024971950_dp, -0.4025426717750242_dp, &
        0.3000976149104752_dp, -0.24970123751468482_dp]
    real(dp), parameter :: y1_special_ref(4) = [ &
        2.197141326031017_dp, 5.429681040794135_dp, &
        8.596005868331169_dp, 11.749154830839881_dp]
    real(dp), parameter :: y1_derivative_ref(4) = [ &
        0.5207864124022674_dp, -0.3403180455234407_dp, &
        0.27145987731153354_dp, -0.23246176601703883_dp]
    real(dp), parameter :: y1p_special_ref(4) = [ &
        3.6830228565851777_dp, 6.9414999536541755_dp, &
        10.123404655436612_dp, 13.285758156782855_dp]
    real(dp), parameter :: y1_value_ref(4) = [ &
        0.4167299281064514_dp, -0.30317374013748943_dp, &
        0.25091253627781257_dp, -0.21897478790909955_dp]

    call y0_zeros(4, yz, companions, status)
    call require(status == bessel_zeros_status_success, 'y0_zeros status')
    call check_array_rel('y0_zeros roots', yz, y0_special_ref, 3.0e-13_dp)
    call check_array_rel('y0_zeros derivatives', companions, y0_derivative_ref, 3.0e-13_dp)

    call y1_zeros(4, yz, companions, status)
    call require(status == bessel_zeros_status_success, 'y1_zeros status')
    call check_array_rel('y1_zeros roots', yz, y1_special_ref, 3.0e-13_dp)
    call check_array_rel('y1_zeros derivatives', companions, y1_derivative_ref, 3.0e-13_dp)

    call y1p_zeros(4, ypz, companions, status)
    call require(status == bessel_zeros_status_success, 'y1p_zeros status')
    call check_array_rel('y1p_zeros roots', ypz, y1p_special_ref, 3.0e-13_dp)
    call check_array_rel('y1p_zeros values', companions, y1_value_ref, 3.0e-13_dp)

    call y0_zeros(141, yz, companions, status)
    call require(status == bessel_zeros_status_success, 'y0_zeros 141 status')
    call check_rel('y0_zeros 141 root', yz(141), 440.60865336367556_dp, 2.0e-12_dp)
    call check_rel('y0_zeros 141 derivative', companions(141), 0.038011384661414446_dp, 2.0e-12_dp)

    call y1_zeros(141, yz, companions, status)
    call require(status == bessel_zeros_status_success, 'y1_zeros 141 status')
    call check_rel('y1_zeros 141 root', yz(141), 442.17831792051356_dp, 2.0e-12_dp)
    call check_rel('y1_zeros 141 derivative', companions(141), 0.037943808813186486_dp, 2.0e-12_dp)

    call y1p_zeros(141, ypz, companions, status)
    call require(status == bessel_zeros_status_success, 'y1p_zeros 141 status')
    call check_rel('y1p_zeros 141 root', ypz(141), 443.74799047555996_dp, 2.0e-12_dp)
    call check_rel('y1p_zeros 141 value', companions(141), 0.03787668826874835_dp, 2.0e-12_dp)

    call y0_zeros(0, yz, companions, status)
    call require(status == bessel_zeros_status_invalid_argument, 'y0_zeros invalid status')
    call require(size(yz) == 0 .and. size(companions) == 0, 'y0_zeros invalid result')
    call y1_zeros(0, yz, companions, status)
    call require(status == bessel_zeros_status_invalid_argument, 'y1_zeros invalid status')
    call require(size(yz) == 0 .and. size(companions) == 0, 'y1_zeros invalid result')
    call y1p_zeros(0, ypz, companions, status)
    call require(status == bessel_zeros_status_invalid_argument, 'y1p_zeros invalid status')
    call require(size(ypz) == 0 .and. size(companions) == 0, 'y1p_zeros invalid result')

    call jnjnp_zeros(20, merged, orders, serials, kinds, status)
    call require(status == bessel_zeros_status_success, 'jnjnp status')
    call check_array_rel('jnjnp zeros', merged, jnjnp_ref, 4.0e-13_dp)
    call require(all(orders == jnjnp_n_ref), 'jnjnp orders')
    call require(all(serials == jnjnp_m_ref), 'jnjnp serials')
    call require(all(kinds == jnjnp_t_ref), 'jnjnp kinds')
    call check_jnjnp_residuals(merged, orders, kinds)

    call jnjnp_zeros(1200, merged, orders, serials, kinds, status)
    call require(status == bessel_zeros_status_success, 'jnjnp 1200 status')
    call require(size(merged) == 1200, 'jnjnp 1200 size')
    call check_rel('jnjnp 1200 last', merged(1200), 68.56313856707126_dp, 3.0e-12_dp)
    call require(orders(1200) == 32, 'jnjnp 1200 last order')
    call require(serials(1200) == 9, 'jnjnp 1200 last serial')
    call require(kinds(1200) == 1, 'jnjnp 1200 last kind')
    call require(maxval(orders) == 65, 'jnjnp 1200 maximum order')
    call require(all(merged(2:) >= merged(:1199)), 'jnjnp 1200 ordering')

    call jnjnp_zeros(0, merged, orders, serials, kinds, status)
    call require(status == bessel_zeros_status_invalid_argument, 'jnjnp invalid zero status')
    call require(size(merged) == 0 .and. size(orders) == 0, 'jnjnp invalid zero result')
    call jnjnp_zeros(1201, merged, orders, serials, kinds, status)
    call require(status == bessel_zeros_status_invalid_argument, 'jnjnp invalid large status')
    call require(size(merged) == 0 .and. size(kinds) == 0, 'jnjnp invalid large result')

    call jnyn_zeros(0, 5, jz, jpz, yz, ypz, status)
    call require(status == bessel_zeros_status_success, 'jnyn n0 status')
    call check_array_rel('jn n0', jz, j0_ref, 2.0e-13_dp)
    call check_array_rel('jnp n0', jpz, jp0_ref, 2.0e-13_dp)
    call check_array_rel('yn n0', yz, y0_ref, 2.0e-13_dp)
    call check_array_rel('ynp n0', ypz, yp0_ref, 2.0e-13_dp)

    call jnyn_zeros(2, 5, jz, jpz, yz, ypz, status)
    call require(status == bessel_zeros_status_success, 'jnyn n2 status')
    call check_array_rel('jn n2', jz, j2_ref, 3.0e-13_dp)
    call check_array_rel('jnp n2', jpz, jp2_ref, 3.0e-13_dp)
    call check_array_rel('yn n2', yz, y2_ref, 3.0e-13_dp)
    call check_array_rel('ynp n2', ypz, yp2_ref, 3.0e-13_dp)

    call jnyn_zeros(100, 3, jz, jpz, yz, ypz, status)
    call require(status == bessel_zeros_status_success, 'jnyn n100 status')
    call check_array_rel('jn n100', jz, j100_ref, 2.0e-12_dp)
    call check_array_rel('jnp n100', jpz, jp100_ref, 2.0e-12_dp)
    call check_array_rel('yn n100', yz, y100_ref, 2.0e-12_dp)
    call check_array_rel('ynp n100', ypz, yp100_ref, 2.0e-12_dp)


    call jnyn_zeros(1000, 3, jz, jpz, yz, ypz, status)
    call require(status == bessel_zeros_status_success, 'jnyn n1000 status')
    call check_array_rel('jn n1000', jz, [ &
        1018.6608809679079_dp, 1032.7618089413058_dp, 1044.3924299671173_dp], 3.0e-12_dp)
    call check_array_rel('jnp n1000', jpz, [ &
        1008.0933633200711_dp, 1025.9763876909700_dp, 1038.6934062914615_dp], 3.0e-12_dp)
    call check_array_rel('yn n1000', yz, [ &
        1009.3418149978422_dp, 1026.1648638987350_dp, 1038.7824519542357_dp], 3.0e-12_dp)
    call check_array_rel('ynp n1000', ypz, [ &
        1018.3049219543557_dp, 1032.6389583061684_dp, 1044.3235911789168_dp], 3.0e-12_dp)

    call jnyn_zeros(0, 1200, jz, jpz, yz, ypz, status)
    call require(status == bessel_zeros_status_success, 'jnyn nt1200 status')
    call check_rel('jn nt1200 last', jz(1200), 3769.125819308542_dp, 2.0e-12_dp)
    call check_rel('jnp nt1200 last', jpz(1200), 3770.696483020029_dp, 2.0e-12_dp)
    call check_rel('yn nt1200 last', yz(1200), 3767.5550229955743_dp, 2.0e-12_dp)
    call check_rel('ynp nt1200 last', ypz(1200), 3769.1256866517874_dp, 2.0e-12_dp)

    call jn_zeros(-3, 4, jz, status)
    call require(status == bessel_zeros_status_success, 'negative order status')
    call check_array_rel('negative order jn zeros', jz, [ &
        6.380161895923984_dp, 9.761023129981670_dp, &
        13.015200721698434_dp, 16.223466160318768_dp], 4.0e-13_dp)

    call jn_zeros(20, 12, jz, status)
    call require(status == bessel_zeros_status_success, 'jn residual status')
    call check_residuals_j(20, jz, .false.)
    call jnp_zeros(20, 12, jpz, status)
    call require(status == bessel_zeros_status_success, 'jnp residual status')
    call check_residuals_j(20, jpz, .true.)
    call yn_zeros(20, 12, yz, status)
    call require(status == bessel_zeros_status_success, 'yn residual status')
    call check_residuals_y(20, yz, .false.)
    call ynp_zeros(20, 12, ypz, status)
    call require(status == bessel_zeros_status_success, 'ynp residual status')
    call check_residuals_y(20, ypz, .true.)

    call jn_zeros(0, 0, jz, status)
    call require(status == bessel_zeros_status_invalid_argument, 'invalid nt status')
    call require(size(jz) == 0, 'invalid nt empty result')

    print *, 'test_special_bessel_zeros: PASS'

contains

    subroutine check_jnjnp_residuals(zeros, n, kind)
        real(dp), intent(in) :: zeros(:) !! merged J/J-prime zeros
        integer, intent(in) :: n(:) !! Bessel order for each zero
        integer, intent(in) :: kind(:) !! zero for J and one for J-prime

        integer :: i
        real(dp) :: residual

        do i = 1, size(zeros)
            if (zeros(i) == 0.0_dp) cycle
            if (kind(i) == 0) then
                residual = jv(real(n(i), dp), zeros(i))
            else
                residual = jvp(real(n(i), dp), zeros(i))
            end if
            call require(abs(residual) <= 2.0e-12_dp, 'jnjnp residual')
        end do
    end subroutine check_jnjnp_residuals

    subroutine check_rel(label, actual, expected, tol)
        character(*), intent(in) :: label !! assertion label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! relative tolerance

        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            write(*, '(a,2(1x,es24.16))') trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_array_rel(label, actual, expected, tol)
        character(*), intent(in) :: label !! assertion label
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tol !! relative tolerance

        integer :: i

        call require(size(actual) == size(expected), trim(label)//' size')
        do i = 1, size(expected)
            if (abs(actual(i) - expected(i)) > tol * max(1.0_dp, abs(expected(i)))) then
                write(*, '(a,1x,i0,2(1x,es24.16))') trim(label), i, actual(i), expected(i)
                error stop 1
            end if
        end do
    end subroutine check_array_rel

    subroutine check_residuals_j(n, zeros, derivative)
        integer, intent(in) :: n !! integer Bessel order
        real(dp), intent(in) :: zeros(:) !! computed positive zeros
        logical, intent(in) :: derivative !! true for derivative roots

        integer :: i
        real(dp) :: residual

        do i = 1, size(zeros)
            if (derivative) then
                residual = jvp(real(n, dp), zeros(i))
            else
                residual = jv(real(n, dp), zeros(i))
            end if
            call require(abs(residual) <= 2.0e-12_dp, 'J residual')
            if (i > 1) call require(zeros(i) > zeros(i - 1), 'J ordering')
        end do
    end subroutine check_residuals_j

    subroutine check_residuals_y(n, zeros, derivative)
        integer, intent(in) :: n !! integer Bessel order
        real(dp), intent(in) :: zeros(:) !! computed positive zeros
        logical, intent(in) :: derivative !! true for derivative roots

        integer :: i
        real(dp) :: residual

        do i = 1, size(zeros)
            if (derivative) then
                residual = yvp(real(n, dp), zeros(i))
            else
                residual = yv(real(n, dp), zeros(i))
            end if
            call require(abs(residual) <= 2.0e-12_dp, 'Y residual')
            if (i > 1) call require(zeros(i) > zeros(i - 1), 'Y ordering')
        end do
    end subroutine check_residuals_y

    subroutine require(condition, message)
        logical, intent(in) :: condition !! assertion condition
        character(*), intent(in) :: message !! failure message

        if (.not. condition) then
            write(*, '(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine require

end program test_special_bessel_zeros
