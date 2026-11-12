! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_orthogonal_roots
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, orthogonal_roots_status_invalid_argument, &
        orthogonal_roots_status_success, roots_chebyc, roots_chebys, roots_chebyt, &
        roots_chebyu, roots_genlaguerre, roots_gegenbauer, roots_hermite, &
        roots_hermitenorm, roots_jacobi, roots_laguerre, roots_legendre, &
        roots_sh_chebyt, roots_sh_chebyu, roots_sh_jacobi, roots_sh_legendre
    implicit none

    real(dp), allocatable :: nodes(:)
    real(dp), allocatable :: weights(:)
    real(dp) :: mu
    integer :: status

    call roots_legendre(5, nodes, weights, mu, status)
    call check_rule('legendre', nodes, weights, &
        [-0.90617984593866396_dp, -0.53846931010568311_dp, 0.0_dp, &
         0.53846931010568311_dp, 0.90617984593866396_dp], &
        [0.23692688505618897_dp, 0.47862867049936653_dp, 0.568888888888889_dp, &
         0.47862867049936653_dp, 0.23692688505618897_dp], mu, 2.0_dp, status)

    call roots_jacobi(5, 0.2_dp, 1.3_dp, nodes, weights, mu, status)
    call check_rule('jacobi', nodes, weights, &
        [-0.77907538376244934_dp, -0.36983639008864427_dp, 0.12827808572192992_dp, &
         0.59205261310780322_dp, 0.90684194458657785_dp], &
        [0.05166385440376671_dp, 0.27657541746946007_dp, 0.57115386572219229_dp, &
         0.62135464977088983_dp, 0.30266465526680003_dp], &
        mu, 1.8234124426331089_dp, status)

    call roots_gegenbauer(5, 0.7_dp, nodes, weights, mu, status)
    call check_rule('gegenbauer', nodes, weights, &
        [-0.8895630676138504_dp, -0.52202235198454749_dp, 0.0_dp, &
         0.52202235198454749_dp, 0.8895630676138504_dp], &
        [0.1815898630680117_dp, 0.43922726238613574_dp, 0.54940949883057233_dp, &
         0.43922726238613574_dp, 0.1815898630680117_dp], &
        mu, 1.7910437497388672_dp, status)

    call roots_laguerre(5, nodes, weights, mu, status)
    call check_rule('laguerre', nodes, weights, &
        [0.26356031971814092_dp, 1.4134030591065168_dp, 3.5964257710407219_dp, &
         7.0858100058588374_dp, 12.640800844275784_dp], &
        [0.52175561058280873_dp, 0.3986668110831757_dp, 0.075942449681707616_dp, &
         0.0036117586799220545_dp, 2.3369972385776238e-05_dp], mu, 1.0_dp, status)

    call roots_genlaguerre(5, 1.2_dp, nodes, weights, mu, status)
    call check_rule('genlaguerre', nodes, weights, &
        [0.69560505674146744_dp, 2.2560123731556181_dp, 4.8130255402966196_dp, &
         8.6579134250376963_dp, 14.577443604768598_dp], &
        [0.35928648938860702_dp, 0.56399906875623596_dp, 0.1676128652107215_dp, &
         0.010815609526096196_dp, 8.8457998052324864e-05_dp], &
        mu, 1.101802490879713_dp, status)

    call roots_hermite(5, nodes, weights, mu, status)
    call check_rule('hermite', nodes, weights, &
        [-2.0201828704560851_dp, -0.95857246461381851_dp, 0.0_dp, &
         0.95857246461381851_dp, 2.0201828704560851_dp], &
        [0.019953242059045882_dp, 0.39361932315224107_dp, 0.94530872048294168_dp, &
         0.39361932315224107_dp, 0.019953242059045882_dp], &
        mu, 1.7724538509055159_dp, status)

    call roots_hermitenorm(5, nodes, weights, mu, status)
    call check_rule('hermitenorm', nodes, weights, &
        [-2.8569700138728056_dp, -1.3556261799742657_dp, 0.0_dp, &
         1.3556261799742657_dp, 2.8569700138728056_dp], &
        [0.028218145533215905_dp, 0.55666178521401732_dp, 1.3368684131365338_dp, &
         0.55666178521401732_dp, 0.028218145533215905_dp], &
        mu, 2.5066282746310002_dp, status)

    call roots_chebyt(5, nodes, weights, mu, status)
    call check_rule('chebyt', nodes, weights, &
        [-0.95105651629515353_dp, -0.58778525229247314_dp, 0.0_dp, &
         0.58778525229247314_dp, 0.95105651629515353_dp], &
        [0.62831853071795862_dp, 0.62831853071795862_dp, 0.62831853071795862_dp, &
         0.62831853071795862_dp, 0.62831853071795862_dp], &
        mu, 3.141592653589793_dp, status)

    call roots_chebyu(5, nodes, weights, mu, status)
    call check_rule('chebyu', nodes, weights, &
        [-0.86602540378443871_dp, -0.49999999999999978_dp, 0.0_dp, &
         0.50000000000000011_dp, 0.86602540378443871_dp], &
        [0.13089969389957468_dp, 0.3926990816987242_dp, 0.52359877559829882_dp, &
         0.39269908169872408_dp, 0.13089969389957468_dp], &
        mu, 1.5707963267948966_dp, status)

    call roots_chebyc(5, nodes, weights, mu, status)
    call check_rule('chebyc', nodes, weights, &
        [-1.9021130325903071_dp, -1.1755705045849463_dp, 0.0_dp, &
         1.1755705045849463_dp, 1.9021130325903071_dp], &
        [1.2566370614359172_dp, 1.2566370614359172_dp, 1.2566370614359172_dp, &
         1.2566370614359172_dp, 1.2566370614359172_dp], &
        mu, 6.283185307179586_dp, status)

    call roots_chebys(5, nodes, weights, mu, status)
    call check_rule('chebys', nodes, weights, &
        [-1.7320508075688774_dp, -0.99999999999999956_dp, 0.0_dp, &
         1.0000000000000002_dp, 1.7320508075688774_dp], &
        [0.26179938779914935_dp, 0.78539816339744839_dp, 1.0471975511965976_dp, &
         0.78539816339744817_dp, 0.26179938779914935_dp], &
        mu, 3.141592653589793_dp, status)

    call roots_sh_legendre(5, nodes, weights, mu, status)
    call check_rule('shifted legendre', nodes, weights, &
        [0.046910077030668018_dp, 0.23076534494715845_dp, 0.5_dp, &
         0.7692346550528415_dp, 0.95308992296933193_dp], &
        [0.11846344252809449_dp, 0.23931433524968326_dp, 0.2844444444444445_dp, &
         0.23931433524968326_dp, 0.11846344252809449_dp], mu, 1.0_dp, status)

    call roots_sh_jacobi(5, 2.5_dp, 1.2_dp, nodes, weights, mu, status)
    call check_rule('shifted jacobi', nodes, weights, &
        [0.046579027706711129_dp, 0.20397369344609834_dp, 0.43586095713903505_dp, &
         0.68491819504432216_dp, 0.88953769188122467_dp], &
        [0.053504057541160555_dp, 0.10984102159369707_dp, 0.10096669288826836_dp, &
         0.048892088300538969_dp, 0.0091329654477844847_dp], &
        mu, 0.32233682577144945_dp, status)

    call roots_sh_chebyt(5, nodes, weights, mu, status)
    call check_rule('shifted chebyt', nodes, weights, &
        [0.024471741852423234_dp, 0.20610737385376343_dp, 0.5_dp, &
         0.79389262614623657_dp, 0.97552825814757682_dp], &
        [0.62831853071795862_dp, 0.62831853071795862_dp, 0.62831853071795862_dp, &
         0.62831853071795862_dp, 0.62831853071795862_dp], &
        mu, 3.141592653589793_dp, status)

    call roots_sh_chebyu(5, nodes, weights, mu, status)
    call check_rule('shifted chebyu', nodes, weights, &
        [0.066987298107780646_dp, 0.25000000000000011_dp, 0.5_dp, &
         0.75_dp, 0.93301270189221941_dp], &
        [0.032724923474893662_dp, 0.098174770424681035_dp, 0.13089969389957468_dp, &
         0.098174770424681007_dp, 0.032724923474893662_dp], &
        mu, 0.3926990816987241_dp, status)

    call roots_jacobi(32, -0.8_dp, 0.2_dp, nodes, weights, mu, status)
    call check_true('moderate-order status', status == orthogonal_roots_status_success)
    call check_true('moderate-order monotone', all(nodes(2:) > nodes(:size(nodes) - 1)))
    call check_true('moderate-order weights', all(weights > 0.0_dp))
    call check_close('moderate-order sum', sum(weights), mu, 2.0e-13_dp)

    call roots_legendre(0, nodes, weights, mu, status)
    call check_true('invalid n status', status == orthogonal_roots_status_invalid_argument)
    call check_true('invalid n arrays', size(nodes) == 0 .and. size(weights) == 0)
    call check_true('invalid n mu', ieee_is_nan(mu))

    call roots_jacobi(3, -1.0_dp, 0.0_dp, nodes, weights, mu, status)
    call check_true('invalid jacobi', status == orthogonal_roots_status_invalid_argument)
    call roots_gegenbauer(3, -0.5_dp, nodes, weights, mu, status)
    call check_true('invalid gegenbauer', status == orthogonal_roots_status_invalid_argument)
    call roots_genlaguerre(3, -1.0_dp, nodes, weights, mu, status)
    call check_true('invalid genlaguerre', status == orthogonal_roots_status_invalid_argument)
    call roots_sh_jacobi(3, 1.0_dp, 0.0_dp, nodes, weights, mu, status)
    call check_true('invalid shifted jacobi', status == orthogonal_roots_status_invalid_argument)

    print *, 'special orthogonal-root quadrature tests passed'

contains

    subroutine check_rule(name, actual_nodes, actual_weights, expected_nodes, &
            expected_weights, actual_mu, expected_mu, actual_status)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual_nodes(:) !! computed quadrature nodes
        real(dp), intent(in) :: actual_weights(:) !! computed quadrature weights
        real(dp), intent(in) :: expected_nodes(:) !! SciPy reference nodes
        real(dp), intent(in) :: expected_weights(:) !! SciPy reference weights
        real(dp), intent(in) :: actual_mu !! computed zeroth moment
        real(dp), intent(in) :: expected_mu !! SciPy reference zeroth moment
        integer, intent(in) :: actual_status !! expected success status

        integer :: i

        call check_true(trim(name)//' status', actual_status == orthogonal_roots_status_success)
        call check_true(trim(name)//' node size', size(actual_nodes) == size(expected_nodes))
        call check_true(trim(name)//' weight size', size(actual_weights) == size(expected_weights))
        do i = 1, size(expected_nodes)
            call check_close(trim(name)//' node', actual_nodes(i), expected_nodes(i), 2.0e-12_dp)
            call check_close(trim(name)//' weight', actual_weights(i), expected_weights(i), 2.0e-12_dp)
        end do
        call check_close(trim(name)//' mu', actual_mu, expected_mu, 5.0e-14_dp)
        call check_close(trim(name)//' sum', sum(actual_weights), actual_mu, 5.0e-13_dp)
    end subroutine check_rule

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_orthogonal_roots
