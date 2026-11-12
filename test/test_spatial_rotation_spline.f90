! SPDX-License-Identifier: MIT
program test_spatial_rotation_spline
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_rotvec, rotation_from_quat, &
        spatial_rotation_success, spatial_rotation_invalid_input
    use scifort_spatial_rotation_spline, only : rotation_spline_type, rotation_spline, &
        spatial_rotation_no_convergence
    implicit none

    type(rotation_type) :: keys(5), samples(2), invalid, orientation, badkeys(12)
    type(rotation_spline_type) :: spline, failed
    type(rotation_type), allocatable :: path(:)
    real(dp), allocatable :: rates(:,:), accelerations(:,:)
    real(dp) :: t(5), pi, expected(3), delta, a(3), b(3), nan
    real(dp) :: t_bad(12), q_bad(4,12)
    integer :: i, status

    pi = acos(-1.0_dp)
    nan = ieee_value(0.0_dp,ieee_quiet_nan)

    ! Two keyframes are the piecewise geodesic special case: constant rate.
    samples(1) = rotation_from_rotvec([0.0_dp,0.0_dp,0.0_dp])
    samples(2) = rotation_from_rotvec([0.0_dp,0.0_dp,pi/2.0_dp])
    spline = rotation_spline([1.0_dp,3.0_dp], samples, status)
    if (status /= spatial_rotation_success .or. .not. spline%valid()) error stop 'two keys'
    orientation = spline%at(2.0_dp,status=status)
    if (status /= 0) error stop 'at query'
    call same_rotation(orientation,rotation_from_rotvec([0.0_dp,0.0_dp,pi/4.0_dp]),2.0e-14_dp)
    expected = [0.0_dp,0.0_dp,pi/4.0_dp]
    call close_vector(spline%angular_rate(2.0_dp),expected,2.0e-14_dp)
    call close_vector(spline%angular_acceleration(2.0_dp),[0.0_dp,0.0_dp,0.0_dp],2.0e-14_dp)
    call same_rotation(spline%at(-1.0_dp),rotation_from_rotvec([0.0_dp,0.0_dp,-pi/2.0_dp]),2.0e-14_dp)
    call same_rotation(spline%at(5.0_dp),rotation_from_rotvec([0.0_dp,0.0_dp,pi]),2.0e-14_dp)

    t = [0.0_dp,0.37_dp,1.1_dp,2.5_dp,3.3_dp]
    keys(1) = rotation_from_rotvec([0.0_dp,0.1_dp,0.2_dp])
    keys(2) = rotation_from_rotvec([0.4_dp,-0.7_dp,0.5_dp])
    keys(3) = rotation_from_rotvec([-0.2_dp,0.3_dp,0.9_dp])
    keys(4) = rotation_from_rotvec([0.2_dp,0.0_dp,-0.1_dp])
    keys(5) = rotation_from_rotvec([0.7_dp,0.5_dp,-0.3_dp])
    spline = rotation_spline(t,keys,status=status)
    if (status /= 0 .or. .not. spline%valid()) error stop 'five keys'
    do i=1,5
        call same_rotation(spline%at(t(i)),keys(i),3.0e-13_dp)
    end do
    ! SciPy's smooth attitude interpolant has continuous angular rate and
    ! acceleration at interior keyframes, even with nonuniform times.
    delta = 1.0e-7_dp
    do i=2,4
        a = spline%angular_rate(t(i)-delta)
        b = spline%angular_rate(t(i)+delta)
        if (maxval(abs(a-b)) > 5.0e-5_dp) error stop 'angular rate discontinuity'
        a = spline%angular_acceleration(t(i)-delta)
        b = spline%angular_acceleration(t(i)+delta)
        if (maxval(abs(a-b)) > 5.0e-4_dp) error stop 'angular acceleration discontinuity'
    end do
    call spline%evaluate_path([3.4_dp,0.37_dp,-0.2_dp,0.37_dp],path,rates,accelerations,status)
    if (status /= 0) error stop 'batched evaluations'
    if (size(path) /= 4 .or. any(shape(rates) /= [4,3]) .or. &
        any(shape(accelerations) /= [4,3])) error stop 'batch shapes'
    do i=1,4
        if (.not. path(i)%valid()) error stop 'batch rotation'
    end do
    call same_rotation(path(2),keys(2),3.0e-13_dp)
    call same_rotation(path(4),keys(2),3.0e-13_dp)
    call close_vector(rates(2,:),rates(4,:),0.0_dp)
    call close_vector(accelerations(2,:),accelerations(4,:),0.0_dp)
    call spline%evaluate_path([real(dp)::],path,status=status)
    if (status /= 0 .or. size(path) /= 0) error stop 'zero query path'

    ! Validation rejects inconsistent keyframes, invalid timestamps and
    ! invalid quaternions and returns invalid objects/NaNs, not crashes.
    failed = rotation_spline([0.0_dp],keys(1:1),status=status)
    if (status /= spatial_rotation_invalid_input .or. failed%valid()) error stop 'single key accepted'
    failed = rotation_spline(t(1:4),keys,status=status)
    if (status /= 1 .or. failed%valid()) error stop 'mismatched lengths accepted'
    failed = rotation_spline([0.0_dp,1.0_dp,1.0_dp,3.0_dp,4.0_dp],keys,status=status)
    if (status /= 1 .or. failed%valid()) error stop 'duplicate time accepted'
    failed = rotation_spline([0.0_dp,1.0_dp,nan,3.0_dp,4.0_dp],keys,status=status)
    if (status /= 1 .or. failed%valid()) error stop 'nan time accepted'
    invalid = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    keys(3) = invalid
    failed = rotation_spline(t,keys,status=status)
    if (status /= 1 .or. failed%valid()) error stop 'invalid key accepted'
    orientation = failed%at(1.0_dp,status=status)
    if (status /= 1 .or. orientation%valid()) error stop 'invalid spline evaluation accepted'
    a = failed%angular_rate(1.0_dp,status=status)
    if (status /= 1 .or. all(ieee_is_finite(a))) error stop 'invalid rate accepted'
    a = failed%angular_acceleration(1.0_dp,status=status)
    if (status /= 1 .or. all(ieee_is_finite(a))) error stop 'invalid acceleration accepted'
    call failed%evaluate_path([0.0_dp],path,rates,accelerations,status)
    if (status /= 1 .or. size(path) /= 0) error stop 'invalid batch accepted'
    orientation = spline%at(nan,status=status)
    if (status /= 1 .or. orientation%valid()) error stop 'nan query accepted'
    call spline%evaluate_path([0.0_dp,nan],path,rates,accelerations,status)
    if (status /= 1 .or. size(path) /= 0) error stop 'nan path accepted'

    ! A nearly antipodal, mixed-scale sequence makes SciPy 1.17's fixed
    ! ten-iteration nonlinear solve diverge. SciFort rejects interpolants
    ! whose cubic coefficients no longer reconstruct their own keyframes.
    t_bad = [0.00000000000000000e+00_dp, &
        1.84049197836532064e+00_dp, &
        2.02892384600011821e+00_dp, &
        5.30822596585675655e+00_dp, &
        8.48353898083114011e+00_dp, &
        1.17116595414112119e+01_dp, &
        1.19211711686024504e+01_dp, &
        1.47224892833094412e+01_dp, &
        1.70810182782083082e+01_dp, &
        1.85011072845271514e+01_dp, &
        2.02104464252284153e+01_dp, &
        2.14307496464621963e+01_dp]
    q_bad = reshape([-7.58742126493651914e-01_dp, &
        -1.37810239600344886e-01_dp, &
        -2.56043165597019118e-01_dp, &
        -5.82889887282532992e-01_dp, &
        4.29547230347342157e-01_dp, &
        -1.32842778254623206e-01_dp, &
        2.08783239868978038e-01_dp, &
        -8.68476558069550486e-01_dp, &
        7.89555021135290191e-01_dp, &
        -2.50424006665246168e-01_dp, &
        3.08905474563373583e-01_dp, &
        4.67405705218232603e-01_dp, &
        1.23933133551265853e-01_dp, &
        9.13868570916582246e-01_dp, &
        3.71729283797244181e-01_dp, &
        1.06311584818569324e-01_dp, &
        4.03887937053387269e-01_dp, &
        1.15582076396904543e-01_dp, &
        -8.76896116928863045e-03_dp, &
        -9.07435079351989193e-01_dp, &
        -1.06583784395323650e-01_dp, &
        -4.67636268642985686e-01_dp, &
        -8.51030295240550960e-01_dp, &
        -2.13784128822551489e-01_dp, &
        -4.43488693547958307e-02_dp, &
        -1.47981408658673763e-01_dp, &
        2.29274141978415280e-01_dp, &
        -9.61024478511555591e-01_dp, &
        8.16037579544758307e-02_dp, &
        9.03678647545582248e-01_dp, &
        3.67033571567794592e-01_dp, &
        -2.04919706226812631e-01_dp, &
        5.84797445380307090e-01_dp, &
        -1.38392061382638737e-01_dp, &
        -1.88869146750226036e-01_dp, &
        -7.76651807844913389e-01_dp, &
        -7.44850221129155887e-01_dp, &
        -3.64595869213521717e-01_dp, &
        1.67315835952144276e-02_dp, &
        -5.58558908573375645e-01_dp, &
        4.90006004716781129e-01_dp, &
        -7.70140521457281624e-01_dp, &
        1.96290517338343851e-01_dp, &
        -3.58116915760830101e-01_dp, &
        -6.29551127720916925e-01_dp, &
        -6.09316908658422363e-01_dp, &
        -4.39344636761165253e-01_dp, &
        -1.98430271272972553e-01_dp], [4,12])
    do i=1,12
        badkeys(i) = rotation_from_quat(q_bad(:,i))
    end do
    failed = rotation_spline(t_bad,badkeys,status=status)
    if (status /= spatial_rotation_no_convergence .or. failed%valid()) &
        error stop 'divergent attitude spline accepted'

    print '(a)', 'test_spatial_rotation_spline PASS'
contains
    subroutine close_vector(actual,expected,tol)
        real(dp), intent(in) :: actual(:) !! result vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! absolute error limit
        if (size(actual) /= size(expected)) error stop 'vector dimensions'
        if (any(abs(actual-expected) > tol)) then
            print *, 'actual = ',actual,' expected = ',expected
            error stop 'vector mismatch'
        end if
    end subroutine close_vector

    subroutine same_rotation(first,second,tol)
        type(rotation_type), intent(in) :: first !! computed rotation
        type(rotation_type), intent(in) :: second !! target rotation
        real(dp), intent(in) :: tol !! quaternion tolerance, modulo sign
        real(dp) :: p(4),q(4)
        if (.not. first%valid() .or. .not. second%valid()) error stop 'rotation invalid'
        p = first%as_quat()
        q = second%as_quat()
        if (min(maxval(abs(p-q)),maxval(abs(p+q))) > tol) then
            print *, 'actual = ',p,' expected = ',q
            error stop 'rotation mismatch'
        end if
    end subroutine same_rotation
end program test_spatial_rotation_spline
