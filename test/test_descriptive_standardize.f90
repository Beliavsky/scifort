program test_descriptive_standardize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : gstd, gstd_all, gstd_keepdims, gzscore, gzscore_all, &
        zmap, zmap_all, zscore, zscore_all
    implicit none

    integer :: i
    real(dp), parameter :: tol = 8.0e-13_dp
    real(dp) :: x1(3)
    real(dp) :: x2(3, 1)
    real(dp) :: x3(3, 1, 1)
    real(dp) :: x4(3, 1, 1, 1)
    real(dp) :: x5(3, 1, 1, 1, 1)
    real(dp) :: x6(3, 1, 1, 1, 1, 1)
    real(dp) :: x7(3, 1, 1, 1, 1, 1, 1)
    real(dp) :: x8(3, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x9(3, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x10(3, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x11(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x12(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x13(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x14(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x15(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: empty15(0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: matrix(3, 2)
    real(dp) :: matrix_nan(3, 2)
    real(dp) :: cube(2, 3, 4)
    real(dp) :: cube_nan(2, 3, 4)
    real(dp) :: zmap_compare(4, 2)
    real(dp) :: zmap_compare_nan(4, 2)
    real(dp) :: zmap_scores(3, 1)
    real(dp) :: zmap_vector_compare(5)
    real(dp) :: zmap_vector_scores(4)
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: values2(:, :)
    real(dp), allocatable :: values3(:, :, :)

    x1 = [1.0_dp, 2.0_dp, 4.0_dp]
    x2 = reshape(x1, shape(x2))
    x3 = reshape(x1, shape(x3))
    x4 = reshape(x1, shape(x4))
    x5 = reshape(x1, shape(x5))
    x6 = reshape(x1, shape(x6))
    x7 = reshape(x1, shape(x7))
    x8 = reshape(x1, shape(x8))
    x9 = reshape(x1, shape(x9))
    x10 = reshape(x1, shape(x10))
    x11 = reshape(x1, shape(x11))
    x12 = reshape(x1, shape(x12))
    x13 = reshape(x1, shape(x13))
    x14 = reshape(x1, shape(x14))
    x15 = reshape(x1, shape(x15))
    matrix = reshape([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp], shape(matrix))
    matrix_nan = matrix
    cube = reshape([(real(i, dp), i=1,24)], shape(cube))
    cube_nan = cube
    cube_nan(2, 2, 3) = ieee_value(0.0_dp, ieee_quiet_nan)
    matrix_nan(2, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    zmap_vector_scores = [0.5_dp, 2.0_dp, 2.5_dp, 3.0_dp]
    zmap_vector_compare = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    zmap_scores(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
    zmap_compare = reshape([0.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, &
        10.0_dp, 12.0_dp, 14.0_dp, 16.0_dp], shape(zmap_compare))
    zmap_compare_nan = zmap_compare
    zmap_compare_nan(3, 1) = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_array(zmap(zmap_vector_scores, zmap_vector_compare), &
        [-1.0606601717798212_dp, 0.0_dp, 0.35355339059327373_dp, &
        0.7071067811865475_dp], 'zmap vector')
    call check_array(zmap(zmap_vector_scores, zmap_vector_compare, ddof=1), &
        [-0.9486832980505138_dp, 0.0_dp, 0.31622776601683794_dp, &
        0.6324555320336759_dp], 'zmap vector ddof 1')
    values2 = zmap(zmap_scores, zmap_compare, axis=1)
    call check_true(all(shape(values2) == [3, 2]), 'zmap broadcast shape')
    call check_array(reshape(values2, [6]), [-0.8944271909999159_dp, &
        -0.4472135954999579_dp, 0.0_dp, -5.366563145999495_dp, &
        -4.919349550499537_dp, -4.47213595499958_dp], 'zmap broadcast values')
    call check_close(sum(zmap(x15, x15)), 0.0_dp, 'zmap rank 15')

    call check_array(reshape(zscore_all(matrix), [6]), &
        [-0.8773384834252593_dp, -0.7849870641173373_dp, &
        -0.6002842255014932_dp, -0.23087854826980508_dp, &
        0.5079328061935712_dp, 1.9855555151203237_dp], 'zscore all')
    call check_array(reshape(gzscore_all(matrix), [6]), &
        [-1.4638501094227998_dp, -0.87831006565368_dp, &
        -0.29277002188455997_dp, 0.2927700218845598_dp, &
        0.8783100656536798_dp, 1.4638501094227998_dp], 'gzscore all')
    call check_array(reshape(zmap_all(zmap_scores, zmap_compare), [3]), &
        [-1.2780193008453875_dp, -1.0954451150103321_dp, &
        -0.9128709291752769_dp], 'zmap all')
    call check_close(gstd_all(matrix), 3.6574260469011497_dp, 'gstd all')
    call check_close(gstd_all(matrix, ddof=0), 3.2666799787475402_dp, &
        'gstd all ddof 0')
    call check_close(sum(zscore_all(x15)), 0.0_dp, 'zscore all rank 15')
    call check_close(sum(gzscore_all(x15)), 0.0_dp, 'gzscore all rank 15')
    call check_close(gstd_all(x15), 2.0_dp, 'gstd all rank 15')

    call check_array(zscore(x1), [-1.0690449676496976_dp, -0.2672612419124244_dp, &
        1.3363062095621219_dp], 'zscore rank 1')
    call check_array(gzscore(x1), [-1.224744871391589_dp, 0.0_dp, &
        1.224744871391589_dp], 'gzscore rank 1')
    call check_close(gstd(x1), 2.0_dp, 'gstd rank 1')
    call check_close(sum(gstd_keepdims(x1)), 2.0_dp, 'gstd rank 1 keepdims')

    call check_close(sum(zscore(x2)), 0.0_dp, 'zscore rank 2')
    call check_close(sum(gzscore(x2)), 0.0_dp, 'gzscore rank 2')
    call check_close(sum(gstd(x2)), 2.0_dp, 'gstd rank 2')
    call check_close(sum(gstd_keepdims(x2)), 2.0_dp, 'gstd rank 2 keepdims')

    call check_close(sum(zscore(x3)), 0.0_dp, 'zscore rank 3')
    call check_close(sum(gzscore(x3)), 0.0_dp, 'gzscore rank 3')
    call check_close(sum(gstd(x3)), 2.0_dp, 'gstd rank 3')
    call check_close(sum(gstd_keepdims(x3)), 2.0_dp, 'gstd rank 3 keepdims')

    call check_close(sum(zscore(x4)), 0.0_dp, 'zscore rank 4')
    call check_close(sum(gzscore(x4)), 0.0_dp, 'gzscore rank 4')
    call check_close(sum(gstd(x4)), 2.0_dp, 'gstd rank 4')
    call check_close(sum(gstd_keepdims(x4)), 2.0_dp, 'gstd rank 4 keepdims')

    call check_close(sum(zscore(x5)), 0.0_dp, 'zscore rank 5')
    call check_close(sum(gzscore(x5)), 0.0_dp, 'gzscore rank 5')
    call check_close(sum(gstd(x5)), 2.0_dp, 'gstd rank 5')
    call check_close(sum(gstd_keepdims(x5)), 2.0_dp, 'gstd rank 5 keepdims')

    call check_close(sum(zscore(x6)), 0.0_dp, 'zscore rank 6')
    call check_close(sum(gzscore(x6)), 0.0_dp, 'gzscore rank 6')
    call check_close(sum(gstd(x6)), 2.0_dp, 'gstd rank 6')
    call check_close(sum(gstd_keepdims(x6)), 2.0_dp, 'gstd rank 6 keepdims')

    call check_close(sum(zscore(x7)), 0.0_dp, 'zscore rank 7')
    call check_close(sum(gzscore(x7)), 0.0_dp, 'gzscore rank 7')
    call check_close(sum(gstd(x7)), 2.0_dp, 'gstd rank 7')
    call check_close(sum(gstd_keepdims(x7)), 2.0_dp, 'gstd rank 7 keepdims')

    call check_close(sum(zscore(x8)), 0.0_dp, 'zscore rank 8')
    call check_close(sum(gzscore(x8)), 0.0_dp, 'gzscore rank 8')
    call check_close(sum(gstd(x8)), 2.0_dp, 'gstd rank 8')
    call check_close(sum(gstd_keepdims(x8)), 2.0_dp, 'gstd rank 8 keepdims')

    call check_close(sum(zscore(x9)), 0.0_dp, 'zscore rank 9')
    call check_close(sum(gzscore(x9)), 0.0_dp, 'gzscore rank 9')
    call check_close(sum(gstd(x9)), 2.0_dp, 'gstd rank 9')
    call check_close(sum(gstd_keepdims(x9)), 2.0_dp, 'gstd rank 9 keepdims')

    call check_close(sum(zscore(x10)), 0.0_dp, 'zscore rank 10')
    call check_close(sum(gzscore(x10)), 0.0_dp, 'gzscore rank 10')
    call check_close(sum(gstd(x10)), 2.0_dp, 'gstd rank 10')
    call check_close(sum(gstd_keepdims(x10)), 2.0_dp, 'gstd rank 10 keepdims')

    call check_close(sum(zscore(x11)), 0.0_dp, 'zscore rank 11')
    call check_close(sum(gzscore(x11)), 0.0_dp, 'gzscore rank 11')
    call check_close(sum(gstd(x11)), 2.0_dp, 'gstd rank 11')
    call check_close(sum(gstd_keepdims(x11)), 2.0_dp, 'gstd rank 11 keepdims')

    call check_close(sum(zscore(x12)), 0.0_dp, 'zscore rank 12')
    call check_close(sum(gzscore(x12)), 0.0_dp, 'gzscore rank 12')
    call check_close(sum(gstd(x12)), 2.0_dp, 'gstd rank 12')
    call check_close(sum(gstd_keepdims(x12)), 2.0_dp, 'gstd rank 12 keepdims')

    call check_close(sum(zscore(x13)), 0.0_dp, 'zscore rank 13')
    call check_close(sum(gzscore(x13)), 0.0_dp, 'gzscore rank 13')
    call check_close(sum(gstd(x13)), 2.0_dp, 'gstd rank 13')
    call check_close(sum(gstd_keepdims(x13)), 2.0_dp, 'gstd rank 13 keepdims')

    call check_close(sum(zscore(x14)), 0.0_dp, 'zscore rank 14')
    call check_close(sum(gzscore(x14)), 0.0_dp, 'gzscore rank 14')
    call check_close(sum(gstd(x14)), 2.0_dp, 'gstd rank 14')
    call check_close(sum(gstd_keepdims(x14)), 2.0_dp, 'gstd rank 14 keepdims')

    call check_close(sum(zscore(x15)), 0.0_dp, 'zscore rank 15')
    call check_close(sum(gzscore(x15)), 0.0_dp, 'gzscore rank 15')
    call check_close(sum(gstd(x15)), 2.0_dp, 'gstd rank 15')
    call check_close(sum(gstd_keepdims(x15)), 2.0_dp, 'gstd rank 15 keepdims')

    values2 = zscore(matrix, axis=1)
    call check_array(reshape(values2, [6]), [-1.0690449676496976_dp, &
        -0.2672612419124244_dp, 1.3363062095621219_dp, &
        -1.0690449676496976_dp, -0.2672612419124244_dp, &
        1.3363062095621219_dp], 'matrix zscore axis 1')
    values2 = zscore(matrix, axis=2)
    call check_array(reshape(values2, [6]), [-1.0_dp, -1.0_dp, -1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], 'matrix zscore axis 2')
    values2 = zscore(matrix, axis=1, ddof=1)
    call check_array(reshape(values2, [6]), [-0.8728715609439694_dp, &
        -0.21821789023599236_dp, 1.0910894511799618_dp, &
        -0.8728715609439694_dp, -0.21821789023599236_dp, &
        1.0910894511799618_dp], 'matrix zscore ddof 1')
    values2 = gzscore(matrix, axis=1)
    call check_array(reshape(values2, [6]), [-1.224744871391589_dp, 0.0_dp, &
        1.224744871391589_dp, -1.224744871391589_dp, 0.0_dp, &
        1.224744871391589_dp], 'matrix gzscore axis 1')
    values = gstd(matrix, axis=1)
    call check_array(values, [2.0_dp, 2.0_dp], 'matrix gstd axis 1')
    values = gstd(matrix, axis=2)
    call check_array(values, [4.350919557124710_dp, 4.350919557124710_dp, &
        4.350919557124710_dp], 'matrix gstd axis 2')
    values = gstd(matrix, axis=1, ddof=0)
    call check_array(values, [1.7611241081939726_dp, 1.7611241081939726_dp], &
        'matrix gstd ddof 0')
    values2 = gstd_keepdims(matrix, axis=2)
    call check_array(reshape(values2, [3]), [4.350919557124710_dp, &
        4.350919557124710_dp, 4.350919557124710_dp], 'matrix gstd keepdims')

    values3 = gstd_keepdims(cube, [1, 3])
    call check_true(all(shape(values3) == [1, 3, 1]), 'gstd tuple shape')
    call check_array(reshape(values3, [3]), [2.9760387117251224_dp, &
        2.100065404942243_dp, 1.8002289874255522_dp], 'gstd tuple values')
    values3 = gstd_keepdims(cube, [3, 1])
    call check_array(reshape(values3, [3]), [2.9760387117251224_dp, &
        2.100065404942243_dp, 1.8002289874255522_dp], 'gstd tuple order')
    values3 = gstd_keepdims(cube_nan, [1, 3], nan_policy='omit')
    call check_array(reshape(values3, [3]), [2.9760387117251224_dp, &
        2.175499629928427_dp, 1.8002289874255522_dp], 'gstd tuple omit')
    values3 = gstd_keepdims(cube_nan, [1, 3], nan_policy='propagate')
    call check_close(values3(1, 1, 1), 2.9760387117251224_dp, &
        'gstd tuple propagate finite')
    call check_nan(values3(1, 2, 1), 'gstd tuple propagate nan')
    values3 = gstd_keepdims(cube, [1, 2, 3])
    call check_close(values3(1, 1, 1), 2.2944076136018947_dp, &
        'gstd tuple all axes')
    call check_true(size(gstd_keepdims(cube, [1, 1])) == 0, &
        'gstd tuple duplicate axis')
    call check_true(size(gstd_keepdims(cube, [4])) == 0, &
        'gstd tuple invalid axis')
    call check_close(sum(gstd_keepdims(x15, [1, 15])), 2.0_dp, &
        'gstd tuple rank 15')

    values2 = zscore(matrix_nan, axis=1, nan_policy='omit')
    call check_nan(values2(2, 1), 'zscore omit preserves nan')
    call check_close(values2(1, 1), -1.0_dp, 'zscore omit first finite')
    call check_close(values2(3, 1), 1.0_dp, 'zscore omit second finite')
    call check_close(values2(1, 2), -1.0690449676496976_dp, &
        'zscore omit unaffected slice')
    values2 = zscore(matrix_nan, axis=1, nan_policy='propagate')
    call check_true(all([(ieee_is_nan(values2(i, 1)), i=1,3)]), &
        'zscore propagate nan slice')
    call check_close(values2(1, 2), -1.0690449676496976_dp, &
        'zscore propagate unaffected slice')
    values2 = gzscore(matrix_nan, axis=1, nan_policy='omit')
    call check_nan(values2(2, 1), 'gzscore omit preserves nan')
    call check_close(values2(1, 1), -1.0_dp, 'gzscore omit first finite')
    call check_close(values2(3, 1), 1.0_dp, 'gzscore omit second finite')
    values = gstd(matrix_nan, axis=1, nan_policy='omit')
    call check_array(values, [2.665144142690225_dp, 2.0_dp], 'gstd omit')
    values = gstd(matrix_nan, axis=1, nan_policy='propagate')
    call check_nan(values(1), 'gstd propagate nan slice')
    call check_close(values(2), 2.0_dp, 'gstd propagate finite slice')

    values2 = zscore_all(matrix_nan, nan_policy='propagate')
    call check_true(all([(ieee_is_nan(values2(i, 1)), i=1,3)]) .and. &
        all([(ieee_is_nan(values2(i, 2)), i=1,3)]), 'zscore all propagate')
    values2 = zscore_all(matrix_nan, nan_policy='omit')
    call check_array(reshape(values2, [6]), [-1.0083953684918676_dp, &
        ieee_value(0.0_dp, ieee_quiet_nan), -0.7382894662172602_dp, &
        -0.37814826318445033_dp, 0.3421341428811695_dp, &
        1.782698955012409_dp], 'zscore all omit', allow_nan=.true.)
    values2 = gzscore_all(matrix_nan, nan_policy='omit')
    call check_array(reshape(values2, [6]), [-1.6274669424134698_dp, &
        ieee_value(0.0_dp, ieee_quiet_nan), -0.46499055497527697_dp, &
        0.1162476387438193_dp, 0.6974858324629157_dp, &
        1.2787240261820123_dp], 'gzscore all omit', allow_nan=.true.)
    call check_nan(gstd_all(matrix_nan, nan_policy='propagate'), &
        'gstd all propagate')
    call check_close(gstd_all(matrix_nan, nan_policy='omit'), &
        3.7935233192378846_dp, 'gstd all omit')
    values2 = zmap_all(zmap_scores, zmap_compare_nan, nan_policy='propagate')
    call check_true(all([(ieee_is_nan(values2(i, 1)), i=1,3)]), &
        'zmap all propagate')
    values2 = zmap_all(zmap_scores, zmap_compare_nan, nan_policy='omit')
    call check_array(reshape(values2, [3]), [-1.3453336688020452_dp, &
        -1.1676480899036619_dp, -0.9899625110052785_dp], 'zmap all omit')

    call check_array(gzscore([1.0_dp, -2.0_dp, 4.0_dp, 8.0_dp], &
        nan_policy='omit'), [-1.3363062095621219_dp, &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.2672612419124244_dp, &
        1.0690449676496976_dp], 'gzscore negative omit', allow_nan=.true.)
    call check_array(zscore(x1, ddof=-1), [-1.2344267996967353_dp, &
        -0.3086066999241838_dp, 1.543033499620919_dp], 'zscore negative ddof')
    call check_close(gstd(x1, ddof=-1), 1.6325269194381529_dp, &
        'gstd negative ddof')

    values = zscore([2.0_dp, 2.0_dp, 2.0_dp])
    call check_true(all([(ieee_is_nan(values(i)), i=1,3)]), 'zscore constant')
    values = zscore([1.0_dp, ieee_value(0.0_dp, ieee_positive_inf), 3.0_dp])
    call check_true(all([(ieee_is_nan(values(i)), i=1,3)]), 'zscore infinity')
    values = gzscore([1.0_dp, 0.0_dp, 4.0_dp])
    call check_true(all([(ieee_is_nan(values(i)), i=1,3)]), 'gzscore zero')
    call check_close(gstd([2.0_dp, 2.0_dp, 2.0_dp]), 1.0_dp, 'gstd constant')
    call check_nan(gstd([1.0_dp, 0.0_dp, 4.0_dp]), 'gstd zero outside domain')
    values = zmap([0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 2.0_dp], &
        [0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 2.0_dp, 4.0_dp], &
        nan_policy='omit')
    call check_array(values, [-1.224744871391589_dp, &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp], 'zmap omit', allow_nan=.true.)
    values = zmap([0.0_dp, 1.0_dp, 2.0_dp], &
        [0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 2.0_dp, 4.0_dp], &
        nan_policy='propagate')
    call check_true(all([(ieee_is_nan(values(i)), i=1,3)]), 'zmap propagate')
    values = zmap([4.0_dp, 5.0_dp, 6.0_dp], [5.0_dp, 5.0_dp, 5.0_dp])
    call check_true(.not. ieee_is_finite(values(1)) .and. values(1) < 0.0_dp, &
        'zmap constant lower')
    call check_nan(values(2), 'zmap constant equal')
    call check_true(.not. ieee_is_finite(values(3)) .and. values(3) > 0.0_dp, &
        'zmap constant upper')

    values = zscore(x1, nan_policy='bad')
    call check_true(all([(ieee_is_nan(values(i)), i=1,3)]), 'zscore invalid policy')
    values = zmap(zmap_vector_scores, zmap_vector_compare, nan_policy='bad')
    call check_true(all([(ieee_is_nan(values(i)), i=1,4)]), 'zmap invalid policy')
    call check_nan(gstd(x1, nan_policy='bad'), 'gstd invalid policy')

    call check_true(size(zscore(matrix, axis=3)) == 0, 'zscore invalid axis')
    call check_true(size(zmap(matrix, zmap_compare, axis=3)) == 0, 'zmap invalid axis')
    call check_true(size(zmap(matrix, reshape([(real(i, dp), i=1,12)], [4, 3]), &
        axis=1)) == 0, 'zmap invalid broadcast')
    call check_true(size(gstd(matrix, axis=3)) == 0, 'gstd invalid axis')
    call check_true(size(zscore(empty15, axis=1)) == 0, 'zscore rank15 empty')
    values = reshape(gstd(empty15, axis=1), [2])
    call check_true(all([(ieee_is_nan(values(i)), i=1,2)]), 'gstd rank15 empty')

    print *, 'descriptive standardization tests passed'

contains

    subroutine check_close(actual, expected, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! SciPy reference value
        character(len=*), intent(in) :: label !! test label

        if (ieee_is_nan(actual) .or. abs(actual - expected) > tol) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_nan(actual, label)
        real(dp), intent(in) :: actual !! computed value expected to be NaN
        character(len=*), intent(in) :: label !! test label

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', label
            print *, ' actual = ', actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_true(condition, label)
        logical, intent(in) :: condition !! condition expected to be true
        character(len=*), intent(in) :: label !! test label

        if (.not. condition) then
            print *, 'FAIL ', label
            error stop 1
        end if
    end subroutine check_true

    subroutine check_array(actual, expected, label, allow_nan)
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! SciPy reference values
        character(len=*), intent(in) :: label !! test label
        logical, intent(in), optional :: allow_nan !! compare NaNs positionally when true

        logical :: compare_nan
        integer :: i

        compare_nan = .false.
        if (present(allow_nan)) compare_nan = allow_nan
        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', label, ' size'
            error stop 1
        end if
        do i = 1, size(actual)
            if (compare_nan .and. ieee_is_nan(expected(i))) then
                if (.not. ieee_is_nan(actual(i))) then
                    print *, 'FAIL ', label, ' at ', i
                    error stop 1
                end if
            else if (ieee_is_nan(actual(i)) .or. abs(actual(i) - expected(i)) > tol) then
                print *, 'FAIL ', label, ' at ', i
                print *, ' actual  = ', actual(i)
                print *, ' expected= ', expected(i)
                error stop 1
            end if
        end do
    end subroutine check_array

end program test_descriptive_standardize
