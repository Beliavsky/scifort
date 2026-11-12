! SPDX-License-Identifier: MIT
! Regression of an incremental 3D hull that formerly left unmatched edges.
! All 26 sites are extreme despite multiple almost-coplanar hull facets.
program test_spatial_convex_hull_nearcoplanar
    use scifort_kinds, only : dp
    use scifort_spatial, only : convex_hull, convex_hull_type, spatial_distance_success
    implicit none
    integer, parameter :: n = 26
    real(dp), parameter :: cloud(n,3) = transpose(reshape([ &
        4.87740523931295300e-01_dp, -1.51620654549295547e-02_dp, 7.32606995285916551e-01_dp, &
        3.45570121096165250e-01_dp, -2.32558835393959440e-01_dp, 6.37545875070491874e-01_dp, &
        -3.32767515896732313e-01_dp, 6.65753440278324282e-01_dp, -3.69405143787377543e-01_dp, &
        -6.03068439452710603e-01_dp, -5.51623436453276494e-01_dp, 2.76370572201077684e-02_dp, &
        -1.64176367521642580e-01_dp, 7.24982909534344078e-01_dp, -4.37735257471461692e-01_dp, &
        2.06701995041655762e-01_dp, 9.33377564326547504e-01_dp, -5.54848641140887819e-02_dp, &
        2.82090872116577962e-01_dp, -6.35690720757650873e-01_dp, 2.55138237772534714e-02_dp, &
        -3.58183654789487549e-01_dp, -6.74877561576116114e-01_dp, -1.22564453670561652e-01_dp, &
        -2.12402954952178119e-01_dp, -6.20096270449605247e-01_dp, -2.33876971411467349e-01_dp, &
        7.07131318431322264e-01_dp, -3.82473212895967196e-01_dp, -5.94709626549475012e-01_dp, &
        6.98531841460448777e-01_dp, -4.19995594293962327e-01_dp, -4.79699658992087608e-01_dp, &
        4.48306164707624799e-01_dp, 2.63248846291697869e-01_dp, 7.13540695498002253e-01_dp, &
        -7.56843345060721262e-01_dp, 5.38898387038615811e-01_dp, -4.66076460162361816e-02_dp, &
        -7.29795685695283858e-01_dp, 3.27626957560259668e-01_dp, -3.39680888393652652e-02_dp, &
        -7.50919089337365953e-01_dp, 5.40029651914884168e-01_dp, -5.54873561026002077e-02_dp, &
        6.51093166391799150e-01_dp, 2.04506852072552474e-02_dp, -6.30645169168178032e-01_dp, &
        5.78024023236730522e-01_dp, 5.06377966111028255e-01_dp, -6.29643263761819671e-01_dp, &
        3.98033467004800912e-01_dp, 5.39864363810346193e-01_dp, -6.12499272326791355e-01_dp, &
        4.88610033238233821e-01_dp, 7.78254830242713469e-01_dp, -2.36898354871166283e-01_dp, &
        2.57334002786996696e-01_dp, 9.34061145059486342e-01_dp, -4.50983360580979459e-02_dp, &
        4.40615844151576785e-01_dp, 5.29973726679257462e-01_dp, 4.52127247842074809e-01_dp, &
        2.63740702890435885e-01_dp, 7.44007782305461940e-01_dp, 3.79307092812650015e-01_dp, &
        2.12890433030052117e-01_dp, 5.93582814833882422e-01_dp, 5.71012078024836778e-01_dp, &
        2.06340851548060034e-01_dp, 6.95273970627813598e-01_dp, 4.74248249464632288e-01_dp, &
        -4.68894894919322125e-01_dp, 6.84476908004677753e-01_dp, 1.39573655604611324e-01_dp, &
        -2.29279706588870230e-01_dp, 7.87017903660879092e-01_dp, 1.67654567134147531e-01_dp &
    ], [3,n]))
    type(convex_hull_type) :: hull
    real(dp) :: points(n,3), max_violation
    integer :: code, trial, i, f, j, order(n), stride

    do trial=1,12
        ! Strides coprime with 26 probe a range of insertion orders.
        stride = 2*trial-1
        if (stride==13) stride=25
        do i=1,n
            order(i) = 1+mod((i-1)*stride+trial-1,n)
        end do
        points = cloud(order,:)
        hull = convex_hull(points,status=code)
        if (code /= spatial_distance_success .or. .not. hull%valid()) &
            error stop 'near-coplanar 3d hull rejected'
        if (size(hull%vertices) /= n .or. size(hull%simplices,2) /= 2*n-4) &
            error stop 'near-coplanar extreme sites or facet count missing'
        if (any(hull%vertices /= [(i,i=1,n)])) &
            error stop 'near-coplanar hull vertex indices'
        if (abs(hull%area-6.4304977220236195_dp) > 1.0e-10_dp) &
            error stop 'near-coplanar hull surface area'
        if (abs(hull%volume-1.1081563594763306_dp) > 1.0e-10_dp) &
            error stop 'near-coplanar hull volume'
        if (any(hull%neighbors < 1) .or. any(hull%neighbors > 2*n-4)) &
            error stop 'near-coplanar hull adjacency'
        do f=1,size(hull%simplices,2)
            do j=1,n
                max_violation = dot_product(points(j,:),hull%equations(1:3,f)) + hull%equations(4,f)
                if (max_violation > 1.0e-11_dp) error stop 'hull not enclosing an input point'
            end do
        end do
    end do
    print '(a)', 'test_spatial_convex_hull_nearcoplanar PASS'
end program test_spatial_convex_hull_nearcoplanar
