! SPDX-License-Identifier: MIT
! Independent deterministic checks for irregular-sampling spectral analysis.
program test_signal_irregular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : lombscargle, lombscargle_amplitude, &
        vectorstrength, vector_strength_result, irregular_status_success, &
        irregular_status_invalid_input, irregular_status_singular
    implicit none

    real(dp), parameter :: x(10) = [0.0_dp,0.4_dp,1.1_dp,1.8_dp,3.0_dp, &
        4.1_dp,4.8_dp,6.3_dp,8.2_dp,9.7_dp]
    real(dp), parameter :: y(10) = [1.5_dp,1.5474837745103478_dp,0.42407453482975876_dp, &
        -0.8552567229791213_dp,0.29066595009171003_dp,1.5449622493622244_dp, &
        0.41525306323449723_dp,-0.5178002465787668_dp, &
        1.0253830840676308_dp,-0.9037385572276493_dp]
    real(dp), parameter :: wt(10) = [1.0_dp,2.0_dp,0.25_dp,1.0_dp,0.5_dp, &
        2.0_dp,1.0_dp,4.0_dp,0.7_dp,0.9_dp]
    real(dp), parameter :: omega(5) = [0.6_dp,1.1_dp,1.7_dp,2.25_dp,3.6_dp]
    real(dp), parameter :: power0(5) = [0.93086927992831248_dp,1.1563803435929578_dp, &
        5.316057755996594_dp,0.38147454483615706_dp,1.6805934051221199_dp]
    real(dp), parameter :: power1(5) = [1.0135370552913028_dp,0.82290061129188263_dp, &
        4.9514207923803752_dp,0.5549391200243553_dp,2.1321582864246_dp]
    real(dp), parameter :: normalized1(5) = [0.20469620696568785_dp, &
        0.16619484503482829_dp,1.0_dp,0.11207674388699469_dp,0.43061544874265778_dp]
    complex(dp), parameter :: amp1(5) = [ &
        cmplx(0.26448256670050319_dp,0.61749373540760588_dp,dp), &
        cmplx(0.30491375985448921_dp,-0.62736531902150983_dp,dp), &
        cmplx(1.2_dp,0.5_dp,dp), &
        cmplx(0.4796912068784196_dp,-0.22859668385063148_dp,dp), &
        cmplx(0.21371018182344881_dp,0.82001413910906262_dp,dp)]
    real(dp), allocatable :: pg(:)
    complex(dp), allocatable :: ag(:)
    type(vector_strength_result) :: vs
    integer :: st

    pg = lombscargle(x, y, omega, weights=wt, status=st)
    if (st /= irregular_status_success) error stop 'default Lomb status'
    call check(pg,power0,2.0e-12_dp)

    pg = lombscargle(x,y,omega,weights=wt,floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'floating mean Lomb status'
    call check(pg,power1,2.0e-12_dp)

    pg = lombscargle(x,y,omega,normalize='normalize', &
        weights=wt,floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'normalized Lomb status'
    call check(pg,normalized1,2.0e-12_dp)

    ag = lombscargle_amplitude(x,y,omega,weights=wt,floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'Lomb amplitude status'
    if (maxval(abs(ag-amp1)) > 2.0e-12_dp) error stop 'Lomb amplitude values'

    ! Shifting the baseline does not change a floating-mean periodogram.
    call check(lombscargle(x,y+123.0_dp,omega,weights=wt,floating_mean=.true.), &
        power1,5.0e-12_dp)
    ! Rescaling all observation weights leaves all three outputs unchanged.
    call check(lombscargle(x,y,omega,weights=8.0_dp*wt,floating_mean=.true.), &
        power1,2.0e-12_dp)
    ! Precentering subtracts the arithmetic mean (not the weighted mean).
    call check(lombscargle(x,y,omega,precenter=.true.,weights=wt), &
        lombscargle(x,y-sum(y)/real(size(y),dp),omega,weights=wt),2.0e-13_dp)

    vs = vectorstrength(x,[1.0_dp,2.0_dp,3.0_dp])
    if (vs%status /= irregular_status_success) error stop 'vector strength status'
    call check(vs%strength,[0.32221058348927883_dp,0.27620020927934968_dp, &
        0.17225974868605831_dp],1.0e-12_dp)
    call check(vs%phase,[0.25485177476248072_dp,0.59094644480749881_dp, &
        1.2630149058025806_dp],1.0e-12_dp)
    vs = vectorstrength(x, 2.0_dp)
    if (vs%status /= irregular_status_success) error stop 'scalar vector strength status'
    if (size(vs%strength) /= 1) error stop 'scalar vector strength shape'
    call check(vs%strength,[0.27620020927934968_dp],1.0e-12_dp)

    ! One constant observation has zero model variance for floating mean.
    pg = lombscargle(x,spread(3.0_dp,1,size(x)),omega, &
        floating_mean=.true., normalize='normalize',status=st)
    if (st /= irregular_status_singular) error stop 'zero-variance Lomb status'
    if (.not. all(ieee_is_nan(pg))) error stop 'zero-variance Lomb output'

    pg = lombscargle(x,y,omega,weights=[real(dp) ::],status=st)
    call reject(pg,st)
    pg = lombscargle(x,y,omega,weights=-wt,status=st)
    call reject(pg,st)
    pg = lombscargle(x,y,omega,weights=0.0_dp*wt,status=st)
    call reject(pg,st)
    pg = lombscargle(x,y,omega,normalize='invalid',status=st)
    call reject(pg,st)
    pg = lombscargle(x,y(:5),omega,status=st)
    call reject(pg,st)
    pg = lombscargle(x,y,[real(dp) ::],status=st)
    if (st /= irregular_status_invalid_input .or. size(pg) /= 0) error stop 'empty frequencies'
    ag = lombscargle_amplitude(x,y(:5),omega,status=st)
    if (st /= irregular_status_invalid_input .or. &
        .not. all(ieee_is_nan(real(ag,dp)))) error stop 'invalid amplitude'
    vs = vectorstrength(x,[-1.0_dp,2.0_dp])
    if (vs%status /= irregular_status_invalid_input) error stop 'invalid vector strength status'
    if (.not. all(ieee_is_nan(vs%strength))) error stop 'invalid vector strength values'
    vs = vectorstrength([real(dp) ::], 1.0_dp)
    if (vs%status /= irregular_status_invalid_input) error stop 'empty events status'

    print '(a)', 'test_signal_irregular passed'
contains
    subroutine check(a,b,tol)
        real(dp), intent(in) :: a(:) !! actual results
        real(dp), intent(in) :: b(:) !! expected results
        real(dp), intent(in) :: tol !! maximum absolute error
        if (size(a) /= size(b)) error stop 'irregular array size mismatch'
        if (any(abs(a-b) > tol)) then
            print *, 'actual',a
            print *, 'reference',b
            error stop 'irregular numerical mismatch'
        end if
    end subroutine check
    subroutine reject(a,code)
        real(dp), intent(in) :: a(:) !! result from invalid arguments
        integer, intent(in) :: code !! status to test
        if (code /= irregular_status_invalid_input) error stop 'expected invalid input'
        if (.not. all(ieee_is_nan(a))) error stop 'expected quiet NaN values'
    end subroutine reject
end program test_signal_irregular
