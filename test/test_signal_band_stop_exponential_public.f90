! SPDX-License-Identifier: MIT
! Compile and exercise the public scifort_signal facade, not just submodules.
program test_signal_band_stop_exponential_public
    use scifort_signal, only : dp, band_stop_obj, band_stop_success, &
        band_stop_invalid_input, get_window
    implicit none
    real(dp), allocatable :: w(:)
    real(dp) :: v
    integer :: st

    v=band_stop_obj(0.5_dp,0,[0.5_dp,4.0_dp],[1.5_dp,2.5_dp],1.0_dp,40.0_dp,'butter',st)
    if (st/=band_stop_success .or. v<0.0_dp) error stop 'public butter objective'
    v=band_stop_obj(0.5_dp,0,[0.5_dp,4.0_dp],[1.5_dp,2.5_dp],1.0_dp,40.0_dp,'ellip',st)
    if (st/=band_stop_success .or. v<0.0_dp) error stop 'public elliptic objective'
    v=band_stop_obj(0.5_dp,10,[0.5_dp,4.0_dp],[1.5_dp,2.5_dp],1.0_dp,40.0_dp,'butter',st)
    if (st/=band_stop_invalid_input) error stop 'public status'
    w=get_window('exponential',8,parameter=2.0_dp,parameter2=0.0_dp,status=st)
    if (st/=0 .or. size(w)/=8) error stop 'public exponential center status'
    if (abs(w(1)-1.0_dp)>1.0e-14_dp) error stop 'public exponential center value'
    print *, 'PASS: public band-stop and explicit-center exponential APIs'
end program test_signal_band_stop_exponential_public
