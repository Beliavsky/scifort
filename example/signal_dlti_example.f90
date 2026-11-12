! SPDX-License-Identifier: MIT
! Demonstrate discrete-time transfer functions, composition, and response.
program signal_dlti_example
    use scifort_signal, only : dp, dlti_tf, dlti_zpk, dlti_system, dlti_gramian_result, &
        discrete_simulation_result, discrete_response_result, frequency_response, dlti_success
    implicit none
    type(dlti_system) :: lowpass, first_order, cascade, closed_loop, minimal
    type(dlti_gramian_result) :: gramians
    type(discrete_simulation_result) :: output
    type(discrete_response_result) :: step_response
    type(frequency_response) :: spectrum
    real(dp) :: input(24), w(4)
    integer :: i

    lowpass = dlti_tf([0.25_dp,0.25_dp],[1.0_dp,-0.5_dp],dt=0.01_dp)
    first_order = dlti_zpk([0.2_dp],[0.6_dp],1.0_dp,dt=0.01_dp)
    cascade = lowpass%series(first_order)
    closed_loop = lowpass%feedback(first_order)
    if (closed_loop%status /= dlti_success) error stop 'invalid feedback'
    minimal = closed_loop%minreal()
    if (minimal%status /= dlti_success) error stop 'minimal realization failed'
    gramians = minimal%gramians()
    if (gramians%status /= dlti_success) error stop 'Gramian computation failed'
    if (cascade%status /= dlti_success) error stop 'invalid series connection'

    input = 0.0_dp
    input(1:12) = 1.0_dp
    output = cascade%output(input)
    step_response = cascade%step(n=size(input))
    w = [0.0_dp,0.25_dp,0.5_dp,1.0_dp]
    spectrum = cascade%freqresp(w=w)
    if (output%status /= 0 .or. step_response%status /= 0 .or. spectrum%status /= 0) then
        error stop 'discrete-time system analysis failed'
    end if
    print '(a,es12.4)', 'sample interval: ',cascade%dt
    print '(a,l1)', 'stable: ',cascade%stable()
    print '(a,es12.4)', 'closed-loop H2 norm: ',gramians%h2_norm
    print '(a,i0)', 'closed-loop minimal state count: ',minimal%nstates()
    do i=1,4
        print '(a,i0,a,2(es12.4,1x))', 'H at frequency ',i,': ', &
            real(spectrum%response(i),dp),aimag(spectrum%response(i))
    end do
    print '(a,es12.4)', 'final output value: ', real(output%y(size(input),1),dp)
end program signal_dlti_example
