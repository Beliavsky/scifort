! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style signal processing compatibility layer.
module scifort_signal
    use scifort_kinds, only : dp
    use scifort_signal_convolve, only : convolve, correlate, correlation_lags, choose_conv_method, &
        convolution_status_success, convolution_status_invalid_input
    use scifort_signal_convolve2d, only : convolve2d, correlate2d
    use scifort_signal_convolve_nd, only : convolve_nd, correlate_nd, &
        convolution_shape_nd, convolve_nd_success, convolve_nd_invalid_input
    use scifort_signal_fftconvolve, only : fftconvolve, fftconvolve_nd, oaconvolve, oaconvolve_nd
    use scifort_signal_detrend, only : detrend, detrend_status_success, &
        detrend_status_invalid_input
    use scifort_signal_spline, only : cspline1d, qspline1d, cspline1d_eval, qspline1d_eval, &
        spline_status_success, spline_status_invalid_input
    use scifort_signal_spline2d, only : cspline2d, qspline2d, spline_filter
    use scifort_signal_spline2d_eval, only : cspline2d_eval, qspline2d_eval, &
        cspline2d_grid, qspline2d_grid
    use scifort_signal_spline_nd, only : spline_filter_nd, spline_filter1d, &
        map_coordinates_nd, map_coordinates_grid_nd, &
        spline_interpolation_plan, prepare_spline_plan_nd, &
        evaluate_spline_plan_nd, evaluate_spline_plan_grid_nd, &
        evaluate_spline_plan_with_gradient_nd, &
        evaluate_spline_plan_grid_with_gradient_nd, &
        evaluate_spline_plan_affine_with_gradient_nd, &
        evaluate_spline_plan_with_hessian_nd, evaluate_spline_plan_grid_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_gradient_nd, &
        evaluate_spline_plan_affine_with_hessian_nd, &
        evaluate_spline_plan_displacement_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_hessian_nd, &
        spline_axis_grid, evaluate_spline_plan_tensor_nd, &
        evaluate_spline_plan_affine_nd, evaluate_spline_plan_geometric_nd, &
        evaluate_spline_plan_geometric_batch_nd, &
        evaluate_spline_plan_displacement_nd, evaluate_spline_plan_displacement_grid_nd, &
        evaluate_spline_plan_displacement_with_gradient_nd, &
        ndimage_coordinate_batch_map, cspline_nd, qspline_nd, &
        spline_nd_status_too_expensive, &
        affine_transform_nd, shift_nd, zoom_nd, rotate_nd, &
        zoom_shape_nd, rotate_shape_nd
    use scifort_signal_denoise, only : medfilt, medfilt2d, wiener, &
        denoise_status_success, denoise_status_invalid_input
    use scifort_signal_spatial, only : order_filter, sepfir2d, &
        spatial_status_success, spatial_status_invalid_input
    use scifort_signal_deconvolve, only : deconvolve, &
        deconvolution_status_success, deconvolution_status_invalid_input
    use scifort_signal_envelope, only : envelope, envelope_into, &
        envelope_status_success, envelope_status_invalid_input
    use scifort_signal_peak_excursion, only : peak_detect, &
        excursion_status_success, excursion_status_invalid_input
    use scifort_signal_relative_extrema, only : argrelmax, argrelmin, argrelextrema, &
        extrema_comparator, extrema_status_success, extrema_status_invalid_input, &
        argrelextrema_mask_into, argrelmax_mask_into, argrelmin_mask_into
    use scifort_signal_place_poles, only : place_poles, pole_placement_result, &
        pole_status_success, pole_status_invalid_input, pole_status_uncontrollable, &
        pole_status_unsupported, pole_status_numerical_failure
    use scifort_signal_cwt_peaks, only : find_peaks_cwt, cwt_ricker, cwt_wavelet, &
        cwt_wavelet_callback, ricker_wavelet, &
        cwt_complex_wavelet, cwt_complex_wavelet_callback, cwt_morlet2, morlet2_wavelet, &
        find_peaks_cwt_complex, find_peaks_cwt_magnitude, &
        cwt_peak_success, cwt_peak_invalid_input
    use scifort_signal_peaks, only : find_peaks, peak_prominences, peak_widths, &
        peak_result, peak_status_success, peak_status_invalid_input
    use scifort_signal_savgol, only : savgol_coeffs, savgol_filter, savgol_filter_nd, savgol_irregular, &
        savgol_status_success, savgol_status_invalid_input, savgol_status_singular
    use scifort_signal_sos, only : sosfilt, sosfilt_nd, sosfilt_zi, &
        sos_status_success, sos_status_invalid_input, sos_status_no_steady_state
    use scifort_signal_filtfilt_nd, only : filtfilt_nd, sosfiltfilt_nd
    use scifort_signal_filtfilt_gust, only : filtfilt_gust, sosfiltfilt_gust, &
        gust_status_success, gust_status_invalid_input, gust_status_numerical_failure
    use scifort_signal_filtfilt, only : filtfilt, sosfiltfilt, &
        filtfilt_status_success, filtfilt_status_invalid_input, filtfilt_status_no_steady_state, &
        filtfilt_status_numerical_failure
    use scifort_signal_lfilter, only : lfilter, lfilter_nd, lfilter_zi, lfiltic, &
        lfilter_status_success, lfilter_status_invalid_input, lfilter_status_no_steady_state
    use scifort_signal_spectral, only : get_window, general_cosine, chebwin, taylor, gauss_spline, &
        kaiser_bessel_derived, &
        periodogram_real => periodogram, &
        welch_real => welch, csd_real => csd, coherence_real => coherence, &
        stft_real => stft, spectrogram_real => spectrogram, &
        spectral_result, cross_spectral_result, &
        stft_result, spectrogram_result, spectral_status_success, spectral_status_invalid_input, &
        fft_plan, forward_fft, forward_fft_planned, spectral_real_detrend_callback
    use scifort_signal_timefreq_nd, only : stft_into, spectrogram_into
    use scifort_signal_istft_nd, only : istft_into
    use scifort_signal_istft_stream, only : istft_stream, istft_stream_init, &
        istft_stream_push, istft_stream_finish
    use scifort_signal_stft_stream, only : stft_stream, stft_stream_init, &
        stft_stream_push, stft_stream_finish, stft_stream_frequencies, &
        spectrogram_stream_init, spectrogram_stream_push, spectrogram_stream_finish
    use scifort_signal_stream_multi, only : stft_stream_multi, stft_stream_multi_init, &
        stft_stream_multi_push, stft_stream_multi_finish, stft_stream_multi_frequencies, &
        spectrogram_stream_multi_init, spectrogram_stream_multi_push, &
        spectrogram_stream_multi_finish, istft_stream_multi, istft_stream_multi_init, &
        istft_stream_multi_push, istft_stream_multi_finish
    use scifort_signal_spectral_nd, only : &
        periodogram_into, welch_into, csd_into, coherence_into
    use scifort_signal_symiir, only : symiirorder1, symiirorder2, &
        symiir_status_success, symiir_status_invalid_input, symiir_status_no_convergence
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_time_fft_result_3d, short_time_fft_result_4d, &
        short_time_fft_result_5d, short_time_fft_result_6d, &
        short_time_fft_result_7d, &
        short_time_fft_result_8d, &
        short_time_fft_result_9d, &
        short_time_fft_result_10d, &
        short_time_fft_result_11d, &
        short_time_fft_result_12d, &
        short_time_fft_result_13d, &
        short_time_fft_result_14d, &
        short_time_fft_result_15d, &
        short_time_spectrogram_result, short_time_spectrogram_result_3d, &
        short_time_spectrogram_result_4d, short_time_spectrogram_result_5d, &
        short_time_spectrogram_result_7d, &
        short_time_spectrogram_result_8d, &
        short_time_spectrogram_result_9d, &
        short_time_spectrogram_result_10d, &
        short_time_spectrogram_result_11d, &
        short_time_spectrogram_result_12d, &
        short_time_spectrogram_result_13d, &
        short_time_spectrogram_result_14d, &
        short_time_spectrogram_result_15d, &
        short_time_spectrogram_result_6d, &
        stft_real_detrend_callback, stft_complex_detrend_callback, &
        short_fft_success, short_fft_invalid_input, short_fft_noninvertible
    use scifort_signal_stft_dual, only : closest_stft_dual_window, &
        stft_dual_status_success, stft_dual_status_invalid_input, stft_dual_status_noninvertible
    use scifort_signal_istft, only : istft, check_cola, check_nola, &
        istft_status_success, istft_status_invalid_input, istft_status_noninvertible
    use scifort_signal_dpss, only : dpss, dpss_result, &
        dpss_status_success, dpss_status_invalid_input, dpss_status_no_convergence
    use scifort_signal_multitaper, only : multitaper_psd, multitaper_csd, &
        multitaper_coherence, multitaper_spectral_matrix, multitaper_effective_tapers, &
        multitaper_matrix_result, multitaper_channels_result, &
        multitaper_psd_ci, multitaper_confidence_result, &
        multitaper_spectrogram, multitaper_spectrogram_result, &
        multitaper_spectrogram_channels_result, multitaper_harmonic_ftest, multitaper_harmonic_result, &
        multitaper_harmonic_scan, multitaper_harmonic_scan_result, &
        multitaper_noise_floor, multitaper_noise_floor_result, multitaper_noise_floor_channels_result, &
        multitaper_noise_floor_spectrogram, multitaper_noise_floor_spectrogram_result, &
        multitaper_noise_floor_spectrogram_channels_result, &
        multitaper_spectral_matrix_spectrogram, multitaper_spectral_matrix_spectrogram_result, &
        multitaper_cross_spectrogram, multitaper_cross_spectrogram_result, &
        multitaper_tf_coherence, multitaper_tf_coherence_result, multitaper_stream
    use scifort_signal_fixed_taper_spectrogram, only : &
        fixed_taper_spectral_matrix_spectrogram
    use scifort_signal_multiple_coherence, only : multiple_coherence_matrix, &
        multiple_coherence_tf_matrix, multitaper_multiple_coherence, &
        multitaper_tf_multiple_coherence, multiple_coherence_result, multiple_coherence_tf_result
    use scifort_signal_group_coherence, only : group_coherence_matrix, group_coherence_tf_matrix, &
        band_group_coherence_matrix, band_group_coherence_tf_matrix, &
        multitaper_group_coherence, multitaper_tf_group_coherence, &
        multitaper_band_group_coherence, multitaper_tf_band_group_coherence, &
        group_coherence_result, group_coherence_tf_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_matrix, &
        canonical_coherence_tf_matrix, band_canonical_coherence_matrix, &
        band_canonical_coherence_tf_matrix, multitaper_canonical_coherence, &
        multitaper_tf_canonical_coherence, multitaper_band_canonical_coherence, &
        multitaper_tf_band_canonical_coherence, canonical_coherence_result, &
        canonical_coherence_tf_result
    use scifort_signal_canonical_jackknife, only : canonical_window_jackknife, &
        multitaper_canonical_window_jackknife, canonical_window_jackknife_result, &
        canonical_block_jackknife, multitaper_canonical_block_jackknife
    use scifort_signal_canonical_block_selection, only : select_canonical_block_length, &
        select_canonical_influence_block_length, canonical_block_length_result
    use scifort_signal_var1_taper_spectrum, only : var1_expected_taper_spectrum, &
        var1_taper_spectrum_result, var1_spectrum_success, &
        var1_spectrum_invalid_input, var1_spectrum_nonstationary
    use scifort_signal_varp_taper_spectrum, only : varp_expected_taper_spectrum, &
        varp_expected_band_taper_spectrum, varp_taper_spectrum_result, &
        varp_spectrum_success, varp_spectrum_invalid_input, varp_spectrum_nonstationary
    use scifort_signal_varp_canonical_target, only : varp_expected_canonical_coherence, &
        varp_expected_band_canonical_coherence
    use scifort_signal_varp_multitaper_target, only : varp_expected_multitaper_spectrum, &
        varp_expected_band_multitaper_spectrum, &
        varp_expected_multitaper_canonical_coherence, &
        varp_expected_band_multitaper_canonical_coherence
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap, &
        multitaper_canonical_moving_block_bootstrap, canonical_moving_block_bootstrap_result
    use scifort_signal_canonical_tracking, only : track_canonical_modes, &
        canonical_mode_tracks_result, canonical_mode_tracks_tf_result, &
        canonical_tracking_success, canonical_tracking_invalid_input
    use scifort_signal_canonical_subspace_tracking, only : track_canonical_subspaces, &
        canonical_subspace_tracks_result, canonical_subspace_tracks_tf_result, &
        canonical_subspace_success, canonical_subspace_invalid_input
    use scifort_signal_band_coherence, only : integrate_spectral_bands, &
        integrate_tf_spectral_bands, band_multiple_coherence_matrix, &
        band_multiple_coherence_tf_matrix, band_partial_coherence_matrix, &
        band_partial_coherence_tf_matrix, multitaper_band_multiple_coherence, &
        multitaper_tf_band_multiple_coherence, multitaper_band_partial_coherence, &
        multitaper_tf_band_partial_coherence
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, &
        multitaper_partial_coherence, partial_coherence_result, &
        partial_coherence_tf_matrix, multitaper_tf_partial_coherence, partial_coherence_tf_result, &
        partial_coherence_success, partial_coherence_invalid_input, partial_coherence_singular
    use scifort_signal_spectral_complex, only : periodogram_complex, welch_complex, &
        csd_complex, coherence_complex, stft_complex, spectrogram_complex, &
        csd_real_complex, csd_complex_real, coherence_real_complex, coherence_complex_real, &
        spectral_complex_detrend_callback
    use scifort_signal_design, only : firwin, kaiser_beta, kaiser_atten, kaiserord, &
        kaiser_order_result, freqz, freqz_zpk, sosfreqz, sos2tf, group_delay, group_delay_sos, &
        group_delay_result, group_delay_zpk, butter, cheby1, &
        frequency_response, butter_result, design_status_success, &
        design_status_invalid_input, design_status_singular
    ! SciPy 1.17 renamed sosfreqz to freqz_sos; retain both spellings.
    use scifort_signal_design, only : freqz_sos => sosfreqz
    use scifort_signal_fir2d, only : firwin_2d
    use scifort_signal_fir_advanced, only : firwin2, firls, minimum_phase, &
        fir_advanced_status_success, fir_advanced_status_invalid_input, fir_advanced_status_singular
    use scifort_signal_unique_roots, only : unique_roots, &
        unique_roots_status_success, unique_roots_status_invalid_input
    use scifort_signal_czt_nd, only : czt_into, czt_plan_into, zoom_fft_into, zoom_fft_plan_into
    use scifort_signal_czt, only : czt, zoom_fft, czt_points, czt_plan, zoom_fft_plan, czt_status_success, &
        czt_status_invalid_input, czt_status_nonfinite_result
    use scifort_signal_remez, only : remez, remez_status_success, &
        remez_status_invalid_input, remez_status_no_convergence
    use scifort_signal_analog, only : freqs, freqs_zpk, findfreqs
    use scifort_signal_ellipap, only : ellipap
    use scifort_signal_bilinear, only : bilinear, bilinear_zpk
    use scifort_signal_analog_prototypes, only : buttap, cheb1ap, cheb2ap, besselap, &
        analog_prototype_result, prototype_success, prototype_invalid_input, &
        prototype_no_convergence, prototype_too_expensive
    use scifort_signal_frequency_transforms, only : lp2lp, lp2hp, lp2bp, lp2bs, &
        lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, &
        analog_tf_transform_result, analog_zpk_transform_result, &
        frequency_transform_success, frequency_transform_invalid_input, &
        frequency_transform_too_expensive
    use scifort_signal_zpk_conversion, only : tf2sos, tf2zpk, zpk2tf, zpk2sos, sos2zpk, &
        filter_zpk_result, filter_tf_result, filter_sos_result, &
        zpk_conversion_success, zpk_conversion_invalid_input, zpk_conversion_numerical_failure
    use scifort_signal_abcd_normalize, only : abcd_normalize
    use scifort_signal_tf_normalize, only : normalize, normalized_tf_result, &
        normalize_status_success, normalize_status_invalid_input
    use scifort_signal_whittaker, only : whittaker_henderson, whittaker_result, &
        whittaker_status_success, whittaker_status_invalid_input, whittaker_status_singular
    use scifort_signal_state_space, only : tf2ss, ss2tf, zpk2ss, ss2zpk, &
        state_space_result, state_space_tf_result, &
        state_space_success, state_space_invalid_input, state_space_numerical_failure
    use scifort_signal_dlti, only : dlti, dlti_tf, dlti_zpk, dlti_ss, dlti_system, dlti_gramian_result, &
        dlti_balanced_result, &
        dlti_success, dlti_invalid_input, dlti_numerical_failure, &
        dlti_transfer_function, dlti_zeros_poles_gain, dlti_state_space
    use scifort_signal_clti, only : clti, clti_tf, clti_zpk, clti_ss, clti_system, &
        clti_gramian_result, clti_balanced_result, &
        clti_success, clti_invalid_input, clti_numerical_failure
    use scifort_signal_discrete_simulation, only : dlsim, dstep, dimpulse, dlsim_sos, dlsim_sos_bank, &
        dstep_sos, dimpulse_sos, dlsim_zpk, dstep_zpk, dimpulse_zpk, &
        discrete_simulation_result, discrete_sos_bank_result, discrete_response_result, &
        discrete_sim_success, discrete_sim_invalid_input, discrete_sim_numerical_failure
    use scifort_signal_continuous_simulation, only : lsim, step, impulse, &
        lsim_zpk, step_zpk, impulse_zpk, continuous_simulation_result, continuous_response_result, &
        continuous_sim_success, continuous_sim_invalid_input, continuous_sim_numerical_failure
    use scifort_signal_lti_frequency, only : freqresp, bode, dfreqresp, dbode, &
        freqresp_zpk, bode_zpk, dfreqresp_zpk, dbode_zpk, bode_result, &
        freqresp_mimo, dfreqresp_mimo, bode_mimo, dbode_mimo, &
        mimo_frequency_response, mimo_bode_result, mimo_frequency_plan, &
        prepare_mimo_frequency
    use scifort_signal_cont2discrete, only : cont2discrete, cont2discrete_zpk, &
        discrete_state_space_result, discrete_tf_result, discrete_zpk_result, &
        cont2discrete_success, cont2discrete_invalid_input, cont2discrete_numerical_failure
    use scifort_signal_partial_fraction, only : residue, residuez, invres, invresz, &
        partial_fraction_result, partial_fraction_tf_result, partial_fraction_success, &
        partial_fraction_invalid_input, partial_fraction_singular
    use scifort_signal_analog_sos, only : iirfilter_zpk, iirfilter_sos, &
        zpk2sos_analog, sosfreqs, analog_zpk_result, analog_sos_result, &
        analog_output_status_success, analog_output_status_invalid_input
    use scifort_signal_bessel_design, only : bessel
    use scifort_signal_band_stop_obj, only : band_stop_obj, band_stop_success, band_stop_invalid_input
    use scifort_signal_iir_extra, only : cheby2, ellip, iirfilter, iirdesign, &
        buttord, cheb1ord, cheb2ord, ellipord, iir_order_result, iirnotch, iirpeak, iircomb, gammatone
    use scifort_signal_hilbert_nd, only : hilbert_into, hilbert2_into
    use scifort_signal_hilbert, only : hilbert, hilbert2, &
        hilbert_status_success, hilbert_status_invalid_input
    use scifort_signal_waveforms, only : chirp, chirp_complex, sweep_poly, &
        gausspulse, gausspulse_cutoff, waveform_status_success, waveform_status_invalid_input
    use scifort_signal_periodic, only : sawtooth, square, unit_impulse, &
        periodic_status_success, periodic_status_invalid_input
    use scifort_signal_sequences, only : max_len_seq, max_len_seq_result, &
        sequence_status_success, sequence_status_invalid_input, sequence_status_allocation_failed
    use scifort_signal_irregular_nd, only : lombscargle_into, lombscargle_amplitude_into, &
        lombscargle_plan, lombscargle_plan_into, lombscargle_plan_amplitude_into
    use scifort_signal_irregular, only : lombscargle, lombscargle_amplitude, &
        lombscargle_multiharmonic, lombscargle_harmonic_fit, lombscargle_harmonic_predict, &
        lombscargle_harmonic_predict_std, &
        lombscargle_harmonic_prepare, lombscargle_harmonic_fit_prepared, &
        lombscargle_harmonic_fit_batch, harmonic_fit_plan, &
        harmonic_fit_result, vectorstrength, vector_strength_result, irregular_status_success, &
        irregular_status_invalid_input, irregular_status_singular
    use scifort_signal_fourier_resample, only : resample, &
        fourier_resample_status_success, fourier_resample_status_invalid_input, &
        real_frequency_window_callback, complex_frequency_window_callback
    use scifort_signal_multirate_nd_into, only : &
        resample_into, resample_poly_into, upfirdn_into, decimate_into
    use scifort_signal_resample, only : upfirdn, resample_poly, decimate, &
        resample_status_success, resample_status_invalid_input
    implicit none
    private
    interface periodogram
        procedure :: periodogram_real, periodogram_complex
    end interface
    interface welch
        procedure :: welch_real, welch_complex
    end interface
    interface csd
        procedure :: csd_real, csd_complex, csd_real_complex, csd_complex_real
    end interface
    interface coherence
        procedure :: coherence_real, coherence_complex, coherence_real_complex, coherence_complex_real
    end interface
    interface stft
        procedure :: stft_real, stft_complex
    end interface
    interface spectrogram
        procedure :: spectrogram_real, spectrogram_complex
    end interface
    public :: real_frequency_window_callback, complex_frequency_window_callback
    public :: upfirdn, resample_poly, decimate, resample
    public :: resample_into, resample_poly_into, upfirdn_into, decimate_into
    public :: stft_into, spectrogram_into, istft_into
    public :: istft_stream, istft_stream_init, istft_stream_push, istft_stream_finish
    public :: stft_stream, stft_stream_init, stft_stream_push, stft_stream_finish
    public :: stft_stream_frequencies, spectrogram_stream_init, &
        spectrogram_stream_push, spectrogram_stream_finish
    public :: stft_stream_multi, stft_stream_multi_init, stft_stream_multi_push, &
        stft_stream_multi_finish, stft_stream_multi_frequencies
    public :: spectrogram_stream_multi_init, spectrogram_stream_multi_push, &
        spectrogram_stream_multi_finish
    public :: istft_stream_multi, istft_stream_multi_init, istft_stream_multi_push, &
        istft_stream_multi_finish
    public :: hilbert, hilbert2, hilbert_into, hilbert2_into, hilbert_status_success, hilbert_status_invalid_input
    public :: chirp, chirp_complex, sweep_poly, gausspulse, gausspulse_cutoff
    public :: sawtooth, square, unit_impulse, max_len_seq, max_len_seq_result
    public :: lombscargle_into, lombscargle_amplitude_into
    public :: lombscargle_plan, lombscargle_plan_into, lombscargle_plan_amplitude_into
    public :: lombscargle, lombscargle_amplitude, lombscargle_multiharmonic, vectorstrength, vector_strength_result
    public :: lombscargle_harmonic_fit, lombscargle_harmonic_predict, &
        lombscargle_harmonic_predict_std, harmonic_fit_result, harmonic_fit_plan, &
        lombscargle_harmonic_prepare, lombscargle_harmonic_fit_prepared, lombscargle_harmonic_fit_batch
    public :: irregular_status_success, irregular_status_invalid_input, irregular_status_singular
    public :: multiple_coherence_matrix, multiple_coherence_tf_matrix
    public :: multitaper_multiple_coherence, multitaper_tf_multiple_coherence
    public :: multiple_coherence_result, multiple_coherence_tf_result
    public :: var1_expected_taper_spectrum, var1_taper_spectrum_result
    public :: var1_spectrum_success, var1_spectrum_invalid_input, var1_spectrum_nonstationary
    public :: varp_expected_taper_spectrum, varp_expected_band_taper_spectrum
    public :: varp_expected_canonical_coherence, varp_expected_band_canonical_coherence
    public :: varp_expected_multitaper_spectrum, varp_expected_band_multitaper_spectrum
    public :: varp_expected_multitaper_canonical_coherence, &
        varp_expected_band_multitaper_canonical_coherence
    public :: varp_taper_spectrum_result, varp_spectrum_success, &
        varp_spectrum_invalid_input, varp_spectrum_nonstationary
    public :: integrate_spectral_bands, integrate_tf_spectral_bands
    public :: band_multiple_coherence_matrix, band_multiple_coherence_tf_matrix
    public :: band_partial_coherence_matrix, band_partial_coherence_tf_matrix
    public :: multitaper_band_multiple_coherence, multitaper_tf_band_multiple_coherence
    public :: multitaper_band_partial_coherence, multitaper_tf_band_partial_coherence
    public :: partial_coherence_matrix, multitaper_partial_coherence, partial_coherence_result
    public :: multitaper_effective_tapers
    public :: partial_coherence_tf_matrix, multitaper_tf_partial_coherence, partial_coherence_tf_result
    public :: partial_coherence_success, partial_coherence_invalid_input, partial_coherence_singular
    public :: periodic_status_success, periodic_status_invalid_input
    public :: sequence_status_success, sequence_status_invalid_input, sequence_status_allocation_failed
    public :: waveform_status_success, waveform_status_invalid_input
    public :: fourier_resample_status_success, fourier_resample_status_invalid_input
    public :: resample_status_success, resample_status_invalid_input
    public :: firwin, firwin_2d, firwin2, firls, minimum_phase, remez
    public :: kaiser_beta, kaiser_atten, kaiserord, kaiser_order_result
    public :: unique_roots, unique_roots_status_success, unique_roots_status_invalid_input
    public :: symiirorder1, symiirorder2, symiir_status_success, &
        symiir_status_invalid_input, symiir_status_no_convergence
    public :: short_time_fft, short_time_fft_result, short_time_fft_result_3d, &
        short_time_fft_result_4d, short_time_fft_result_5d, &
        short_time_fft_result_6d, &
        short_time_fft_result_7d, &
        short_time_fft_result_8d, &
        short_time_fft_result_9d, &
        short_time_fft_result_10d, &
        short_time_fft_result_11d, &
        short_time_fft_result_12d, &
        short_time_fft_result_13d, &
        short_time_fft_result_14d, &
        short_time_fft_result_15d, &
        short_time_spectrogram_result, &
        short_time_spectrogram_result_3d, short_time_spectrogram_result_4d, &
        short_time_spectrogram_result_5d, short_time_spectrogram_result_6d, &
        short_time_spectrogram_result_7d, &
        short_time_spectrogram_result_8d, &
        short_time_spectrogram_result_9d, &
        short_time_spectrogram_result_10d, &
        short_time_spectrogram_result_11d, &
        short_time_spectrogram_result_12d, &
        short_time_spectrogram_result_13d, &
        short_time_spectrogram_result_14d, &
        short_time_spectrogram_result_15d, &
        short_fft_success, short_fft_invalid_input, short_fft_noninvertible, &
        stft_real_detrend_callback, stft_complex_detrend_callback
    public :: closest_stft_dual_window, stft_dual_status_success, &
        stft_dual_status_invalid_input, stft_dual_status_noninvertible
    public :: czt, zoom_fft, czt_points, czt_plan, zoom_fft_plan, &
        czt_into, czt_plan_into, zoom_fft_into, zoom_fft_plan_into, &
        czt_status_success, czt_status_invalid_input, czt_status_nonfinite_result
    public :: remez_status_success, remez_status_invalid_input, remez_status_no_convergence
    public :: fir_advanced_status_success, fir_advanced_status_invalid_input, fir_advanced_status_singular
    public :: freqz, freqz_zpk, freqs, freqs_zpk, findfreqs, sosfreqz, freqz_sos, sos2tf, &
        group_delay, group_delay_sos, group_delay_zpk, group_delay_result
    public :: band_stop_obj, band_stop_success, band_stop_invalid_input
    public :: buttap, cheb1ap, cheb2ap, besselap, ellipap, analog_prototype_result
    public :: bilinear, bilinear_zpk
    public :: prototype_success, prototype_invalid_input, prototype_no_convergence, prototype_too_expensive
    public :: lp2lp, lp2hp, lp2bp, lp2bs
    public :: lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk
    public :: analog_tf_transform_result, analog_zpk_transform_result
    public :: frequency_transform_success, frequency_transform_invalid_input
    public :: frequency_transform_too_expensive
    public :: tf2sos, tf2zpk, zpk2tf, zpk2sos, sos2zpk, filter_zpk_result, filter_tf_result, filter_sos_result
    public :: zpk_conversion_success, zpk_conversion_invalid_input, zpk_conversion_numerical_failure
    public :: tf2ss, ss2tf, zpk2ss, ss2zpk, state_space_result, state_space_tf_result
    public :: state_space_success, state_space_invalid_input, state_space_numerical_failure
    public :: freqresp, bode, dfreqresp, dbode, bode_result
    public :: freqresp_mimo, dfreqresp_mimo, bode_mimo, dbode_mimo
    public :: mimo_frequency_response, mimo_bode_result, mimo_frequency_plan
    public :: prepare_mimo_frequency
    public :: freqresp_zpk, bode_zpk, dfreqresp_zpk, dbode_zpk
    public :: dlti, dlti_tf, dlti_zpk, dlti_ss, dlti_system, dlti_gramian_result
    public :: dlti_balanced_result
    public :: clti, clti_tf, clti_zpk, clti_ss, clti_system
    public :: clti_gramian_result, clti_balanced_result
    public :: clti_success, clti_invalid_input, clti_numerical_failure
    public :: dlti_success, dlti_invalid_input, dlti_numerical_failure
    public :: dlti_transfer_function, dlti_zeros_poles_gain, dlti_state_space
    public :: dlsim, dstep, dimpulse, dlsim_zpk, dstep_zpk, dimpulse_zpk
    public :: dlsim_sos, dlsim_sos_bank, dstep_sos, dimpulse_sos
    public :: lsim, step, impulse, lsim_zpk, step_zpk, impulse_zpk
    public :: continuous_simulation_result, continuous_response_result
    public :: continuous_sim_success, continuous_sim_invalid_input, continuous_sim_numerical_failure
    public :: discrete_simulation_result, discrete_sos_bank_result, discrete_response_result
    public :: discrete_sim_success, discrete_sim_invalid_input, discrete_sim_numerical_failure
    public :: cont2discrete, cont2discrete_zpk, discrete_state_space_result, discrete_tf_result, discrete_zpk_result
    public :: cont2discrete_success, cont2discrete_invalid_input, cont2discrete_numerical_failure
    public :: residue, residuez, invres, invresz, partial_fraction_result, partial_fraction_tf_result
    public :: partial_fraction_success, partial_fraction_invalid_input, partial_fraction_singular
    public :: iirfilter_zpk, iirfilter_sos, zpk2sos_analog, sosfreqs
    public :: analog_zpk_result, analog_sos_result
    public :: analog_output_status_success, analog_output_status_invalid_input
    public :: butter, cheby1, cheby2, ellip, bessel, iirfilter, &
        iirnotch, iirpeak, iircomb, gammatone, iirdesign, buttord, cheb1ord, cheb2ord, ellipord, &
        iir_order_result, frequency_response, butter_result
    public :: design_status_success, design_status_invalid_input, design_status_singular
    public :: savgol_coeffs, savgol_filter, savgol_filter_nd, savgol_irregular
    public :: savgol_status_success, savgol_status_invalid_input, savgol_status_singular
    public :: lfilter, lfilter_nd, lfilter_zi, lfiltic, sosfilt, sosfilt_nd, sosfilt_zi, filtfilt, sosfiltfilt
    public :: filtfilt_gust, sosfiltfilt_gust, gust_status_success, gust_status_invalid_input
    public :: gust_status_numerical_failure, filtfilt_status_numerical_failure
    public :: filtfilt_nd, sosfiltfilt_nd
    public :: sos_status_success, sos_status_invalid_input, sos_status_no_steady_state
    public :: filtfilt_status_success, filtfilt_status_invalid_input, filtfilt_status_no_steady_state
    public :: lfilter_status_success, lfilter_status_invalid_input, lfilter_status_no_steady_state
    public :: dp
    public :: periodogram_into, welch_into, csd_into, coherence_into
    public :: get_window, general_cosine, chebwin, taylor, kaiser_bessel_derived, gauss_spline, &
        periodogram, welch, csd, coherence, stft, spectrogram
    public :: dpss, dpss_result, dpss_status_success, dpss_status_invalid_input, &
        dpss_status_no_convergence, multitaper_psd, multitaper_csd, &
        multitaper_coherence, multitaper_spectral_matrix, &
        multitaper_matrix_result, multitaper_channels_result, &
        multitaper_psd_ci, multitaper_confidence_result, &
        multitaper_spectrogram, multitaper_spectrogram_result, &
        multitaper_spectrogram_channels_result, multitaper_harmonic_ftest, multitaper_harmonic_result, &
        multitaper_harmonic_scan, multitaper_harmonic_scan_result, &
        multitaper_noise_floor, multitaper_noise_floor_result, multitaper_noise_floor_channels_result, &
        multitaper_noise_floor_spectrogram, multitaper_noise_floor_spectrogram_result, &
        multitaper_noise_floor_spectrogram_channels_result, &
        multitaper_spectral_matrix_spectrogram, multitaper_spectral_matrix_spectrogram_result, &
        multitaper_cross_spectrogram, multitaper_cross_spectrogram_result, &
        multitaper_tf_coherence, multitaper_tf_coherence_result, multitaper_stream
    public :: fixed_taper_spectral_matrix_spectrogram
    public :: group_coherence_matrix, group_coherence_tf_matrix, &
        band_group_coherence_matrix, band_group_coherence_tf_matrix, &
        multitaper_group_coherence, multitaper_tf_group_coherence, &
        multitaper_band_group_coherence, multitaper_tf_band_group_coherence, &
        group_coherence_result, group_coherence_tf_result
    public :: canonical_coherence_matrix, canonical_coherence_tf_matrix, &
        band_canonical_coherence_matrix, band_canonical_coherence_tf_matrix, &
        multitaper_canonical_coherence, multitaper_tf_canonical_coherence, &
        multitaper_band_canonical_coherence, multitaper_tf_band_canonical_coherence, &
        canonical_coherence_result, canonical_coherence_tf_result
    public :: canonical_window_jackknife, multitaper_canonical_window_jackknife, &
        canonical_window_jackknife_result, canonical_block_jackknife, &
        multitaper_canonical_block_jackknife
    public :: select_canonical_block_length, select_canonical_influence_block_length, &
        canonical_block_length_result
    public :: canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap, &
        canonical_moving_block_bootstrap_result
    public :: track_canonical_modes, canonical_mode_tracks_result, canonical_mode_tracks_tf_result
    public :: canonical_tracking_success, canonical_tracking_invalid_input
    public :: track_canonical_subspaces, canonical_subspace_tracks_result, &
        canonical_subspace_tracks_tf_result
    public :: canonical_subspace_success, canonical_subspace_invalid_input
    public :: istft, check_cola, check_nola
    public :: istft_status_success, istft_status_invalid_input, istft_status_noninvertible
    public :: spectral_result, cross_spectral_result, stft_result, spectrogram_result
    public :: spectral_status_success, spectral_status_invalid_input
    public :: spectral_real_detrend_callback, spectral_complex_detrend_callback
    public :: fft_plan, forward_fft, forward_fft_planned
    public :: convolve, correlate, correlation_lags, choose_conv_method
    public :: convolve2d, correlate2d, fftconvolve, fftconvolve_nd, oaconvolve, oaconvolve_nd
    public :: convolve_nd, correlate_nd, convolution_shape_nd
    public :: convolve_nd_success, convolve_nd_invalid_input
    public :: normalize, normalized_tf_result, normalize_status_success, normalize_status_invalid_input
    public :: abcd_normalize, whittaker_henderson, whittaker_result
    public :: whittaker_status_success, whittaker_status_invalid_input, whittaker_status_singular
    public :: convolution_status_success, convolution_status_invalid_input
    public :: place_poles, pole_placement_result
    public :: pole_status_success, pole_status_invalid_input, pole_status_uncontrollable
    public :: pole_status_unsupported, pole_status_numerical_failure
    public :: find_peaks_cwt, cwt_ricker, cwt_wavelet, cwt_wavelet_callback, ricker_wavelet
    public :: cwt_complex_wavelet, cwt_complex_wavelet_callback, cwt_morlet2, morlet2_wavelet
    public :: find_peaks_cwt_complex, find_peaks_cwt_magnitude
    public :: cwt_peak_success, cwt_peak_invalid_input
    public :: find_peaks, peak_prominences, peak_widths, peak_result
    public :: peak_detect, excursion_status_success, excursion_status_invalid_input
    public :: argrelextrema, argrelmax, argrelmin, extrema_comparator
    public :: argrelextrema_mask_into, argrelmax_mask_into, argrelmin_mask_into
    public :: extrema_status_success, extrema_status_invalid_input
    public :: peak_status_success, peak_status_invalid_input
    public :: detrend, detrend_status_success, detrend_status_invalid_input
    public :: cspline1d, qspline1d, cspline1d_eval, qspline1d_eval
    public :: cspline2d, qspline2d, spline_filter
    public :: cspline2d_eval, qspline2d_eval, cspline2d_grid, qspline2d_grid
    public :: spline_filter_nd, spline_filter1d, map_coordinates_nd
    public :: map_coordinates_grid_nd, cspline_nd, qspline_nd
    public :: spline_interpolation_plan, prepare_spline_plan_nd
    public :: evaluate_spline_plan_nd, evaluate_spline_plan_grid_nd
    public :: evaluate_spline_plan_with_gradient_nd
    public :: evaluate_spline_plan_grid_with_gradient_nd
    public :: evaluate_spline_plan_affine_with_gradient_nd
    public :: evaluate_spline_plan_with_hessian_nd, evaluate_spline_plan_grid_with_hessian_nd
    public :: evaluate_spline_plan_displacement_grid_with_gradient_nd
    public :: evaluate_spline_plan_affine_with_hessian_nd
    public :: evaluate_spline_plan_displacement_with_hessian_nd
    public :: evaluate_spline_plan_displacement_grid_with_hessian_nd
    public :: spline_axis_grid, evaluate_spline_plan_tensor_nd
    public :: evaluate_spline_plan_affine_nd, evaluate_spline_plan_geometric_nd
    public :: evaluate_spline_plan_geometric_batch_nd, ndimage_coordinate_batch_map
    public :: evaluate_spline_plan_displacement_nd, evaluate_spline_plan_displacement_grid_nd
    public :: evaluate_spline_plan_displacement_with_gradient_nd
    public :: affine_transform_nd, shift_nd, zoom_nd, rotate_nd
    public :: zoom_shape_nd, rotate_shape_nd
    public :: spline_nd_status_too_expensive
    public :: spline_status_success, spline_status_invalid_input
    public :: medfilt, medfilt2d, wiener, denoise_status_success, denoise_status_invalid_input
    public :: order_filter, sepfir2d, spatial_status_success, spatial_status_invalid_input
    public :: deconvolve, deconvolution_status_success, deconvolution_status_invalid_input
    public :: envelope, envelope_into, envelope_status_success, envelope_status_invalid_input
end module scifort_signal
