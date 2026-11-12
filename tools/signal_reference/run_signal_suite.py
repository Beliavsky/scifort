#!/usr/bin/env python3
"""Build/run standalone signal-module tests without FPM (GNU Fortran)."""
from __future__ import annotations

from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
MODULES = [
    'kinds', 'constants', 'math', 'special/bessel_i', 'special/bessel_k',
    'special/bessel_jy', 'special/bessel_order', 'special/bessel_derivatives',
    'special/hankel', 'special/bessel_complex', 'linalg',
    # FFT is a dependency of the direct/FFT convolution dispatcher.
    'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_fftconvolve',
    'signal_denoise', 'signal_spatial',
    'signal_sequences', 'signal_convolve', 'signal_convolve2d', 'signal_convolve_nd', 'signal_detrend', 'signal_peaks',
    'signal_savgol', 'signal_spline', 'signal_spline_nd', 'signal_spline2d', 'signal_spline2d_eval', 'signal_symiir', 'signal_lfilter', 'signal_sos', 'signal_filtfilt_gust', 'signal_filtfilt', 'signal_filtfilt_nd',
    'signal_short_time_fft', 'signal_multitaper', 'signal_stft_dual', 'signal_istft',
    'signal_istft_stream', 'signal_stft_stream', 'signal_stream_multi',
    'signal_spectral_complex', 'signal_spectral_nd', 'signal_timefreq_nd', 'signal_istft_nd', 'signal_fir_advanced',
    'signal_remez', 'signal_czt', 'signal_czt_nd', 'signal_unique_roots', 'signal_analog_design', 'signal_elliptic', 'signal_analog_prototypes', 'signal_frequency_transforms', 'signal_bilinear', 'signal_ellipap', 'signal_design',
    'signal_fir2d', 'signal_analog', 'signal_zpk_conversion', 'signal_state_space', 'signal_cont2discrete', 'signal_discrete_simulation', 'signal_continuous_simulation', 'signal_lti_frequency', 'signal_dlti', 'signal_clti', 'signal_partial_fraction', 'signal_analog_sos', 'signal_bessel_design', 'signal_iir_extra',
    'signal_fourier_resample', 'signal_hilbert', 'signal_hilbert_nd', 'signal_waveforms',
    'signal_periodic', 'signal_irregular', 'signal_irregular_nd', 'signal_resample', 'signal_multirate_nd_into',
    'signal_deconvolve', 'signal_envelope', 'signal_peak_excursion',
    'signal_relative_extrema', 'signal_cwt_peaks', 'signal_place_poles_mimo',
    'signal_place_poles', 'signal_abcd_normalize', 'signal_tf_normalize', 'signal_whittaker',
    'ndimage_distance', 'ndimage_chamfer', 'ndimage_measurements',
    'ndimage_watershed', 'ndimage_morphology', 'ndimage', 'signal',
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized', action='store_true', help='build with -O3')
    parser.add_argument('--match', default='test_signal_*', help='test file stem glob')
    parser.add_argument('--jobs', type=int, default=1, help='number of parallel test processes')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    flags = ['-std=f2018', '-Wimplicit-interface', '-Werror=implicit-interface']
    flags += ['-O3'] if args.optimized else ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    tests = sorted((ROOT / 'test').glob(f'{args.match}.f90'))
    if not tests:
        parser.error('no tests matched --match')
    with tempfile.TemporaryDirectory(prefix='scifort-signal-') as temp:
        work = Path(temp)
        modules = []
        for name in MODULES:
            path = ROOT / 'src' / f'{name}.f90'
            obj = work / f"{name.replace('/', '_')}.o"
            subprocess.run(['gfortran', *flags, '-I',str(work), '-J',str(work),
                '-c', str(path), '-o', str(obj)],check=True,cwd=work)
            modules.append(str(obj))
        def test_one(source):
            exe = work / source.stem
            subprocess.run(['gfortran', *flags, '-I',str(work), '-J',str(work),
                str(source), *modules, '-o',str(exe)],check=True,cwd=work)
            subprocess.run([str(exe)],check=True,cwd=work,stdout=subprocess.DEVNULL)

        if args.jobs == 1:
            for source in tests:
                test_one(source)
        else:
            with ThreadPoolExecutor(max_workers=args.jobs) as executor:
                list(executor.map(test_one, tests))
    print(f'PASS: {len(tests)}/{len(tests)} signal tests, '
          f'{"optimized -O3" if args.optimized else "runtime checked -O0"}')


if __name__ == '__main__':
    main()
