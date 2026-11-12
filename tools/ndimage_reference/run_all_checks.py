#!/usr/bin/env python3
"""Run focused ndimage SciPy comparisons and Fortran regression binaries."""
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=[
    'compare_ndimage_scipy.py',
    'compare_morphology_scipy.py',
    'compare_extended_filters_scipy.py',
    'compare_extended_morphology_scipy.py',
    'compare_distance_scipy.py',
    'compare_distance_random_scipy.py',
    'compare_footprint_feature_scipy.py',
    'compare_feature_random_scipy.py',
    'compare_advanced_gray_scipy.py',
    'compare_binary_origin_scipy.py',
    'compare_measurements_scipy.py',
    'compare_watershed_scipy.py',
    'compare_regions_masks_rank_scipy.py',
    'compare_extrema_histogram_scipy.py',
    'compare_chamfer_random_scipy.py',
    'compare_callbacks_weighted_scipy.py',
    'compare_weighted_stencil_random.py',
    'compare_mixed_complex_scipy.py',
    'compare_fourier_scipy.py',
    'compare_ellipsoid_order_line_scipy.py',
    'compare_bruteforce_generic_scipy.py',
    'compare_structure_axes_scipy.py',
    'compare_axis_extrema_scipy.py',
    'compare_selected_axes_scipy.py',
    'compare_gaussian_axes_scipy.py',
    'compare_order_axes_running_scipy.py',
    'compare_gray_generic_axes_scipy.py',
    'compare_complex_line_binary_axes_scipy.py',
    'compare_value_indices_scipy.py',
    'compare_complex_derivatives_labels_scipy.py',
    'compare_propagation_axes_scipy.py',
    'compare_geometric_comprehension_scipy.py',
    'compare_label_axes_scipy.py',
    'compare_vectorized_filter_scipy.py',
    'compare_vectorized_multi_scipy.py',
    'compare_vectorized_mixed_scipy.py',
    '../signal_reference/compare_spline_tensor_scipy.py',
    '../signal_reference/compare_spline_plan_transforms_scipy.py',
    '../signal_reference/compare_spline_plan_batch_mapping_scipy.py',
    '../signal_reference/compare_spline_plan_displacements_scipy.py',
    '../signal_reference/compare_spline_plan_gradients_scipy.py',
    '../signal_reference/compare_spline_grid_affine_gradients_scipy.py',
    '../signal_reference/compare_spline_hessians_scipy.py',
    '../signal_reference/compare_spline_parameter_hessians_scipy.py',
]
FILES=[
    'src/kinds.f90',
    'src/signal_spline.f90',
    'src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90',
    'src/ndimage_distance.f90',
    'src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90',
    'src/ndimage_watershed.f90',
    'src/ndimage.f90',
]
TESTS=[
    'test/test_ndimage_filters.f90',
    'test/test_ndimage_morphology.f90',
    'test/test_ndimage_extended.f90',
    'test/test_ndimage_footprints.f90',
    'test/test_ndimage_advanced.f90',
    'test/test_ndimage_regions_masks.f90',
    'test/test_ndimage_chamfer_regions.f90',
    'test/test_ndimage_callbacks_weighted.f90',
    'test/test_ndimage_mixed_fourier.f90',
    'test/test_ndimage_ellipsoid_order_line.f90',
    'test/test_ndimage_bruteforce_generic.f90',
    'test/test_ndimage_structure_axes.f90',
    'test/test_ndimage_axis_extrema.f90',
    'test/test_ndimage_selected_axes_deque.f90',
    'test/test_ndimage_order_axes_running.f90',
    'test/test_ndimage_gray_generic_axes.f90',
    'test/test_ndimage_complex_line_binary_axes.f90',
    'test/test_ndimage_value_indices.f90',
    'test/test_ndimage_complex_derivatives_label_lookup.f90',
    'test/test_ndimage_propagation_axes.f90',
    'test/test_ndimage_geometric_comprehension.f90',
    'test/test_ndimage_label_selected_axes.f90',
    'test/test_ndimage_vectorized.f90',
    'test/test_ndimage_vectorized_multi.f90',
    'test/test_ndimage_vectorized_mixed.f90',
    'test/test_signal_spline_tensor.f90',
    'test/test_signal_spline_plan_transforms.f90',
    'test/test_signal_spline_plan_batch_mapping.f90',
    'test/test_signal_spline_plan_displacements.f90',
    'test/test_signal_spline_plan_gradients.f90',
    'test/test_signal_spline_displacement_gradients.f90',
    'test/test_signal_spline_grid_affine_gradients.f90',
    'test/test_signal_spline_hessians.f90',
    'test/test_signal_spline_parameter_hessians.f90',
    'test/test_ndimage_spline_hessian_exports.f90',
]


def main():
    for script in SCRIPTS:
        print(f'Running {script}',flush=True)
        subprocess.run(['python',str(ROOT/'tools/ndimage_reference'/script)],check=True)
    for profile,flags in (
        ('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
        ('optimized',['-O1','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
    ):
        with tempfile.TemporaryDirectory(prefix=f'scifort-ndimage-{profile}-') as tmp:
            tmp=Path(tmp)
            cmd=['gfortran','-std=f2018','-Werror=implicit-interface',
                 '-J',str(tmp),'-I',str(tmp),*flags]
            objects=[]
            for path in FILES:
                obj=tmp/(Path(path).stem+'.o')
                subprocess.run([*cmd,'-c',str(ROOT/path),'-o',str(obj)],cwd=tmp,check=True)
                objects.append(str(obj))
            for test in TESTS:
                exe=tmp/Path(test).stem
                subprocess.run([*cmd,*objects,str(ROOT/test),'-o',str(exe)],cwd=tmp,check=True)
                subprocess.run([str(exe)],cwd=tmp,check=True)
        print(f'PASS: {profile} ndimage Fortran regressions',flush=True)
    print('PASS: all focused ndimage reference and Fortran tests',flush=True)


if __name__=='__main__':
    main()
