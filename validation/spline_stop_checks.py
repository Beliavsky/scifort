# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Residual-aware comparisons of convex derivative-jump smoothing fits."""
import numpy as np


def check_stopping_difference(design, rhs, actual, reference):
    """Bound coefficient differences attributable to different residual stops.

    For two exact minimizers of ||D*c||**2 at different residual constraints,
    convexity gives ||A*(c1-c2)||**2 <= abs(fp1-fp2). With full-column-rank A,
    ||c1-c2|| <= sqrt(abs(fp1-fp2))/sigma_min(A). This comparison supplements,
    rather than replaces, each checker's target and KKT checks. A roundoff
    budget for residual evaluation is included before taking the square root.
    """
    if actual.ndim == 1:
        actual = actual[:, None]
        reference = reference[:, None]
        rhs = rhs[:, None]
    singular = np.linalg.svd(design, compute_uv=False)
    assert singular[-1] > max(design.shape)*np.finfo(float).eps*singular[0]
    fp = np.sum((design@actual-rhs)**2)
    ref_fp = np.sum((design@reference-rhs)**2)
    scale = np.linalg.norm(rhs)+singular[0]*max(np.linalg.norm(actual), np.linalg.norm(reference))
    roundoff = 256*np.finfo(float).eps*scale**2
    bound = np.sqrt(abs(fp-ref_fp)+roundoff)/singular[-1]
    error = np.linalg.norm(actual-reference)
    assert error <= bound, ("residual stopping bound", error, bound, fp, ref_fp)
    return error/bound if bound else 0.0
