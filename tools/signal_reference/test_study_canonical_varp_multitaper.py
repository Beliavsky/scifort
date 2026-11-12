#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent algebra, band-integral and DPSS study setup checks."""
from __future__ import annotations

import numpy as np

from study_canonical_varp_multitaper import (
    BANDS, FREQUENCIES, Setting, models, population_spectra,
    sample_spectra, simulate_windows, taper_bank,
)
from verify_varp_taper_spectrum import reference
from verify_varp_canonical_target import canonical


def run():
    rng=np.random.default_rng(20261102)
    for setting in (Setting(2.0,1,'equal'),Setting(2.5,3,'eigen'),
                    Setting(3.0,4,'equal'),Setting(3.0,4,'eigen')):
        tapers,weights=taper_bank(15,setting)
        if tapers.shape != (15,setting.tapers):
            raise AssertionError('DPSS shape mismatch')
        if not np.all((np.linalg.norm(tapers,axis=0) > 0.85) &
                      (np.linalg.norm(tapers,axis=0) <= 1.000001)):
            raise AssertionError('bad periodic DPSS energy')
        for model in models():
            windows=simulate_windows(rng,model,3,15)
            estimated=sample_spectra(windows,tapers,weights)
            if estimated.shape!=(3,3,5,5):
                raise AssertionError('wrong spectra dimensions')
            if not np.allclose(estimated,estimated.swapaxes(-2,-1).conj(),atol=2e-12):
                raise AssertionError('sample spectral matrix not Hermitian')
            if np.min(np.linalg.eigvalsh(estimated)) < -5e-12:
                raise AssertionError('sample spectral matrix not PSD')
            if not np.allclose(estimated,sample_spectra(windows,tapers,weights*1e4),
                               rtol=2e-13,atol=2e-13):
                raise AssertionError('relative weight scale changed estimate')
            if not np.allclose(estimated,sample_spectra(windows,tapers[:,::-1],weights[::-1]),
                               rtol=2e-13,atol=2e-13):
                raise AssertionError('taper ordering changed spectrum')
            if not np.allclose(estimated,sample_spectra(windows,-tapers,weights),
                               rtol=2e-13,atol=2e-13):
                raise AssertionError('taper sign changed spectrum')
            population=population_spectra(model,tapers,weights)
            roots=canonical(population)
            if not np.all((roots>=-1e-12)&(roots<=1+1e-12)):
                raise AssertionError('invalid canonical population targets')
            if not np.all(roots[:,1]>=roots[:,2]-1e-12):
                raise AssertionError('canonical roots out of order')
            if not np.allclose(1-(1-roots[:,1])*(1-roots[:,2]),
                               roots[:,0],atol=5e-13):
                raise AssertionError('aggregate canonical identity failed')
            if setting.tapers==1:
                _,p,b=reference(model.a,model.sigma,tapers[:,0],FREQUENCIES,BANDS,
                                normalized=True,average=True)
                if not np.allclose(population,np.concatenate((p,b)),atol=1e-12):
                    raise AssertionError('one-taper population identity')
            # Different construction: integrate the raw tapered Fourier
            # transforms at 96 Gauss-Legendre nodes (not the sinc kernel).
            nodes,quadrature=np.polynomial.legendre.leggauss(96)
            relative=weights/weights.sum()
            energy=np.sum(tapers**2,axis=0)
            time=np.arange(15)
            for band_idx,(lo,hi) in enumerate(BANDS):
                exact=np.zeros((3,5,5),complex)
                for z,w in zip(nodes,quadrature):
                    f=(hi-lo)*z/2+(hi+lo)/2
                    phase=np.exp(-2j*np.pi*f*time)
                    transformed=np.einsum('wti,tk,t->wki',windows,tapers,phase,
                                         optimize=True)
                    exact+=w/2*np.einsum('k,wki,wkj->wij',relative/energy,
                                         transformed,transformed.conj(),optimize=True)
                if not np.allclose(estimated[1+band_idx],exact,atol=5e-12,rtol=5e-12):
                    raise AssertionError(('exact DPSS band quadrature mismatch',setting,model.name))
    try:
        Setting.parse('2.5:3:adaptive')
    except ValueError:
        pass
    else:
        raise AssertionError('adaptive/data-dependent weights must not be accepted')
    try:
        taper_bank(16,Setting(9.0,2,'equal'))
    except ValueError:
        pass
    else:
        raise AssertionError('invalid DPSS NW accepted')
    print('PASS DPSS fixed-weight population, sample kernel, quadrature and invariances')


if __name__=='__main__':
    run()
