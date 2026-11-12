#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact-band fixed-weight DPSS multitaper VAR(p) bootstrap coverage study.

The Fortran driver performs the entire finite-window spectral estimation and
canonical moving-block bootstrap from *raw simulated samples*. An independent
NumPy finite Fourier/kernel estimator and SciPy Lyapunov population reference
check both the estimator and the estimand. Coverage is descriptive only.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import subprocess

import numpy as np
from scipy.signal.windows import dpss

from study_canonical_varp_conditional import (
    BANDS, FREQUENCIES, LABELS, METHODS, STATISTICS, models,
    parse_output, simulate_windows, wilson,
)
from verify_varp_canonical_target import canonical
from verify_varp_taper_spectrum import reference


@dataclass(frozen=True)
class Setting:
    bandwidth: float
    tapers: int
    weighting: str

    @classmethod
    def parse(cls, text: str) -> 'Setting':
        """Read 'NW:K:equal' or 'NW:K:eigen'."""
        try:
            nw, k, method = text.split(':')
            out = cls(float(nw), int(k), method)
        except (ValueError, TypeError) as exc:
            raise ValueError(f'invalid taper setting {text!r}') from exc
        if (not np.isfinite(out.bandwidth) or out.bandwidth <= 0 or
                out.tapers < 1 or out.weighting not in ('equal', 'eigen')):
            raise ValueError(f'invalid taper setting {text!r}')
        return out

    def label(self) -> str:
        return f'NW={self.bandwidth:g}, K={self.tapers}, {self.weighting}'


DEFAULT_SETTINGS = (Setting(2.5, 1, 'equal'), Setting(2.5, 3, 'eigen'),
                    Setting(3.5, 5, 'equal'), Setting(3.5, 5, 'eigen'))


def taper_bank(length: int, setting: Setting) -> tuple[np.ndarray, np.ndarray]:
    """SciPy periodic DPSS tapers, with deterministic, frequency-fixed weights."""
    if setting.bandwidth >= length/2 or setting.tapers > length:
        raise ValueError('time-bandwidth product or taper count exceeds length')
    taper, ratios = dpss(length, NW=setting.bandwidth, Kmax=setting.tapers,
                         sym=False, norm=2, return_ratios=True)
    tapers = taper.T.copy()
    weights = (np.ones(setting.tapers) if setting.weighting == 'equal'
               else np.asarray(ratios).copy())
    if weights.min() < 0 or weights.sum() <= 0:
        raise AssertionError('invalid DPSS concentration ratios')
    return tapers, weights


def sample_spectra(windows: np.ndarray, tapers: np.ndarray,
                   weights: np.ndarray) -> np.ndarray:
    """Complex multichannel periodogram/band averages from exact sinc kernels."""
    nw, length, nc = windows.shape
    if tapers.shape[0] != length or tapers.shape[1] != len(weights):
        raise ValueError('incompatible windows and weights')
    if np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError('invalid weights')
    energy = np.sum(tapers**2,axis=0)
    if np.any((weights>0) & (energy<=0)):
        raise ValueError('zero-energy active taper')
    relative = (weights/weights.sum())/np.where(energy>0,energy,1)
    tapered = windows[:,None,:,:]*tapers.T[None,:,:,None]
    time = np.arange(length)
    f = FREQUENCIES[0]
    z = np.einsum('wkti,t->wki',tapered,np.exp(-2j*np.pi*f*time))
    point = np.einsum('k,wki,wkj->wij',relative,z,z.conj(),optimize=True)
    outputs = [point]
    lags = time[:,None]-time[None,:]
    for low,high in BANDS:
        width,midpoint = high-low,(high+low)/2
        kernel = np.sinc(width*lags)*np.exp(-2j*np.pi*midpoint*lags)
        exact = np.einsum('k,wkti,ts,wksj->wij',relative,tapered,kernel,
                          tapered.conj(),optimize=True)
        outputs.append(exact)
    result = np.stack(outputs)
    if result.shape != (1+len(BANDS),nw,nc,nc):
        raise AssertionError('unexpected spectral shape')
    return result


def population_spectra(model, tapers: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Independent time-pair expectation from SciPy Lyapunov covariance."""
    effective = weights/weights.sum()
    result = np.zeros((1+len(BANDS),model.a.shape[0],model.a.shape[0]),complex)
    for k, weight in enumerate(effective):
        if weight <= 0:
            continue
        _, frequency_spectra, bands = reference(model.a,model.sigma,tapers[:,k],
                                                FREQUENCIES,BANDS,normalized=True,
                                                average=True)
        result += weight*np.concatenate((frequency_spectra,bands),axis=0)
    return result


def encode_raw(windows: np.ndarray, tapers: np.ndarray, weights: np.ndarray,
               block: int, n_resamples: int, seed: int) -> list[str]:
    """Pass raw time samples (not Python spectra) to the native Fortran driver."""
    nw,length,nc = windows.shape
    data=windows.reshape(nw*length,nc)
    result=[f'{nw} {nc} {length} {tapers.shape[1]} {len(FREQUENCIES)} '
            f'{len(BANDS)} {block} {n_resamples} {seed}']
    for j in range(nc):
        result.append(' '.join(f'{x:.17g}' for x in data[:,j]))
    for k in range(tapers.shape[1]):
        result.append(' '.join(f'{x:.17g}' for x in tapers[:,k]))
    result.append(' '.join(f'{x:.17g}' for x in weights))
    result.append(' '.join(f'{x:.17g}' for x in FREQUENCIES))
    for low,high in BANDS:
        result.append(f'{low:.17g} {high:.17g}')
    return result


def conduct_study(executable: str | Path, *, trials: int = 2, nw: int = 40,
                  length: int = 16, resamples: int = 39,
                  blocks: tuple[int,...] = (2,4,8),
                  settings: tuple[Setting,...] = DEFAULT_SETTINGS,
                  seed: int = 20261101, report: bool = True,
                  csv_path: str | Path | None = None) -> dict:
    """Paired dataset study varying NW, taper count/weighting, and block length."""
    if trials < 1 or nw < 8 or length < 3 or resamples < 20:
        raise ValueError('invalid experiment dimensions')
    if not settings or not blocks or any(b < 1 or b >= nw for b in blocks):
        raise ValueError('invalid taper settings or block lengths')
    if len(set(blocks)) != len(blocks) or len(set(settings)) != len(settings):
        raise ValueError('repeated block length or taper setting')
    all_models=models()
    rng=np.random.default_rng(seed)
    # Hold every underlying dataset fixed across taper and block choices.
    observed_windows = [[simulate_windows(rng,model,nw,length)
                         for _ in range(trials)] for model in all_models]
    banks={setting:taper_bank(length,setting) for setting in settings}
    targets={(setting,i):canonical(population_spectra(model,*banks[setting]))[:,[1,2,0]]
             for setting in settings for i,model in enumerate(all_models)}
    requests=[]
    for setting in settings:
        taper,weights=banks[setting]
        for block in blocks:
            for m,model in enumerate(all_models):
                for t in range(trials):
                    requests.append((setting,block,m,t))
    lines=[str(len(requests))]
    predictions=[]
    for setting,block,m,t in requests:
        taper,weights=banks[setting]
        windows=observed_windows[m][t]
        lines.extend(encode_raw(windows,taper,weights,block,resamples,
                                seed=1000+m*trials+t))
        predictions.append(canonical(sample_spectra(windows,taper,weights).mean(axis=1))[:,[1,2,0]])
    process=subprocess.run([str(executable)],input='\n'.join(lines)+'\n',
                           capture_output=True,text=True,check=True)
    output=parse_output(process.stdout,len(requests))
    prediction=np.stack(predictions)
    max_error=float(np.max(abs(output[...,0]-prediction)))
    if max_error>5e-9:
        raise AssertionError(f'native Fortran/independent NumPy canonical disagreement: {max_error}')
    rows=[]
    if report:
        print(f'DPSS VAR(p) conditional coherence: {len(all_models)} models, '
              f'{trials} datasets/model, {len(settings)} settings, {len(blocks)} blocks, '
              f'{len(requests)} native Fortran bootstrap runs')
        print(f'Max independent finite-spectrum/canonical error {max_error:.3e}')
        print('Coverage is descriptive; near-tied roots lack reliable rank identification.')
    for setting in settings:
        for block in blocks:
            if report:
                print(f'  {setting.label()}, block={block}:')
            for m,model in enumerate(all_models):
                indices=[i for i,request in enumerate(requests)
                         if request[0]==setting and request[1]==block and request[2]==m]
                obs=output[indices]
                target=targets[setting,m]
                for f,frequency_label in enumerate(LABELS):
                    for s,statistic in enumerate(STATISTICS):
                        bias=float(np.mean(obs[:,f,s,0]-target[f,s]))
                        rmse=float(np.sqrt(np.mean((obs[:,f,s,0]-target[f,s])**2)))
                        for method,method_name in enumerate(METHODS):
                            lower,upper,flags=(obs[:,f,s,1+3*method+k] for k in range(3))
                            valid=(flags==1)&np.isfinite(lower)&np.isfinite(upper)
                            count=int(valid.sum())
                            hits=int(np.count_nonzero((lower[valid]<=target[f,s]) &
                                                      (target[f,s]<=upper[valid])))
                            wilson_bounds=wilson(hits,count)
                            width=float(np.mean(upper[valid]-lower[valid])) if count else float('nan')
                            rows.append(dict(model=model.name,
                                time_bandwidth=setting.bandwidth,taper_count=setting.tapers,
                                weighting=setting.weighting,block_length=block,
                                frequency=frequency_label,statistic=statistic,method=method_name,
                                population_target=float(target[f,s]),
                                population_root_gap=float(target[f,0]-target[f,1]),
                                datasets=trials,valid_intervals=count,covered=hits,
                                coverage=hits/count if count else float('nan'),
                                wilson_lower=wilson_bounds[0],wilson_upper=wilson_bounds[1],
                                mean_interval_width=width,full_estimator_bias=bias,
                                full_estimator_rmse=rmse))
                if report:
                    subset=[r for r in rows if r['model']==model.name and
                            r['time_bandwidth']==setting.bandwidth and
                            r['taper_count']==setting.tapers and
                            r['weighting']==setting.weighting and
                            r['block_length']==block and
                            r['frequency']==LABELS[0] and r['statistic']=='mode 1' and
                            r['method']=='studentized']
                    r=subset[-1]
                    print(f'    {model.name:25s} leading target={target[0,0]:.4f}, '
                          f'studentized hits={r["covered"]}/{r["valid_intervals"]}, '
                          f'mean bias={r["full_estimator_bias"]:+.4f}')
    if csv_path is not None:
        destination=Path(csv_path)
        destination.parent.mkdir(parents=True,exist_ok=True)
        with destination.open('w',newline='',encoding='utf-8') as file:
            writer=csv.DictWriter(file,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        if report:
            print(f'Wrote {len(rows)} coverage/bias records to {destination}')
    return {'rows':rows,'data':output,'targets':targets,'max_point_error':max_error,
            'requests':requests}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver',type=Path)
    parser.add_argument('--trials',type=int,default=2)
    parser.add_argument('--windows',type=int,default=40)
    parser.add_argument('--samples-per-window',type=int,default=16)
    parser.add_argument('--resamples',type=int,default=39)
    parser.add_argument('--block-lengths',default='2,4,8')
    parser.add_argument('--settings',default=','.join(f'{s.bandwidth}:{s.tapers}:{s.weighting}'
                                                    for s in DEFAULT_SETTINGS))
    parser.add_argument('--seed',type=int,default=20261101)
    parser.add_argument('--csv',type=Path)
    args=parser.parse_args()
    conduct_study(args.driver,trials=args.trials,nw=args.windows,
                  length=args.samples_per_window,resamples=args.resamples,
                  blocks=tuple(int(x) for x in args.block_lengths.split(',')),
                  settings=tuple(Setting.parse(x) for x in args.settings.split(',')),
                  seed=args.seed,csv_path=args.csv)
