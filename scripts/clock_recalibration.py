#!/usr/bin/env python3
"""Rederiva os benchmarks de custo depois da correcao do EPP (17/09/2026).

Reexecuta o MESMO micro-benchmark de g_parametrix_contour_norms.py (item 4:
uma aplicacao matriz-livre completa de B, 18 convolucoes FFT, por grade) e
compara com os valores gravados em build/spectrum/g-parametrix-contour-norms.json,
que foram medidos com a maquina a ~1/3 do clock (EPP=power).

DIAGNOSTICO de custo; nao toca em nenhuma alegacao matematica.
Nao sobrescreve o JSON antigo: ele e a referencia "antes".
"""
from __future__ import annotations

import json
import pathlib
import subprocess
import time

import numpy as np
from scipy.signal import fftconvolve

import independent_column_diagnostic as D

ROOT = pathlib.Path(__file__).resolve().parent.parent
GRIDS = ((24, 72), (64, 200), (107, 384))
REPEATS = 5


def cpu_state():
    def read(p):
        try:
            return pathlib.Path(p).read_text().strip()
        except OSError:
            return None
    mhz = [float(l.split(':')[1]) for l in pathlib.Path('/proc/cpuinfo').read_text().splitlines()
           if l.startswith('cpu MHz')]
    return dict(
        epp=sorted({read(f'/sys/devices/system/cpu/cpu{i}/cpufreq/energy_performance_preference')
                    for i in range(len(mhz))}),
        governor=sorted({read(f'/sys/devices/system/cpu/cpu{i}/cpufreq/scaling_governor')
                         for i in range(len(mhz))}),
        mhz_idle=mhz,
    )


def bench():
    coefs = D.coef_arrays()
    kern = [a for arrs in coefs.values() for a in arrs if np.any(a)]
    out = {}
    for M, N in GRIDS:
        x = np.random.rand(2*M+1, 2*N+1)+1j*np.random.rand(2*M+1, 2*N+1)
        runs = []
        for _ in range(REPEATS):
            t = time.perf_counter()
            for k in kern:
                fftconvolve(k, x, mode='full')
            runs.append(time.perf_counter()-t)
        out[f'{M}x{N}'] = dict(best=min(runs), median=float(np.median(runs)), runs=runs)
    return out, len(kern)


def main():
    before = json.loads((ROOT/'build/spectrum/g-parametrix-contour-norms.json').read_text())
    old = before['benchmark_seconds_per_full_application_1core']
    state_idle = cpu_state()
    new, nkern = bench()
    state_load = cpu_state()
    speedup = {g: old[g]/new[g]['best'] for g in old}
    rep = dict(
        note='DIAGNOSTICO de custo apos corrigir CPU_ENERGY_PERF_POLICY_ON_AC=performance',
        date='2026-09-17',
        convolutions_per_application=nkern,
        cpu_before_description='EPP=power, 900-1000 MHz sob carga (medida original)',
        cpu_now=dict(idle=state_idle, after_bench=state_load),
        seconds_per_full_application_1core_before=old,
        seconds_per_full_application_1core_after={g: v['best'] for g, v in new.items()},
        raw_after=new,
        speedup=speedup,
        speedup_median=float(np.median(list(speedup.values()))),
    )
    p = ROOT/'build/spectrum/clock-recalibration-2026-09-17.json'
    p.write_text(json.dumps(rep, indent=2)+'\n')
    print(json.dumps({k: rep[k] for k in ('seconds_per_full_application_1core_before',
                                          'seconds_per_full_application_1core_after',
                                          'speedup', 'speedup_median')}, indent=2))
    print('escrito:', p)


if __name__ == '__main__':
    main()
