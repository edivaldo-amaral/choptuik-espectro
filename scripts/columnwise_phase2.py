#!/usr/bin/env python3
"""Fase 2 de docs/COLUMNWISE_EXTERIOR.md: faixa radial coluna a coluna, EXATA.

Para cada coluna e=e_{d,p,m,n} com N_lo<=n<N_far e m<M (todas), calcula um
majorante racional EXATO
    U(e) = ||B_ref Q e|| + (s_max+e_B)||Q e|| + e_mu + resto  >=  ||(H_true(s)-I) e||,
uniforme em |s|<=s_max. Ingredientes:
  * coluna de Q (RT (22)) arredondada a 2^-P com raio em modulo propagado;
  * convolucao EXATA com os campos de 2 S Xi Gamma2(RefA,.) (ou Gamma2 sharp) por
    substituicao de Kronecker (produto de inteiros grandes), mais o termo mu R;
  * profundidade D: graus j<n-D majorados por grau com beta_glob;
  * n-D>=101: imagens em n desdobradas (tau=+ basta, linha N=0 com peso 1/2);
    m>=41: imagens em m desdobradas (sigma=+ basta) e p=1 identico a p=0.
Checkpoint: uma linha JSON por (sistema, n) em build/spectrum/columnwise-phase2/.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import json
import os
from pathlib import Path

import columnwise_exterior as CE

P_BITS = 64
ROOT = CE.ROOT
CKPT = ROOT/'build/spectrum/columnwise-phase2'


def column(sysname, d, m, n, depth, pbits=None):
    """mids inteiros (2^-P) e raio em modulo, j=n..n-depth (passo k)."""
    kind = CE.SYSTEMS[sysname]['kinds'][d]
    k = 1 if kind == 'J' else 2
    mu, b = CE.MU, Q(m, 2)
    pb = pbits or P_BITS
    grid = Q(1, 1 << pb)
    sc = 1 << pb

    def rnd(x):
        return (2*x.numerator*sc+x.denominator)//(2*x.denominator)

    def inv(re, im):
        den = re*re+im*im
        return re/den, -im/den
    out = {}
    re, im = inv(mu*(n+1), b)
    out[n] = (rnd(re), rnd(im), grid)
    if n-k < 0:
        return out, k
    a1, b1 = inv(mu*(n+1), b)
    a2, b2 = inv(mu*(n+1-k), b)
    pr, pi = a1*a2-b1*b2, a1*b2+b1*a2
    xr, xi = rnd(-2*mu*n*pr), rnd(-2*mu*n*pi)
    rad = 2*grid
    j = n-k
    while True:
        out[j] = (xr, xi, rad)
        if j-k < 0 or n-(j-k) > depth:
            break
        jj = j-k
        nr, ni = (mu*jj, -b) if k == 1 else (mu*(jj+1), -b)
        dr, di = inv(mu*(jj+1), b)
        fr, fi = nr*dr-ni*di, nr*di+ni*dr
        mr, mi = Q(xr, sc), Q(xi, sc)
        xr, xi = rnd(-(fr*mr-fi*mi)), rnd(-(fr*mi+fi*mr))
        rad += grid
        j = jj
    return out, k


_PACK = {}


def pack_signed(entries, nslots, nb):
    """entries: (slot, valor inteiro) -> inteiro sum v 2^(8 nb slot), via bytes (linear)."""
    pos = bytearray(nslots*nb)
    neg = bytearray(nslots*nb)
    for slot, v in entries:
        v = int(v)
        if v:
            buf = pos if v > 0 else neg
            buf[slot*nb:(slot+1)*nb] = abs(v).to_bytes(nb, 'little')
    return int.from_bytes(pos, 'little')-int.from_bytes(neg, 'little')


def packed(sysname, W, Bw):
    key = (sysname, W, Bw)
    if key in _PACK:
        return _PACK[key]
    tabs = CE.tables(sysname)
    nb = Bw//8
    res = {}
    for tkey, comps in tabs.items():
        lst = []
        for field in comps:
            slots = 81*W
            A = pack_signed([(((m1+40)*W+(n1+100)), a) for (m1, n1), (a, b) in field.items()], slots, nb)
            B = pack_signed([(((m1+40)*W+(n1+100)), b) for (m1, n1), (a, b) in field.items()], slots, nb)
            lst.append((A, B))
        res[tkey] = lst
    _PACK[key] = res
    return res


_FABS = {}
_RADF = {}


def rad_factor(sysname, d):
    """max sobre chaves de sum_i eta_i/eta_d field_abs_i + parte R (cota de coluna p/ o raio)."""
    if (sysname, d) not in _RADF:
        cfg = CE.SYSTEMS[sysname]
        r = CE.r_of(sysname)
        best = Q(0)
        for j in (101, 102):
            tkey = CE.key_for(sysname, d, j)
            val = sum((Q(cfg['eta'][i], cfg['eta'][d])*field_abs(sysname, tkey, i) for i in range(cfg['ncomp'])), Q(0))
            val += sum((Q(cfg['eta'][c], cfg['eta'][d])*abs(mult)*2*CE.MU*r/(1-r*r)
                        for c, mult in CE.r_terms(sysname, d, j)), Q(0))
            best = max(best, val)
        _RADF[sysname, d] = best
    return _RADF[sysname, d]


def field_abs(sysname, tkey, i):
    key = (sysname, tkey, i)
    if key not in _FABS:
        cfg = CE.SYSTEMS[sysname]
        (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
        field = CE.tables(sysname)[tkey][i]
        _FABS[key] = sum((abs(a)+abs(b))*Q(p1, q1)**abs(m1)*Q(p2, q2)**abs(n1)
                         for (m1, n1), (a, b) in field.items())/(1 << 72)
    return _FABS[key]


def unpack(P, S, nb):
    raw = P.to_bytes((S+1)*nb, 'little', signed=True)
    half, full = 1 << (8*nb-1), 1 << (8*nb)
    out = [0]*S
    carry = 0
    frm = int.from_bytes
    for s in range(S):
        v = frm(raw[s*nb:(s+1)*nb], 'little')+carry
        if v >= half:
            v -= full
            carry = 1
        else:
            carry = 0
        out[s] = v
    return out


def column_upper(args):
    sysname, d, p, m, n, depth = args
    cfg = CE.SYSTEMS[sysname]
    r = CE.r_of(sysname)
    col, k = column(sysname, d, m, n, depth)
    j0 = min(col)
    if j0-100 < 1:
        raise ValueError('faixa exige n-D>=101 (imagens em n desdobradas)')
    W = 201+depth+1
    maxcol = max(max(abs(x), abs(y)) for x, y, _ in col.values()).bit_length()
    Bw = 72+maxcol+(2*(depth+1)).bit_length()+4
    Bw = ((Bw+7)//8)*8
    nb = Bw//8
    pk = packed(sysname, W, Bw)
    sig = (1,) if m >= 41 or m == 0 else (1, -1)
    base_p = p if m < 41 else 0
    # colunas empacotadas por chave (paridade de j para J)
    ent = {}
    for j, (xr, xi, rad) in col.items():
        tkey = CE.key_for(sysname, d, j)
        vr, vi = (xr, xi) if base_p == 0 else (-xi, xr)
        er, ei = ent.setdefault(tkey, ([], []))
        er.append((j-j0, vr))
        ei.append((j-j0, vi))
    cols = {tkey: (pack_signed(er, depth+1, nb), pack_signed(ei, depth+1, nb))
            for tkey, (er, ei) in ent.items()}
    rows = 81+(2*m if len(sig) == 2 else 0)
    S = rows*W
    shift_plus = Bw*W*(2*m) if len(sig) == 2 else 0
    ncomp = cfg['ncomp']
    total = 0
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    assert q2 & (q2-1) == 0
    shq2 = q2.bit_length()-1
    E1 = m+40
    W_ = 201+depth+1
    E2 = j0-100+W_-1   # maior grau de saida possivel
    P1 = [p1**e*q1**(E1-e) for e in range(E1+1)]
    P2 = [p2**e*q2**(E2-e) for e in range(E2+1)]
    for i in range(ncomp):
        PR = 0
        PI = 0
        for tkey, (XR, XI) in cols.items():
            A, B = pk[tkey][i]
            if not A and not B:
                continue
            for sg in sig:
                xr_, xi_ = (XR, XI) if sg == 1 else (XR, -XI)
                pr = A*xr_-B*xi_
                pi = A*xi_+B*xr_
                if sg == 1 and shift_plus:
                    pr <<= shift_plus
                    pi <<= shift_plus
                PR += pr
                PI += pi
        re = unpack(PR, S, nb) if PR else [0]*S
        im = unpack(PI, S, nb) if PI else [0]*S
        # linha M: indice de linha q -> M = q-40-(m se dobrado)
        off = 40+(m if len(sig) == 2 else -m)
        # termo mu R (Gamma1), linhas M=sigma*m, graus N=j-1-2L>=0
        rterms = [(mult, j) for j in col for (c_, mult) in CE.r_terms(sysname, d, j) if c_ == i]
        rtail_contrib = 0
        if rterms:
            mult_by_parity = {}
            for mult, j in rterms:
                mult_by_parity[j % 2] = mult
            for sg in sig:
                # S(N) = c_{N+1} - S(N+2), so com j da coluna e multiplicador por paridade
                Sr = {}
                for N in range(n, -1, -1):
                    j = N+1
                    cr, ci = 0, 0
                    if j in col and (j % 2) in mult_by_parity:
                        xr, xi, _ = col[j]
                        vr, vi = (xr, xi) if base_p == 0 else (-xi, xr)
                        if sg == -1:
                            vi = -vi
                        mlt = mult_by_parity[j % 2]
                        cr, ci = mlt*vr, mlt*vi
                    nr, ni = Sr.get(N+2, (0, 0))
                    Sr[N] = (cr-nr, ci-ni)
                amp = CE.MU_INT72*2
                M = sg*m
                qrow = M+off
                for N, (sr, si) in Sr.items():
                    if not sr and not si:
                        continue
                    vr, vi = amp*sr, amp*si
                    t = N+100-j0
                    if 0 <= t < W and 0 <= qrow < rows:
                        re[qrow*W+t] += vr
                        im[qrow*W+t] += vi
                    else:
                        wN = P2[N]*(1 if N > 0 else Q(1, 2))
                        rtail_contrib += (abs(vr)+abs(vi))*P1[abs(M)]*wN
        acc = 0
        # linha q: sum_t v_t p2^(N) q2^(E2-N), N=t-100+j0, via Horner homogeneo em t
        N0 = j0-100
        lead = p2**N0*q2**(E2-N0-(W-1))
        for q in range(rows):
            base = q*W
            H = 0
            for i_ in range(W):
                t = W-1-i_
                v = abs(re[base+t])+abs(im[base+t])
                H = H*p2+(v << (shq2*i_) if v else 0)
            if H:
                acc += H*lead*P1[abs(q-off)]
        total += cfg['eta'][i]*(acc+rtail_contrib)
    sigcount = 1 if (m >= 41 or m == 0) else 2
    kappa_in = Q(p1**m*p2**n, q1**m*q2**n)
    bq = Q(total, q1**E1*q2**E2*(1 << (72+P_BITS)))/(sigcount*cfg['eta'][d]*kappa_in)
    res = assemble(sysname, d, m, n, col, k, bq, P_BITS)
    res['p'] = p
    return res


def assemble(sysname, d, m, n, col, k, bq, pb):
    """U = bq + raio + (s_max+e_B) massa(Q) + e_mu + resto de profundidade (EXATO)."""
    cfg = CE.SYSTEMS[sysname]
    r = CE.r_of(sysname)
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    shq2 = q2.bit_length()-1
    j0 = min(col)
    radmax = max(rad for _, _, rad in col.values())
    geo = 1/(1-r**k)
    H = 0
    for jj in range(n, j0-1, -1):
        v = abs(col[jj][0])+abs(col[jj][1]) if jj in col else 0
        H = H*p2+(v << (shq2*(n-jj)))
    qnorm = Q(H, p2**(n-j0)*(1 << pb))+2*radmax*geo
    radpart = 2*radmax*geo*rad_factor(sysname, d)
    remainder = Q(0)
    if j0-k >= 0:
        nu = Q(1) if k == 1 else Q(n, n-1)
        remainder = (cfg['beta_glob']+cfg['s_max'])*CE.SQ2*2*nu/(CE.MU*(n+1))*r**(n-j0+k)/(1-r**k)
    if sysname == 'A':
        e_B = CE.EPS_MU*Q(43390, 2061)+6*Q(4096, 243)*CE.EPS_OM
    else:
        e_B = CE.EPS_MU*Q(144, 17)+3*CE.EPS_OM
    e_mu = CE.EPS_MU*2*(1+r)/((1-r)*CE.MU)
    U = bq+radpart+(cfg['s_max']+e_B)*qnorm+e_mu+remainder
    return dict(d=d, m=m, n=n, U=str(U), U_decimal=float(U), BQ=float(bq), q=float(qnorm))


# ------------------------------------------------------------- motor C
KERNEL_SRC = Path(__file__).resolve().parent/'columnwise_kernel.c'
KERNEL_BIN = CKPT/'columnwise_kernel'
C_PBITS = 40


def kernel_binary():
    CKPT.mkdir(parents=True, exist_ok=True)
    if not KERNEL_BIN.exists() or KERNEL_BIN.stat().st_mtime < KERNEL_SRC.stat().st_mtime:
        import subprocess
        subprocess.run(['gcc', '-O2', '-fopenmp', '-o', str(KERNEL_BIN), str(KERNEL_SRC), '-lgmp'], check=True)
    return KERNEL_BIN


_HEADER = {}


def c_header(sysname):
    if sysname in _HEADER:
        return _HEADER[sysname]
    cfg = CE.SYSTEMS[sysname]
    tabs = CE.tables(sysname)
    order = list(tabs.keys())
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    lines = [f"SYS {cfg['ncomp']} {p1} {q1} {p2} {q2} {CE.MU_INT72}",
             'ETA '+' '.join(str(e) for e in cfg['eta']), f'NKEYS {len(order)}']
    for key in order:
        for field in tabs[key]:
            lines.append(f'COMP {len(field)}')
            lines.extend(f'{m1} {n1} {int(a)} {int(b)}' for (m1, n1), (a, b) in field.items())
    _HEADER[sysname] = ('\n'.join(lines)+'\n', {k: i for i, k in enumerate(order)})
    return _HEADER[sysname]


def run_row_c(sysname, jobs, threads):
    import subprocess
    header, kidx = c_header(sysname)
    cfg = CE.SYSTEMS[sysname]
    parts = [header, f'NJOBS {len(jobs)}\n']
    meta = []
    for q, (_, d, p, m, n, depth) in enumerate(jobs):
        col, k = column(sysname, d, m, n, depth, C_PBITS)
        j0 = min(col)
        full = 1 if j0-100 < 1 else 0          # n-D<101: imagens em n dobram
        if full:
            sig = 1 if m == 0 else 2
            base_p = p
        else:
            sig = 1 if (m >= 41 or m == 0) else 2
            base_p = p if m < 41 else 0
        ke = kidx[CE.key_for(sysname, d, 0)]
        ko = kidx[CE.key_for(sysname, d, 1)]
        rts = []
        for par in (0, 1):
            for comp, mult in CE.r_terms(sysname, d, 101+(1-par) if par == 0 else 101):
                pass
        rts = []
        for par, jrep in ((0, 102), (1, 101)):
            for comp, mult in CE.r_terms(sysname, d, jrep):
                rts.append((comp, par, mult))
        parts.append(f'JOB {q} {d} {m} {n} {j0} {sig} {full} {ke} {ko} {len(rts)} '
                     + ' '.join(f'{c} {pa} {mu}' for c, pa, mu in rts)+'\n')
        parts.append(f'NCOL {len(col)}\n')
        for j, (xr, xi, rad) in sorted(col.items()):
            vr, vi = (xr, xi) if base_p == 0 else (-xi, xr)
            parts.append(f'{j} {vr} {vi}\n')
        meta.append((d, p, m, n, col, k, sig, full))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    proc = subprocess.run([str(kernel_binary())], input=''.join(parts), capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr)
    results = []
    out = {}
    for line in proc.stdout.split('\n'):
        if line.strip():
            q, hx = line.split()
            out[int(q)] = int(hx, 16)
    (p1, q1), (p2, q2) = cfg['k1'], cfg['k2']
    for q, (d, p, m, n, col, k, sig, full) in enumerate(meta):
        E1, E2 = m+40, n+100
        kappa_in = Q(p1**m*p2**n, q1**m*q2**n)
        if full:
            n_img = (1 if m == 0 else 2)*(1 if n == 0 else 2)
            bq = Q(out[q], q1**E1*q2**E2*(1 << (72+C_PBITS)))/(n_img*cfg['eta'][d]*kappa_in)
        else:
            bq = Q(out[q], 2*q1**E1*q2**E2*(1 << (72+C_PBITS)))/(sig*cfg['eta'][d]*kappa_in)
        res = assemble(sysname, d, m, n, col, k, bq, C_PBITS)
        res['p'] = p
        results.append(res)
    return results


def jobs_for_row(sysname, n, M, depth, m_lo=0):
    cfg = CE.SYSTEMS[sysname]
    full = (n-depth) < 101
    out = []
    for d in range(cfg['ncomp']):
        if cfg['kinds'][d] == 'K' and n % 2 == 0:
            continue
        for m in range(m_lo, M):
            if m % 2 != cfg['mpar'][d]:
                continue
            ps = (0,) if m == 0 else ((0, 1) if (full or m < 41) else (0,))
            for p in ps:
                out.append((sysname, d, p, m, n, depth))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('system', choices=('A', 'sharp'))
    parser.add_argument('--n-lo', type=int, default=384)
    parser.add_argument('--n-hi', type=int, required=True, help='exclusivo (N_far)')
    parser.add_argument('--M', type=int, default=107)
    parser.add_argument('--m-lo', type=int, default=0)
    parser.add_argument('--depth', type=int, default=None)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--order', choices=('up', 'down'), default='up')
    parser.add_argument('--engine', choices=('python', 'c'), default='c')
    parser.add_argument('--rows-per-call', type=int, default=4)
    args = parser.parse_args()
    depth = args.depth or (100 if args.system == 'A' else 150)
    CKPT.mkdir(parents=True, exist_ok=True)
    ck = CKPT/(f'{args.system}-M{args.M}-D{depth}.jsonl' if args.m_lo == 0
               else f'{args.system}-m{args.m_lo}-M{args.M}-D{depth}.jsonl')
    done = set()
    if ck.exists():
        for line in ck.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)['n'])
    rows = range(args.n_lo, args.n_hi) if args.order == 'up' else range(args.n_hi-1, args.n_lo-1, -1)
    CE.tables(args.system)
    todo = [n for n in rows if n not in done]
    tag = 'C' if args.engine == 'c' else 'PY'

    def record(n, results):
        worst, worst_col = Q(0), None
        for res in results:
            u = Q(res['U'])
            if u > worst:
                worst, worst_col = u, res
        rec = dict(system=args.system, n=n, M=args.M, m_lo=args.m_lo, depth=depth, engine=tag, columns=len(results),
                   max_U=str(worst), max_U_decimal=float(worst),
                   argmax=dict(d=worst_col['d'], p=worst_col['p'], m=worst_col['m']),
                   all_lt_1=worst < 1)
        with ck.open('a') as fh:
            fh.write(json.dumps(rec)+'\n')
            fh.flush()
            os.fsync(fh.fileno())
        print(json.dumps({k: rec[k] for k in ('n', 'columns', 'max_U_decimal', 'argmax')}), flush=True)

    if args.engine == 'c':
        for start in range(0, len(todo), args.rows_per_call):
            chunk = todo[start:start+args.rows_per_call]
            jobs, owner = [], []
            for n in chunk:
                js = jobs_for_row(args.system, n, args.M, depth, args.m_lo)
                jobs.extend(js)
                owner.extend([n]*len(js))
            results = run_row_c(args.system, jobs, args.workers)
            for n in chunk:
                record(n, [r for r, o in zip(results, owner) if o == n])
        return
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for n in todo:
            jobs = jobs_for_row(args.system, n, args.M, depth, args.m_lo)
            record(n, list(ex.map(column_upper, jobs, chunksize=4)))


if __name__ == '__main__':
    main()
