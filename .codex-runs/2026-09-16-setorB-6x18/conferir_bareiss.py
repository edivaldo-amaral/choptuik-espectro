"""Bareiss independente em Fraction, diretamente sobre A + zI."""
import sys, time, hashlib, json
from pathlib import Path
from fractions import Fraction
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
pasta = Path('.codex-runs/2026-09-16-setorB-6x18')
origem = Path('build/spectrum/rt-A-6x18.dat')
sha = hashlib.sha256(origem.read_bytes()).hexdigest()
linhas = origem.read_text().splitlines()
n = int(linhas[1].split()[1])
matriz = [[Fraction(0) for _ in range(n)] for _ in range(n)]
for linha in linhas[linhas.index('BEGIN_ENTRIES')+1:]:
    campos = list(map(int, linha.split()))
    i,j,c = campos[:3]
    e = campos[3] if len(campos)==4 else 0
    matriz[i][j] = Fraction(c)*Fraction(2)**e
resultados = []
for z in (Fraction(1,32), Fraction(1,10)):
    inicio = time.time()
    a = [[matriz[i][j] + (z if i==j else 0) for j in range(n)] for i in range(n)]
    sinal = 1
    anterior = Fraction(1)
    print(f'Bareiss Fraction: início em z={z}', flush=True)
    for k in range(n-1):
        if a[k][k] == 0:
            candidato = next((i for i in range(k+1,n) if a[i][k]),None)
            if candidato is None:
                raise RuntimeError('Pivô nulo: determinante zero')
            a[k],a[candidato] = a[candidato],a[k]
            sinal = -sinal
        pivo = a[k][k]
        fonte = a[k]
        for i in range(k+1,n):
            alvo = a[i]
            fator = alvo[k]
            for j in range(k+1,n):
                alvo[j] = (pivo*alvo[j]-fator*fonte[j])/anterior
            alvo[k] = Fraction(0)
        anterior = pivo
        if (k+1)%10==0:
            print(f'z={z}: pivô {k+1}/{n-1}, {time.time()-inicio:.1f}s', flush=True)
    valor = sinal*a[-1][-1]
    resultados.append({'z':str(z), 'determinante':str(valor), 'segundos':time.time()-inicio})
    with (pasta / f'bareiss-ponto-{len(resultados)}.json').open('x') as f:
        json.dump({'sha256':sha, **resultados[-1]},f,indent=2)
    print(f'z={z}: determinante exato salvo, {time.time()-inicio:.1f}s', flush=True)
assert hashlib.sha256(origem.read_bytes()).hexdigest() == sha
print('Dois determinantes Bareiss Fraction concluídos.', flush=True)
