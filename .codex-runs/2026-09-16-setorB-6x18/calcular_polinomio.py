"""Calcula do zero; não lê cache anterior."""
import sys, time, hashlib, json
from pathlib import Path
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
sys.path.insert(0, str(Path.cwd() / 'scripts'))
from build_contour import read_spectral_matrix, charpoly_exact
pasta = Path('.codex-runs/2026-09-16-setorB-6x18')
origem = Path('build/spectrum/rt-A-6x18.dat')
cache = pasta / 'polinomio-novo.txt'
assert not cache.exists()
inicio = time.time()
sha = hashlib.sha256(origem.read_bytes()).hexdigest()
with (pasta / 'origem-polinomio.json').open('x') as f:
    json.dump({'matriz': str(origem), 'sha256': sha, 'inicio': inicio}, f, indent=2)
n, matriz, expoente = read_spectral_matrix(origem)
print(f'Matriz n={n}, expoente={expoente}, sha256={sha}', flush=True)
coef = charpoly_exact(matriz, n)
assert hashlib.sha256(origem.read_bytes()).hexdigest() == sha
with cache.open('x') as f:
    f.write('\n'.join(map(str,coef)) + '\n')
print(f'Polinômio novo concluído em {time.time()-inicio:.1f}s', flush=True)
