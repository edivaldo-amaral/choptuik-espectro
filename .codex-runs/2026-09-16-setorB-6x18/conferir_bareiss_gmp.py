"""Bareiss racional com Fraction e operações inteiras aceleradas por GMP."""
import sys,time,hashlib,json,math,ctypes
from pathlib import Path
from fractions import Fraction
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
pasta=Path('.codex-runs/2026-09-16-setorB-6x18')
lib=ctypes.CDLL(str((pasta/'bareiss_gmp.so').resolve()))
lib.det_bareiss.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_char_p)]
lib.det_bareiss.restype=ctypes.c_void_p
lib.liberar_determinante.argtypes=[ctypes.c_void_p]
def bareiss_racional(m):
    n=len(m)
    den=math.lcm(*(v.denominator for linha in m for v in linha))
    entradas=(ctypes.c_char_p*(n*n))(*(str(v.numerator*(den//v.denominator)).encode('ascii') for linha in m for v in linha))
    p=lib.det_bareiss(n,entradas)
    assert p, 'Falha de alocação'
    try: inteiro=int(ctypes.string_at(p))
    finally: lib.liberar_determinante(p)
    return Fraction(inteiro,den**n)
# Testes independentes: expansão por permutações, com troca de pivô e singularidade.
from itertools import permutations
for m in [ [[0,2,1],[3,4,5],[6,7,9]], [[1,2,3],[2,4,6],[3,7,8]], [[Fraction(1,3),Fraction(2,7)],[Fraction(-4,5),Fraction(3,11)]], [[7]] ]:
    m=[[Fraction(v) for v in linha] for linha in m]
    esperado=Fraction(0)
    for perm in permutations(range(len(m))):
        inv=sum(perm[i]>perm[j] for i in range(len(m)) for j in range(i+1,len(m)))
        esperado+=(-1)**inv*math.prod(m[i][perm[i]] for i in range(len(m)))
    assert bareiss_racional(m)==esperado
print('Testes Bareiss × expansão racional: OK',flush=True)
if '--testar' in sys.argv: sys.exit(0)
origem=Path('build/spectrum/rt-A-6x18.dat')
sha=hashlib.sha256(origem.read_bytes()).hexdigest()
linhas=origem.read_text().splitlines()
n=int(linhas[1].split()[1])
matriz=[[Fraction(0) for _ in range(n)] for _ in range(n)]
for linha in linhas[linhas.index('BEGIN_ENTRIES')+1:]:
    campos=list(map(int,linha.split())); i,j,c=campos[:3]
    e=campos[3] if len(campos)==4 else 0
    matriz[i][j]=Fraction(c)*Fraction(2)**e
for indice,z in enumerate((Fraction(1,32),Fraction(1,10)),start=1):
    inicio=time.time()
    print(f'Bareiss Fraction/GMP: z={z}',flush=True)
    racional=[[matriz[i][j]+(z if i==j else 0) for j in range(n)] for i in range(n)]
    valor=bareiss_racional(racional)
    with (pasta/f'bareiss-ponto-{indice}.json').open('x') as f:
        json.dump({'sha256':sha,'z':str(z),'determinante':str(valor),'segundos':time.time()-inicio,'metodo':'Fraction para leitura, eliminação de denominadores e resultado; Bareiss inteiro GMP com divisibilidade exata conferida em todas as operações'},f,indent=2)
    print(f'z={z}: determinante Fraction exato salvo, {time.time()-inicio:.1f}s',flush=True)
assert hashlib.sha256(origem.read_bytes()).hexdigest()==sha
print('Dois determinantes Bareiss exatos concluídos.',flush=True)
