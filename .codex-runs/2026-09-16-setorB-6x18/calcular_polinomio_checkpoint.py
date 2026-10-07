"""Execução 2: funções modulares existentes, checkpoints e CRT incremental."""
import os,sys,time,json,hashlib,resource,uuid
from pathlib import Path
from fractions import Fraction
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
raiz=Path.cwd()
pasta=raiz/'.codex-runs/2026-09-16-setorB-6x18'
sys.path.insert(0,str(raiz/'scripts'))
from build_contour import read_spectral_matrix,_coefficient_bound,_primes_above,_charpoly_mod
matriz_path=raiz/'build/spectrum/rt-A-6x18.dat'
sha=hashlib.sha256(matriz_path.read_bytes()).hexdigest()
sha_pipeline=hashlib.sha256((raiz/'scripts/build_contour.py').read_bytes()).hexdigest()
assert sha=='a3335d9550f9567db287c1ae105674a5914446f9d3772d3955e3628a747c5846'
origem_path=pasta/'origem-polinomio-execucao2.json'
if origem_path.exists():
    origem=json.loads(origem_path.read_text())
    assert origem['sha256']==sha and origem['sha256_pipeline']==sha_pipeline
else:
    origem={'matriz':str(matriz_path),'sha256':sha,'sha256_pipeline':sha_pipeline,'inicio':time.time(),'limite_memoria_virtual_bytes':1536*1024**2,'cache':'polinomio-novo.txt'}
    with origem_path.open('x') as f: json.dump(origem,f,indent=2)
n,matriz,expoente=read_spectral_matrix(matriz_path)
assert n==333 and expoente==-72
cota=_coefficient_bound(matriz,n)
necessario=2*cota+1
primos=[]; modulo_total=1; reservatorio=_primes_above(1<<30,8); indice=0
while modulo_total<=necessario:
    if indice>=len(reservatorio): reservatorio.extend(_primes_above(reservatorio[-1]+2,32))
    primos.append(reservatorio[indice]); modulo_total*=reservatorio[indice]; indice+=1
print(f'Matriz: n={n}; sha256={sha}; {len(primos)} primos; cota={cota.bit_length()} bits',flush=True)
checkpoints=pasta/'checkpoints-execucao2'
checkpoints.mkdir(exist_ok=True)
coef=[0]*(n+1); modulo=1; concluidos=0

def incorporar(primo,residuo):
    global modulo,concluidos
    assert primo==primos[concluidos] and len(residuo)==n+1
    assert residuo[-1]==1 and all(type(x) is int and 0<=x<primo for x in residuo)
    inverso=pow(modulo%primo,primo-2,primo)
    for grau,alvo in enumerate(residuo):
        coef[grau]+=modulo*((alvo-coef[grau])%primo*inverso%primo)
    modulo*=primo; concluidos+=1

# Apenas arquivos desta execução, completos e com hash confirmado, são retomados.
for marcador in sorted(checkpoints.glob('*.ok')):
    arquivo=marcador.with_suffix('.json')
    conteudo=arquivo.read_bytes()
    assert hashlib.sha256(conteudo).hexdigest()==marcador.read_text().strip()
    dados=json.loads(conteudo)
    assert dados['sha256_matriz']==sha and dados['sha256_pipeline']==sha_pipeline
    assert dados['inicio']==concluidos
    for primo,residuo in dados['residuos']: incorporar(primo,residuo)
if concluidos: print(f'Retomados {concluidos} resíduos completos desta execução.',flush=True)
inicio_atual=time.time(); lote=[]; inicio_lote=concluidos
for primo in primos[concluidos:]:
    residuo=_charpoly_mod(matriz,n,primo)
    incorporar(primo,residuo); lote.append([primo,residuo])
    if concluidos%50==0 or concluidos==len(primos):
        assert hashlib.sha256(matriz_path.read_bytes()).hexdigest()==sha
        dados={'sha256_matriz':sha,'sha256_pipeline':sha_pipeline,'inicio':inicio_lote,'residuos':lote}
        conteudo=json.dumps(dados,separators=(',',':')).encode()
        arquivo=checkpoints/f'lote-{inicio_lote:04d}-{concluidos:04d}-{uuid.uuid4().hex}.json'
        with arquivo.open('xb') as f:
            f.write(conteudo); f.flush(); os.fsync(f.fileno())
        with arquivo.with_suffix('.ok').open('x') as f:
            f.write(hashlib.sha256(conteudo).hexdigest()+'\n'); f.flush(); os.fsync(f.fileno())
        pico=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        print(f'Checkpoint {concluidos}/{len(primos)}; {time.time()-inicio_atual:.1f}s nesta sessão; pico RSS={pico:.1f} MiB',flush=True)
        lote=[]; inicio_lote=concluidos
assert modulo==modulo_total and modulo>necessario
coef=[c-modulo if c>modulo//2 else c for c in coef]
assert coef[-1]==1 and all(abs(c)<=cota for c in coef)
igualdades=[]
for indice in (1,2):
    dado=json.loads((pasta/f'bareiss-ponto-{indice}.json').read_text())
    assert dado['sha256']==sha
    z=Fraction(dado['z']); y=z*2**72; valor=Fraction(0)
    for c in reversed(coef): valor=valor*y+c
    valor*=Fraction(1,2**(72*n))
    direto=Fraction(dado['determinante'])
    assert valor==direto,f'Bareiss e polinômio divergem em {z}'
    igualdades.append({'z':str(z),'bareiss':str(direto),'polinomio':str(valor),'iguais':True})
    print(f'Bareiss salvo × polinômio novo em z={z}: igualdade racional EXATA',flush=True)
cache=pasta/'polinomio-novo.txt'
texto='\n'.join(map(str,coef))+'\n'
if cache.exists(): assert cache.read_text()==texto
else:
    with cache.open('x') as f: f.write(texto)
registro=pasta/f'conclusao-polinomio-execucao2-{uuid.uuid4().hex}.json'
with registro.open('x') as f:
    json.dump({'sha256_matriz':sha,'sha256_cache':hashlib.sha256(cache.read_bytes()).hexdigest(),'pontos':igualdades,'segundos_sessao':time.time()-inicio_atual,'fim':time.time(),'pico_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},f,indent=2)
print(f'Polinômio concluído e validado: {cache}',flush=True)
