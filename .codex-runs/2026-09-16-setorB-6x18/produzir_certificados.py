"""Orquestra o pipeline inalterado, após conferir a proveniência do cache novo."""
import os, sys, json, time, hashlib, subprocess, shlex
from pathlib import Path
from fractions import Fraction
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
raiz = Path.cwd()
pasta = raiz / '.codex-runs/2026-09-16-setorB-6x18'
python = str(raiz / '.venv/bin/python')
cache = pasta / 'polinomio-novo.txt'
matriz = raiz / 'build/spectrum/rt-A-6x18.dat'
# Só a presença dos dois determinantes concluídos autoriza a comparação.
while not all(p.exists() for p in (cache,pasta/'bareiss-ponto-1.json',pasta/'bareiss-ponto-2.json')):
    time.sleep(20)
origem = json.loads((pasta/'origem-polinomio.json').read_text())
assert hashlib.sha256(matriz.read_bytes()).hexdigest() == origem['sha256']
coef = [int(x) for x in cache.read_text().split()]
assert len(coef)==334 and coef[-1]==1
conferencias = []
for indice in (1,2):
    dado = json.loads((pasta/f'bareiss-ponto-{indice}.json').read_text())
    assert dado['sha256']==origem['sha256']
    z = Fraction(dado['z'])
    y = z*2**72
    p = Fraction(0)
    for c in reversed(coef):
        p = p*y+c
    p *= Fraction(1,2**(72*333))
    d = Fraction(dado['determinante'])
    assert d==p, f'Divergência em {z}'
    conferencias.append({'z':str(z),'bareiss':str(d),'polinomio':str(p),'iguais':d==p})
    print(f'Bareiss Fraction × polinômio: z={z}, igualdade EXATA',flush=True)
with (pasta/'igualdades-exatas.json').open('x') as f:
    json.dump({'sha256_matriz':origem['sha256'],'sha256_cache':hashlib.sha256(cache.read_bytes()).hexdigest(),'pontos':conferencias},f,indent=2)
comandos = []
def executar(args, log):
    cmd = shlex.join(args)
    comandos.append(cmd)
    print(cmd,flush=True)
    with (pasta/log).open('x') as f:
        retorno = subprocess.run(args,stdout=f,stderr=subprocess.STDOUT).returncode
    print(f'{log}: código {retorno}',flush=True)
    return retorno
resultados=[]
for nome,parametros,h,v in [
    ('contour-A-6x18-reduced-B-sigma1_64',['--sigma0','1/64','--radius','1','--height','1/4','--imag-center','1/2'],640,333),
    ('contour-A-6x18-B-root-box',['--sigma0','1/32','--radius','1/10','--height','1/50','--imag-center','12/25'],128,97),
]:
    for tentativa in range(7):
        bruto=pasta/f'{nome}-tentativa-{tentativa}-nodes.txt'
        log=f'{nome}-tentativa-{tentativa}.log'
        args=[python,'scripts/build_contour.py',str(matriz),*parametros,'--horizontal',str(h),'--vertical',str(v),'--bits','64','--name',('setorB6x18Grande' if 'sigma' in nome else 'setorB6x18Raiz'),'--label','6x18','--cache',str(cache),'--out',str(bruto)]
        codigo=executar(args,log)
        if codigo==0:
            break
        assert codigo==2, (pasta/log).read_text()
        h*=2
        v=2*v+1
    else:
        raise RuntimeError('Não fechou após sete refinamentos; consultar logs.')
    conteudo=bruto.read_text()
    assert '\ntarget 1\n' in conteudo
    nos=Path('examples')/f'{nome}-nodes.txt'
    with nos.open('x') as f:
        f.write(conteudo)
    inc=Path('examples')/f'{nome}-increments.txt'
    assert not inc.exists()
    assert executar([python,'scripts/make_winding_certificate.py',str(nos),'--coarsen','--standalone','--quiet','--lean',str(pasta/f'{nome}.lean'),'--increments',str(inc)],f'{nome}-emissao.log')==0
    numero=sum(1 for linha in conteudo.splitlines() if linha.startswith('node '))
    arestas=sum(1 for linha in inc.read_text().splitlines() if linha.strip() and not linha.startswith('#'))
    resultados.append({'nome':nome,'nos':numero,'arestas_originais':numero,'arestas_coarsen':arestas,'none':0,'winding':1,'horizontal':h,'vertical':v,'log':log})
with (pasta/'comandos-pipeline.txt').open('x') as f:
    f.write('\n'.join(comandos)+'\n')
with (pasta/'resultados.json').open('x') as f:
    json.dump(resultados,f,indent=2)
codigo=executar(['bash','scripts/verify_reduced_b_6x18.sh'],'verificacao-final.log')
with (pasta/'estado-final.json').open('x') as f:
    json.dump({'codigo_verificacao':codigo,'segundos_totais':time.time()-origem['inicio']},f,indent=2)
assert codigo==0
print('CONCLUÍDO: ambos os certificados aceitos pelo kernel.',flush=True)
