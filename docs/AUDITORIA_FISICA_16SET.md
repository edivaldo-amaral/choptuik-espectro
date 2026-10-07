# Auditoria física de 16/09/2026: constraints, gauge e setor B

Sub-agente FÍSICO. Escopo: a cadeia de constraints, gauge e setor B (objetivos
2 e 4), escrita pelo ChatGPT entre 14/09 16:40 e 16/09 13:27. Nenhum arquivo
existente foi editado. O ledger e a RETOMADA serviram só como mapa das alegações.

**Regra.** Tudo o que aparece como "numérico" abaixo é diagnóstico em ponto
flutuante. Os únicos resultados rigorosos desta auditoria são as identidades
algébricas/diferenciais verificadas simbolicamente (§0.2) e as contas racionais
exatas indicadas. Nenhuma contagem foi provada.

## 0. O que foi feito e como reproduzir

### 0.1 Arquivos criados

| arquivo | conteúdo | tempo |
|---|---|---|
| `scripts/independent_rt_symbolic.py` | verificação simbólica **universal** e exata de 12 identidades, transcritas direto do `choptuik.tex` de RT, sem usar os macros C nem `sharp_propagation.py`; também compara as tabelas dos macros com o TeX | 1,5 s |
| `scripts/test_independent_rt_symbolic.py` | regressão: exige que as identidades valham **e** que 5 mutações plausíveis falhem (prova de que o teste não é vácuo) | 2 s |
| `scripts/independent_numeric_audit.py` | reimplementação float64 própria dos operadores Fourier–Chebyshev de RT (a partir do TeX), com parsers próprios; subcomandos `validate`, `gauge`, `sharp`, `census`, `exact`, `identity` (e `m13`, exploratório e fora das conclusões) | 1–5 min por comando |
| scratchpad `fisico/{gauge,sharp,census,validate}.txt`, `axis.py`, `cosdigits.py`, `symbol.py` | saídas brutas | — |

```bash
.venv/bin/python scripts/independent_rt_symbolic.py                  # identidades universais
.venv/bin/python -m unittest scripts/test_independent_rt_symbolic.py # (a partir de scripts/)
export OPENBLAS_NUM_THREADS=2
.venv/bin/python scripts/independent_numeric_audit.py validate 4x12 --all-columns
.venv/bin/python scripts/independent_numeric_audit.py validate 6x18 12x36
.venv/bin/python scripts/independent_numeric_audit.py gauge 4x12 6x18 8x24 10x30 12x36 --decompose
.venv/bin/python scripts/independent_numeric_audit.py sharp 4x12 6x18 8x24 12x36
.venv/bin/python scripts/independent_numeric_audit.py sharp 6x18 8x24 12x36 --targets "0.0414+0.5j"
.venv/bin/python scripts/independent_numeric_audit.py census 4x12 6x18 8x24 10x30 12x36
.venv/bin/python scripts/independent_numeric_audit.py identity
.venv/bin/python scripts/independent_numeric_audit.py exact
make verify-reduced-b
```

### 0.2 Método da verificação simbólica (por que ela é uma prova das identidades)

Cada campo é um polinômio em `(tau, xi)` com coeficientes **simbólicos**; `s` e `mu`
são variáveis. Uma identidade de grau `d` nos campos, com derivadas de ordem `<= k`,
só depende dos k-jatos dos campos em `(tau0, +xi0)` e `(tau0, -xi0)`, porque a reflexão
`P` é o único termo não local. A base `tau^a xi^b` (com `b <= 5-2a` para k=2,
`b <= 3-2a` para k=1 e `b <= 1` para k=0; só `b` ímpar na componente ímpar) realiza
por interpolação de Hermite **qualquer** k-jato nesses dois pontos. As identidades
não têm `tau` explícito. Logo, anular o polinômio para todos os coeficientes prova
a identidade para todo campo C^k, todo `s` e todo `mu`. Cada grau é verificado com
a base de jatos adequada. As divisões por `xi` só são feitas depois de provar a
divisibilidade (senão o código levanta exceção). As mutações do teste confirmam que
o verificador detecta erros de coeficiente, de sinal, de fator (1 contra 2) e de
expoente (`s` genérico contra `s=mu`).

### 0.3 O que eu rederivei, recomputei, ou só li

- **Rederivado em papel:**
  - ação geométrica do gerador nulo (δζ, δW, δQ, δF, δφ e os quatro ω);
  - comutador `[D^sigma, e^{mu tau}Z] = e^{mu tau}D^sigma`;
  - classificação `s_n = mu(1-n)`;
  - álgebra `KC=ML ⇒ M1=C1, M0=C0+BC1-C1A, BC0=M0A, BD=DA` e `D h_j`;
  - lema de seção do quociente;
  - inversa de `T_s` e transporte `T_c(s)`, com ressonâncias e condição de cone;
  - redução do setor B e realidade;
  - cota 9/8 de `C1` (forte→fraco);
  - matriz 5×5 dos resíduos geométricos (determinante 128, exato).
- **Recomputado independentemente:**
  - as 12 identidades simbólicas;
  - as matrizes exportadas (selecionada, C0, C1, sharp);
  - todos os números de B1, B3 e B6;
  - transversalidade (em `Fraction`);
  - `make verify-reduced-b`.
- **Só li, sem verificar:**
  - `SHARP_DOMAIN.md` (fechos D_min=D_max e perdas de raio);
  - `RADIUS_RECOVERY.md`;
  - `SHARP_EXTERIOR_BOUNDS.md`, `SHARP_INTERMEDIATE_SHELL.md`, `ALGEBRAIC_TAIL.md` (exclusão em 414);
  - a fórmula RT `dfkhdjhsdshkfd` (Ricci em termos de Ω, que RT não publicam derivada);
  - o produtor de caixas `build_contour.py`.

## 1. Tabela de vereditos

| # | alegação | arquivo | veredito | evidência |
|---|---|---|---|---|
| B1.1 | `rt_sharp_matrix.c` exporta `K_N(0)=Π♯[O♯_mu + mu R_x Γ1♯ + Γ2♯(ω_RefA,·)]Π♯`, com fator 1 em Γ2♯ e fundo recortado em 2×caixa | `rt_sharp_matrix.c`, `SHARP_SPECTRAL_LOCATOR.md` | **CONFIRMADO** | reimplementação a partir do TeX (eq. 6383hi3h39z e 3i3hj4kh4) reproduz a matriz: erro 4,4e-16 (4×12, todas as colunas); 8,9e-16 (6×18 e 12×36, 40 colunas aleatórias) |
| B1.2 | a matriz selecionada é `Π L Π` (sem truncar intermediários) e `.constraints` é `Π♯(C0+sC1)Π` com `C1h=(0,-ξh2)` | patch, `analyze_sharp.py` | **CONFIRMADO** | erro 1,3e-15 / 2,7e-15 / 7,1e-15 (4×12, 6×18, 12×36). A variante "truncar antes de R_x" erra por O(1), o que discrimina as duas hipóteses. C0 até 1,8e-15; C1 até 2,2e-16 |
| B1.3 | `K(s)=K(0)+sI` é a linearização de Floquet da identidade sharp de RT | `SHARP_SPECTRAL_LOCATOR.md` | **CONFIRMADO** | §0.2: a identidade homogênea de RT vale off-shell (universal); `G(w,Ic)=Γ2♯(w,c)` (universal) |
| B1.4 | raiz sharp 0,40102473 contra selecionada 0,40102441 (3,2e-7); resíduo de propagação 0,0592 | `sharp-comparison-*.json` | **CONFIRMADO** (números) | recomputado: 3,20e-7; 5,921e-2; sequência do resíduo relativo 0,734 → 0,668 → 0,441 → 0,0592 |
| B1.5 | "reforça a violação de constraints; não prova" | idem | **CONFIRMADO, e a evidência é mais forte do que a nota diz** | ver §2: `C(s)h` fica colinear ao autovetor sharp com `1-|cos|` = 1,3e-4 (8×24) e 6,4e-8 (12×36), enquanto `‖Ch‖/‖h‖_2 = 0,324` fica estável |
| B1.6 | "há raízes sharp complexas com Re<1/8 e aliases" | idem | **CONFIRMADO COM CORREÇÃO (omissão grave)** | a raiz sharp `0,041419±0,5i` é também raiz de `L_A` (`0,041421+0,500008i` em 12×36), com a mesma assinatura de violação. É um **candidato do setor B com s_B≈0,0414 real**, fora do retângulo certificado. Ver §7 |
| B1.7 | testemunho `‖C_true(s)h‖ >= c_ref - ε_apl - c_op ε_v - (9/8)ε_s(‖v‖+ε_v)` | `SHARP_SPECTRAL_LOCATOR.md`, `constraint_witness.py` | **CONFIRMADO COM CORREÇÃO** | desigualdade triangular e cota 9/8 rederivadas (o máximo é atingido em n=0). **Falta** a hipótese de multiplicidade algébrica 1 da raiz de L; sem ela, excluir um vetor não exclui o espaço de raízes (é preciso `D|E` injetivo) |
| B2.1 | identidade off-shell completa `T_0(w)e(w)=0`, com `e=-PΩ^+` | `SHARP_LINEARIZED_IDENTITY.md` | **CONFIRMADO** | universal, graus 1–3 nos campos; as tabelas dos macros coincidem com o TeX (0 divergências) |
| B2.2 | reconstrução `e = Ic + RF` e divisibilidade das componentes selecionadas por ξ | idem | **CONFIRMADO** | papel + divisibilidade exigida no simbólico |
| B2.3 | `K(s)C(s)h = M(s)L(s)h - G(h,e(w))`, com as fórmulas explícitas de K e M | idem | **CONFIRMADO** | universal (graus 1–3, `s` e `mu` simbólicos); float com RefA: resíduo 2,8e-15, e 1,2e-4 **sem** o termo de correção. `M1=C1` confirmado |
| B2.4 | os testes de `test_sharp_propagation.py` não são tautológicos | testes | **CONFIRMADO, com ressalva** | testam uma identidade não trivial com tabelas RT, mas só em exemplos polinomiais. A verificação universal (§0.2) os supera |
| B2.5 | extensão aos domínios fechados com perda de raio | `SHARP_DOMAIN.md` | **NÃO VERIFICADO** | escopo analítico do Matemático; custo estimado ~1 dia de leitura crítica |
| B3.1 | `h=(1+mu^-1∂τ+ξ∂ξ)ω` é a ação da translação nula conjunta `u±→u±+ε` | `GAUGE_RT_ACTION.md` | **CONFIRMADO** | rederivado: métrica radial ∝ e^{-2ζ}(u_++u_-)^{-2}; δζ=e^{μτ}(Zζ-1); δW=e^{μτ}ZW; δQ; δF; δφ; `[D^σ,X]=e^{μτ}D^σ`; as quatro componentes dão (Z+1)ω_i |
| B3.2 | identidade off-shell `DΩ(w)[e^{μτ}(Z+1)w]=e^{μτ}(Z+1)Ω(w)`; `L(μ)(Z+1)w=(Z+2)F`, `C(μ)(Z+1)w=(Z+1)c` | idem | **CONFIRMADO** | universal; a mutação `s≠μ` falha |
| B3.3 | h ≠ 0 porque `(1+n+im/(2μ))a_mn=0` | idem | **CONFIRMADO** | papel |
| B3.4 | classificação `s_n=mu(1-n)` das reparametrizações nulas analíticas que preservam o centro; só n=0 é positiva e move o cone | `GAUGE_DOMAIN_CLASSIFICATION.md` | **CONFIRMADO** (dentro da classe declarada) | papel: `f(e^{-2πμ}u)=e^{2π(s-μ)}f(u)` força monômio único. Os testes do pushforward são não triviais |
| B3.5 | `|b(∂τω*)| > 117/1000` | idem | **CONFIRMADO** | `Fraction`: Im RefA4_(1,0)=138672388040959354547/2^69; ε=2^-25+2^-277 (L2+R de RT); peso 65/32; cota 0,1174600742 |
| B3.6 | projetar a fase em s≠0 não preserva o kernel | idem | **CONFIRMADO** (trivial) | `L(s)Πh=-s p b(h)/b(p)`. O teste `test_neutral_phase_projection_does_not_preserve_positive_kernel` é **vácuo**: compara duas tuplas literais e nunca aplica L |
| B3.7 | **consistência entre o modo gauge exato em s=μ e a raiz numérica estagnada a 2,21e-6** | `MODE_DISCRIMINATION.md` §10 vs `GAUGE_RT_ACTION.md` | **CONSISTENTE; a estagnação é truncamento** (diagnóstico) | §4: identidade exata `s_N-μ = -y^H r/(y^H v)`, com parcela do truncamento +2,2126e-6 e parcela do defeito do fundo +1,1e-10. `|cos|` entre (Z+1)ω_RefA e o autovetor vai a 1-O(1e-6). `|μ_true-μ_RefA|≈3,9e-11` não explica nada |
| B3.8 | o gerador está no domínio RT forte (via RADIUS_RECOVERY) | `RADIUS_RECOVERY.md` | **NÃO VERIFICADO** | escopo do Matemático |
| B4.1 | K(s) injetivo ⇒ ker L(s) ⊂ ker C(s); injetividade só no contorno não basta | `CONSTRAINT_GAUGE_EQUIVALENCE.md` | **CONFIRMADO** (inclusive a correção) | papel. Precisa de K(s) injetivo em **cada** raiz de L no interior, no espaço de raio menor |
| B4.2 | cadeias: `C(s)h(s)=O((s-s*)^k)` quando K é invertível numa vizinhança | idem | **CONFIRMADO** | papel |
| B4.3 | lema de seção `Z/imG ≅ Z∩kerB` | idem | **CONFIRMADO** | papel |
| B4.4 | `M1=C1`, `M0=C0+BC1-C1A`, `BC0=M0A`, `BD=DA`, `Dh_j=C(s0)h_j+C1h_{j-1}`, `dim E_c = dim E - rank(D|E)`, quociente `-1` | `CONSTRAINT_ROOT_QUOTIENT.md` | **CONFIRMADO** (álgebra) | papel + `M1=C1` universal. É condicional às hipóteses de domínio (B2.5, B3.8) |
| B4.5 | não há cadeia positiva extra dentro do gauge residual analítico | idem | **CONFIRMADO** | `(-μu∂u)^k u^n = (-μn)^k u^n`. Não exclui cadeias mistas (a nota diz isso) |
| B4.6 | testes de `test_root_constraints.py` e `test_constraint_witness.py` | testes | **CONFIRMADO, mas são brinquedos** | verificam álgebra linear 2×2 e aritmética do checker; nada sobre os operadores RT |
| B5.1 | compatibilidade `-2μξ(T_s b - ∂ξ a) = δΩ_j^- - δΩ_j^+`, com (i,j)=(2,2),(3,3),(4,4); identidade de W | `GEOMETRIC_FLOQUET_RECONSTRUCTION.md` | **CONFIRMADO** | universal (fundo e `s` simbólicos) |
| B5.2 | `T_s^{-1}` em funções 4π-periódicas, sem obstrução de período para Re s>0, nos dois setores | idem | **CONFIRMADO** | símbolo `im/2+s` só se anula em `s=-im/2` (Re s=0). No setor B, m ímpar muda a paridade de z, p e q, e nada mais. A cota `2/Re s` é válida (√2 bastaria) |
| B5.3 | 5 resíduos geométricos forçam os 5 resíduos RT restantes | idem | **CONFIRMADO COM RESSALVA** | reconstruí a matriz a partir de `dfkhdjhsdshkfd`: determinante 128. O teste só inverte uma matriz digitada à mão; `dfkhdjhsdshkfd` é de RT e não foi verificada (NÃO VERIFICADO; custo ~1 dia com CAS, que não há no .venv) |
| B5.4 | equivalência linearizada global **no ansatz** para Re s>0, com jets | idem | **CONFIRMADO** (com o sistema completo L e C) | papel. Não é fixação global de gauge (a nota diz isso) |
| B5.5 | transporte `T_c(s)`, fórmula (*), cotas, unicidade, ressonância `s=μ-im/2`, obstrução de cadeia, cone fixo ⇔ `f_-(τ,1)≡0` | `GLOBAL_NULL_GAUGE_TRANSPORT.md` | **CONFIRMADO** | papel (∂± em (τ,ξ), características, Taylor em ξ=c). A ressonância em s=μ é exatamente a translação nula de B3.1 |
| B6.1 | `M^{-1}L_B(s)M=L_A(s+i/2)`, Lema 5.1, modos reais de B em Im s=1/2 | `SECTOR_B.md` | **CONFIRMADO** | papel |
| B6.2 | H-dof, H-constr (seleção 0 de `IndexDOF`, C afim em s) | idem | **CONFIRMADO** | código RT lido; C1 validado numericamente (B1.2) |
| B6.3 | 4×12: 15 raízes com Re>0, 3 na faixa A, 0 na faixa B | idem §9 | **CONFIRMADO** | census |
| B6.4 | `make verify-reduced-b`: winding 0 sobre `1/8<=Re<=1, 1/4<=Im<=3/4` | Makefile, `examples/` | **CONFIRMADO** (só a regressão finita 4×12) | rodei: 1946→55 arestas, kernel OK, 19,6 s |
| B6.5 | C4 = "winding 0 no setor B" como alvo de `L_A` (Cor. 5.3), e "zero na faixa B, primeiro sinal a favor de C4" | `SECTOR_B.md` §5.3/§9, ledger C4 | **SOBREVENDIDO** (contradito numericamente) | de 6×18 a 12×36, `L_A` tem **uma** raiz convergente na faixa B, `0,063+0,480i → 0,04183+0,50037i → 0,04137+0,5001i → 0,041421+0,500008i`, com `‖Ch‖/‖h‖_2≈0,50` estável e colinearidade sharp `1-|cos|=7,5e-8`. O zero do 4×12 é artefato da malha grossa. O alvo correto é "0 **físico** em B", não "winding 0 de L_A em Re s>0" |
| B6.6 | não há gauge nulo analítico no setor B | (nova, minha) | **derivado em papel** | `f(e^{-2πμ}u)=-e^{2π(s-μ)}f(u)` dá Im s = ±1/2 mod 1, fora da faixa fundamental de B |

## 2. B1: o que é o operador sharp e o que a raiz 0,401 sustenta

**O operador.** RT (§2.3, identidade antes de `kdjskkkkksjsjhd`; §2.4, eq. `6383hi3h39z`)
mostram que `ω♯ = S♯Ω^+ = (Ω1 ímpar, Ω2♯^-)` satisfaz uma identidade homogênea
sempre que `SΩ^+=0`. Ela é o análogo 2D, no formalismo de frames, da identidade de
Bianchi contraída: é a lei de propagação das constraints. Linearizando em `ω♯=0`
com o fator de Floquet, obtém-se `K(s)=O♯_μ+sI+μR_xΓ1♯+Γ2♯(ω*,·)`. O exportador
calcula exatamente o truncamento de Galerkin desse operador com `ω*→RefA` recortado
em 2×caixa (B1.1). Verifiquei a identidade off-shell a partir do TeX, não a sua
derivação a partir de Bianchi. Para o uso aqui, isso é mais forte: é exatamente a
identidade que se consome.

**Números recomputados** (`sharp`, 12×36; norma l2 salvo indicação):

| raiz de L_A | ‖C(s)h‖/‖h‖ | raiz sharp mais próxima (distância) | 1-\|cos\|(Ch, autovetor sharp) | ‖K(s)Ch‖/‖Ch‖ (pond.) |
|---|---|---|---|---|
| 0,733179678 (físico) | 9,17e-5 | 0,401 (0,33) | 0,95 | 2,30 |
| **0,401024411** | **0,3243** | 0,401024731 (**3,2e-7**) | **6,4e-8** | **0,0592** |
| 0,168309292 (μ) | 1,47e-4 | 0,401 (0,23) | 0,95 | 2,21 |
| 4,6e-6 (fase) | 1,55e-4 | — | — | 2,33 |
| **0,041421+0,500008i** (faixa B) | **0,5026** | 0,041419+0,5i (**8,1e-6**) | **7,5e-8** | **0,0143** |

Em 8×24: `1-|cos|` = 1,3e-4 (0,401) e 1,8e-4 (B). Em 6×18: 0,0040 e 0,0068.

**O que isso sustenta (como diagnóstico).** Um modo que viola constraints precisa
satisfazer, no operador infinito, `L(s)h=0`, `C(s)h=c≠0` e `K(s)c=0`, pela identidade
de B2.3 no fundo exato. As três assinaturas aparecem, convergindo juntas:

- a raiz de L coincide com uma raiz de K;
- a constraint do autovetor de L é colinear, a 6e-8, ao autovetor de K;
- `‖Ch‖/‖h‖_2` não cai (0,323 → 0,320 → 0,325 → 0,324).

Os modos físico, gauge e de fase têm a assinatura oposta: constraint pequena, que
cai com o truncamento, e sem alinhamento sharp. A mesma assinatura aparece na raiz
da faixa B.

**O que NÃO sustenta.**

1. Nenhuma prova: é float, truncado e com fundo RefA.
2. Na norma ponderada (129/128, 9/8), que é a que o testemunho usaria,
   `‖Ch‖/‖h‖` ainda **cai**: 0,172 → 0,115 → 0,074 → 0,056. Sem limite certificado,
   não se pode fixar `c_ref`.
3. Nada exclui que a raiz infinita em 0,401 tenha espaço de raízes de dimensão >1
   contendo um vetor que satisfaz as constraints. O testemunho precisa de
   multiplicidade 1 certificada (B1.7).
4. A coincidência sugere que toda raiz sharp em Re>0 se levanta a uma raiz de L
   (vale para as duas raízes sharp da faixa fundamental). Isso é **observação**, não
   teorema: a fórmula "N_L = N_físico + N_sharp" continua sem prova de
   sobrejetividade (CONSTRAINT_ROOT_QUOTIENT tem razão em não usá-la).

## 3. B2: identidade linearizada

Tudo confirmado, universal e exato:

- a identidade homogênea completa de RT, off-shell;
- `Ω^+` Fourier–Chebyshev = `Ω^+` 2D (`f23ljf3lfj`);
- `PΩ^+=-Ω^-`;
- a simetria de Γ2 (necessária ao fator 2 do pencil);
- `G(w,Ic)=Γ2♯`;
- `K(s)C(s)h - M(s)L(s)h + G(h,e(w)) = 0` com as fórmulas **explícitas** de K e M da nota;
- `dC/ds=(0,-ξh2)`.

A derivação da linearização (regra do produto sobre uma identidade válida fora das
soluções; o fator de Floquet entra em todas as derivadas do resíduo linearizado)
está correta. O termo `-G(h,e(w))` é necessário. Numericamente com RefA, omiti-lo
deixa resíduo 1,2e-4 contra 2,8e-15 com ele.

## 4. B3: o modo gauge exato em s=μ explica a estagnação em 2,21e-6

A identidade B3.2 vale para **qualquer** fundo `w`. Com `w_c` = RefA recortado,
`v = Π(Z+1)w_c` e `y` o autovetor à esquerda da raiz `s_N` de `A_N=ΠLΠ`, tem-se
exatamente `s_N-μ = -y^H r/(y^H v)` com `r=(A_N+μ)v`. Como `A_N=ΠLΠ` (B1.2), ainda
`r = Π(Z+2)F_{w_c} - ΠL(1-Π)(Z+1)w_c` (defeito do fundo + truncamento). A
decomposição confere até 1,5e-14.

| caixa | s_N - μ_RefA | parcela do truncamento | parcela do defeito | \|cos\|((Z+1)ω, autovetor) | cond |
|---|---|---|---|---|---|
| 4×12 | +6,349e-2 | +6,361e-2 | -1,25e-4 | 0,9472 | 14 |
| 6×18 | -2,430e-2 | -2,437e-2 | +7,7e-5 | 0,99592 | 105 |
| 8×24 | -1,918e-3 | -1,915e-3 | -2,7e-6 | 0,999941 | 71 |
| 10×30 | +3,808e-6 | +3,824e-6 | -1,6e-8 | 0,999998 | 75 |
| 12×36 | **+2,213e-6** | **+2,2126e-6** | **+1,1e-10** | 1,000000 | 80 |

**Conclusão.** Não há contradição. Três explicações foram testadas:

- **μ_RT ≠ μ exato:** descartada. O pencil usa o próprio μ_RefA, e
  `|μ_true-μ_RefA| ≈ 3,9e-11`.
- **Erro do fundo ε_μ ou defeito de RefA:** descartado. Contribui 1e-10.
- **Truncamento:** explica tudo. A distância a μ é a projeção espectral do
  acoplamento `ΠL(1-Π)` aplicado ao **gerador exato**. Esse gerador tem coeficientes
  que decaem mais devagar que os do fundo, por causa das derivadas `∂τ` e `ξ∂ξ`.

O fator 504 entre 8×24 e 10×30 foi uma coincidência de mudança de sinal
(-1,9e-3 → +3,8e-6), e a previsão `<1e-7` extrapolou esse passo. O modo de fase se
comporta igual: 4,6e-6 em 12×36, com cond 602. A leitura 2 da §10 de
MODE_DISCRIMINATION ("o limite não é μ_RT") fica **desfavorecida**. O que ainda
não se prova: o valor do limite, a multiplicidade algébrica em μ, e a admissibilidade
do vetor no domínio físico (cone móvel).

**Leitura física (papel, não rigorosa).** A translação nula `u±→u±+ε` desloca o
ponto singular (0,0): é o representante, no gauge analítico de RT, da variação do
instante de acumulação. Seu expoente é `λ = 2πμ/K ≈ 0,6139`, e não `λ=1`. No centro
vale `dt_próprio = e^{-ζ}dτ` e `ζ(τ+2π)=ζ+K`, logo `|t| ∝ |T|^{K/(2πμ)}`, com
`T=(u_++u_-)/2`. O representante `λ=1` do gauge de tempo próprio central corresponde
a `f(u) ∝ |u|^{1-K/(2πμ)}`, **não analítico no cone**, e por isso não aparece no
pencil de RT. Isso corrige a expectativa de `SPECTRAL_PROBLEM.md` §4 e de
`MODE_DISCRIMINATION.md` §9.6.1 (`s=K/2π=0,274`).

## 5. B4: lemas condicionais

Álgebra correta, incluindo a correção "invertibilidade sharp só na fronteira não
elimina violações interiores". Dois pontos:

- **(a)** Nos dois pontos onde violações aparecem (0,401 e 0,0414+i/2), K é
  numericamente singular. O lema de injetividade não se aplica ali. O caminho
  obrigatório é o testemunho, com multiplicidade 1.
- **(b)** A linha gauge entra em `E_c` porque `Dg=C(μ)g=0` exatamente no fundo exato.
  Isso está confirmado universalmente (B3.2). A admissibilidade de g no domínio
  (RADIUS_RECOVERY) não foi verificada.

## 6. B5: reconstrução e transporte

Confirmado (B5.1–B5.5). Normalização de período: séries `e^{imτ/2}`, cilindro de
período 4π.

- `T_s=∂τ+s` só é singular em `s=-im/2`.
- `T_c(s)` só é singular em `s=μ-im/2`.
- No setor B, os m ímpares entram em h1..h3, e z, p, q herdam a paridade. Para
  Re s>0 nada muda na inversão de `T_s`.
- As ressonâncias do transporte em `Re s=μ` afetam os dois setores, conforme m par
  ou ímpar. Pela B6.6, em B elas caem em Im s = ±1/2 mod 1, fora da faixa
  fundamental de B.

## 7. B6: setor B, o que está coberto e o que falta

**Coberto hoje, e só no nível finito 4×12:** winding 0 do determinante da matriz
138×138 sobre o bordo de `R = [1/8, 1] × [1/4, 3/4]` (coordenadas A). Condicionado à
correção do produtor de caixas. **Nada** do operador infinito está coberto.

**Candidato novo.** Em todas as malhas de 6×18 a 12×36, `L_A` tem exatamente uma raiz
na faixa B com Re s > -0,18: `s ≈ 0,04142+0,50001i`. Em coordenadas do setor B, é um
modo **real** `s_B≈0,0414` (`λ≈0,151`), com assinatura de violação de constraints
(§2). Ele está em `0<Re s<1/8`, fora de R. Em 4×12, essa raiz ainda não tinha entrado
na faixa. Nenhuma raiz com Re s>1 aparece em nenhuma malha, nas duas faixas.

**Regiões ainda descobertas** para o fechamento de C4 (faixa `1/4 < Im s <= 3/4`,
coordenadas A):

1. **R e seu bordo, para o operador infinito:** falta todo o transporte
   finito→infinito (ε_N/C1, P1–P3), igual ao de Γ_A. A faixa B é o mesmo `L_A`, então
   os bounds de cauda precisam ser uniformes em `-1/4 <= Im s <= 3/4`, e não só
   em `|Im s| <= 1/4`.
2. **`0 < Re s < 1/8`:** totalmente descoberta, e **contém** o candidato
   0,0414+0,5i. O símbolo livre tem parte real `μ(n+1)+Re s > 0` até Re s > -μ, então
   não há degenerescência Fredholm no eixo. A faixa pode ser coberta por contorno,
   mas o winding de `L_A` ali é numericamente **1**, não 0.
3. **`Re s = 0`:** ausência de espectro sobre o eixo neutro na faixa B, ou tratamento
   explícito. Numericamente, a raiz mais próxima está em Re = +0,041. As demais estão
   em Re <= -0,18 (8×24), -0,215 (10×30) e -0,265 (12×36).
4. **`1 < Re s < 414`:** descoberta. A exclusão em `Re s >= 414` (ALGEBRAIC_TAIL) é
   condicional aos bounds RT e não foi auditada por mim.
5. **Bordos horizontais `Im s = 1/4` e `3/4`** fora de `[1/8, 1]`, isto é, em
   `(0, 1/8)` e em `(1, 414)`: é preciso certificar a ausência de espectro sobre eles.
   Pela simetria `s ↦ conj(s)+i`, os dois bordos têm o mesmo conteúdo espectral e
   basta tratar um.
6. **S3 no setor B:** obrigatório por causa do item 2. A contagem física 0 em B
   exige um testemunho de constraint (com multiplicidade 1) para a raiz
   0,0414+0,5i, ou uma exclusão equivalente. **S4 não é necessário em B** (B6.6).
7. **Simetria de conjugação:** usar a metade `1/4 < Im <= 1/2` exigiria tratar a reta
   `Im=1/2`, onde vivem os modos reais de B, **inclusive o candidato**. Não é
   recomendável; a nota já desaconselha pelo verificador.

## 8. Texto exato proposto para o ledger

**C4 (substituir a célula de estado/evidência):**
> C4 (setor B): regressão FINITA 4×12, winding 0 sobre ∂([1/8,1]×[1/4,3/4]), kernel OK
> (`make verify-reduced-b`); nada transferido ao operador infinito. Diagnóstico
> (AUDITORIA_FISICA_16SET §7): de 6×18 a 12×36, L_A tem uma raiz convergente na faixa
> B, s≈0,04142+0,50001i (s_B≈0,0414 real), com assinatura de violação de constraints
> (‖Ch‖/‖h‖≈0,50 estável, colinear ao autovetor sharp em 0,041419+0,5i a 7,5e-8). Logo
> o alvo "winding 0 de L_A na faixa B em Re s>0" é FALSO numericamente. O alvo correto
> é "multiplicidade FÍSICA 0 em B": winding 0 em Re s>=σ0 com σ0>0,0415, MAIS a faixa
> 0<Re s<σ0 (winding numérico 1) com testemunho S3 para essa raiz. Descobertos: 0<Re
> s<1/8, Re s=0, 1<Re s<414, bordos Im=1/4,3/4 fora de [1/8,1], e toda a transferência
> infinita.

**S3 (acrescentar):**
> Diagnóstico independente: o exportador sharp calcula Π♯K(0)Π♯ de RT
> (eq. 6383hi3h39z), reproduzido a 1e-15. Na raiz 0,401024 (A) e em 0,041421+0,500008i
> (B), C(s)h é colinear ao autovetor sharp (1-|cos|≈7e-8 em 12×36), com
> ‖Ch‖/‖h‖_2≈0,32 e ≈0,50 estáveis. Na norma (129/128,9/8) a razão da raiz 0,401 ainda
> cai (0,172→0,056). O testemunho exige enclosure certificado do autovetor E
> multiplicidade algébrica 1 da raiz.

**S4 (acrescentar):**
> A identidade DΩ(w)[e^{μτ}(Z+1)w]=e^{μτ}(Z+1)Ω(w) foi verificada universalmente.
> Diagnóstico: com w=RefA, s_N-μ_RefA=-y^H r/(y^H v) exatamente. Em 12×36 a parcela de
> truncamento é +2,2126e-6 e a do defeito do fundo +1,1e-10; |cos|((Z+1)ω,
> autovetor)→1. A estagnação em 2,21e-6 é truncamento, não evidência contra s=μ.
> Leitura física: s=μ (λ=2πμ/K≈0,614) é o representante no gauge analítico RT da
> translação do ponto singular; o λ=1 do tempo próprio central não é analítico no cone
> nessas coordenadas.

**SECTOR_B.md §9, leitura 2 (corrigir):**
> O zero na faixa da imagem de B em 4×12 é artefato da malha: de 6×18 em diante há uma
> raiz em 0,0414+0,5i, fora do retângulo 1/8<=Re<=1, com assinatura de violação de
> constraints.

**MODE_DISCRIMINATION §10 e README ("a identificação com mu_RT fica enfraquecida"):**
> Revisto em 16/09: a derivação do gerador (GAUGE_RT_ACTION) e a decomposição exata
> do deslocamento mostram que a distância 2,21e-6 é truncamento do gerador exato (a
> contribuição do defeito de RefA é 1e-10). A identificação numérica com s=μ fica
> fortalecida; a subtração continua dependendo de S4 (domínio, cone móvel,
> multiplicidade).

**constraint_witness.py / SHARP_SPECTRAL_LOCATOR (acrescentar hipótese):**
> O testemunho exclui só o vetor encerrado. Para excluir a raiz da contagem é preciso
> também certificar multiplicidade algébrica 1 (winding 1 num disco que contenha o
> enclosure de s), de modo que D|E seja injetivo.

**Testes vácuos ou fracos (registrar):**
> `test_uniform_domain_gauge.test_neutral_phase_projection_does_not_preserve_positive_kernel`
> não aplica operador algum (compara tuplas literais). `test_geometric_residual_combinations_are_invertible`
> inverte uma matriz digitada, sem derivá-la de `dfkhdjhsdshkfd`.
> `test_root_constraints` e `test_constraint_witness` testam só álgebra de brinquedo e
> aritmética. A verificação universal das identidades está em
> `independent_rt_symbolic.py`.

## 9. Próximo passo de maior alavancagem

### Objetivo 2 (equivalência gauge/constraints): testemunho certificado de violação na raiz 0,401

- **Por quê este:** das três raízes de Γ_A, é a única cuja exclusão não depende da
  decisão geométrica cone fixo/móvel. A margem numérica é grande (‖Ch‖/‖h‖_2≈0,32).
  A mesma máquina serve sem alteração à raiz B 0,0414+0,5i (objetivo 4).
- **Arquivos:**
  - `scripts/certify_constraint_witness.py` (novo);
  - relatório `build/spectrum/constraint-witness-0401.json`;
  - verificação escalar em Lean (Std, sem `sorry`/`axiom`/`native_decide`) em `formal/ChoptuikFormal/ConstraintWitness.lean`.
- **Enunciado:** existe um disco `|s-s0| <= ε_s` com winding 1 de L (operador
  infinito, na realização RT forte, incluindo cauda e erro da bola RT), e um vetor
  racional v com `‖h-v‖_{Y+} <= ε_v` para o autovetor h, tais que o lado direito
  `c_ref - ε_apl - c_op ε_v - (9/8)ε_s(‖v‖+ε_v)` é > 0, com `c_op` e `ε_apl` incluindo
  `C_true-C_RefA`.
- **Aceitação:**
  - `lower > 0` exato em `Fraction` e aceito pelo kernel;
  - todos os ε rastreáveis a relatórios com hash;
  - multiplicidade 1 do disco certificada pelo mesmo pipeline de tiles/cauda usado em
    Γ_A;
  - teste de mutação: com `ε_v` 10× maior, `lower` passa a ser <= 0.

### Objetivo 4 (fechamento analítico do setor B): reformular C4 e cobrir `0 < Re s < 1/8`

- **Por quê este:** o alvo atual está numericamente errado para o pencil, e a região
  descoberta contém um candidato.
- **Arquivos:**
  - `docs/SECTOR_B_TARGET.md`: enunciado corrigido, "multiplicidade física 0 =
    winding de L_A em Γ_B(σ0) menos raízes com testemunho S3";
  - contorno exato `examples/contour-A-6x18-reduced-B-sigma1_64-*.txt`, com
    `Γ_B' = ∂([1/64, 1] × [1/4, 3/4])` em 6×18. Dimensão 333; custo estimado ~40 min
    de Hessenberg modular, contra 166 s em 4×12;
  - um disco exato em torno da raiz B do 6×18 (0,0635+0,4795i).
- **Aceitação:**
  - kernel aceita winding **1** em Γ_B' e winding 1 no disco, com a mesma raiz no
    interior de ambos;
  - `make verify-reduced-b-sigma` reproduz;
  - o documento lista as regiões 3–5 da §7 como pendências com enunciado.
- **Depois:** o fechamento analítico exige (i) a transferência infinita com bounds
  uniformes em `-1/4 <= Im s <= 3/4` (obra do Matemático) e (ii) o testemunho do
  Objetivo 2 aplicado a essa raiz.

## 10. Não verificado, com custo

| item | por quê | custo estimado |
|---|---|---|
| `SHARP_DOMAIN.md` (D_min=D_max, perda de raio 8385/1440) | análise funcional, escopo do Matemático | 0,5–1 dia |
| `RADIUS_RECOVERY.md` (gerador no domínio forte, Fredholm índice 0) | idem | 1 dia |
| `dfkhdjhsdshkfd` de RT (Ricci em termos de Ω) | RT não publicam a conta; sem CAS no .venv | ~1 dia com motor de frames próprio |
| exclusão em Re s >= 414 | ALGEBRAIC_TAIL, escopo do Matemático | — |
| correção do produtor de caixas do winding B | `build_contour.py`, escopo do Matemático | — |
| completude do gauge físico global (toda perturbação → ansatz) | aberto nas próprias notas; não é uma verificação, é pesquisa | semanas |

## 11. Revisão cruzada do majorante exterior (rodada 3)

Alvo: `docs/COLUMNWISE_EXTERIOR.md`, Marcos 1–7. Os números são d*_A(G)=0,99820 e
d*_♯(G)=0,99903 em G=107×384, com margens de Schur da ordem de 10⁻³.

Minhas ferramentas não importam `scripts/columnwise_*.py` nem
`independent_component_beta.py`:

- `scripts/independent_columnwise_crosscheck.py`, com o motor TeX de
  `independent_numeric_audit.py`, validado a 1e-15 contra as matrizes RT
  exportadas (§1, B1.1–B1.2);
- a coluna `Q e` vem de uma resolução triangular própria de `O_μ c = e_n`, construída
  a partir da definição do operador e não da fórmula (22) de RT.

É diagnóstico em ponto flutuante (double) com margem declarada.

**Armadilha encontrada no meu próprio motor.** A convolução por FFT tem ruído
absoluto ~1e-17 por entrada. Os pesos `(5/4)^n` amplificam esse ruído até ~1e-6 nas
normas β do setor A; no sharp o efeito fica em ~1e-12. Troquei por soma direta sobre
as entradas não nulas da coluna. Com isso o ruído fica no nível relativo de ~1e-15.

### Vereditos

| item | veredito | evidência |
|---|---|---|
| (a1) β independe de j para j>=101 (por paridade de j), com cauda do termo μR | **CONFIRMADO** | papel + checagem em pontos que o Matemático não usou (ver abaixo) |
| (a1) β independe de m para m>=41; p=1 ≡ p=0 nesse regime | **CONFIRMADO** | igualdade a <=1,8e-15 em m=41/53/87 (A, d=3), 42/88 (A, d=2 e d=0; sharp, d=1), 42/64 (sharp, d=0) |
| (a1) desdobramento das imagens em n (n-D>=101, só τ=+, linha N=0 com peso 1/2) e em m (m>=41) | **CONFIRMADO** | papel. Empírico: meu U, sem desdobramento, com 4 imagens e coluna inteira, fica ABAIXO do U deles por 1e-6 a 1,2e-5 em todas as colunas da faixa, e coincide a 1e-17 nas colunas de n pequeno calculadas por eles com 4 imagens |
| (a1) valores absolutos de β_R | **CONFIRMADO** | max_m β(j=101/102) por componente. Meus: A 8,074493930 / 5,771748140 / 6,755589794 / 6,512887989; sharp 2,501013009 / 3,517539542. Deles (com cauda): 8,07449393 / 5,77174814 / 6,755589794 / 6,512887989; 2,501022795 / 3,517539542. A diferença no sharp d=0 (9,8e-6) é exatamente a cauda R, com meu R completo contra a janela + cota deles |
| (a2) U(e) recomputado em 12 colunas | **CONFIRMADO** | tabela abaixo; U_meu <= U_deles em todas |
| (a3) montagem d*(G) = max(faixa, d_R(606), d_F(107)) cobre todo o exterior | **CONFIRMADO** (montagem e contagem) | abaixo |
| lemas analíticos da Fase 1 (cotas de coluna de Q com √2, (1+√2)/2, 2/√3; resto de graus baixos com β_glob; constantes e_B, e_μ) | **NÃO VERIFICADO** (só amostrado) | minhas colunas além do corte estão muito abaixo de d_R(606) e d_F(107) (tabela), mas amostrar não prova um supremo. Custo de uma rederivação em papel: ~0,5 dia |
| amostra do critério 1 | **CORREÇÃO menor** | `columnwise-criteria.json` inclui duas colunas sharp inadmissíveis, `sharp d=1 p=1 m=7 n=700` e `sharp d=0 p=1 m=107 n=1`: m ímpar num espaço TP (m par). É inofensivo para d*, mas essas colunas não pertencem ao espaço |

### (a1) Invariância em j: diferenças medidas contra a cauda

Pontos: j=101, 137, 211 (ímpares) e 102, 138, 210 (pares).

| sistema, d | max \|β(j)-β(j')\| | cauda do termo R |
|---|---:|---:|
| A, d=1 (J, j par, sem termo R) | 0 | 0 |
| A, d=3 (J, j ímpar, mult 2) | 2,50e-10 | 3,04e-10 |
| A, d=0 (K) | 5,92e-10 | 7,2e-10 |
| sharp, d=1 | 0 | 0 |
| sharp, d=0 | 9,79e-6 | 1,09e-5 (margem ~10%) |

**Papel.**

- **Em j.** RefA tem suporte |m1|<=40 e |n1|<=100. Para j>=101 as imagens ±j da
  coluna produzem saídas em intervalos de N disjuntos, e nenhuma cai em N=0. Os pesos
  fatoram: κ2^{|±j+n1|} = κ2^j κ2^{±n1}. A razão só depende de j pela paridade, via P
  nas partes (h±Ph)/2 e na projeção (1-P)/2 da saída.
- **O termo μR_x Γ1.** É a soma alternada 2μ(-1)^L nas saídas N=j-1-2L. A janela
  L<=49 é idêntica para todo j>=101. As parcelas L>=50 somam no máximo
  2μ|mult|(η_c/η_d) r^{101}/(1-r²).
- **Em m.** Com m>=41 e |m1|<=40 não há interação entre as imagens ±m, nem saída em
  M=0.
- **p=1.** ‖i z‖_C = ‖z‖_C, e as imagens não interagem.

### (a2) U(e) recomputado (float, com a coluna Q e inteira)

Fórmula: `U_meu = ‖B_ref Qe‖ + (s_max+e_B)‖Qe‖ + e_μ`. As constantes e_B e e_μ são
as do documento deles (não rederivadas); contribuem <1e-5.

| sistema d p m n | U_meu | U_deles | U_deles - U_meu |
|---|---:|---:|---:|
| A 0 1 18 385 (**argmax A**) | 0,5049280066 | 0,5049292133 | +1,2e-6 |
| sharp 1 0 104 385 (**argmax sharp**) | 0,4629923907 | 0,4630044971 | +1,2e-5 |
| A 3 0 105 384 | 0,4565857303 | 0,4565883388 | +2,6e-6 |
| A 1 0 0 385 | 0,2871047009 | 0,2871076387 | +2,9e-6 |
| A 2 1 40 500 | 0,2903566344 | 0,2903588710 | +2,2e-6 |
| sharp 0 1 4 605 | 0,2246587705 | 0,2246625907 | +3,8e-6 |
| A 3 0 57 700 (n>=606) | 0,2220132732 | 0,2220148323 | +1,6e-6 |
| A 2 0 108 384 (m>=107) | 0,4211643177 | 0,4211669636 | +2,6e-6 |
| A 1 1 108 1 (m>=107, n pequeno) | 0,05780609978277352 | 0,05780609978277352 | 0 (<1e-17) |
| sharp 1 0 108 2 (m>=107) | 0,07126240410476363 | 0,07126240410476364 | ~1e-17 |
| A 1 1 18 606 (logo além de N_far) | 0,20785 | — (coberta por d_R(606)=0,96656) | — |
| sharp 1 1 4 606 (logo além de N_far) | 0,21828 | — (coberta por d_R(606)=0,99903) | — |

Margem do float: todas as somas são de termos positivos com erro relativo <=1e-13.
Nas colunas da faixa, a diferença (>=1,2e-6) é 10⁷ vezes maior que o erro. Nas duas
colunas de n pequeno a concordância está no nível do double; ali a desigualdade é
decidida pelo cálculo EXATO deles, e o meu só o confirma.

### (a3) Cobertura

O exterior de G={m<=106, n<=383} é a união disjunta de três partes:

- {m>=107} → d_F(107), qualquer n;
- {m<=106, 384<=n<=605} → faixa;
- {m<=106, n>=606} → d_R(606), qualquer m.

Não há buraco nas fronteiras: m=106 está na faixa, m=107 em d_F; n=605 está na faixa,
n=606 em d_R.

Conferi nos checkpoints (leitura direta do JSONL, sem o código deles):

- **A-M107-D60:** 222 linhas contíguas n=384..605, sem duplicatas, M=107,
  profundidade 60 (n-D>=324>=101), `all_lt_1` em todas, máximo 0,504929.
  Colunas: 57276 = 111·295 (n ímpar) + 111·221 (n par).
  - n ímpar: K com m par <=106 (m=0 com p=0; m=2..40 com p=0,1; m=42..106 com p=0) =
    74; o mesmo para d=1 e d=2 (74+74); d=3 com m ímpar (1..39 com p=0,1; 41..105 com
    p=0) = 73. Total 295.
  - n par: sem K (XiOdd), 221.
- **sharp-M107-D100:** 222 linhas, 24642 = 111·148 + 111·74, máximo 0,463004.

A omissão de p=1 só ocorre em m=0 (não existe) e em m>=41 (justificada por (a1)).
s_max=5/4 (A) cobre a faixa A e o retângulo B (|1+0,75i|=5/4). s_max=4/5 (sharp)
cobre o disco |s-0,733|<=1/32.

**O que esta revisão não cobre.**

- A validade dos supremos analíticos d_R e d_F; como são eles que dominam d*, essa
  é a parte com margem de 10⁻³ que continua sem segunda verificação independente.
- Os escalares do Schur (0,998763 e 0,999033).

Nenhum erro foi encontrado; por isso a Parte B foi executada
(`docs/CONSTRAINT_WITNESS_NORMS.md`).

## 12. Verificação das cotas analíticas e do modo completo (rodada 5b)

Fecha o elo que a §11 deixou como **NÃO VERIFICADO**: as cotas d_R(N) e d_F(M) da
Fase 1, que são o que fixa d*(G). Scripts novos:
`scripts/independent_phase1_bounds.py` (não importa `columnwise_*`; os β vêm do meu
motor) e as novas colunas de `scripts/independent_columnwise_crosscheck.py`.

### 12.1 Rederivação em papel dos três lemas de coluna

Notação: coluna `e` no grau n, modo m; `b=m/2`, `x=μ(n+1)`, `D_j=ib+μ(j+1)`,
`r=1/κ2`, `k=1` (J) ou `2` (K), `R_k=r^k/(1-r^k)`. Da fórmula (22) de RT:
`c_n=1/D_n`, `c_{n-k}=-2μn/(D_n D_{n-k})`, `c_j=-f_j c_{j+k}` com `|f_j|<=1`
(em K vale `|f_j|=1`). Logo `|c_{n-ik}| <= 2μn/(|D_n||D_{n-k}|)` para **todo** i>=1.

- **L1.** `‖c_n‖_C = (x+b)/(x²+b²)`. Maximizando em x com b fixo: máximo em
  `x=b(√2-1)`, valor `(1+√2)/(2b)`. **Confirma** `‖c_n‖_C <= G0/b`, `G0=(1+√2)/2`.
- **L2.** `|c_{n-ik}| <= 2y/√(((y+μ)²+b²)((y-δ)²+b²))` com `y=μn`. Em k=1
  (δ=0, com (y+μ)²>=y²) dá `<= 2y/(y²+b²) <= 1/b`. Em k=2 (δ=μ) o denominador é
  `√((y²+μ²+b²)²-4y²μ²)`; a derivada zera em `y²=μ²+b²` e o valor é exatamente
  `1/b`. Portanto `|c_{n-ik}| <= 1/b` nos dois casos e `‖c_{n-ik}‖_C <= √2/b`.
  O autor usa `ν_F=2/√3` em K: é **conservador** (vale, com folga de 15%).
- **L3 (a peça central).** `(n+1)‖Qe‖ <= (Cx²+xb)/(μ(x²+b²))` com
  `C=1+2√2 ν R_k`. Prova: o termo j=n dá `x(x+b)/(μ(x²+b²))` por L1; para i>=1,
  `(x/μ)·√2·2y R_k/(|D_n||D_{n-k}|)`, e
  `y√(x²+b²)/(x√(z²+b²)) <= y/z` (pois `z²b²<=x²b²`), com `z=x-kμ`. Em k=1,
  `y/z=1`, logo `ν=1`; em k=2, `y/z=n/(n-1)<=N/(N-1)`, logo `ν=N/(N-1)`.
  Maximizando `(C+u)/(1+u²)` em `u=b/x`: máximo `(C+√(C²+1))/2 <= C+1/(4C)`.
  **Confirma** `S(N)=(C+1/(4C))/μ` e as duas definições de C do autor.
- **Cauda de graus baixos.** De L3, `‖c_j‖_C <= √2·2ν/(μ(n+1))`, e
  `Σ_{j<=100} r^{n-j} <= r^{N-100}/(1-r^k)`: **confirma** o termo `low(N)`.

### 12.2 Recomputação independente dos quatro valores

Com os **meus** β (motor próprio) e os lemas acima; β_glob, e_B e e_μ são
constantes externas do Matemático, não auditadas por mim.

| valor | autor | meu | diferença |
|---|---:|---:|---:|
| d_R(900), A | 0,65116 | 0,651134 | -2,6e-5 |
| d_R(900), sharp | 0,67304 | 0,673007 | -3,3e-5 |
| d_F(160), A | 0,66755 | 0,667519 | -3,1e-5 |
| d_F(160), sharp | 0,66784 | 0,667808 | -3,2e-5 |
| d_R(606), A / sharp | 0,96656 / 0,99903 | 0,966511 / 0,998978 | -4,5e-5 / -4,9e-5 |
| d_F(107), A / sharp | 0,99820 / 0,99864 | 0,998160 / 0,998591 | -4,2e-5 / -4,6e-5 |
| d_R(384), A / sharp | 1,52389 / 1,57509 | 1,523824 / 1,575012 | -7e-5 / -7,7e-5 |

**Veredito: CONFIRMADO.** Meus valores ficam sistematicamente ~5e-5 ABAIXO dos dele,
com diferença relativa constante de 5,1e-5, que é exatamente a folga dos majorantes
racionais que ele usa: `99/70 >= √2` (razão 1+5,1e-5), `1207107/1000000 >= (1+√2)/2`
e `1155/1000 >= 2/√3`. Os valores dele são cotas superiores válidas da mesma
expressão. As componentes que mandam são d=2 (J) no setor A e d=1 (J) no sharp;
por isso a escolha `ν_K=2/√3` não afeta nenhum dos quatro números (com ν_K=1 os
valores não mudam).

**Forma funcional.**

- `d_F(M) ∝ 1/M`: **CONFIRMADO COM CORREÇÃO.** O exato é
  `d_F(M) = Φ/(M/2) + e_μ`, ou seja `(d_F(M)-e_μ)·M` é constante — verifiquei
  `(d_F(160)-e_μ)·160 = (d_F(107)-e_μ)·107` com diferença relativa <=1,3e-16 nos dois
  sistemas. Há o termo aditivo `e_μ` (4,2e-9 em A, 7,9e-9 no sharp), pequeno mas não
  nulo: a afirmação "d_F é exatamente Φ/b+e_μ, logo ∝1/M" do Marco 8 é correta na
  primeira metade e imprecisa na segunda.
- `d_R(N) ≈ 1/N`: **CONFIRMADO.** 0,966511/0,651134 = 1,4844 contra 900/606 = 1,4851;
  o desvio vem de ν=N/(N-1) e do termo `low`.

**Os cortes se sustentam:** as duas caudas ficam em 0,6511 e 0,6730 (A) e 0,6675 e
0,6678 (sharp), abaixo da meta 0,7. Amostrei ainda colunas reais nessas regiões com o
meu motor: A (d=2,m=0,n=900) 0,1411; A (d=2,m=40,n=901) 0,1530; sharp (d=1,m=0,n=900)
0,1296; A (d=1,m=160,n=5) 0,0611; A (d=1,m=200,n=50) 0,0869; sharp (d=1,m=160,n=3)
0,0462 — todas 4 a 14 vezes abaixo das cotas, como esperado de um supremo analítico.

### 12.3 Modo "completo" do núcleo C (região R3)

Recomputei U(e) sem `columnwise_*`, com a coluna inteira e as quatro imagens.

| coluna | autor (modo completo) | meu | diferença |
|---|---:|---:|---:|
| sharp d=1 p=0 m=108 n=0 | 0,0903737368185 | 0,0903737367282 | +9,0e-11 |
| sharp d=1 p=0 m=108 n=1 | 0,0498376271134 | 0,0497962017120 | +4,1e-5 |
| sharp d=1 p=0 m=108 n=2 | 0,0713367185678 | 0,0712624041048 | +7,4e-5 |
| sharp d=1 p=0 m=108 n=3 | 0,0689835063946 | 0,0688829524716 | +1,0e-4 |
| sharp d=1 p=0 m=108 n=25 | 0,1205106101863 | 0,1204323271524 | +7,8e-5 |
| A d=3 p=0 m=107 n=3 (via Python 4 imagens) | 0,1240099 | 0,1240099407614 | 0 |
| A d=0 p=1 m=150 n=5 (idem) | 0,1095716 | 0,1095715834938 | 0 |
| A d=1 p=0 m=110 n=25 | — (ainda não calculada) | 0,1247448142876 | — |
| sharp d=1 p=0 m=108 n=220 (modo desdobrado) | — (linha ainda não calculada) | 0,3719449938290 | — |

**Veredito: CONFIRMADO COM RESSALVA.** Em todas as colunas comparáveis o valor dele é
**maior ou igual** ao meu, como deve ser num majorante. Mas há um viés sistemático:

- em n=0 (onde as duas imagens em n coincidem) a concordância é de 9e-11, o nível do
  arredondamento 2^-40 do núcleo C;
- em n>=1 o modo completo fica 4e-5 a 1e-4 acima (8e-4 a 1,5e-3 em termos relativos),
  ordem de grandeza **cinco a seis** vezes maior que o arredondamento;
- o caminho Python de 4 imagens, usado no critério 1 da rodada 2, coincide comigo até
  o último dígito.

Não ameaça nada (os U de R3 ficam em ~0,12, contra a meta 0,7), e a direção do viés é
segura. Mas a validação declarada do modo completo contra o Python não deve ser lida
como concordância ao nível do arredondamento: alguma perda de cancelamento entre as
imagens em n (que se sobrepõem quando n<101) explicaria o sinal e a ordem do desvio.
Sugiro ao autor comparar, numa coluna com n>=1, as parcelas `bq` das duas
implementações.

**Equivalência p=1 ≡ p=0 no regime de n pequeno:** verificada diretamente
(sharp d=1 m=108 n=2: 0,07126240410476363 com p=0 e 0,07126240410476364 com p=1).
O autor **não** usa esse atalho em R3 (calcula os dois p), o que é conservador.

### 12.4 Cobertura das regiões

**Veredito: CONFIRMADO** (em papel e na estrutura dos checkpoints).

Uma coluna fora de G=107×384 tem m>=107 ou n>=384. Se m>=160, cai na cauda de Fourier
d_F(160), válida para qualquer n. Se n>=900, cai na cauda radial d_R(900), válida para
qualquer m. Caso contrário m<160 e n<900, e então: n<384 força m>=107 (senão estaria em
G) e cai em R3; 384<=n<606 cai em R0 (m<107) ou R1 (107<=m<160); 606<=n<900 cai em R2
(m<160). As três fronteiras críticas ficam do lado certo: m=106 em R0/R2, m=107 em
R1/R3, m=159 em R1/R2/R3, m=160 na cauda; n=383 em R3, n=384 em R0/R1, n=605 em R0/R1,
n=606 em R2, n=899 em R2, n=900 na cauda. Não há buraco; as sobreposições (as caudas
cobrem regiões já cobertas) só podem reduzir o máximo, nunca escondê-lo.

Conferi a contagem de colunas por linha nos checkpoints novos, que é onde um gerador
com paridade errada apareceria:

| arquivo | linhas | colunas por linha | esperado |
|---|---|---|---|
| A-m107-M160-D60, R1 (n>=384) | 222 | 105 (n ímpar) / 79 (n par) | (3 comp × 26 m pares + 27 m ímpares) e sem K em n par, p=0 só ✔ |
| A-m107-M160-D60, R3 (n=0..39) | 40 | 210 / 158 | o mesmo, com p=0 e p=1 ✔ |
| sharp-m107-M160-D100, R1 | 222 | 52 / 26 | 2 comp × 26 m pares, p=0 só ✔ |
| sharp-m107-M160-D100, R3 (n=0..95) | 96 | 104 / 52 | o mesmo, com os dois p ✔ |

K só aparece em n ímpar nos quatro casos, as paridades de m batem por componente, e o
atalho p=1≡p=0 só é usado no regime desdobrado, onde a §11 já o tinha confirmado.

**Pendências registradas:** R2 ([606,900) × [0,160)) ainda não tinha checkpoint quando
escrevi isto, e R1/R3 estavam incompletas (A até n=591 em R1 e n=39 em R3; sharp até
n=95 em R3). A montagem final tem de conferir a contiguidade das três regiões e que
todas as linhas têm U<1 — o `holes_count` do autor faz isso, e o que verifiquei aqui é
a regra de cobertura que ele implementa, não a execução completa.
