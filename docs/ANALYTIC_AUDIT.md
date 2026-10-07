# Auditoria da ponte entre o pencil e os certificados

## 1. Três matrizes diferentes

Escreva `P0 = L(z0)`, `Π` para a projeção de Galerkin, e `ι` para
a inclusão do espaço finito. Mesmo no caso `J=I`, são objetos diferentes:

| Objeto | Definição |
|---|---|
| Matriz espectral exportada | `A_N = Π L(0) ι` |
| Resolvente finita | `R_N = (A_N + z0 I)^-1` |
| Compressão da resolvente | `M_N = Π P0^-1 ι` |

O determinante calculado por `build_contour.py` corresponde a
`det(A_N+sI)`. Sua normalização em `z0` é exatamente
`det(I+(s-z0)R_N)`. Até esta auditoria, `build_tile.py` construía
`I+(s-z0)A_N`: certificava uma família diferente. Agora calcula `R_N`
em aritmética racional; `--base-point` fixa `z0` (padrão 1).
`--raw-coefficient` mantém explicitamente o antigo benchmark.

Os tiles gerados continuam sendo finitos e, neste gerador, têm centros e
deslocamentos reais. Sua aceitação não cobre por si só um contorno complexo.
Submatrizes escolhidas com `--dim` são apenas experiências de custo.
Os exemplos Lean antigos permanecem testes aritméticos; não se tornam
certificados do pencil por esta correção.

## 2. Inverter e comprimir não comutam

Para `P0 = [[2,1],[1,2]]`, sua inversa é
`[[2/3,-1/3],[-1/3,2/3]]`. Projetando sobre a primeira coordenada:
`M_N = [2/3]`, mas `R_N = [1/2]`. A diferença é `1/6`.
As duas identidades de inversa e a diferença são verificadas pelo kernel em
`formal/ChoptuikFormal/ResolventCompression.lean`.

Logo a identificação de `M_N` com a matriz calculada, na §2 de
`TAIL_BOUND.md`, requer correção. O lema de Sylvester permanece válido para
`ΠC` e `ΠCΠ`, mas não identifica nenhuma delas com a inversa de `ΠP0ι`.

Há uma identidade útil para construir o bound que falta. Se `R_N` e
`R0=P0^-1` existem e `ι V_N` está no domínio de `P0`, defina

```
E_N = P0 ι R_N - ι,
R0 ι - ι R_N = -R0 E_N,
M_N - R_N = -Π R0 E_N.
```

A segunda igualdade resulta de aplicar `R0` à definição de `E_N`; a
terceira usa `Πι=I`. Portanto
`||M_N-R_N|| <= ||Π|| ||R0|| ||E_N||`, nas normas compatíveis.
O resíduo `E_N` inclui componentes fora do truncamento. Calcular somente
`Π E_N = 0` não fornece esse bound.

Para usar a aproximação bilateral implementável `ι R_N Π`, a desigualdade
triangular fornece, se as projeções têm norma no máximo 1,

```
||C - ι R_N Π|| <= ||C - ι M_N Π|| + ||M_N-R_N||.
```

Isso não resolve a escolha unilateral: `ΠC` depende também de dados fora
de `V_N`. A estimativa unilateral antiga não deve ser simplesmente aplicada
a `ι R_N Π`. É necessário justificar a cauda bilateral ou construir outra
aproximação com seu próprio majorante. Esta é a obrigação nova C1-G.

## 3. P2 não é consequência de coeficientes analíticos

Considere o modelo periódico `P=1+∂τ` no espaço de Wiener sem pesos
`Y=ℓ¹(Z)`, com domínio
`X={u: sum_n (1+|n|)|u_n| < infinity}`.
Para `f_n(τ)=exp(inτ)`, temos

```
||f_n||_Y = 1,
R f_n = exp(inτ)/(1+in),
||R f_n||_{σ} = exp(σ n)/sqrt(1+n²).
```

A inversa é limitada `Y→X` pois
`(1+|n|)/sqrt(1+n²) <= sqrt(2)`, mas para todo `σ>0`
a última expressão diverge quando `n→infinity`.
Assim, P1 vale neste modelo e P2 não: a composição P1/P2 exigiria
`||R f_n||_σ <= A0 M0` uniformemente.

O mesmo modelo tem resolvente compacta em `Y`: sua cauda é exatamente
`||(I-Π_N)R|| = 1/sqrt(1+(N+1)²)`, que tende a zero.
Logo compactação não exige o ganho exponencial proposto em P2.

Mais geralmente, se um operador diferencial de ordem finita tem modos
admissíveis `u_n` com `||u_n||_σ >= exp(σ n)` e
`||P0 u_n||_Y <= c(1+n)^r`, uma inversa `Y→𝒜_σ` limitada é impossível:
aplicar o bound a `f=P0u_n` exigiria que uma exponencial fosse majorada por
um polinômio. A aplicação ao domínio exato de Choptuik ainda exige fixar esse
domínio e seus modos admissíveis; o exemplo não é uma refutação do teorema
físico. Ele refuta a inferência de P2 a partir da mera analiticidade do fundo.

Uma alternativa concreta é buscar ganho finito de regularidade com cauda
algébrica, ou trabalhar desde o início nos espaços de funções suaves com
pesos apropriados. Reiterer constrói resolventes compactas em espaços de
funções infinitamente diferenciáveis sob hipóteses específicas; a aplicação
ao linearizado de Choptuik é explicitamente deixada em aberto (p. 4 do PDF,
numeração impressa 4 após a introdução motivacional).
Fonte primária: [Reiterer, 1510.05310](https://arxiv.org/pdf/1510.05310),
introdução e §§3–4. Não há nesse resultado um bound P2 para `Y=ℓ¹` sem pesos.

## 4. O que o Lean atual conclui

`CertifiedSpectralObligations` é um ledger condicional. Seus campos de tipo
`Prop` são escolhidos pelo produtor; não são definições dos espaços físicos
ou do operador Einstein–escalar. O campo `countTransfer` já supõe a igualdade
da multiplicidade física com as contagens. Não há erro lógico no teorema
de contabilidade, mas preencher campos com proposições triviais não provaria
o enunciado físico. A ponte PDE precisa ser definida e demonstrada.

Não há, neste estado, certificado completo de Choptuik. As correções acima
impedem que a validação das contas finitas seja confundida com essa ponte.
