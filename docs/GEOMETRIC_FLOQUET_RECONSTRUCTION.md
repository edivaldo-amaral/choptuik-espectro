# Reconstrução geométrica global no ansatz, para Re s>0

## Atualização 03/10/2026: revisada (a volta), enunciado preciso e ligação com L e C

**Revisão independente** (`.codex-runs/2026-10-03-revisao-V/`): **não quebrou.**
- As compatibilidades e a identidade de log W foram provadas em sympy com funções genéricas, fora das
  soluções (`v123_identidades.py`).
- L(s) e C(s) do código coincidem ponto a ponto, a 10⁻¹⁵, com SδΩ⁺ e S♯δΩ⁺ do TeX (`v1_L_C_pontual.py`).
- As fórmulas de curvatura de RT batem com o Ricci 4D calculado direto da métrica (`v6_ricci_rt.py`).
- As cadeias foram conferidas em `v7_cadeias.py`.
- Os testes de `test_geometric_reconstruction.py` são de regressão, numa instância polinomial; a
  verificação genérica está nos scripts da revisão (L5).

**Enunciado (volta).** Seja Re s > 0 e h = (h1, ..., h4) tal que:
- h é 4π-periódico, com h1 ímpar;
- h ∈ D(Y+), setor A, B ou soma; para os autovetores de T2 isso vem do lema L-real;
- L(s)h = 0 e C(s)h = 0.

Então existem q, z, p únicos, 4π-periódicos, pares e real-analíticos em ℝ × (−41/40, 41/40), tais que
(δQ, δζ, δφ) = e^{sτ}(q, z, p), no ansatz, tem estas propriedades:
1. as variáveis de RT dela são e^{sτ}h;
2. (δg, δφ) é real-analítica em M̂, em coordenadas cartesianas, com eixo e cone (ξ = 1) incluídos;
3. resolve as equações de Einstein–campo escalar linearizadas;
4. (Θ²)*δg = e^{4πs}e^{−4K}δg e (Θ²)*δφ = e^{4πs}δφ.

Isto é, (δg, δφ) é um modo físico de Floquet no sentido de `HREC_IDA.md`. A aplicação h ↦ δ é linear e
injetiva. O semieixo da elipse de Bernstein do peso 5/4 é (5/4 + 4/5)/2 = 41/40, o domínio M de RT.

**Lema V1 (L e C juntos dão os dez δΩ nulos; lacuna L1).** Definições:
- S é Ξ⁻¹ aplicado a ½(1 + P) na componente 1 e às componentes 2, 4 e 5;
- S♯ é ½(1 − P) na componente 1 e −P na componente 3;
- Ξ⁻¹f = (f − f|_{ξ=0})/ξ;
- L(s)h = SδΩ⁺[h] e C(s)h = S♯δΩ⁺[h].

Prova:
1. As componentes selecionadas se anulam em ξ = 0 para todo h com h1 ímpar, porque os termos principal e
   quadrático têm fator ξ. Logo nada se perde em Ξ⁻¹.
2. L(s)h = 0 dá ½(1 + P)δΩ1⁺ = δΩ2⁺ = δΩ3⁺ = δΩ4⁺ = 0. C(s)h = 0 dá ½(1 − P)δΩ1⁺ = 0 e δΩ2♯⁺ = 0.
   Logo δΩ⁺ = 0.
3. PδΩ_i^σ = −δΩ_i^{−σ}, logo δΩ⁻ = 0.

Os dez δΩ^σ são nulos. É a Remark de RT após 2DprobSeries, linearizada, e `SHARP_LINEARIZED_IDENTITY.md`
(e = Ic + RF). **C só é necessário para recuperar h2** pela identidade de log W; as compatibilidades
(2,2), (3,3) e (4,4) já saem de L.

**De onde vem C(s)h = 0 em T2.** Da injetividade de K e de KC = ML. Os certificados de K estão no
toro♯(129/128, 9/8), e C(s)h cai nesse espaço pela perda de raio de (65/64, 5/4) para (129/128, 9/8)
(`SHARP_DOMAIN.md`, fechamento de 03/10). Era a observação O1 da revisão.

**Cadeias (L4).** Para uma cadeia, exigem-se as constraints de **jato**: C(s0)h0 = 0 e
C(s0)h_j + C1h_{j−1} = 0, e não C(s0)h_j = 0 em cada vetor, como em `CONSTRAINT_ROOT_QUOTIENT.md`. T2 não
depende disso, porque as raízes são simples.

**Cotas (L6).** As constantes de norma abaixo não são usadas em T2. Elas vêm de |1/(s + im/2)| <= 1/Re s
e de ‖·‖_C <= √2|·|; b é a cota de ‖ξ‖ no espaço.

**A "pendência geométrica" do fim deste documento** é a ida, a fixação do gauge, provada em
`HREC_IDA.md` (Lema R, revisado). A escolha de cone fixo ou móvel é a modelagem H-cone.

Continuação: GLOBAL_NULL_GAUGE_TRANSPORT.md resolve o transporte para
remover as componentes radiais diagonais numa classe analítica explícita
no cilindro inteiro. A inversa usa o traço no cone separado da integral
de características; nas ressonâncias há uma condição de compatibilidade.
A identificação com o domínio espectral RT e todas as condições físicas
de admissibilidade continuam pendentes.

Fixe o fundo exato RT, mu fixo, W=mu+xi²Q>0 no domínio real e as
normalizações/paridades RT. Esta nota estabelece a reconstrução no
CILINDRO inteiro e sua cobertura temporal, não a fixação global do gauge
de uma perturbação geométrica arbitrária. O sistema linearizado COMPLETO
significa todos os resíduos delta Omega^sigma nulos, não só a seleção.

## Potenciais sem obstrução de períodos

Para h_i^-=h_i, h_i^+=-P h_i (i=2,3,4), defina

```
a_i=[(1-xi)h_i^+-(1+xi)h_i^-]/2,
b_i=(h_i^++h_i^-)/(2mu),
T_s=partial_tau+s.
```

As equações D_s^sigma z=h_i^sigma equivalem a T_s z=a_i e
partial_xi z=b_i. Em funções 4pi-periódicas, para Re s>0,

```
T_s^-1 a(tau,xi)= integral_0^infty exp(-s t) a(tau-t,xi) dt.
```

Essa inversa preserva as paridades temporal e radial, com norma no máximo
2/(Re s) nos DOFs reais complexificados ponderados (as translações em
tau têm norma no máximo 2). É holomorfa em s no semiplano positivo.
Se T_s b_i=partial_xi a_i, então z=T_s^-1 a_i também satisfaz z_xi=b_i.
Não há condição de período adicional: o multiplicador exp(4pi s) difere
de 1. Equivalentemente, a primitiva na cobertura tem um único ajuste de
constante que satisfaz a condição Floquet. O argumento NÃO vale em s=0.

Os resíduos completos fornecem a compatibilidade:

```
-2mu xi (T_s b_i-partial_xi a_i)
    = delta Omega_j^- - delta Omega_j^+,
(i,j)=(3,3),(4,4),(2,2).
```

As três identidades foram verificadas diretamente com polinômios
racionais e as tabelas RT. Os termos algébricos cancelam na subtração
dos sinais; não se supõe o fundo exato para essas identidades de diferença.

## Reconstrução explícita

Os perfis das variações são

```
q = (h1-h2^- -h2^+)/(2xi),
z = T_s^-1 a_3,
p = T_s^-1 a_4.
```

O numerador de q é ímpar; a divisão no centro é removível. z,p são pares.
As perturbações geométricas são exp(s tau)(q,z,p) de (Q,zeta,phi).
O inverso em tau preserva os pesos; por exemplo,
||z|| <= 2(1+b)/(Re s)||h3||. Para q usa-se o bound da divisão
regularizada por xi. Não foi necessária uma integral radial atravessando
o cone característico xi=1.

Falta verificar a definição de h2, que contém log W. A identidade
ALGÉBRICA, válida fora das soluções, é

```
D^sigma W - W(omega2^sigma-omega3^sigma)
   = (Omega1^sigma-Omega2^sigma-Omega2sharp^sigma)/2.
```

Linearizando no fundo exato, com delta W=xi²q, obtemos

```
D_s^sigma(xi²q/W)-(h2^sigma-h3^sigma)
   = (delta Omega1^sigma-delta Omega2^sigma-delta Omega2sharp^sigma)/(2W).
```

Logo h2^sigma=D_s^sigma(z+xi²q/W), h3^sigma=D_s^sigma z e
h4^sigma=D_s^sigma p. A fórmula de q dá também a definição de h1.
As identidades off-shell e sua linearização sem omitir o resíduo foram
testadas em test_geometric_reconstruction.py. Em fundo aproximado, ao
linearizar a identidade DIVIDIDA por W, há ainda o termo
`-(D^sigma W-W(omega2^sigma-omega3^sigma))*delta W/W²`.

## Centro, cone e métrica

Por paridade, q,z,p são funções analíticas das coordenadas cartesianas
radiais. O centro não é uma singularidade da reconstrução. No domínio
real W>0, e as fórmulas são analíticas também em xi=1. Não presumimos
que 1/W pertença ao mesmo espaço de coeficientes na elipse complexa:
a conclusão geométrica usa a não degenerescência REAL e analiticidade
em vizinhanças dos compactos reais.

No frame RT, delta e0=0 e o perfil de delta e_i é
`q sum_k x^k(x^k partial_i-x^i partial_k)`. Assim

```
delta g^-1 = exp(s tau) exp(2zeta)
  [2z(-e0 tensor e0+sum_i e_i tensor e_i)
   +sum_i(delta e_i tensor e_i+e_i tensor delta e_i)],
delta g = -g (delta g^-1) g,      delta phi=exp(s tau)p.
```

Aqui delta e_i dentro dos colchetes designa somente seu perfil, sem
repetir exp(s tau). Essas fórmulas são globais na cobertura temporal
do domínio RT. Não afirmam pequenez uniforme de uma perturbação instável
quando tau tende a infinito, nem existência de uma família não linear.

## Equações geométricas e unicidade dentro do ansatz

As fórmulas RT `dfkhdjhsdshkfd` expressam Ric-2dphi tensor dphi e Box phi
linearmente nos resíduos Omega. No fundo exato, sua derivada aplicada
a delta Omega=0 prova as equações geométricas linearizadas. Não se usa
integrabilidade de h numa família de soluções não lineares.

Reciprocamente, perfis geométricos no ansatz fornecem, por suas definições,
Omega1^sigma=Omega2^sigma+Omega2sharp^sigma e igualdade entre os sinais
para Omega2, Omega3, Omega4, todas entendidas LINEARIZADAMENTE. Substitua
essas relações nas fórmulas de curvatura: os cinco resíduos geométricos
independentes forçam os cinco resíduos RT restantes a zero fora do
centro; analiticidade estende a igualdade ao centro.

O núcleo da reconstrução é trivial para Re s>0: h=0 implica q=0 e
T_s z=T_s p=0, logo z=p=0. As constantes conformes/de phi invisíveis
a omega não possuem perfis periódicos com expoente real positivo.
Assim há uma equivalência linearizada GLOBAL NO ANSATZ, em Re s>0.
Por holomorfia em s, a reconstrução também se aplica a jets/cadeias.

Explicitamente, para uma cadeia h_j no expoente s0, reconstrua
`T_s0 z_j=a_3(h_j)-z_(j-1)` e `T_s0 p_j=a_4(h_j)-p_(j-1)`, com
z_-1=p_-1=0; q_j continua algébrico. Reconstruir cada h_j isoladamente
com T_s0^-1 perderia os termos da cadeia. A recorrência foi testada
racionalmente, assim como a invertibilidade das cinco combinações de
resíduos geométricos usadas no argumento recíproco.

## Pendência geométrica precisa

Ainda é necessário provar que uma perturbação geométrica arbitrária da
classe física pode ser levada a esse ansatz com as condições globais de
centro, cone e Floquet, por um gauge admissível e com estimativas. A
equivalência acima não prova isso. A escolha de cone fixo ou móvel também
não foi feita por conveniência espectral. A contagem física permanece aberta.

Há uma condição necessária explícita para essa fixação de gauge. Em
coordenadas nulas, g_++=g_--=0 e g_+- não é zero. Para remover os termos
radiais diagonais de uma perturbação k por k+Lie_X g, é preciso resolver

```
partial_+ X^- = -k_++/(2g_+-),
partial_- X^+ = -k_--/(2g_+-).
```

Se o gauge deve preservar o cone u_-=0, então X^- vale zero ao longo
dele e a primeira equação exige k_++|_(u_-=0)=0. Não se pode resolver
isso para dados arbitrários mantendo essa condição de fronteira.
No domínio de cone móvel, o valor de X^- no cone é uma incógnita de
deslocamento, não zero imposto a priori. Restam as condições de centro,
Floquet e a solvabilidade global desses transportes. Essa distinção
explica por que a escolha do domínio físico não é uma mera normalização.
