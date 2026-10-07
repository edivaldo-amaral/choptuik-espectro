# Domínio sharp por fechamento e perda controlada de raio

## Realizações mínima e máxima coincidem para a parte livre

Fixe pesos a>1,b>1 e o espaço Y(a,b) de coeficientes RT somáveis.
No sistema sharp a parte livre é Jsharp=diag(K_mu,J_mu), com primeira
componente ímpar. Defina D_min como o fechamento do gráfico sobre
polinômios Fourier–Chebyshev admissíveis. Defina

```
D_max = {u em Y(a,b): Jsharp u, calculado diferencialmente, pertence a Y(a,b)}.
```

Temos D_min=D_max. Eis a demonstração, usando a inversa livre já construída:

1. Para z>0 real, Qsharp(z) é limitado em Y(a,b), preserva polinômios e
   satisfaz as identidades de inversa nesses polinômios. Por aproximação
   dos dados, Qsharp(z) leva Y(a,b) para D_min e é inversa de Jsharp+z ali.
2. D_min está contido em D_max: convergência em Y implica convergência
   local uniforme, inclusive de derivadas em compactos interiores; portanto
   o limite do gráfico coincide com o operador diferencial.
3. Para u em D_max, ponha v=Qsharp(z)(Jsharp+z)u. Então v está em D_min
   e w=u-v satisfaz (Jsharp+z)w=0, como função analítica.
4. Em cada modo Fourier assinado m, a equação escalar J tem solução local
   proporcional a (1+xi)^(-1-(z+im/2)/mu), e K tem solução proporcional
   a xi^(-1-(z+im/2)/mu). Como z>0, nenhuma solução não nula é analítica
   em xi=-1, respectivamente xi=0. Esses pontos estão no interior da
   elipse de Bernstein de parâmetro b>1. Logo cada coeficiente é zero.
   Assim u=v pertence a D_min.

O argumento usa a realização complexa de Fourier; a conversão para os
DOFs reais complexificados é uma isomorfia limitada. Vale também para
o principal selecionado de quatro componentes.

Como Bsharp=mu R_x Gamma1sharp+Gamma2sharp(w,.) é limitado no espaço
ponderado (w pertence a ele), K(s)=Jsharp+Bsharp+sI é fechado nesse
mesmo domínio, independente de s. Isso dá uma família afim no domínio
fixo, sem exigir que derivadas separadas sejam limitadas em Y.

## A identidade passa aos domínios com perda de raio

Escolha 1<a'<a e 1<b'<b. As derivadas são limitadas de Y(a,b) para
Y(a',b'). Bounds suficientes, não otimizados, são

```
||partial_tau|| <= q_a/[2(1-q_a)^2],       q_a=a'/a,
||partial_xi||  <= 2 q_b/[(b'-1)(1-q_b)^2], q_b=b'/b.
```

A primeira estimativa usa a diagonal m/2. Para a segunda, T'_n=n U_(n-1)
e a soma geométrica dos pesos dão 2n(b')^n/(b'-1); a razão com b^n
é majorada usando sum n q_b^n=q_b/(1-q_b)^2.
Reflexão, multiplicação por xi e R_x são limitados nos espaços envolvidos.

Por exemplo, de (65/64,5/4) para (129/128,9/8), os bounds são 8385 e
1440. São grandes, mas finitos: servem ao domínio, não a uma contração.

Escreva D_L para o domínio selecionado em Y(a,b), com norma de gráfico.
As fórmulas explícitas de SHARP_LINEARIZED_IDENTITY mostram que
C(s):Y(a,b)->Ysharp(a',b') e M(s):Y(a,b)->Ysharp(a',b') são limitados,
localmente uniformemente em s. Para h em D_L, L(s)h pertence a Y(a,b).
A identidade diferencial no fundo exato dá

```
Jsharp C(s)h = M(s)L(s)h - (Bsharp+sI)C(s)h em Ysharp(a',b').
```

Logo C(s)h está em D_max sharp no raio menor e, pelo lema anterior,
em D_min sharp. Além disso essa fórmula limita sua norma de gráfico
pela norma de gráfico de h. Obtivemos, portanto, uma aplicação limitada
`C(s):D_L(a,b)->D_K(a',b')`, e a identidade K(s)C(s)=M(s)L(s)
nesse sentido de operadores fechados.

Em fundo aproximado, acrescente -G(h,e(w)). Para obter a mesma conclusão
é necessário que esse termo pertença ao espaço de saída; isso vale com
perda de raio quando w é analítico com os pesos maiores e o resíduo é
formado pelas fórmulas diferenciais. O termo não pode ser omitido.

## Consequência e limite

Para h no kernel selecionado, C(s)h pertence ao domínio sharp NO RAIO MENOR
e está no kernel de K(s). Assim a injetividade sharp nesse espaço menor
elimina a constraint. Isso fecha a questão de domínio necessária para
essa implicação, sob as hipóteses de pesos e inversa livre explicitadas.

Não provamos C(s):D_L(a,b)->D_K(a,b), sem perda de raio. Também não
provamos injetividade sharp nos novos pesos nem equivalência com o domínio
geométrico físico. Bounds sharp certificados apenas nos pesos originais
não podem ser reaproveitados automaticamente: o espaço menor em raio é
MAIOR em funções. Por exemplo, o bound livre conservador muda de
18/(mu+z) para 34/(mu+z) ao trocar o peso radial 5/4 por 9/8.

Esta é uma demonstração funcional em papel, não uma formalização Lean.

## Fechamento (03/10/2026): a injetividade certificada já está nos pesos menores

A revisão independente da volta (`.codex-runs/2026-10-03-revisao-V/`, observação O1) apontou que "não
provamos injetividade sharp nos novos pesos" poderia deixar a implicação K injetivo ⇒ C(s)h = 0 sem
apoio. Não deixa, porque **todos os certificados de injetividade de K estão nos pesos menores**:
- S3a: tiles de Perron, `signed_operator.py` / `fourier_tail_bound.py`;
- S3b: Rouché de K, `rouche_K.py`;
- R♯ <= 1,765: símbolo, `certify_symbol_numerical_range.py`.

Todos estão no L² do toro com sinal de raios **(129/128, 9/8)**.

A cadeia é a seguinte:
1. h ∈ D_L(65/64, 5/4), com L(s)h = 0. Para os autovetores de T2, isso vem do lema L-real.
2. Pelo exemplo acima, com (a, b) = (65/64, 5/4) e (a', b') = (129/128, 9/8), cujas constantes são 8385 e
   1440, C(s)h ∈ D_K(129/128, 9/8) ⊂ ℓ¹(129/128, 9/8), e K(s)C(s)h = M(s)L(s)h = 0.
3. Inclusão ℓ¹(129/128, 9/8) ⊂ toro♯(129/128, 9/8), com norma <= 1. Nos índices com sinal,
   κ1^m <= κ1^{|m|} e √(κ2^{2j} + κ2^{−2j}) <= 2κ2^j, que é o peso de ℓ¹ do par ±j.
   (`S3_CONDICIONAMENTO_SHARP.md` §9.1.)
4. **Domínio.** Seja u = C(s)h. Então J♯u = −(B♯ + s)u, que está em ℓ¹ ⊂ toro. Logo u está no domínio
   máximo da realização no toro. O lema D_min = D_max vale também no toro, com o mesmo argumento: Q♯(z)
   é limitado ali, e as soluções homogêneas são singulares em ξ = 0 ou ξ = −1, dentro da elipse. Assim
   u está no domínio em que os três tipos de certificado afirmam a injetividade de K(s):
   - Perron e Rouché: H = I + (...)Q0 invertível ou sem zeros, com Q0 = (J♯ + s0)⁻¹;
   - imagem numérica: ‖K(s)u‖ >= (Re s − R♯)‖u‖ no domínio.

   Como K(s)u = 0, vem u = 0.
5. Logo C(s)h = 0 em toda s da faixa A com Re s >= 0, exceto s_K.

Para cadeias, o mesmo vale ordem a ordem, derivando KC = ML. A perda de raio de (65/64, 5/4) para
(129/128, 9/8) cabe exatamente no intervalo entre os espaços de L e de K. Não é preciso recuperar o raio
para K.

