# Formulação do problema espectral

## 1. Fundo que está realmente certificado

Usamos a normalização de RT,

\[
  \operatorname{Ric}_g=2\,d\phi\otimes d\phi,\qquad \Box_g\phi=0.
\]

O Teorema 1 de RT fornece uma solução esfericamente simétrica, real-analítica
e não plana no domínio

\[
  \mathcal M=\{(u_-,u_+):u_+<0,\ u_+<u_-<-u_+/81\}.
\]

Ela atravessa o cone de autossimilaridade passado `u_-=0`; não inclui o ponto
singular `(0,0)`, o horizonte de Cauchy futuro ou o infinito assintoticamente
plano. Para a transformação fundamental `Theta`,

\[
 \phi\circ\Theta=-\phi,\quad \Theta^*g=e^{-2K}g,
 \quad K\simeq1.7227262011139.
\]

O período escalar usual é `Delta=2K ≈ 3.4454524022278`. O parâmetro
`mu_RT ≈ 0.168307078963` controla a reescala dos nulos `u_±`; ele não é o
período de eco.

## 2. Pencil direto nas variáveis RT

O fundo é representado por `omega*=(omega_1,...,omega_4)` em séries de Fourier
em `tau_R` e Chebyshev em `xi`. RT escrevem o sistema selecionado como

\[
 S\Omega^+(\mu,\omega)=
 S\{\Xi H_\mu\omega+\mu\Gamma_1\omega+
       \Xi\Gamma_2(\omega,\omega)\}=0.
\]

Para uma perturbação de Floquet

\[
  \delta\omega(\tau_R,\xi)=e^{s\tau_R}h(\tau_R,\xi),
\]

substitui-se `partial_tau` por `partial_tau+s` apenas onde a derivada atua em
`h`. Como `H_mu` é afim nessa derivada, resulta o pencil

\[
 \boxed{L_A(s)h=L_A(0)h+s h},
\]

\[
 L_A(0)h=O_\mu h+\mu S\Gamma_1h+
                 2S\Xi\Gamma_2(\omega_*,h).
\]

Esta é exatamente a parte reutilizada do código RT pela extensão em
[`patches/rt-spectral.patch`](../patches/rt-spectral.patch). Se `q` é um
autovalor da matriz de Galerkin de `L_A(0)`, o localizador retorna `s=-q`.
A conversão para a convenção física `tau=-log(-t/L)` é

\[
  \lambda=\frac{2\pi}{K}s.
\]

Assim, `lambda ≈ 2.674` corresponde a `s ≈ 0.733`.

## 3. Dois setores, não um

A involução “meio eco + troca de sinal de `phi`” decompõe toda perturbação de
período completo em dois setores invariantes:

- `A`, mesma simetria do fundo: `h_1,h_2,h_3` periódicos em `2pi` e `h_4`
  antiperiódico;
- `B`, simetria oposta: paridades trocadas.

Os tipos `V,W` do código RT implementam somente `A`. Multiplicar o perfil por
`exp(i tau_R/2)` troca `A` por `B` e desloca `s` por `i/2`. Uma contagem física
usa ambos, sem dupla contagem.

**Isto deixou de ser uma pendência de implementação.** A conjugação
`M^{-1} L_B(s) M = L_A(s+i/2)` foi verificada e transfere multiplicidades por
Gohberg--Sigal, de modo que `Sigma_B = Sigma_A - i/2` como conjunto. Contar em
`B` numa região é idêntico a contar em `A` na região transladada: **nada precisa
ser complexificado**, e a obrigação S2 saiu da lista de dependências do teorema
final. A derivação, as hipóteses que ela consome e a consequência operacional
exata para o contorno estão em [`SECTOR_B.md`](SECTOR_B.md).

Um detalhe que a redução expõe e que não é cosmético: os modos **reais** do setor
`B`, transportados para `A`, vivem sobre `Im(s) = 1/2`. A janela simétrica
`[-1/2, 1/2]` os colocaria exatamente **sobre** o contorno de contagem; daí a
faixa deslocada `(-1/4, 3/4]`.

O diagnóstico em `scripts/analyze_spectrum.py --show-b` apresenta as raízes
`s_B = s_A - i/2` e conserva o mesmo resíduo de constraints.

## 4. Constraints e gauge

O resíduo que não pode ser ignorado é

\[
 C_A(s)h=S^\sharp\{\Xi H_{\mu,s}h+\mu\Gamma_1h+
                    2\Xi\Gamma_2(\omega_*,h)\}.
\]

RT provam, para o problema não linear periódico, que as equações selecionadas
implicam as constraints por invertibilidade de um sistema homogêneo. Para o
pencil é necessário repetir a identidade com o deslocamento `s` e excluir
kernel do sistema `sharp` nos pontos característicos interiores relevantes
(não basta verificar apenas o contorno). Até isso ser provado,
raízes de `L_A` são apenas candidatas.

Veja `CONSTRAINT_GAUGE_EQUIVALENCE.md` para os lemas de kernels, cadeias e
quociente por gauge, com suas hipóteses ainda não verificadas.

Também é necessário classificar difeomorfismos que preservam o domínio. O
modo neutro `s=0`, tangente à fase DSS, pode ser removido por uma condição de
fase. A variação do tempo de acumulação corresponde geometricamente a
`lambda=1`, mas move o cone/domínio e pode surgir como polo em vez de autovetor
do BVP. Ela não deve ser subtraída da contagem por decreto.

## 5. Teorema-alvo mínimo

Fixado o fundo no ball certificado de RT, considere perturbações esféricas
analíticas no centro e através do cone passado, satisfazendo as constraints,
na realização gauge-fixada que se prove equivalente ao quociente físico.
Módulo a equivalência de Floquet, o pencil tem multiplicidade algébrica total
um em `Re(lambda)>0`; a única classe é real e simples, e não há espectro físico
em `Re(lambda)=0` após a condição de fase.

Este teorema é deliberadamente sobre o **espectro analítico esférico**. Ele não
implica por si só estabilidade de semigrupo, variedade estável não linear,
universalidade do limiar ou separação global entre dispersão e buraco negro.

## 6. Evidência numérica reproduzida

O localizador dyádico produziu:

| Fourier x Chebyshev | dimensão real | `s` candidato | `lambda` candidato |
|---:|---:|---:|---:|
| 4 x 12 | 138 | 0.7096266582 | 2.5881743655 |
| 6 x 18 | 333 | 0.7311839087 | 2.6667986990 |
| 8 x 24 | 612 | 0.7330938195 | 2.6737645905 |
| 10 x 30 | 975 | 0.7331778411 | 2.6740710368 |
| 12 x 36 | 1422 | 0.7331796783 | 2.6740777376 |

Na extensão atual também calculamos o resíduo `l2` das constraints para o
autovetor truncado: `1.14e-1`, `1.94e-2`, `3.61e-3`, `5.61e-4` e `9.17e-5`,
respectivamente. Essa queda é um teste de consistência do discretizado; continua
sem ser um bound para a solução infinita.

O `10 x 30` cai em `lambda = 2.67407`, dentro da barra de Gundlach
(`2.674 ± 0.009`) e a `7·10^-5` do valor central. Isso é concordância
numérica, não prova: o objeto calculado continua sendo o determinante de um
truncamento.

**Não é o único modo instável do truncamento.** As mesmas matrizes têm outras
duas raízes reais com `Re s > 0` dentro do contorno de contagem, com assinaturas
distintas — uma com resíduo de constraints estagnado, outra convergindo para
`mu_RT`. Ver [`MODE_DISCRIMINATION.md`](MODE_DISCRIMINATION.md).

O valor publicado por Gundlach é `lambda=2.674±0.009` na convenção crescente
(o artigo usa a convenção de sinal oposta). A convergência acima é evidência
de que o pencil foi identificado corretamente, mas não é um enclosure nem uma
contagem.
