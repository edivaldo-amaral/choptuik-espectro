# O que é nosso, o que é conhecido

> **Atualização 03/10/2026:** documento de 22/09, anterior a T2. T2 está provado e revisado; ver `PANORAMA.md` e `T2_ENUNCIADO.md`. A §4 ("Não contamos os modos instáveis") está superada.

Levantamento de 22/09/2026, feito a pedido, antes de qualquer conclusão sobre
originalidade. Verificação contra a literatura, com fontes.

## 1. Técnicas que usamos e que são padrão

Nenhuma das ferramentas que renderam os resultados desta semana é nova. Todas
são clássicas ou bem estabelecidas; o que fizemos foi **aplicá-las a este
operador**, onde não tinham sido aplicadas.

| o que usamos | estado na literatura |
|---|---|
| decomposição da multiplicação de Chebyshev em **Toeplitz + Hankel + posto 1** (§5.19) | **conhecida**. Olver & Townsend, *A fast and well-conditioned spectral method* (arXiv:1202.1347), afirma exatamente isso. **Redescobrimos**, correção de posto 1 inclusa |
| símbolo na **elipse de Bernstein** ↔ decaimento geométrico dos coeficientes de Chebyshev | clássico (Trefethen, *ATAP*) |
| `spec(T) ⊂ fecho(W(T))`, campo de valores para enclausurar espectro | clássico; há linha ativa em *Pseudo Numerical Ranges and Spectral Enclosures* (arXiv:2204.03584) e em enclosures/exclosures para problemas não auto-adjuntos (LMS J. Comput. Math.) |
| semisseparável + diagonal ⇒ `LDL^T` em recursão escalar | teoria madura (Vandebril, Van Barel, Mastronardi) |
| pencil generalizado para evitar irracionais | álgebra linear numérica padrão |
| convergência da seção finita | teoria de operadores limite (Rabinovich–Roch–Silbermann) |
| Combes–Thomas / decaimento exponencial de resolvente | clássico em teoria espectral |
| piso de Perron para reponderação diagonal | clássico |

**Consequência honesta: não temos técnica nova.** Temos aplicação de técnica
conhecida a um problema onde ela não tinha sido aplicada — o que é contribuição
legítima, mas de outra natureza, e deve ser apresentada assim.

## 2. O que o projeto herda, e que não é nosso

- A **existência** da solução crítica: Reiterer & Trubowitz, *Choptuik's
  critical spacetime exists*, CMP 368 (2019) 143 (arXiv:1203.3766). Prova
  assistida por computador, aritmética racional, **normas de operador em `l^1`**
  nos argumentos de contração.
- A compacidade da inversa do operador linear relevante, na expansão
  Fourier–Chebyshev — também deles.
- O fundo `RefA.dat` com ~80 dígitos e cotas rigorosas — deles.
- Os bounds `mu ∈ [1/6, 17/100]`, `||B|| <= 23`, e a decomposição
  `1,511 + 21 + eps` — de `ALGEBRAIC_TAIL.md`, derivados dos bounds publicados.

**O `l^1` não é escolha nossa: é herança.** RT montaram tudo em `l^1` porque é
onde a prova de existência fecha. Nossa passagem a `l^2` é a mudança de
realização que o problema espectral pedia, e a inclusão gratuita
`l^1(w) ⊂ l^2(w)` é o que a torna legítima na direção da exclusão.

## 3. O que parece ser nosso

Não a técnica — os **números para este operador**, e o **mapa**:

| resultado | natureza |
|---|---|
| `R <= 6,66` contra os **414** em uso, fator 62 | constante nova para este operador |
| `R_J = mu · f(kappa2)`, com `f` monótona, `f -> -1`, zero em `kappa2 = 1,3204632636`, dimensão efetiva 18 | caracterização completa de uma constante do problema |
| razão dos pivôs `-> 1 - 1/kappa2 = 1/5` | consequência da estrutura semisseparável, específica deste `J_mu` |
| `sup/l^1 = 0,149` para o fundo publicado | medida sobre `RefA` |
| `||M_fundo||_2 <= 2,53` contra `12,40` da norma `l^1` | constante nova |
| dimensionamento da rota da resolvente: `~9800×`, e `366×` com todas as alavancas nos pisos | mapa de custo |
| sete alavancas fechadas com piso medido | resultado negativo, com constantes |
| contagem finita estável em `n0 ~ 16`, contra `n > 61` do certificado | medida da distância entre verdade e certificado |

## 4. O que NÃO reivindicamos

- **Não contamos os modos instáveis.** `T2` segue aberta, e depende de S3, S4,
  F1–F3 e C1–C4.
- **F3 tem agora `R <= 2,93`, certificado** (24/09, pelo símbolo pontual; antes 10,92 em 23/09 e a versão 6,05 de
  22/09 misturava convenções, `ALTERNATIVE_ROUTES.md` §5.25). O texto abaixo é de 22/09:
  temos uma cadeia com um elo certificado (`R_J`) e três
  somas finitas por refazer em racionais, mais três ressalvas de identificação.
- **Não provamos nada em Lean.** Tudo está em diagnóstico; a transcrição é
  trabalho não começado.
- **Não temos técnica nova**, como a §1 registra.

## 5. Contexto recente que vale acompanhar

- *Discovery of Unstable Singularities* (arXiv:2509.14185): descoberta
  sistemática de singularidades instáveis com precisão de máquina, com a
  observação empírica de que a `n`-ésima solução instável tem exatamente `n`
  modos instáveis. É outra família de problemas, mas a metodologia — precisão
  alta para viabilizar prova assistida depois — é a mesma aposta deste projeto.
- Linha ativa de enclosures/exclosures validados para problemas não
  auto-adjuntos, que é exatamente a caixa de ferramentas de que a contagem
  precisaria.

## 6. Como isto deveria ser apresentado

Não como "novo método", e sim como:

> **um estudo de viabilidade quantitativo do problema espectral de Choptuik na
> realização de Reiterer–Trubowitz**: sete rotas de barateamento fechadas com
> piso medido, um dimensionamento honesto da rota certificada, e uma melhora de
> 62× na fronteira de exclusão de F3 obtida ao trocar a realização `l^1` por
> `l^2` e usar o campo de valores.

O valor está nos negativos com constante e no mapa — não em ferramenta nova.
E o público que mais se beneficia é quem fosse tentar as mesmas sete rotas.
