# Revisão independente do certificado S3a

Escopo: leitura dos §§10.4–10.12 e dos arquivos indicados, com conferência algébrica das fórmulas e inspeções em caixas pequenas. Não executei a cobertura inteira nem alterei arquivos do projeto.

## Resumo

| pergunta | veredicto | evidência (arquivo:linha) |
|---|---|---|
| 1. Estrutura (A2) | OK | `scripts/tile_perron.py:1–22`, `459–486`; `signed_operator.py:69–83` |
| 2. Lema de Perron (A3) | OK | `scripts/tile_perron.py:22–27`, `174–185` |
| 3. Blocos (A4) | OK | `scripts/tile_perron.py:187–443`, `459–486` |
| 4. Rump (D1) | OK | `scripts/verified_norms.py:1–158`; Rump, Teorema 2.3 e Corolário 2.4 |
| 5. Erros (D3) | OK | `scripts/tile_perron.py:93–116`, `223–443`, `459–486` |
| 6. Forma fechada (C2) | OK | `scripts/fourier_tail_bound.py:119–181`; `scripts/tile_perron.py:164–171` |
| 7. Cobertura (E1) | OK | `scripts/check_cover.py:23–46`, `78–129`; `scripts/tile_cover.py:45–57` |

## 1. Estrutura (A2) — OK

`Jrad` é diagonal em `m` e triangular superior em `n` (`signed_operator.py:44–63`), enquanto o fundo só desloca `m` por uma banda finita. Portanto as caixas são preservadas pelo operador livre e a ordenação `[Z,T1,T2,far]` dá uma matriz triangular superior por blocos.

Em cada coluna, `Q0(s_l)` é o inverso do bloco livre na caixa correspondente. Como os blocos diagonais são inversíveis, a soma de colunas `P = Σ Q0(s_l)P_l` é bijetiva.

A expansão usada no código é exatamente
`(K_livre(s_l)+(s-s_l)I+B)Q0(s_l)P_l = P_l +(B+s-s_l)Q0(s_l)P_l`. Não há erro de sinal ou de ordem nessa identidade.

## 2. Lema de Perron (A3) — OK

É necessário, além de `v_i>0` e `Nv <= theta v`, que `N` seja não negativa e que cada entrada majorada seja uma norma de bloco subordinada. Essas hipóteses são respeitadas: `N` é montada somente com cotas não negativas (`tile_perron.py:459–486`) e os blocos são medidos na norma espectral ponderada.

Para `x=(x_i)`, definindo `||x||_v=max_i ||x_i||/v_i`, a desigualdade de blocos dá `||Ex||_v <= theta||x||_v`. Como `theta<1`, `I-E=AH` é injetivo; como `A` é invertível, `H`, e então `K`, também são injetivos. A condição é suficiente, não exige que `v` seja exatamente autovetor de Perron.

## 3. Completude e correção dos blocos (A4) — OK

As linhas de `parte_Z` incluem o inverso denso `V` quando necessário (`tile_perron.py:224–293`); as linhas `b1,b2` são cotas dos blocos que saem de Z. `parte_T1` trata a casca finita, inclusive `T1<-Z` e `T2<-T1` (`295–361`). `parte_T2` acumula os Grams das saídas em Z e T1 e trata `T2<-T2` por deslocamento de Fourier (`362–443`).

O ponto potencialmente perigoso, `far`, está de fato presente em todas as linhas relevantes de `montar`: `far<-Z` e `far<-T1` recebem `dB`; `Z<-far` recebe `V(nB dZ+dB Qf)`; e as colunas `T2,far` recebem `Qf`/`beta` (`tile_perron.py:478–485`). Assim, a perturbação global não é indevidamente tratada como banda local.

Os termos `(s-s_l)` aparecem precisamente nas três matrizes de variação `NZm,N1m,N2m` (`474–476`). A inflação `ARRED` também é aplicada às quatro matrizes antes do JSON (`481–486`), e `check_cover.py:31–33` repete uma inflação conservadora ao reverificar.

Não encontrei bloco omitido: os blocos que parecem zeros por triangularidade são zeros do operador livre; os termos de resto global são recolocados por `dB`.

## 4. Critério de Rump (D1) — OK

O artigo citado permite matriz real simétrica ou complexa hermitiana. O código realifica a hermitiana (`verified_norms.py:72–78`); a realificação usada pode trocar o sinal da parte imaginária conforme a convenção, mas continua semelhante à realificação padrão e preserva positividade.

Em `cota_gram`, a parte hermitiana é formada antes do Cholesky (`verified_norms.py:125–139`). O erro declarado de formar `Y*Y`, a simetrização, o arredondamento de `t²-G_ii` e o termo de underflow entram em `extra` e em `_margem_rump` (`50–68`, `140–145`). A dimensão da realificação é usada porque `m=2*n` (`121–124`), e `dpotrf` recebe uma matriz simétrica real (`62–69`).

Não confundi a condição `c >= ||Delta(A)||` com a margem de erro do produto: ambas são somadas. A cota I do artigo é precisamente `gamma_{n+1}/(1-gamma_{n+1}) tr(A)+n M eta`, que é a expressão em `_margem_rump`.

## 5. Contabilidade de erros (D3) — OK

`err_mm` usa a cota Frobenius para o erro de produto (`tile_perron.py:93–96`), `err_gram` usa a dimensão de soma (`106–108`), e `Rw_cert` aplica a identidade de inverso por resíduo `||R|| ||E||/(1-||E||)` (`110–120`).

Nas partes Z/T1/T2, os erros são propagados por produtos, Grams e pelos erros de `Q0`; o sanduíche `V S V*` inclui tanto os erros de formação como o erro do Gram (`montar`, `465–475`). As linhas que contêm `V` multiplicam explicitamente as cotas correspondentes; não há uso acidental de uma cota sem o fator `V`.

`FATOR_PESO` é aplicado em `N2` e `L2` (`98–104`), enquanto `dB` contém resto de banda, arredondamento dos campos e erro de RT (`Ctx.__init__`, `187–215`). A montagem inclui `dB` também nos blocos envolvendo `far`, como descrito acima.

## 6. Cotas da forma fechada (C2) — OK

Para `J` triangular, com diagonal `d_n=mu(n+1)+sigma` e entrada superior `v_n=2mu n`, a entrada da inversa para `j<c` é `(1/d_j) Π_{k=j}^{c-1}(-v_k/d_k)`.

`cauda_Rinv` implementa exatamente esse produto, incluindo o fator de pesos `sqrt(omega(j)/omega(c))`, e usa soma máxima de colunas e linhas para Schur (`fourier_tail_bound.py:137–171`). O cálculo em log evita o antigo `inf/inf`.

Para `Re(s)>=0`, `|1-v_k/d_k|<=1`; para Fourier distante, `|d_k| >= |m|/2-|Im(s_0)|` e `v_k/|d_k|<2`. Isso produz a fórmula de `Rinv_uniforme` (`175–181`). A fórmula de `descida` usa a mesma entrada fechada, com soma harmônica nas colunas e geométrica nas linhas (`tile_perron.py:164–171`). Os fatores `1.0001` são aplicados depois das cotas analíticas e as hipóteses `Re(s)>=0`/`h>0` são verificadas ou afirmadas pelos chamadores.

## 7. Cobertura (E1) — OK

`check_cover.py` não aceita simplesmente a marca de sucesso: reconstrói `N=N0+rZ NZ+r1 N1+r2 N2`, verifica `v>0`, `Nv< v` em racionais e as duas contenções de centros (`23–46`). Depois testa os quatro vértices de cada célula. Como discos e células quadradas são convexos, esse teste cobre a célula inteira.

O disco de S3b é tratado como aberto (`87–88`); portanto pontos da sua circunferência não são descartados. A quadtree cobre a região fechada e os tiles, incluindo a margem dos raios. A conjugação está corretamente condicionada ao espaço simétrico: a norma de sinal isolada não é invariante, mas o subespaço simétrico tem os pesos pareados e fica controlado pelo certificado no semiplano superior; conjugação leva esse espaço ao semiplano inferior.

Os limites usados pelo verificador são exatamente `Re s in [0,1.765]` e `Im s in [0,0.25]` (`check_cover.py:70–79`). A extensão ao semiplano inferior é por conjugação, não por uma suposição de cobertura geométrica.

## Lista final

Não encontrei ERROS nem LACUNAS que possam ser sustentados por contraexemplo ou por uma cota violada. Os pontos que permanecem como dependências externas do certificado — a cota do símbolo, o erro de RT e a escolha final do disco de S3b — estão fora das sete alegações examinadas ou são entradas declaradas delas; não foram silenciosamente reprovados como se fossem demonstrados por estes scripts.
