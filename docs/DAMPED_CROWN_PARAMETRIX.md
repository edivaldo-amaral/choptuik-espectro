# Parametriz amortecida da coroa e auditoria por tiles

## Objeto exato: não confundir coroa com exterior

Aqui P é a caixa 4×12, G a caixa 6×18 e S=G-P, de dimensão 195.
A matriz estudada é H_SS(s), o bloco de coroa de
`H(s)=(A+sI)Q_ref`, com Q_ref=J_(mu_RefA)^-1. Os pesos das componentes
são (916,4096,243,1233), além dos pesos RT. Não estamos invertendo
THT com T=I-G, nem ignorando os acoplamentos P/S/T.

## Por que amortecer

A parametriz polinomial usa

```
E_j=I-omega H_SS(s_j),
V_j=omega sum_(k=0..31) E_j^k.
```

Em aritmética exata, I-V_j H_SS(s_j)=E_j^32. Na implementação, V_j é
proposto numericamente, arredondado para diádicos e seu defeito é
recalculado racionalmente: não confiamos nos autovalores numéricos.

O diagnóstico prévio mostrou que omega=1 não serve perto de s=1:
o raio espectral NUMÉRICO de E era 1.221 e a norma da potência 32
aproximadamente 14765. Com omega=1/2, a mesma norma numérica cai para
1.081e-5. No ponto 1/8+i/4, a potência 32 fica aproximadamente 0.0633.
Esses números motivaram a proposta; não são os bounds certificados finais.

## Como a aritmética é validada

1. Arredondamos uma inversa candidata de J_ref para a grade 2^-48.
   O resíduo I-J_ref Q_hat é calculado exatamente em inteiros/racionais.
   A norma ponderada confere ||Q_hat||/(1-||resíduo||)<=108, de modo que
   ||Q_ref-Q_hat||<=108||resíduo||. Isso é uma validação matricial finita,
   independente de uma hipótese de EDP para o número 108.
2. A Q_hat e Q_hat são conjugados pelos pesos e arredondados para a grade.
   Cada erro de arredondamento real por coluna é no máximo 195/(2*2^48).
   O erro de H(s) inclui (||A||+|s|)||Q_ref-Q_hat||.
3. V_j é proposto pelo polinômio amortecido e arredondado. O produto
   I-V_j H_hat(s_j) é recalculado exatamente, com partes real/imaginária
   inteiras. Somar |Re|+|Im| fornece uma cota superior para a norma complexa.
4. Para um V_j FIXO no tile, calculamos bounds racionais a_j e v_j de
   ||I-V_j H_SS(s_j)|| e ||V_j (Q_ref)_SS||. No disco |s-s_j|<=r_j,
   o defeito é no máximo a_j+r_j v_j. O produtor escolhe r_j diádico
   estritamente dentro da margem de contração.

O produto inteiro rápido usa limbs assinados de 15 bits. Cada produto
matricial int64 intermediário é limitado por n(2^15-1)^2<2^63, condição
verificada antes da chamada. A recombinação usa inteiros arbitrários.
Não há arredondamento em ponto flutuante na verificação dos produtos.

O percurso liga pontos consecutivos por segmentos de comprimento no máximo
r_j. Ele cobre a borda superior e as metades superiores das bordas verticais
do retângulo 1/8<=Re s<=1, |Im s|<=1/4. Como A e Q_ref são reais, os tiles
conjugados cobrem a metade inferior. O relatório só marca cobertura após
alcançar todos os extremos previstos sem falha.

## Erro do fundo, incluindo mu

O precondicionador Q_ref fica FIXO ao passar ao fundo exato. Com
epsilon_mu=|mu_RefAplusB-mu_RefA|+2^-277, o cabeçalho dos arquivos dá
epsilon_mu≈3.88112e-11 <2^-34. Além disso epsilon_omega=2^-25+2^-277.

Como J_ref=partial_tau+mu_ref Jslash, segue
`||Jslash Q_ref||<=(1+||partial_tau Q_ref||)/mu_ref<=114`:
a fatorização livre dá ||partial_tau Q_ref||<=18 e mu_ref>=1/6.
Esse operador relativo não é tratado como uma multiplicação pelo fundo.
Com g1_eta=||S Gamma1|| majorado pela matriz de componentes,

```
||H_true(s)-H_ref(s)|| <=
  114 epsilon_mu +108[g1_eta epsilon_mu +6(4096/243)epsilon_omega]
  < 0.0003257.
```

O bound é independente de s. `reference_error_budget.py` conserva os
racionais completos; sua validade depende da bola publicada de RT.
O termo relativo de mu não decai com o corte exterior: não deve ser
absorvido artificialmente numa cauda que tende a zero.

O exportador selecionado usa RefA cortado em 2G. Isso não altera o bloco
G H_ref G: Q_ref preserva G, e uma multiplicação Fourier–Chebyshev entre
entrada e saída abaixo de G só pode usar coeficientes do fundo com índices
abaixo de 2G. S Xi Gamma2 reduz-se à multiplicação e à projeção de paridade,
sem uma divisão que traria graus arbitrariamente altos para baixo.
Os termos Gamma1 e principal não dependem desse corte do fundo.
Essa observação justifica aplicar o erro acima ao bloco S exportado.

Com --include-background, o produtor acrescenta
`||V_j|| * ||H_true-H_ref||` ao defeito do centro. A derivada em s continua
Q_ref, portanto seu bound não muda. A ponte está habilitada somente para
a matriz 6×18 auditada, identificada por hash.

## Alcance e reprodução

Mesmo uma cobertura aceita prova somente a invertibilidade uniforme do
bloco finito H_SS. Restam os acoplamentos com P e com o exterior infinito,
a parametriz exterior, as constraints e o quociente por gauge. Não há
contagem física ou transferência de winding por essa etapa isolada.
Os cálculos são auditorias racionais Python, não certificados Lean.

`check_crown_cover.py` verifica separadamente as desigualdades escalares e
o percurso: continuidade, orientação, ausência de saltos nos cantos,
comprimento de cada segmento e conjugação declarada. Não recalcula os
resíduos matriciais nem prova que os dados representam um operador real.
A cobertura de referência concluída tem 51 tiles e defeito máximo
17525650539575/35184372088832 < 0.499.

A execução COM erro do fundo também terminou: 51 tiles, defeito máximo
17674845549527/35184372088832 < 0.503. Em cada tile,
`||H_SS^-1|| <= ||V_j||/(1-a_j-r_j v_j)`; o máximo racional é

```
kappa = 669414234798400/18478486838353 < 36.227.
```

O máximo dessa cota ocorre no tile centrado em 1/8+i/4. A validação
finita de Q_ref forneceu norma menor que 5.942 (o produtor usou o
majorante conservador 108 para o erro). O relatório completo é
`build/spectrum/crown-contour-6x18-background.json`; a versão sem
incerteza é `build/spectrum/crown-contour-6x18.json`.
A conclusão para o fundo exato é condicional à bola RT e à ponte
analítica descrita acima, não uma recertificação da existência do fundo.

```bash
.venv/bin/python scripts/check_crown_cover.py \
  build/spectrum/crown-contour-6x18-background.json
```

## Uso correto na eliminação de Schur

Escreva R=P+T, mantendo S como a coroa FINITA acima. Se
`kappa >= sup_contorno ||H_SS^-1||`, o operador reduzido sobre R é

```
Hred_XY = H_XY - H_XS H_SS^-1 H_SY,    X,Y em {P,T}.
```

Assim, na mesma norma e no mesmo contorno,

```
||Hred_TT-I|| <= ||H_TT-I|| + kappa ||H_TS|| ||H_ST||,
||Hred_PT||   <= ||H_PT||   + kappa ||H_PS|| ||H_ST||,
||Hred_TP||   <= ||H_TP||   + kappa ||H_TS|| ||H_SP||.
```

Para uma parametriz V_P, a correção no defeito PP é no máximo
`kappa ||V_P H_PS|| ||H_SP||`. Esses bounds mostram exatamente onde
entra o novo resultado, mas não são números já certificados para T.
O bound antigo TT>1 não é melhorado por somar correções positivas:
será necessário capturar cancelamentos no complemento ou construir
uma parametriz exterior mais forte.

Além disso, invertibilidade de H_SS SOMENTE NO CONTORNO não implica
winding zero de det H_SS. Na fatorização de Schur, sua contribuição
espectral deve ser contabilizada separadamente; descartá-la requer
um certificado adicional. A ausência de raízes observada numericamente
no pencil da coroa não fornece esse certificado.

```bash
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/certify_crown_tiles.py \
  build/spectrum/rt-A-6x18.dat --cover-contour --include-background \
  --out build/spectrum/crown-contour-6x18-background.json
```
