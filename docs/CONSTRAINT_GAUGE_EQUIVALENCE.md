# Equivalência física: lema e hipóteses ainda abertas

Atualização: `SHARP_LINEARIZED_IDENTITY.md` agora constrói a identidade
off-shell e os operadores K(s), M(s), incluindo o defeito do fundo
aproximado. A lacuna restante dessa etapa é a extensão aos domínios
fechados e o controle espectral de K, não mais a fórmula diferencial.

Atualização de domínio: `SHARP_DOMAIN.md` prova a passagem
C(s):D_L(a,b)->D_K(a',b') para raios estritamente menores. O controle
espectral de K deve agora ser feito NESSES pesos. `GAUGE_DOMAIN_CLASSIFICATION.md`
classifica uma classe de gauge residual e verifica quantitativamente a
seção de fase neutra; a equivalência física completa permanece aberta.

## 1. O que a identidade de propagação permite provar

Fonte primária: [Reiterer–Trubowitz, equações (16), (24)–(26)](https://arxiv.org/pdf/1203.3766).
O artigo fornece o sistema homogêneo sharp quando as equações selecionadas
se anulam. Isso não é, sozinho, uma identidade off-shell linearizada.
O código `MultiField_SharpIdentities_complete_alloc` é um ponto de partida
para extrair os termos de equações selecionadas da identidade completa.

Denote por F as equações selecionadas e por C as constraints, no fundo
exato u*. Uma demonstração suficiente deve estabelecer, ANTES de impor F=0,

```
K(u) C(u) = M(u) F(u).
```

Com F(u*)=C(u*)=0, a regra do produto elimina os termos que derivam K e M:
`K(u*) DC(u*)h = M(u*) DF(u*)h`.
Conjugando todos os operadores diferenciais pelo mesmo fator de Floquet,
obtém-se `K(s) C(s) = M(s) L(s)`. Esta é uma derivação condicional exata;
a extração de M é feita na nota acima; seus domínios continuam abertos.
Não é legítimo diferenciar apenas ao longo de soluções e concluir a
identidade para todo vetor do kernel linearizado: nem todo vetor é
necessariamente tangente a uma família não linear.

Se a identidade vale e K(s) é injetivo, então
`L(s)h=0 => C(s)h=0`. A recíproca relevante é imediata: uma solução
das equações completas satisfaz as selecionadas. Logo os kernels coincidem
nesse s. Para identificar TODOS os modos dentro de um contorno, K precisa
ser injetivo nos pontos característicos interiores, não apenas no contorno.
Invertibilidade de K na fronteira não exclui modos violadores dentro dela.

Para multiplicidades algébricas, suponha K(s) holomorfo e invertível numa
vizinhança do ponto. Se um polinômio de cadeia h(s) satisfaz
`L(s)h(s)=O((s-s*)^k)`, a identidade e K(s)^-1 implicam
`C(s)h(s)=O((s-s*)^k)`. Assim toda a cadeia satisfaz as constraints à ordem
necessária. Isso fecha a etapa algébrica sob as hipóteses especificadas.
Não autoriza subtrair o winding de K do de L: tal fórmula exigiria ainda
uma decomposição exata, com sobrejetividade e levantamentos apropriados.

## 2. Lema do quociente por gauge

Fixe um s e seja Z=ker L(s) intersectado com ker C(s), no domínio físico.
Seja G o operador das transformações infinitesimais admissíveis, com
imagem contida em Z. Para uma condição de gauge linear B, suponha que
`BG: A -> W` seja bijetivo, após remover estabilizadores dos parâmetros A.
Então

```
Pi h = h - G (BG)^-1 B h
```

pertence a Z intersectado com ker B, é constante nas classes módulo im G,
e tem a mesma classe que h. Se dois representantes nessa seção diferem
por Ga, então BGa=0 e a=0. Isso prova a bijeção
`Z/im G ≅ Z intersectado com ker B`.
Em espaços de Banach, exigir operadores limitados e inversa limitada;
para cadeias espectrais, exigir versões holomorfas locais compatíveis com
o pencil. A bijeção pontual não basta para transportar multiplicidades.

No sistema atual ainda faltam: construir G em variáveis RT, provar que
preserva regularidade no centro/cone e paridades, classificar seu kernel,
e verificar a condição BG no domínio escolhido. O argumento geométrico
da seção 1.4 do artigo é explicitamente motivacional, não um teorema de
equivalência para perturbações de Floquet.

A translação de fase produz o candidato neutro derivada_tau u*, pela
invariância das equações, mas não identifica o modo positivo próximo de mu.
Mover o tempo de acumulação pode mover o cone: não se pode declará-lo
gauge admissível de um domínio fixo sem construir essa identificação.

## 3. Estado verificável

Os dois lemas acima estão demonstrados algebricamente sob hipóteses
explícitas. As hipóteses analíticas RT/Floquet e a classificação de gauge
não foram demonstradas nem formalizadas em Lean. Portanto a equivalência
física e a contagem de um único modo físico continuam abertas.
