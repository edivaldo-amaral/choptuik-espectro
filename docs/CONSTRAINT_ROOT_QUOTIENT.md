# Constraints no espaço de raízes e quociente gauge

Esta nota resolve a álgebra das multiplicidades, não a equivalência
geométrica global. Fixe uma raiz isolada s0 na realização RT e escreva
L(s)=A+sI, K(s)=B+sI, C(s)=C0+sC1, M(s)=M0+sM1.
As identidades de SHARP_LINEARIZED_IDENTITY têm exatamente essa forma.

## Um mapa de constraints independente do parâmetro espectral

Comparando potências de s em KC=ML, obtemos

```
M1=C1,
M0=C0+B C1-C1 A,
B C0=M0 A.
```

Sobre um espaço de raízes generalizadas E de A, defina
`D=C0-C1 A`. Então `B D=D A`. A identidade é entendida nas realizações
com perda de raio: E é finito-dimensional, seus vetores estão no domínio
de todas as potências necessárias, e SHARP_DOMAIN dá o domínio de saída.
Assim ker(D|E) é A-invariante. Para uma cadeia de L,

```
A h_j=-s0 h_j-h_(j-1),
D h_j=C(s0)h_j+C1 h_(j-1).
```

Isso reproduz precisamente as constraints de jets da nota de exclusão
sharp. Cortar apenas ker C(s0) em todos os vetores de uma cadeia pode
dar uma dimensão errada. A matriz testada em test_root_constraints.py
exibe esse erro já em dimensão dois.

O teorema posto-nulidade dá

```
dim E_constraints = dim E - rank(D|E).
```

Não se substitui rank(D|E) pela multiplicidade sharp inteira: o mapa
pode não ser sobrejetivo no espaço de raízes sharp.

## Gauge positivo: linha invariante e multiplicidade

Pela ação RT e pela recuperação de raio, o gerador não nulo
g=(1+mu^-1 partial_tau+xi partial_xi)omega* está no domínio forte,
Ag=-mu g e Dg=0. Logo G=span(g) é uma linha invariante em E_constraints
na raiz s0=mu. Se esse gauge for admitido pelo domínio físico (cone móvel),
a álgebra do quociente é inequívoca:

```
dim(E_constraints/G) = dim E - rank(D|E) - 1.
```

Essa fórmula NÃO exige que g tenha complemento A-invariante.
Uma cadeia mista que começa em g pode sobreviver com comprimento menor
no quociente. Por exemplo, A=[[-mu,1],[0,-mu]], D=0, G=span(e1):
o quociente tem dimensão um, não zero. Portanto não se deve excluir
uma cadeia inteira simplesmente porque seu primeiro vetor é gauge.

Dentro do gauge residual analítico de coordenadas nulas, o gerador
espectral em f(u) é mu(1-u partial_u). Em s=mu, sua k-ésima potência
deslocada atua nos coeficientes por (-mu n)^k. Seu kernel, para todo
k>=1, é somente o termo constante. Não há cadeia generalizada positiva
adicional nesse gauge residual analítico. Logaritmos necessários para
tal cadeia não seriam germes analíticos em u=0.

## Hipóteses ainda abertas

A identificação física exige que toda perturbação admissível possa ser
posta no ansatz RT, com reconstrução e normalizações corretas, e que
o gauge global permitido seja o residual aqui classificado. Para cone
fixo, g move a fronteira e NÃO se subtrai essa linha automaticamente.
Além disso dim E e rank(D|E) ainda não foram certificados para o operador
infinito nos discos relevantes. Esta nota não atribui valores físicos
a essas dimensões nem valida a subtração de uma raiz finita observada.
