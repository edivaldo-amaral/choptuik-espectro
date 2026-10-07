# Identidade completa de propagação linearizada

Continuação: `SHARP_DOMAIN.md` demonstra sua extensão aos domínios fechados
com perda de raio analítico; não se afirma a extensão nos mesmos pesos.

## Convenções e reconstrução do resíduo

P nesta nota é a reflexão espacial `Pf(tau,xi)=f(tau,-xi)`, NÃO uma caixa
de frequências. Ponha O=(1-P)/2, E=(1+P)/2, X=multiplicação por xi e
R_x f=(f-f|_{xi=0})/xi, a divisão regularizada de RT.
Trabalhamos com mu fixo e uma perturbação `exp(s*tau) h`, h1 ímpar em xi.

Se q=OmegaPlus(mu,w) tem componentes (1,2,2sharp,3,4), defina
`F=S q`, `c=Ssharp q=(O q1,-P q2sharp)` e e=-Pq, componente a componente.
Para resíduos das equações RT (não para um vetor arbitrário de cinco
funções), as componentes selecionadas são divisíveis por xi. Portanto

```
e = I c + R F,
I c = (c1,0,c2,0,0),
R f = (X P f1, X P f2,0,X P f3,X P f4).
```

Prova: q1=X F1+c1, q2=X F2, q3=X F3, q4=X F4,
q2sharp=-P c2. Aplicar -P e usar PX=-XP e Pc1=-c1 dá a fórmula.
A divisibilidade segue da forma de OmegaPlus: os termos principal e
quadrático têm fator X, e as componentes selecionadas de Gamma1 são
ímpares. A mesma demonstração vale para DOmegaPlus[h] com o deslocamento s.
Isso resolve, nessa classe de resíduos, a possível perda do valor em xi=0
pela divisão regularizada. Não afirmamos R_x invertível no espaço inteiro.

## Operador completo, antes de impor as equações

Defina J_s=partial_tau+s+mu*((1+xi)*partial_xi+1).
Denote por G(w,e) a aplicação bilinear de cinco resíduos para duas saídas
implementada por `GAMMA2_SHARP_macro.c`, versão NÃO restrita.
Ela é a soma das dez parcelas com as partes par/ímpar de cada resíduo;
`sharp_propagation.py` lê essas vinte entradas diretamente do fonte RT.

O operador completo é

```
T_s(w)e = X ( O J_s e1, J_s e2sharp - P J_s e2 )
          + mu ( E e1, e1+e2+P e2-2 P e3 ) + X G(w,e).
```

A identidade não linear é `T_0(w) e(w)=0`, inclusive quando e(w) não zera.
É a identidade off-shell da seção 2.3 de
[Reiterer–Trubowitz](https://arxiv.org/pdf/1203.3766), implementada em
`MultiField_SharpIdentities_complete_alloc`; o teste original de RT passa
`-P OmegaPlus`, e não `OmegaPlus`, a essa rotina.

## Linearização, incluindo o defeito do fundo aproximado

Escreva e'=exp(-s*tau) De(w)[exp(s*tau)h]. A única dependência de T em w
está em G. Logo a regra do produto dá a identidade completa

```
T_s(w)e' + X G(h,e(w)) = 0.
```

O fator s entra em TODAS as derivadas de e', não nas do fundo w. Essa
fórmula é válida em fundo aproximado: omitir G(h,e(w)) ali seria um erro.
No fundo exato w*, e(w*)=0, e o termo adicional desaparece.

Agora ponha `F'=L(s)h`, `c'=C(s)h` e use a reconstrução anterior. Defina

```
K(s)c = R_x T_s(w) I c,
M(s)f = -R_x T_s(w) R f.
```

Como R_x X=I, obtemos

```
K(s) C(s)h = M(s) L(s)h - G(h,e(w)).
```

Em particular, **no fundo exato**:
`K(s) C(s) = M(s) L(s)`.
Não se diferenciou apenas uma família de soluções: a identidade usada
antes da diferenciação vale fora do conjunto de soluções.

## Fórmulas explícitas de K e da fonte M

Escreva G_R=G(w,Rf), G_I=G(w,Ic). Usando c1 ímpar e f1 ímpar:

```
(K(s)c)1 = [partial_tau+s+mu*(xi*partial_xi+1)] c1 + (G_I)1,
(K(s)c)2 = J_s c2 + mu R_x c1 + (G_I)2,

(M(s)f)1 = -mu*(xi*partial_xi+2) P f1 - (G_R)1,
(M(s)f)2 = P J_s(X P f2)
            - mu*(P f1+P f2-f2+2 f3) - (G_R)2.
```

K tem exatamente a parte principal diag(K_mu,J_mu)+sI do sistema sharp,
e G_I é a versão restrita Gamma2sharp. O primeiro componente de M não
contém s explicitamente; o segundo contém, via J_s. Assim os termos de
fonte estão identificados, e não escondidos sob a hipótese F=0.

## Verificação e alcance

`test_sharp_propagation.py` usa polinômios racionais em tau, xi e s.
Verifica a identidade não linear em fundo NÃO solução, a reconstrução,
a identidade linearizada com s simbólico e correção não nula, o sinal de M
e o caso de fundo zero exato. As tabelas bilineares são as dos macros RT;
derivadas, reflexão e divisão são implementadas em monômios, independentemente
da representação Fourier–Chebyshev do código C. Esses testes são regressões
exatas em exemplos, não prova simbólica universal nem certificado Lean.

As manipulações acima provam a identidade diferencial local a partir da
identidade completa RT. Sua aplicação a funções analíticas com paridades
admissíveis é legítima localmente; o fator de Floquet não altera a reflexão
espacial. Falta provar a extensão aos domínios fechados ponderados usados
na contagem: M contém derivadas, e não foi provado limitado no espaço
ambiente sem norma de gráfico. Também faltam a exclusão de kernel de K
nos pontos característicos relevantes e a classificação de gauge.
Logo a identidade algébrica está construída; a equivalência espectral
física NÃO está concluída.
