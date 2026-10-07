# Índice dos documentos

**Comece por [`PANORAMA.md`](PANORAMA.md)** (03/10/2026): é a visão de conjunto, com o resultado T2, o inventário, o que falta e os caminhos.
Os documentos `ESTADO_DA_ARTE.md`, `RELATORIO_SETEMBRO_2026.md` e `PROBLEMAS_ABERTOS.md` são anteriores a T2 e estão parcialmente desatualizados.

Quarenta e quatro documentos, agrupados por função. Os títulos são os dos
próprios arquivos. Se você está chegando agora, leia os quatro primeiros da
primeira seção e pare — o resto é consulta.

---

## Entrada

| documento | conteúdo |
|---|---|
| [RELATORIO_SETEMBRO_2026](RELATORIO_SETEMBRO_2026.md) | **Comece aqui.** O problema, o resultado, o mapa e o método, numa narrativa só |
| [PROBLEMAS_ABERTOS](PROBLEMAS_ABERTOS.md) | **Para atacar.** Sete aberturas, cada uma com o que bloqueia e o que conta como progresso |
| [NAO_FAZER](NAO_FAZER.md) | **Para não perder tempo.** Rotas fechadas: o que não tentar, e por qual constante |
| [ESTADO_DA_ARTE](ESTADO_DA_ARTE.md) | O que é nosso, o que é conhecido |
| [METODO](METODO.md) | O que funcionou e o que não funcionou, como método |
| [PRIORIDADES](PRIORIDADES.md) | Sequência de prioridades |
| [RETOMADA](RETOMADA.md) | Ledger cronológico, com as retratações preservadas |

## O problema e o que falta provar

| documento | conteúdo |
|---|---|
| [SPECTRAL_PROBLEM](SPECTRAL_PROBLEM.md) | Formulação do problema espectral — o enunciado preciso |
| [PROOF_OBLIGATIONS](PROOF_OBLIGATIONS.md) | Ledger das obrigações de prova (S*, C*, F*, T*) |
| [KERNEL_TRUST](KERNEL_TRUST.md) | Ledger de confiança do verificador formal |
| [WINDING_DATA](WINDING_DATA.md) | O que um produtor de certificados de winding precisa entregar |

## Auditorias

| documento | conteúdo |
|---|---|
| [ANALYTIC_AUDIT](ANALYTIC_AUDIT.md) | Auditoria da ponte entre o pencil e os certificados |
| [AUDITORIA_MATEMATICA_16SET](AUDITORIA_MATEMATICA_16SET.md) | Cauda, contração e Schur |
| [AUDITORIA_FISICA_16SET](AUDITORIA_FISICA_16SET.md) | Constraints, gauge e setor B |

## Rotas, custo e dimensionamento

| documento | conteúdo |
|---|---|
| [ALTERNATIVE_ROUTES](ALTERNATIVE_ROUTES.md) | Levantamento de rotas alternativas — o documento mais longo, crescido por acreção |
| [TT_CEILING](TT_CEILING.md) | Análise de teto do defeito TT; o dimensionamento de 9794× |
| [MAJORANT_ROUTE](MAJORANT_ROUTE.md) | A rota do majorante: estado de fechamento |
| [SHELL_PRECONDITIONER](SHELL_PRECONDITIONER.md) | Pré-condicionador de casca: rota fechada |
| [G_PARAMETRIX_FEASIBILITY](G_PARAMETRIX_FEASIBILITY.md) | Viabilidade de uma parametriz certificada na caixa G |
| [DAMPED_CROWN_PARAMETRIX](DAMPED_CROWN_PARAMETRIX.md) | Parametriz amortecida da coroa e auditoria por tiles |
| [UNIFORM_COUPLED_CONTRACTION](UNIFORM_COUPLED_CONTRACTION.md) | Critério acoplado uniforme, com homotopia explícita |
| [GUARDED_INFINITE_COUPLING](GUARDED_INFINITE_COUPLING.md) | Acoplamento com o exterior infinito usando faixa de separação |

## Cauda infinita e estimativas de exterior

| documento | conteúdo |
|---|---|
| [TAIL_BOUND](TAIL_BOUND.md) | O bound de cauda da obrigação C1 |
| [ALGEBRAIC_TAIL](ALGEBRAIC_TAIL.md) | Rota de cauda com pesos fixos — origem do `R = 414` que F3 substituiu |
| [COLUMNWISE_EXTERIOR](COLUMNWISE_EXTERIOR.md) | Majorante exterior por colunas |
| [COMPONENT_EXTERIOR](COMPONENT_EXTERIOR.md) | Exterior com estrutura entre componentes |
| [STRUCTURED_EXTERIOR](STRUCTURED_EXTERIOR.md) | Estrutura das colunas da inversa livre |
| [FREQUENCY_BLOCKS](FREQUENCY_BLOCKS.md) | Estimativas por blocos e o fechamento físico |
| [RADIUS_RECOVERY](RADIUS_RECOVERY.md) | Recuperação do raio nos espaços de raízes selecionados |

## A família *sharp*

| documento | conteúdo |
|---|---|
| [SHARP_SPECTRAL_LOCATOR](SHARP_SPECTRAL_LOCATOR.md) | Localizador espectral sharp e testemunho de violação |
| [SHARP_LINEARIZED_IDENTITY](SHARP_LINEARIZED_IDENTITY.md) | Identidade completa de propagação linearizada |
| [SHARP_EXTERIOR_BOUNDS](SHARP_EXTERIOR_BOUNDS.md) | Bounds sharp próprios, nos pesos de propagação |
| [SHARP_DOMAIN](SHARP_DOMAIN.md) | Domínio sharp por fechamento e perda controlada de raio |
| [SHARP_INTERMEDIATE_SHELL](SHARP_INTERMEDIATE_SHELL.md) | Faixa sharp incorporada e redução a um Schur finito |
| [SHARP_LOCAL_EXCLUSION](SHARP_LOCAL_EXCLUSION.md) | Exclusão sharp local: componente finita e extensão necessária |

## Constraints e gauge — as obrigações S3 e S4

| documento | conteúdo |
|---|---|
| [MODE_DISCRIMINATION](MODE_DISCRIMINATION.md) | As constraints discriminam as três raízes instáveis? |
| [CONSTRAINT_GAUGE_EQUIVALENCE](CONSTRAINT_GAUGE_EQUIVALENCE.md) | Equivalência física: lema e hipóteses ainda abertas |
| [CONSTRAINT_ROOT_QUOTIENT](CONSTRAINT_ROOT_QUOTIENT.md) | Constraints no espaço de raízes e quociente gauge |
| [CONSTRAINT_WITNESS_NORMS](CONSTRAINT_WITNESS_NORMS.md) | Normas do testemunho de constraints: hipótese H3 |
| [GAUGE_RT_ACTION](GAUGE_RT_ACTION.md) | Ação RT da translação nula e completude residual local |
| [GAUGE_DOMAIN_CLASSIFICATION](GAUGE_DOMAIN_CLASSIFICATION.md) | Gauge residual, cone fixo e seção de fase |
| [GLOBAL_NULL_GAUGE_TRANSPORT](GLOBAL_NULL_GAUGE_TRANSPORT.md) | Transporte global do gauge nulo com condição no cone |

## Geometria e o segundo setor

| documento | conteúdo |
|---|---|
| [SECTOR_B](SECTOR_B.md) | Setor B por redução ao setor A (obrigações S2/C4) |
| [GEOMETRIC_FLOQUET_RECONSTRUCTION](GEOMETRIC_FLOQUET_RECONSTRUCTION.md) | Reconstrução geométrica global no ansatz, para `Re s > 0` |

---

## Uma nota sobre a forma destes documentos

Vários cresceram por acreção, registrando o que se descobriu na ordem em que se
descobriu — `ALTERNATIVE_ROUTES.md` tem mais de vinte seções assim. Isso é
deliberado: reescrever para a conclusão final apagaria as **ressalvas**, e as
ressalvas são o que tornou possível pegar cinco erros nossos em cinco dias. O
`RELATORIO_SETEMBRO_2026.md` existe para dar a versão organizada sem destruir a
versão cronológica.
- [`AUDITORIA_S3A_TILES_24SET.md`](AUDITORIA_S3A_TILES_24SET.md) — auditoria do modo rigoroso dos tiles de S3a (Perron), 24/09
