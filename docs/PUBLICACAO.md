# Publicação: checklist (preparado em 03/10/2026, nada publicado)

**Proposta:**
- código, documentos e certificados pequenos num repositório GitHub;
- dados pesados num arquivo do Zenodo, com DOI citável no artigo;
- começar **privado** e abrir quando o artigo for ao arXiv.

**Publicar é uma ação externa:** só com a decisão explícita do autor.

## Decisões do autor (07/10/2026)

1. **GitHub privado** agora; **público quando o artigo for ao arXiv**.
2. **Histórico limpo:** repositório novo a partir de `scripts/exporta_publico.py`.
3. **Licenças:** MIT para o código (`LICENSE`); CC BY 4.0 para documentos, certificados e dados
   (`LICENSE-DOCS.md`). Citação em `CITATION.cff`.
4. **Zenodo** para os dados pesados, como rascunho com DOI pré-reservado até o artigo. Rascunho criado em
   07/10: id 23222363, DOI **10.5281/zenodo.23222363** (`build/reproducao/zenodo.json`).
5. **Anonimização:** como sugerido ("PC 3" no texto; os rótulos `joao` dentro dos certificados ficam).
6. **Trilha de revisões:** como sugerido (TAREFA, RELATORIO e scripts; os logs brutos e os snapshots vão
   para o Zenodo).

## O que já está pronto

- o inventário alegação → artefato (`docs/REPRODUCAO_MAPA.md`, `scripts/inventario_certificados.py`; o
  host remoto vem de `CHOPTUIK_HOST_REMOTO`, não fica fixo no código);
- o verificador de nível 1 (`scripts/verifica_T2.py`) e o guia (`docs/REPRODUCAO.md`);
- os certificados pequenos todos no git; nenhuma credencial no repositório (auditado em 03/10).

## O que falta antes de publicar

- [ ] Rodada completa do nível 1, com os 22 discos (em curso no laboratorio2, 03/10).
- [ ] Janelas no nível 1, com uma rodada longa.
- [ ] README de entrada (o que é, T2 em uma página, como verificar), CITATION.cff e LICENSE.
- [ ] Script de exportação: copia a árvore limpa, aplica a redação de IPs, nomes e caminhos, e confere com
      um grep que nada sensível sobrou.
- [ ] Pacote do Zenodo:
  - A.npy, os 13 F_*.npy, T.npy e U.npy, com manifesto sha256;
  - os logs brutos e os snapshots das revisões.
- [ ] Atualizar os documentos desatualizados (`ESTADO_DA_ARTE`, `RELATORIO_SETEMBRO`, `PROBLEMAS_ABERTOS`)
      ou movê-los para `docs/historico/`.
