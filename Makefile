LEAN_LAKE ?= <home>/.elan/bin/lake

.PHONY: setup formal test rt-audit spectrum-small verify-winding verify-reduced-b verify-reduced-b-6x18 clean-results

setup:
	python3 -m venv .venv
	.venv/bin/pip install -r requirements.txt

formal:
	./scripts/check_formal.sh

spectrum-small:
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-1} ./scripts/generate_spectrum.sh 4 12

verify-winding:
	.venv/bin/python scripts/verify_winding_certificate.py examples/winding-regression-A.txt 1
	.venv/bin/python scripts/verify_winding_certificate.py examples/winding-regression-B.txt 0
	.venv/bin/python scripts/make_winding_certificate.py examples/winding-interval-square-plus1-nodes.txt >/dev/null
	.venv/bin/python scripts/make_winding_certificate.py examples/winding-interval-loop-zero-nodes.txt >/dev/null
	.venv/bin/python scripts/verify_winding_certificate.py examples/winding-interval-square-plus1-increments.txt 1
	.venv/bin/python scripts/verify_winding_certificate.py examples/winding-interval-loop-zero-increments.txt 0
	.venv/bin/python scripts/make_winding_certificate.py examples/winding-interval-dip-zero-nodes.txt >/dev/null
	.venv/bin/python scripts/verify_winding_certificate.py examples/winding-interval-dip-zero-increments.txt 0
	@tmp=$$(mktemp -t choptuik-gerado-XXXXXX.lean); \
	{ printf 'import ChoptuikFormal.WindingInterval\nset_option maxRecDepth 100000\n'; \
	  for base in square-plus1 loop-zero dip-zero; do \
	    .venv/bin/python scripts/make_winding_certificate.py \
	      examples/winding-interval-$$base-nodes.txt; \
	  done; } > $$tmp; \
	(cd formal && $(LEAN_LAKE) env lean $$tmp); status=$$?; rm -f $$tmp; \
	if [ $$status -ne 0 ]; then \
	  printf 'ERRO: o Lean emitido pelo produtor nao typechecka contra o verificador.\n' >&2; \
	  exit 1; \
	fi; \
	printf 'OK: o Lean emitido pelo produtor typechecka contra o verificador.\n'

rt-audit:
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-1} ./scripts/reproduce_rt.sh

verify-reduced-b:
	.venv/bin/python scripts/verify_winding_certificate.py examples/contour-A-4x12-reduced-B-increments.txt 0
	@tmp=$$(mktemp -t choptuik-setor-b-XXXXXX.lean); \
	trap 'rm -f "$$tmp"' EXIT; \
	.venv/bin/python scripts/make_winding_certificate.py \
	  examples/contour-A-4x12-reduced-B-nodes.txt --coarsen --standalone --quiet \
	  --lean "$$tmp" && (cd formal && $(LEAN_LAKE) env lean "$$tmp")

verify-reduced-b-6x18:
	bash scripts/verify_reduced_b_6x18.sh

test: formal
	.venv/bin/python -m py_compile scripts/*.py
	.venv/bin/python -m unittest discover -s scripts -p 'test_*.py'
	bash -n scripts/*.sh

clean-results:
	@printf '%s\n' 'Remova manualmente build/ se desejar apagar os resultados reproduzíveis.'
