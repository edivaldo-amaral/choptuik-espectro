"""Continua a produção após a validação do cache, com memória limitada."""
import sys,resource,runpy
sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
runpy.run_path('.codex-runs/2026-09-16-setorB-6x18/produzir_certificados.py',run_name='__main__')
