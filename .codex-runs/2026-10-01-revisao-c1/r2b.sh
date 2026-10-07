cd ~/revisao_c1
export CHOPTUIK_EPS_FUNDO_L=4.7e-08 OPENBLAS_NUM_THREADS=4
PY=~/Choptuik/.venv/bin/python
$PY r2_lu_densa.py --dir ~/Choptuik/build/rouche_L --c 0j --F F_+0_0000+0_0000i.npy --rig rig_+0_0000+0_0000i_rig.json --disco disco_0.json --ip 6,10 --extra '6:0j' --saida r2b_resultado.json > r2b_saida.txt 2>&1
$PY r2_lu_densa.py --dir ~/Choptuik/build/rouche_L --rig rig_+0_0000+0_2500i_rig.json --disco disco_0_25.json --ip 8 --saida r2c_resultado.json > r2c_saida.txt 2>&1
echo fim > r2_fim.txt
