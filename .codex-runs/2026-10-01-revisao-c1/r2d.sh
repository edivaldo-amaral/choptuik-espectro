cd ~/revisao_c1
until [ -f r2_fim.txt ]; do sleep 20; done
mkdir -p d10 dm25
ln -sf ~/Choptuik/build/rouche_L/A.npy d10/A.npy; ln -sf ~/Choptuik/build/certF_extra/F_+1_0000+0_0000i.npy d10/F_+1_0000+0_0000i.npy
ln -sf ~/Choptuik/build/rouche_L/A.npy dm25/A.npy; ln -sf ~/Choptuik/build/rouche_L/F_+0_0000-0_2500i.npy dm25/F_+0_0000-0_2500i.npy
export CHOPTUIK_EPS_FUNDO_L=4.7e-08 OPENBLAS_NUM_THREADS=4
PY=~/Choptuik/.venv/bin/python
$PY r2_lu_densa.py --dir d10 --c 1 --F F_+1_0000+0_0000i.npy --rig rig_+1_0000+0_0000i_rig.json --disco disco_1_0.json --ip 0 --saida r2d_resultado.json > r2d_saida.txt 2>&1
$PY r2_lu_densa.py --dir dm25 --c=-0.25j --F F_+0_0000-0_2500i.npy --rig rig_+0_0000-0_2500i_rig.json --disco disco_-0_25.json --ip 0 --saida r2e_resultado.json > r2e_saida.txt 2>&1
echo fim > r2d_fim.txt
