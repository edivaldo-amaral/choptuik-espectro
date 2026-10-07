import sys, time
from fractions import Fraction as Q
sys.path.insert(0,'scripts'); sys.path.insert(0,'<scratch>')
from sweep_V import norma_V
print("||V|| em s=1 (longe de qualquer raiz): cresce com a malha ou satura?")
print(f"{'malha':>8} {'dim':>6} {'||V||':>10}")
for m in ['4x12','6x18','8x24']:
    t=time.time(); dim,V,eta = norma_V(m, Q(1))
    print(f"{m:>8} {dim:>6} {float(V):>10.2f}  ({time.time()-t:.0f}s)")
