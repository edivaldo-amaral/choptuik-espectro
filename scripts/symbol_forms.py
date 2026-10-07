"""Formas pontuais do fundo no L^2 do toro, para K (sistema sharp) e L (selecionado).

Convencao: valores no par (z, -z). A primeira componente e impar
(x1(-z) = -x1(z) = -a); as demais tem valores independentes b_j, b_j'.
Devolve S (forma Hermitiana em u) e M (mapa u -> saidas em z e -z).
Usada em ponto flutuante (prototipo, folga) e em Arb (certificado).
"""


def _consts(like):
    try:
        from flint import acb, arb, fmpq
        if isinstance(like, acb):
            return acb(arb(fmpq(1, 2))), acb(arb(fmpq(1, 4))), acb(2)
    except ImportError:
        pass
    return 0.5, 0.25, 2.0


def forma_L(vz, vm, mu_xi):
    """B_L = mu S Gamma1 + 2 S Xi Gamma2(omega, .), RT (skdfhkfdjhdkhjdgfdjh33hk).
    u = (a, b2, b2', b3, b3', b4, b4').  mu_xi = mu/xi(z)."""
    h, q, dois = _consts(mu_xi)

    def As(v, Pv):
        v1, v2, v3, v4 = v; P1, P2, P3, P4 = Pv
        z0 = 0*v1
        return [(P2+v2-P3-v3, v1-v2+P2, -v1+v2-P2, z0),
                (v1+P3-v3, 2*v2+2*P2, -v2-P2, P4+v4),
                (z0, -v1+P2-v2, v1, P4-v4),
                (-v1, z0, z0, z0),
                (-v2-P2, z0, z0, z0),
                (v4-P4, z0, P4+v4, P2+v2),
                (v4+P4, z0, P4-v4, P2-v2)]
    Az, Am = As(vz, vm), As(vm, vz)

    def saida(A, sgn):
        # coeficientes de y_i no ponto (sgn=+1: z; -1: -z) em u = (a, b2, b2', b3, b3', b4, b4')
        out = []
        for i in range(4):
            c = [dois*q*A[0][i]*sgn]                        # 2*(1/4) A1 x1, x1(-z) = -a
            for j, (e, o) in enumerate(((1, 2), (3, 4), (5, 6))):
                Ae, Ao = A[e][i], A[o][i]                   # 2*(1/2) A E x + 2*(1/2) A O x
                if sgn > 0:   # E=(b+b')/2, O=(b-b')/2 em z
                    c += [dois*h*(Ae+Ao)*h, dois*h*(Ae-Ao)*h]
                else:         # em -z: E=(b+b')/2, O=(b'-b)/2
                    c += [dois*h*(Ae-Ao)*h, dois*h*(Ae+Ao)*h]
            out.append(c)
        return out
    oz, om = saida(Az, +1), saida(Am, -1)
    # mu S Gamma1: y2 += mu (x1 + 2 O x2)/xi ; y3 += -mu x1/xi ; y4 += mu 2 O x4 / xi
    # em z: x1=a, O x2 = (b2-b2')/2 ; em -z: x1/xi(-z) = a/xi(z), O x2(-z)/xi(-z) = O x2(z)/xi(z)
    for o in (oz, om):
        o[1][0] = o[1][0] + mu_xi
        o[1][1] = o[1][1] + mu_xi
        o[1][2] = o[1][2] - mu_xi
        o[2][0] = o[2][0] - mu_xi
        o[3][5] = o[3][5] + mu_xi
        o[3][6] = o[3][6] - mu_xi
    S = [[oz[0][j] - om[0][j] for j in range(7)],
         oz[1][:], om[1][:], oz[2][:], om[2][:], oz[3][:], om[3][:]]
    M = [oz[0][:], om[0][:], oz[1][:], om[1][:], oz[2][:], om[2][:], oz[3][:], om[3][:]]
    return S, M


def forma_K(vz, vm, mu_xi):
    """vz, vm: valores (v1..v4) em z e em -z. mu_xi = mu/xi(z).
    Devolve S (3x3, forma em (a,b,b')) e M (4x3, saidas out1(z),out1(-z),out2(z),out2(-z))."""
    def Fs(v, Pv):
        v1, v2, v3, v4 = v; P1, P2, P3, P4 = Pv
        return ((v2+P2-P3-v3, v1-v2+P2), (-v1+P3-v3, -v1+v2+3*P2), (v2+P2+P3+v3, -v1+v2+3*P2))
    F1z, F2z, F3z = Fs(vz, vm)
    F1m, F2m, F3m = Fs(vm, vz)
    h, q, _ = _consts(mu_xi)
    # out_i(z)  = h F1_i(z) a + q (F2_i+F3_i)(z) b + q (F2_i-F3_i)(z) b'
    # out_i(-z) = -h F1_i(-z) a + q (F2_i-F3_i)(-z) b + q (F2_i+F3_i)(-z) b'
    oz = [[h*F1z[i], q*(F2z[i]+F3z[i]), q*(F2z[i]-F3z[i])] for i in range(2)]
    om = [[-h*F1m[i], q*(F2m[i]-F3m[i]), q*(F2m[i]+F3m[i])] for i in range(2)]
    oz[1][0] = oz[1][0] + mu_xi
    om[1][0] = om[1][0] + mu_xi                     # x1(-z)/xi(-z) = x1(z)/xi(z)
    S = [[oz[0][j] - om[0][j] for j in range(3)], oz[1][:], om[1][:]]
    M = [oz[0][:], om[0][:], oz[1][:], om[1][:]]
    return S, M
