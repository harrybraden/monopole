from numpy import complex, complex64, mat, dot, trace, pi, sqrt, mat
from mpmath import ellipk, ellipe
from cmath import exp


def dmus(zeta, x, k):

    K = complex64(ellipk(k**2))

    E = complex64(ellipe(k**2))

    cm= (2*E-K)/K

    k1 = sqrt(1-k**2)

    xp = x[0]+complex(0,1)*x[1]
    xm = x[0]-complex(0,1)*x[1]
    S =  sqrt(K**2-4*xp*xm)
    SP = sqrt(K**2-4*xp**2)
    SM = sqrt(K**2-4*xm**2)
    SPM = sqrt(-k1**2*(K**2*k**2-4*xm*xp)+(xm-xp)**2)
    R = 2*K**2*k1**2-S**2-8*x[2]**2
    RM = complex(0,1)*SM**2*(xm*(2*k1**2-1)+xp)-(16*complex(0,1))*xm*x[2]**2
    RP = complex(0,1)*SM**2*(xp*(2*k1**2-1)+xm)+(16*complex(0,1))*xp*x[2]**2
    RMBAR=-complex(0,1)*SP**2*( xp*(2*k1**2-1)+xm ) +16*complex(0,1)*xp*x[2]**2
    RPBAR=-complex(0,1)*SP**2*( xm*(2*k1**2-1)+xp ) -16*complex(0,1)*xm*x[2]**2
    r=sqrt(x[0]**2+x[1]**2+x[2]**2)

    DM = [None,None,None,None]

    DM[0] = [None]*3
    DM[1] = [None]*3
    DM[2] = [None]*3
    DM[3] = [None]*3

    DM[0][0] = complex(0, 1) * (E * K * zeta[0] ** 2 - K ** 2 * zeta[0] ** 2 + 4 * zeta[0] ** 2 * x[0] ** 2 + complex(0, -4) * zeta[0] ** 2 * x[1] * x[0] + E * K - K ** 2 + 4 * x[0] ** 2 + complex(0, -8) * zeta[0] * x[2] * x[0] + complex(0, 4) * x[0] * x[1]) / (4 * zeta[0] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[0] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[0] - K ** 2 * zeta[0] ** 3 + 4 * x[0] ** 2 * zeta[0] + complex(0, -12) * zeta[0] ** 2 * x[0] * x[2] - 4 * zeta[0] ** 3 * x[1] ** 2 - 12 * zeta[0] ** 2 * x[1] * x[2] - K ** 2 * zeta[0] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[0] - 8 * zeta[0] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[1][0] = complex(0, 1) * (E * K * zeta[1] ** 2 - K ** 2 * zeta[1] ** 2 + 4 * zeta[1] ** 2 * x[0] ** 2 + complex(0, -4) * zeta[1] ** 2 * x[1] * x[0] + E * K - K ** 2 + 4 * x[0] ** 2 + complex(0, -8) * zeta[1] * x[2] * x[0] + complex(0, 4) * x[0] * x[1]) / (4 * zeta[1] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[1] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[1] - K ** 2 * zeta[1] ** 3 + 4 * x[0] ** 2 * zeta[1] + complex(0, -12) * zeta[1] ** 2 * x[0] * x[2] - 4 * zeta[1] ** 3 * x[1] ** 2 - 12 * zeta[1] ** 2 * x[1] * x[2] - K ** 2 * zeta[1] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[1] - 8 * zeta[1] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[2][0] = complex(0, 1) * (E * K * zeta[2] ** 2 - K ** 2 * zeta[2] ** 2 + 4 * zeta[2] ** 2 * x[0] ** 2 + complex(0, -4) * zeta[2] ** 2 * x[1] * x[0] + E * K - K ** 2 + 4 * x[0] ** 2 + complex(0, -8) * zeta[2] * x[2] * x[0] + complex(0, 4) * x[0] * x[1]) / (4 * zeta[2] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[2] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[2] - K ** 2 * zeta[2] ** 3 + 4 * x[0] ** 2 * zeta[2] + complex(0, -12) * zeta[2] ** 2 * x[0] * x[2] - 4 * zeta[2] ** 3 * x[1] ** 2 - 12 * zeta[2] ** 2 * x[1] * x[2] - K ** 2 * zeta[2] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[2] - 8 * zeta[2] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[3][0] = complex(0, 1) * (E * K * zeta[3] ** 2 - K ** 2 * zeta[3] ** 2 + 4 * zeta[3] ** 2 * x[0] ** 2 + complex(0, -4) * zeta[3] ** 2 * x[1] * x[0] + E * K - K ** 2 + 4 * x[0] ** 2 + complex(0, -8) * zeta[3] * x[2] * x[0] + complex(0, 4) * x[0] * x[1]) / (4 * zeta[3] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[3] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[3] - K ** 2 * zeta[3] ** 3 + 4 * x[0] ** 2 * zeta[3] + complex(0, -12) * zeta[3] ** 2 * x[0] * x[2] - 4 * zeta[3] ** 3 * x[1] ** 2 - 12 * zeta[3] ** 2 * x[1] * x[2] - K ** 2 * zeta[3] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[3] - 8 * zeta[3] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[0][1] = (complex(0, 4) * zeta[0] ** 2 * x[1] * x[0] + E * K * zeta[0] ** 2 + 4 * zeta[0] ** 2 * x[1] ** 2 + complex(0, 4) * x[0] * x[1] + 8 * zeta[0] * x[1] * x[2] - E * K - 4 * x[1] ** 2) / (4 * zeta[0] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[0] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[0] - K ** 2 * zeta[0] ** 3 + 4 * x[0] ** 2 * zeta[0] + complex(0, -12) * zeta[0] ** 2 * x[0] * x[2] - 4 * zeta[0] ** 3 * x[1] ** 2 - 12 * zeta[0] ** 2 * x[1] * x[2] - K ** 2 * zeta[0] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[0] - 8 * zeta[0] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[1][1] = (complex(0, 4) * zeta[1] ** 2 * x[1] * x[0] + E * K * zeta[1] ** 2 + 4 * zeta[1] ** 2 * x[1] ** 2 + complex(0, 4) * x[0] * x[1] + 8 * zeta[1] * x[1] * x[2] - E * K - 4 * x[1] ** 2) / (4 * zeta[1] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[1] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[1] - K ** 2 * zeta[1] ** 3 + 4 * x[0] ** 2 * zeta[1] + complex(0, -12) * zeta[1] ** 2 * x[0] * x[2] - 4 * zeta[1] ** 3 * x[1] ** 2 - 12 * zeta[1] ** 2 * x[1] * x[2] - K ** 2 * zeta[1] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[1] - 8 * zeta[1] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[2][1] = (complex(0, 4) * zeta[2] ** 2 * x[1] * x[0] + E * K * zeta[2] ** 2 + 4 * zeta[2] ** 2 * x[1] ** 2 + complex(0, 4) * x[0] * x[1] + 8 * zeta[2] * x[1] * x[2] - E * K - 4 * x[1] ** 2) / (4 * zeta[2] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[2] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[2] - K ** 2 * zeta[2] ** 3 + 4 * x[0] ** 2 * zeta[2] + complex(0, -12) * zeta[2] ** 2 * x[0] * x[2] - 4 * zeta[2] ** 3 * x[1] ** 2 - 12 * zeta[2] ** 2 * x[1] * x[2] - K ** 2 * zeta[2] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[2] - 8 * zeta[2] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[3][1] = (complex(0, 4) * zeta[3] ** 2 * x[1] * x[0] + E * K * zeta[3] ** 2 + 4 * zeta[3] ** 2 * x[1] ** 2 + complex(0, 4) * x[0] * x[1] + 8 * zeta[3] * x[1] * x[2] - E * K - 4 * x[1] ** 2) / (4 * zeta[3] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[3] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[3] - K ** 2 * zeta[3] ** 3 + 4 * x[0] ** 2 * zeta[3] + complex(0, -12) * zeta[3] ** 2 * x[0] * x[2] - 4 * zeta[3] ** 3 * x[1] ** 2 - 12 * zeta[3] ** 2 * x[1] * x[2] - K ** 2 * zeta[3] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[3] - 8 * zeta[3] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[0][2] = 2 * (-K ** 2 * k1 ** 2 * zeta[0] + complex(0, 2) * zeta[0] ** 2 * x[0] * x[2] + 2 * zeta[0] ** 2 * x[1] * x[2] + E * K * zeta[0] + complex(0, 2) * x[0] * x[2] + 4 * zeta[0] * x[2] ** 2 - 2 * x[1] * x[2]) / (4 * zeta[0] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[0] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[0] - K ** 2 * zeta[0] ** 3 + 4 * x[0] ** 2 * zeta[0] + complex(0, -12) * zeta[0] ** 2 * x[0] * x[2] - 4 * zeta[0] ** 3 * x[1] ** 2 - 12 * zeta[0] ** 2 * x[1] * x[2] - K ** 2 * zeta[0] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[0] - 8 * zeta[0] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[1][2] = 2 * (-K ** 2 * k1 ** 2 * zeta[1] + complex(0, 2) * zeta[1] ** 2 * x[0] * x[2] + 2 * zeta[1] ** 2 * x[1] * x[2] + E * K * zeta[1] + complex(0, 2) * x[0] * x[2] + 4 * zeta[1] * x[2] ** 2 - 2 * x[1] * x[2]) / (4 * zeta[1] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[1] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[1] - K ** 2 * zeta[1] ** 3 + 4 * x[0] ** 2 * zeta[1] + complex(0, -12) * zeta[1] ** 2 * x[0] * x[2] - 4 * zeta[1] ** 3 * x[1] ** 2 - 12 * zeta[1] ** 2 * x[1] * x[2] - K ** 2 * zeta[1] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[1] - 8 * zeta[1] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[2][2] = 2 * (-K ** 2 * k1 ** 2 * zeta[2] + complex(0, 2) * zeta[2] ** 2 * x[0] * x[2] + 2 * zeta[2] ** 2 * x[1] * x[2] + E * K * zeta[2] + complex(0, 2) * x[0] * x[2] + 4 * zeta[2] * x[2] ** 2 - 2 * x[1] * x[2]) / (4 * zeta[2] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[2] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[2] - K ** 2 * zeta[2] ** 3 + 4 * x[0] ** 2 * zeta[2] + complex(0, -12) * zeta[2] ** 2 * x[0] * x[2] - 4 * zeta[2] ** 3 * x[1] ** 2 - 12 * zeta[2] ** 2 * x[1] * x[2] - K ** 2 * zeta[2] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[2] - 8 * zeta[2] * x[2] ** 2 + 4 * x[1] * x[2])

    DM[3][2] = 2 * (-K ** 2 * k1 ** 2 * zeta[3] + complex(0, 2) * zeta[3] ** 2 * x[0] * x[2] + 2 * zeta[3] ** 2 * x[1] * x[2] + E * K * zeta[3] + complex(0, 2) * x[0] * x[2] + 4 * zeta[3] * x[2] ** 2 - 2 * x[1] * x[2]) / (4 * zeta[3] ** 3 * x[0] ** 2 + complex(0, -8) * zeta[3] ** 3 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 * zeta[3] - K ** 2 * zeta[3] ** 3 + 4 * x[0] ** 2 * zeta[3] + complex(0, -12) * zeta[3] ** 2 * x[0] * x[2] - 4 * zeta[3] ** 3 * x[1] ** 2 - 12 * zeta[3] ** 2 * x[1] * x[2] - K ** 2 * zeta[3] + complex(0, -4) * x[0] * x[2] + 4 * x[1] ** 2 * zeta[3] - 8 * zeta[3] * x[2] ** 2 + 4 * x[1] * x[2])

    return DM