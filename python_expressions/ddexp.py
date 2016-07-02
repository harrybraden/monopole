from numpy import complex, complex64, mat, dot, trace, pi, sqrt, mat
from mpmath import ellipk, ellipe
from cmath import exp


def ddexp(x, k):

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

    DDS = [None]*2
    DDSP = [None]*2
    DDSM = [None]*2
    DDSPM = [None]*2
    DDRP = [None]*3
    DDRM = [None]*3
    DDRPBAR = [None]*3
    DDRMBAR = [None]*3

    DDS[0] = -4 * (K ** 2 - 4 * x[1] ** 2) * (K ** 2 - 4 * x[0] ** 2 - 4 * x[1] ** 2) ** (-0.3e1 / 0.2e1)

    DDS[1] = -4 * (K ** 2 - 4 * x[0] ** 2) * (K ** 2 - 4 * x[0] ** 2 - 4 * x[1] ** 2) ** (-0.3e1 / 0.2e1)

    DDSP[0] = -4 * K ** 2 * (4 * x[1] ** 2 + complex(0, -8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.3e1 / 0.2e1)

    DDSP[1] = 4 * K ** 2 * (4 * x[1] ** 2 + complex(0, -8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.3e1 / 0.2e1)

    DDSM[0] = -4 * K ** 2 * (4 * x[1] ** 2 + complex(0, 8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.3e1 / 0.2e1)

    DDSM[1] = 4 * K ** 2 * (4 * x[1] ** 2 + complex(0, 8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.3e1 / 0.2e1)

    DDSPM[0] = -4 * k1 ** 2 * (K ** 2 * k ** 2 * k1 ** 2 - 4 * k1 ** 2 * x[1] ** 2 + 4 * x[1] ** 2) * (-K ** 2 * k ** 2 * k1 ** 2 + 4 * k1 ** 2 * x[0] ** 2 + 4 * k1 ** 2 * x[1] ** 2 - 4 * x[1] ** 2) ** (-0.3e1 / 0.2e1)

    DDSPM[1] = -4 * (k1 ** 2 - 1) * k1 ** 2 * (K ** 2 * k ** 2 - 4 * x[0] ** 2) * (-K ** 2 * k ** 2 * k1 ** 2 + 4 * k1 ** 2 * x[0] ** 2 + 4 * k1 ** 2 * x[1] ** 2 - 4 * x[1] ** 2) ** (-0.3e1 / 0.2e1)

    DDRP[0] = -16 * k1 ** 2 * x[1] + complex(0, -48) * k1 ** 2 * x[0] - 16 * x[1]

    DDRP[1] = -48 * k1 ** 2 * x[1] + complex(0, -16) * k1 ** 2 * x[0] + 48 * x[1] + complex(0, 32) * x[0]

    DDRP[2] = -32 * x[1] + complex(0, 32) * x[0]

    DDRM[0] = -48 * k1 ** 2 * x[1] + complex(0, -48) * k1 ** 2 * x[0] + 16 * x[1]

    DDRM[1] = 48 * k1 ** 2 * x[1] + complex(0, 48) * k1 ** 2 * x[0] - 48 * x[1] + complex(0, -32) * x[0]

    DDRM[2] = -32 * x[1] + complex(0, -32) * x[0]

    DDRPBAR[0] = -16 * k1 ** 2 * x[1] + complex(0, 48) * k1 ** 2 * x[0] - 16 * x[1]

    DDRPBAR[1] = -48 * k1 ** 2 * x[1] + complex(0, 16) * k1 ** 2 * x[0] + 48 * x[1] + complex(0, -32) * x[0]

    DDRPBAR[2] = -32 * x[1] + complex(0, -32) * x[0]

    DDRMBAR[0] = -48 * k1 ** 2 * x[1] + complex(0, 48) * k1 ** 2 * x[0] + 16 * x[1]

    DDRMBAR[1] = 48 * k1 ** 2 * x[1] + complex(0, -48) * k1 ** 2 * x[0] - 48 * x[1] + complex(0, 32) * x[0]

    DDRMBAR[2] = -32 * x[1] + complex(0, 32) * x[0]

