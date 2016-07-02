from numpy import complex, complex64, mat, dot, trace, pi, sqrt, mat
from mpmath import ellipk, ellipe
from cmath import exp


def dexp(x, k):

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

    DS = [None]*2
    DSP = [None]*2
    DSM = [None]*2
    DSPM = [None]*2
    DRP = [None]*3
    DRM = [None]*3
    DRPBAR = [None]*3
    DRMBAR = [None]*3


    DS[0] = -4 * (K ** 2 - 4 * x[0] ** 2 - 4 * x[1] ** 2) ** (-0.1e1 / 0.2e1) * x[0]

    DS[1] = -4 * (K ** 2 - 4 * x[0] ** 2 - 4 * x[1] ** 2) ** (-0.1e1 / 0.2e1) * x[1]

    DSP[0] = -4 * (x[0] + complex(0, 1) * x[1]) * (4 * x[1] ** 2 + complex(0, -8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.1e1 / 0.2e1)

    DSP[1] = -4 * (complex(0, 1) * x[0] - x[1]) * (4 * x[1] ** 2 + complex(0, -8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.1e1 / 0.2e1)

    DSM[0] = 4 * (complex(0, 1) * x[1] - x[0]) * (4 * x[1] ** 2 + complex(0, 8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.1e1 / 0.2e1)

    DSM[1] = 4 * (x[1] + complex(0, 1) * x[0]) * (4 * x[1] ** 2 + complex(0, 8) * x[0] * x[1] + K ** 2 - 4 * x[0] ** 2) ** (-0.1e1 / 0.2e1)

    DSPM[0] = 4 * (-K ** 2 * k ** 2 * k1 ** 2 + 4 * k1 ** 2 * x[0] ** 2 + 4 * k1 ** 2 * x[1] ** 2 - 4 * x[1] ** 2) ** (-0.1e1 / 0.2e1) * k1 ** 2 * x[0]

    DSPM[1] = 4 * x[1] * (k1 ** 2 - 1) * (-K ** 2 * k ** 2 * k1 ** 2 + 4 * k1 ** 2 * x[0] ** 2 + 4 * k1 ** 2 * x[1] ** 2 - 4 * x[1] ** 2) ** (-0.1e1 / 0.2e1)

    DRP[0] = complex(0, -8) * k1 ** 2 * x[1] ** 2 - 16 * k1 ** 2 * x[0] * x[1] + complex(0, 2) * K ** 2 * k1 ** 2 + complex(0, 16) * x[1] ** 2 + complex(0, -24) * k1 ** 2 * x[0] ** 2 - 16 * x[0] * x[1] + complex(0, 16) * x[2] ** 2

    DRP[1] = -24 * k1 ** 2 * x[1] ** 2 + complex(0, -16) * k1 ** 2 * x[0] * x[1] - 2 * K ** 2 * k1 ** 2 + 24 * x[1] ** 2 - 8 * k1 ** 2 * x[0] ** 2 + complex(0, 32) * x[0] * x[1] + 2 * K ** 2 - 8 * x[0] ** 2 - 16 * x[2] ** 2

    DRP[2] = -32 * x[1] * x[2] + complex(0, 32) * x[0] * x[2]

    DRM[0] = complex(0, 24) * k1 ** 2 * x[1] ** 2 - 48 * k1 ** 2 * x[0] * x[1] + complex(0, 2) * K ** 2 * k1 ** 2 + complex(0, -16) * x[1] ** 2 + complex(0, -24) * k1 ** 2 * x[0] ** 2 + 16 * x[0] * x[1] + complex(0, -16) * x[2] ** 2

    DRM[1] = 24 * k1 ** 2 * x[1] ** 2 + complex(0, 48) * k1 ** 2 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 - 24 * x[1] ** 2 - 24 * k1 ** 2 * x[0] ** 2 + complex(0, -32) * x[0] * x[1] - 2 * K ** 2 + 8 * x[0] ** 2 - 16 * x[2] ** 2

    DRM[2] = -32 * x[1] * x[2] + complex(0, -32) * x[0] * x[2]

    DRPBAR[0] = complex(0, 8) * k1 ** 2 * x[1] ** 2 - 16 * k1 ** 2 * x[0] * x[1] + complex(0, -2) * K ** 2 * k1 ** 2 + complex(0, -16) * x[1] ** 2 + complex(0, 24) * k1 ** 2 * x[0] ** 2 - 16 * x[0] * x[1] + complex(0, -16) * x[2] ** 2

    DRPBAR[1] = -24 * k1 ** 2 * x[1] ** 2 + complex(0, 16) * k1 ** 2 * x[0] * x[1] - 2 * K ** 2 * k1 ** 2 + 24 * x[1] ** 2 - 8 * k1 ** 2 * x[0] ** 2 + complex(0, -32) * x[0] * x[1] + 2 * K ** 2 - 8 * x[0] ** 2 - 16 * x[2] ** 2

    DRPBAR[2] = -32 * x[1] * x[2] + complex(0, -32) * x[0] * x[2]

    DRMBAR[0] = complex(0, -24) * k1 ** 2 * x[1] ** 2 - 48 * k1 ** 2 * x[0] * x[1] + complex(0, -2) * K ** 2 * k1 ** 2 + complex(0, 16) * x[1] ** 2 + complex(0, 24) * k1 ** 2 * x[0] ** 2 + 16 * x[0] * x[1] + complex(0, 16) * x[2] ** 2

    DRMBAR[1] = 24 * k1 ** 2 * x[1] ** 2 + complex(0, -48) * k1 ** 2 * x[0] * x[1] + 2 * K ** 2 * k1 ** 2 - 24 * x[1] ** 2 - 24 * k1 ** 2 * x[0] ** 2 + complex(0, 32) * x[0] * x[1] - 2 * K ** 2 + 8 * x[0] ** 2 - 16 * x[2] ** 2

    DRMBAR[2] = -32 * x[1] * x[2] + complex(0, 32) * x[0] * x[2]

    return DS, DSP, DSM, DSPM, DRP, DRM, DRPBAR, DRMBAR