from numpy import complex, complex64, mat, dot, trace, pi, sqrt, mat
from mpmath import ellipk, ellipe
from cmath import exp


def dzetas(zeta, x, k):

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

    DZ = [None]*4
    DZ[0] = [None]*3
    DZ[1] = [None]*3
    DZ[2] = [None]*3
    DZ[3] = [None]*3

    DZ[0][0] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[0] ** 2 + 2 * x[2] * zeta[0] - x[1] + complex(0, 1) * x[0]) * (complex(0, 1) * zeta[0] ** 2 + complex(0, 1)) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[0] ** 2 + 2 * x[2] * zeta[0] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[0] + 2 * x[2]) + K ** 2 * (4 * zeta[0] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[0]) / 4)

    DZ[1][0] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[1] ** 2 + 2 * zeta[1] * x[2] - x[1] + complex(0, 1) * x[0]) * (complex(0, 1) * zeta[1] ** 2 + complex(0, 1)) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[1] ** 2 + 2 * zeta[1] * x[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[1] + 2 * x[2]) + K ** 2 * (4 * zeta[1] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[1]) / 4)

    DZ[2][0] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[2] ** 2 + 2 * x[2] * zeta[2] - x[1] + complex(0, 1) * x[0]) * (complex(0, 1) * zeta[2] ** 2 + complex(0, 1)) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[2] ** 2 + 2 * x[2] * zeta[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[2] + 2 * x[2]) + K ** 2 * (4 * zeta[2] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[2]) / 4)

    DZ[3][0] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[3] ** 2 + 2 * zeta[3] * x[2] - x[1] + complex(0, 1) * x[0]) * (complex(0, 1) * zeta[3] ** 2 + complex(0, 1)) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[3] ** 2 + 2 * zeta[3] * x[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[3] + 2 * x[2]) + K ** 2 * (4 * zeta[3] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[3]) / 4)

    DZ[0][1] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[0] ** 2 + 2 * x[2] * zeta[0] - x[1] + complex(0, 1) * x[0]) * (zeta[0] ** 2 - 1) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[0] ** 2 + 2 * x[2] * zeta[0] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[0] + 2 * x[2]) + K ** 2 * (4 * zeta[0] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[0]) / 4)

    DZ[1][1] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[1] ** 2 + 2 * zeta[1] * x[2] - x[1] + complex(0, 1) * x[0]) * (zeta[1] ** 2 - 1) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[1] ** 2 + 2 * zeta[1] * x[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[1] + 2 * x[2]) + K ** 2 * (4 * zeta[1] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[1]) / 4)

    DZ[2][1] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[2] ** 2 + 2 * x[2] * zeta[2] - x[1] + complex(0, 1) * x[0]) * (zeta[2] ** 2 - 1) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[2] ** 2 + 2 * x[2] * zeta[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[2] + 2 * x[2]) + K ** 2 * (4 * zeta[2] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[2]) / 4)

    DZ[3][1] = -2 * ((x[1] + complex(0, 1) * x[0]) * zeta[3] ** 2 + 2 * zeta[3] * x[2] - x[1] + complex(0, 1) * x[0]) * (zeta[3] ** 2 - 1) / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[3] ** 2 + 2 * zeta[3] * x[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[3] + 2 * x[2]) + K ** 2 * (4 * zeta[3] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[3]) / 4)

    DZ[0][2] = -4 * ((x[1] + complex(0, 1) * x[0]) * zeta[0] ** 2 + 2 * x[2] * zeta[0] - x[1] + complex(0, 1) * x[0]) * zeta[0] / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[0] ** 2 + 2 * x[2] * zeta[0] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[0] + 2 * x[2]) + K ** 2 * (4 * zeta[0] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[0]) / 4)

    DZ[1][2] = -4 * ((x[1] + complex(0, 1) * x[0]) * zeta[1] ** 2 + 2 * zeta[1] * x[2] - x[1] + complex(0, 1) * x[0]) * zeta[1] / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[1] ** 2 + 2 * zeta[1] * x[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[1] + 2 * x[2]) + K ** 2 * (4 * zeta[1] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[1]) / 4)

    DZ[2][2] = -4 * ((x[1] + complex(0, 1) * x[0]) * zeta[2] ** 2 + 2 * x[2] * zeta[2] - x[1] + complex(0, 1) * x[0]) * zeta[2] / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[2] ** 2 + 2 * x[2] * zeta[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[2] + 2 * x[2]) + K ** 2 * (4 * zeta[2] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[2]) / 4)

    DZ[3][2] = -4 * ((x[1] + complex(0, 1) * x[0]) * zeta[3] ** 2 + 2 * zeta[3] * x[2] - x[1] + complex(0, 1) * x[0]) * zeta[3] / (2 * ((x[1] + complex(0, 1) * x[0]) * zeta[3] ** 2 + 2 * zeta[3] * x[2] - x[1] + complex(0, 1) * x[0]) * (2 * (x[1] + complex(0, 1) * x[0]) * zeta[3] + 2 * x[2]) + K ** 2 * (4 * zeta[3] ** 3 + 4 * (-2 * k1 ** 2 + 1) * zeta[3]) / 4)

    return DZ