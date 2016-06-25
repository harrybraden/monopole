import ddphi111
import ddphi112
import ddphi121
import ddphi122
import ddphi211
import ddphi212
import ddphi221
import ddphi222
import ddphi311
import ddphi312
import ddphi321
import ddphi322
from numpy import mat

def ddphi1(zeta, mu, x, k):

    return mat([[ddphi111(zeta, mu, x, k), ddphi112(zeta, mu, x, k)],[ddphi121(zeta, mu, x, k), ddphi122(zeta, mu, x, k)]])


def ddphi2(zeta, mu, x, k):

    return mat([[ddphi211(zeta, mu, x, k), ddphi212(zeta, mu, x, k)],[ddphi221(zeta, mu, x, k), ddphi222(zeta, mu, x, k)]])


def ddphi3(zeta, mu, x, k):

    return mat([[ddphi311(zeta, mu, x, k), ddphi312(zeta, mu, x, k)],[ddphi321(zeta, mu, x, k), ddphi322(zeta, mu, x, k)]])