from ddphi111 import ddphi111
from ddphi112 import ddphi112
from ddphi121 import ddphi121
from ddphi122 import ddphi122
from ddphi211 import ddphi211
from ddphi212 import ddphi212
from ddphi221 import ddphi221
from ddphi222 import ddphi222
from ddphi311 import ddphi311
from ddphi312 import ddphi312
from ddphi321 import ddphi321
from ddphi322 import ddphi322
from numpy import mat

def ddphi1(zeta, mu, x, k):

    return mat([[ddphi111(zeta, mu, x, k), ddphi112(zeta, mu, x, k)],[ddphi121(zeta, mu, x, k), ddphi122(zeta, mu, x, k)]])


def ddphi2(zeta, mu, x, k):

    return mat([[ddphi211(zeta, mu, x, k), ddphi212(zeta, mu, x, k)],[ddphi221(zeta, mu, x, k), ddphi222(zeta, mu, x, k)]])


def ddphi3(zeta, mu, x, k):

    return mat([[ddphi311(zeta, mu, x, k), ddphi312(zeta, mu, x, k)],[ddphi321(zeta, mu, x, k), ddphi322(zeta, mu, x, k)]])