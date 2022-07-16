from enum import Enum

import numpy as np
import matplotlib.pyplot as plt


# CONSTANTS
AVOGADRO_CONSTANT   = 6.02214086e23
SYSTEM_VOLUME       = 1e-15
MOLECULES_SUBSTRATE = int(5e-7 * AVOGADRO_CONSTANT * SYSTEM_VOLUME)
MOLECULES_ENZYME    = int(2e-7 * AVOGADRO_CONSTANT * SYSTEM_VOLUME)

v = np.array([[-1, 1, 0], [-1, 1, 1], [1, -1, -1], [0, 0, 1]]) # 4x3 stoichiometric matrix

def ssa(tfinal):
    a = np.zeros(3)
    c = np.array([1e6/(AVOGADRO_CONSTANT * SYSTEM_VOLUME), 1e-4, 0.1]) # 1x3

    x = np.zeros((4,1)) # 4x1
    x[0] = MOLECULES_SUBSTRATE
    x[1] = MOLECULES_ENZYME

    t = 0
    while t < tfinal:
        a[0] = c[0] * x[0] * x[1]
        a[1] = c[1] * x[2]
        a[2] = c[2] * x[2]
        asum = np.sum(a)
        rand = np.random.random_sample()
        j = np.amin(np.argwhere(rand < np.cumsum(a / asum))) # lesser element larger than rand
        tau = np.log(1 / rand) / asum
        x = x + v[:, [j]]
        t = t + tau
    return

def cle(tfinal, L):
    tau = tfinal / L

    a = np.zeros(3)
    d = np.zeros(3)
    c = np.array([1e6/(AVOGADRO_CONSTANT * SYSTEM_VOLUME), 1e-4, 0.1]) # 1x3

    y = np.zeros((4,1)) # 4x1
    y[0] = MOLECULES_SUBSTRATE
    y[1] = MOLECULES_ENZYME

    for k in range(L):
        rand = np.random.random_sample()
        a[0] = c[0] * y[0] * y[1]
        a[1] = c[1] * y[2]
        a[2] = c[2] * y[2]
        d[0] = tau * a[0] + np.sqrt(np.absolute(tau * a[0])) * rand
        d[1] = tau * a[1] + np.sqrt(np.absolute(tau * a[1])) * rand
        d[2] = tau * a[2] + np.sqrt(np.absolute(tau * a[2])) * rand
        y = y + d[0]*v[:, [0]] + d[1]*v[:, [1]] + d[2]*v[:, [2]]
    return

class ReactionType(Enum):
    SECOND_ORDER = 'SECOND_ORDER'
    DIMERIZATION = 'DIMERIZATION'
    FIRST_ORDER = 'FIRST_ORDER'


def get_constant_coefficient(reaction_type: str):
    pass


def second_order(X_1, X_2, X_3):
    # substrate (S1)
    # enzyme (S2)
    # complex (S3)
    pass


def dimerization(X_1, X_3):
    # substrate (S1)
    # complex (S3)
    pass


def first_order(X_2, X_3, X_4):
    # enzyme (S2)
    # complex (S3)
    # product (S4)
    pass


def main(reaction_type: ReactionType):
    # get a_k from 1 to M then sum
    #

    ssa(50)
    cle(50, 250)

    if ReactionType.SECOND_ORDER == reaction_type:
        pass
    elif ReactionType.DIMERIZATION == reaction_type:
        pass
    elif ReactionType.FIRST_ORDER == reaction_type:
        pass
    else:
        raise ValueError('Invalid reaction type')
    pass


if __name__ == '__main__':
    from sys import argv
    if len(argv) < 2:
        raise ValueError('No arguments given')

    main(ReactionType(argv[1]))
