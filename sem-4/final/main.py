from enum import Enum

import numpy as np
import matplotlib.pyplot as plt


# CONSTANTS
AVOGADRO_CONSTANT   = 6.02214086e23
SYSTEM_VOLUME       = 1e-15
MOLECULES_SUBSTRATE = int(5e-7 * AVOGADRO_CONSTANT * SYSTEM_VOLUME)
MOLECULES_ENZYME    = int(2e-7 * AVOGADRO_CONSTANT * SYSTEM_VOLUME)

STROICHIMETRIC_MATRIX = np.array([[-1, 1, 0], [-1, 1, 1], [1, -1, -1], [0, 0, 1]]) # 4x3

x = np.zeros((4,1)) # 4x1
x[0] = MOLECULES_SUBSTRATE
x[1] = MOLECULES_SUBSTRATE
c = np.array([1e6/(AVOGADRO_CONSTANT * SYSTEM_VOLUME), 1e-4, 0.1]) # 1x3

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
