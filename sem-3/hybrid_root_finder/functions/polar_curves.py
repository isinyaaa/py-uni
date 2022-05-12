import math
from functions.function_meta import (ContinuousPolarFunction,
                                     ContinuousFunctionMeta)


def butterfly(theta: float) -> float:
    return math.exp(math.sin(theta)) - 2 * math.cos(4 * theta)


def cardioid(theta: float) -> float:
    return 1 - math.sin(theta)


def miscOne(theta: float) -> float:
    return (math.sin(4 * theta)) ** 2 + math.cos(4 * theta)


def miscTwo(theta: float) -> float:
    return 1 + 2 * math.sin(theta / 2)


ONE_X_TWO = ContinuousPolarFunction(
    ContinuousFunctionMeta(
        "Butterfly curve × Cardioid",
        "This is the intersection between the butterfly curve and the cardioid"
    ),
    lambda theta: butterfly(theta) - cardioid(theta),
    2 * math.pi
)


ONE_X_THREE = ContinuousPolarFunction(
    ContinuousFunctionMeta(
        "Butterfly curve × Misc curve #1",
        ("This is the intersection between the butterfly curve and the third"
         "curve of the programming exercise")
    ),
    lambda theta: butterfly(theta) - miscOne(theta),
    2 * math.pi
)


ONE_X_FOUR = ContinuousPolarFunction(
    ContinuousFunctionMeta(
        "Butterfly curve × Misc curve #2",
        ("This is the intersection between the butterfly curve and the fourth"
         "curve of the programming exercise")
    ),
    lambda theta: butterfly(theta) - miscTwo(theta),
    4 * math.pi
)


TWO_X_THREE = ContinuousPolarFunction(
    ContinuousFunctionMeta(
        "Cardioid × Misc curve #1",
        ("This is the intersection between the cardioid curve and the third"
         "curve of the programming exercise")
    ),
    lambda theta: cardioid(theta) - miscOne(theta),
    2 * math.pi
)


TWO_X_FOUR = ContinuousPolarFunction(
    ContinuousFunctionMeta(
        "Cardioid × Misc curve #2",
        ("This is the intersection between the cardioid curve and the fourth"
         "curve of the programming exercise")
    ),
    lambda theta: cardioid(theta) - miscTwo(theta),
    4 * math.pi
)


THREE_X_FOUR = ContinuousPolarFunction(
    ContinuousFunctionMeta(
        "Misc curve #1 × Misc curve #2",
        ("This is the intersection between the third and fourth curves of the"
         "programming exercise")
    ),
    lambda theta: miscOne(theta) - miscTwo(theta),
    4 * math.pi
)
