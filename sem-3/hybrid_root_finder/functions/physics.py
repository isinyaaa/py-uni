import math
from typing import Callable

from functions.function_meta import ContinuousFunction, ContinuousFunctionMeta
from point import Point


def testFunction(mass: float, gravity: float) -> Callable[[float], float]:
    def f(x: float) -> float:
        return gravity*x**2 - mass*math.exp(x)

    return f


TEST_FUNCTION = ContinuousFunction(
    ContinuousFunctionMeta(
        "Test function",
        "Function used for testing under normal conditions"
    ),
    testFunction,
    [1, 10],
    [-1, 0]
)


def _fallingBodySpeed(boundaries: list[Point]):
    def fallingBodySpeed(mass: float, gravity: float):
        def f(x: float) -> float:
            return (mass * gravity - math.exp(-x * boundaries[1].x / mass)
                    * (mass * gravity - boundaries[0].y * x)) / x \
                    - boundaries[1].y

        return f

    return fallingBodySpeed


FALLING_BODY_SPEED_FUNCTION = ContinuousFunction(
    ContinuousFunctionMeta(
        "Speed of falling body",
        "This models the speed of a falling body under air resistance"
    ),
    _fallingBodySpeed,
    [1, 10],
    [-1, 1],
    [
        Point(0, 3),
        Point(2, 20)
    ]
)


def _wireHeight(boundaries: list[Point]):
    def wireHeight(delta: float) -> Callable[[float], float]:
        def f(x: float) -> float:
            return x * math.cosh(boundaries[1].x / x) \
                - x * math.cosh(boundaries[0].x / x) - delta

        return f

    return wireHeight


WIRE_HEIGHT_FUNCTION = ContinuousFunction(
    ContinuousFunctionMeta(
        "Wire height",
        "This function models the height of a wire hanging between two poles"
    ),
    _wireHeight,
    [0.5],
    [99, 101],
    [
        Point(0, 0),
        Point(10, 10)
    ]
)
