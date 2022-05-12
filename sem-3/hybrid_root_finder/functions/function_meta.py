from dataclasses import dataclass, field
from typing import Callable

from point import Point


@dataclass
class ContinuousFunctionMeta:
    '''
    This should store all metadata about the functions
    '''

    name: str = ""
    context: str = ""


@dataclass
class ContinuousFunctionBase:
    '''
    This is the most general continuous function definition
    '''

    meta: ContinuousFunctionMeta
    function: Callable


@dataclass
class ContinuousFunction(ContinuousFunctionBase):
    '''
    This is useful in the context of physics-based functions
    '''

    constants: list[float]
    interval: list[float]
    boundaries: list[Point] = None

    def __post_init__(self):
        if self.boundaries is not None:
            self.function = self.function(self.boundaries)


@dataclass
class ContinuousPolarFunction(ContinuousFunctionBase):
    '''
    This is focused on polar functions
    '''

    period: float
    interval: list[float] = field(init=False)

    def __post_init__(self):
        self.interval = [0, self.period]
