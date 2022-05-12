from typing import NamedTuple


class Point(NamedTuple):
    '''
    Holds coordinates to a point
    '''

    x: float
    y: float

    def getTuple(self):
        return self.x, self.y
