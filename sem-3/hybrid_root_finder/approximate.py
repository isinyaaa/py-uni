import math
from typing import Callable

from point import Point


# we define this so that we can use properties as class-owned stuff
class classproperty(property):
    def __get__(self, cls, owner):
        return classmethod(self.fget).__get__(None, owner)()


class FunctionApproximation:
    '''
    This is responsible for all approximation calculations
    '''

    @staticmethod
    def sign(x: float) -> int:
        return (x > 0) - (x < 0)

    _base_epsilon = None

    @classproperty
    def floatEpsilon(cls):
        if cls._base_epsilon is None:
            eps = 1.0
            while eps + 1 > 1:
                eps /= 2
            cls._base_epsilon = eps
        return cls._base_epsilon

    @classmethod
    def getEpsilon(cls, x: float) -> float:
        '''
        Get the least measurable (float-wise) difference around x
        '''

        mag = math.floor(math.log2(math.fabs(x))) + 1 if x != 0 else 0
        return cls.floatEpsilon * (2 ** mag)

    @classmethod
    def isBelowThreshold(cls, x1: float, x2: float) -> bool:
        '''
        Compare current approximations with float epsilon
        '''

        return math.fabs(x2 - x1) <= cls.getEpsilon(math.fabs(x1))

    @staticmethod
    def midpoint(x1: float, x2: float) -> float:
        return (x1 + x2)/2

    @classmethod
    def secantMethod(cls, point_a: Point, point_b: Point, point_c: Point)\
            -> float:

        a = point_a.x
        b, fb = point_b.getTuple()
        c, fc = point_c.getTuple()
        # recalculating this shouldn't be much overhead
        m = cls.midpoint(a, b)

        # now we get p and q (p/q is the ratio we want to analyse)
        p = (b - c) * fb
        p_sign = cls.sign(p)
        q = p_sign * (fc - fb)
        p = p_sign * p

        # if the numerator is less than the denominator
        # we won't ever divide by zero!
        if p <= cls.getEpsilon(q):
            return b + cls.sign(a - b) * cls.getEpsilon(b)
        elif p <= (m - b) * q:
            # secant
            return b + p/q
        else:
            # midpoint
            return m

    @classmethod
    def dekkerMethod(cls, f: Callable[[float], float], errors: list[float],
                     a: float, b: float) -> float:

        if cls.sign(f(a)) == cls.sign(f(b)):
            raise Exception("Expected antipodal points!")

        if (errors is not None and
                len(errors) == 2 and
                all([e != 0 for e in errors])):
            user_defined_errors = True
        else:
            user_defined_errors = False

        # b is the best zero so far
        point_b = Point(b, f(b))
        # (a, f(a)) should be antipodal to (b, f(b)) all the time
        # so that [a, b] contains the wanted root
        point_a = Point(a, f(a))
        # c is the previous value of b
        point_c = Point(a, f(a))

        # times 2 so that we don't incur into the threshold too soon
        previous_m = 2 * cls.midpoint(point_a.x, point_b.x)
        previous_x = 2 * point_b.x

        # another arbitrary constant
        slow_convergence_limit = 1000

        # initialize m here so that we can get the last epsilon after quitting
        m = 0

        while True:
            # if they're the same sign, a should be c
            if cls.sign(point_a.y) == cls.sign(point_b.y):
                point_a = point_c

            # fb should always be the smallest value
            if math.fabs(point_a.y) < math.fabs(point_b.y):
                point_c = point_b
                point_b = point_a
                point_a = point_c

            m = cls.midpoint(point_a.x, point_b.x)

            if not user_defined_errors:
                if cls.isBelowThreshold(m, point_b.x):
                    break
            else:
                TOL = max(errors[0], math.fabs(point_b.x) * errors[1])
                if math.fabs(point_a.x - point_b.x)/2 < TOL:
                    break

            # test if interval isn't changing
            if ((m == previous_m and point_b.x == previous_x) or
                    slow_convergence_limit <= 0):
                raise Exception("This interval won't converge")

            if (cls.isBelowThreshold(m, previous_m) and
                    cls.isBelowThreshold(point_b.x, previous_x)):
                slow_convergence_limit -= 1

            previous_m = m
            previous_x = point_b.x

            s = cls.secantMethod(point_a, point_b, point_c)

            # record previous point
            point_c = point_b
            # then update the best root
            point_b = Point(s, f(s))

        if not user_defined_errors:
            final_error = cls.getEpsilon(math.fabs(m))
        else:
            final_error = max(errors[0], math.fabs(point_b.x) * errors[1])

        return point_b.x, final_error
