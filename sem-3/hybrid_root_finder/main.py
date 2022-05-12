#!/usr/bin/python
from approximate import FunctionApproximation
from functions import FUNCTIONS
from dataclasses import dataclass


@dataclass
class Params:
    '''
    This should hold all params to be used for Dekker method calls
    '''

    function: int
    errors: list[float] = None
    interval: list[float] = None
    constants: list[float] = None

    def __post_init__(self):
        if self.interval is not None and self.constants is not None:
            return

        f = FUNCTIONS[self.function]
        if self.interval is None:
            self.interval = f.interval
        if self.function < 3 and self.constants is None:
            self.constants = f.constants


def get_args(args: list[str]) -> Params:
    '''
    Processes command line arguments
    '''

    if len(args) < 1:
        function_index = 0
    else:
        function_index = int(args[0])

    # we have to parse all other arguments
    args = [
            float(i)
            for i in args[1:]
        ]

    consts = None
    errors = None
    interval = None

    if len(args) >= 2:
        errors = args[:2]

    if len(args) >= 4:
        interval = args[2:4]

    if len(args) > 4:
        consts = args[4:]

    return Params(function_index, errors, interval, consts)


def main(params: Params) -> None:
    function = FUNCTIONS[params.function]

    print(function.meta.name)
    print()
    print(function.meta.context)
    print()
    print('> Calculating...')

    # if the function isn't polar...
    if params.function < 3:
        root, error = FunctionApproximation.dekkerMethod(
            function.function(*params.constants),
            params.errors,
            *params.interval  # notice that we have to unpack parameters
        )

        print('> Done!')
        print()

        print('For the interval', tuple(params.interval))
        print('Dekker method found a zero at x =', root)
        print('with calculated error of', error)
    else:
        # if the function is polar, we must record all of its roots
        roots = []
        # arbitrary interval amount
        steps = 100
        # the size of each interval
        delta = (params.interval[1] - params.interval[0]) / steps
        # and then we generate them all
        intervals = [(delta * i, delta * (i + 1)) for i in range(steps)]
        # then we iterate through each one of them
        for a, b in intervals:
            # and if they either don't converge or have
            # both endpoints with the same sign
            # dekker method will throw us an error
            try:
                new_root, error = FunctionApproximation.dekkerMethod(
                    function.function,
                    params.errors,
                    a, b
                )
                roots.append([new_root, (a, b), error])

            except Exception as e:
                print('Could not find root for', (a, b), 'as', e)

        print('> Done!')
        print()

        for theta, interval, error in roots:
            print('For the interval', interval)
            print('Dekker method found a zero at theta =', theta)
            print('with calculated error of', error)


if __name__ == '__main__':
    import sys

    main(get_args(sys.argv[1:]))
