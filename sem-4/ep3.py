#!/usr/bin/env python3

def main(p: float, values: [int], iterations):
    from random import random

    probability = 0.0
    pc = 1 - p

    for i in values:
        wins = 0
        for _ in range(iterations):
            trial = True
            for j in range(i - 1):
                if random() > pc:
                    trial = False
                    break
            if random() < p and trial:
                wins += 1

        probability += wins/iterations

    print("%.5f" %probability)


# This tries to find how many iterations we need until we achieve a difference
# of machine-epsilon between subsequent values of p(X)
def range_for_min_error(p: float, range_start: int) -> int:
    # complement of p
    pc = 1 - p

    epsilon = 0.0
    # here is the only part where we use the knowledge of the formula to make
    # calculations more effective (otherwise we could use just p ** range_start
    # to make this calculation, but that wouldn't be as optimized)
    init_value = pc ** (range_start - 1) * p
    tmp = init_value
    i = 0

    while init_value + tmp > init_value:
        epsilon = tmp
        tmp *= pc
        i += 1

    return i + range_start


def parse_input(argv: [str]) -> (float, [int]):
    if len(argv) == 1:
        raise Exception(("Usage: ./ep3.py (probability) (ranges/values)\n"
                         "where the range is given by\n"
                         "num1 - num2\n"
                         "you can also do\n"
                         "num1 - num2 num3 num4 -\n"
                         "where the last range goes to infinity"));

    try:
        p = float(argv[1])
    except Exception:
        raise Exception("Could not parse", argv[1])

    if p > 1.0:
        raise Exception(p, "will not converge!")

    argv = argv[2:]
    values = []
    range_start = -1
    skip = False
    for i, val in enumerate(argv):
        if skip:
            skip = False
            continue

        if val != '-':
            if range_start != -1:
                values.extend([range_start])
            try:
                range_start = int(val)
            except Exception:
                raise Exception("Could not parse", i)

        else:
            if range_start == -1:
                raise Exception("You must specify where the range begins");

            if len(argv) == i + 1:
                values.extend(list(range(range_start, range_for_min_error(p, range_start))))
                break

            try:
                values.extend(list(range(range_start, int(argv[i + 1]) + 1)))
            except Exception:
                raise Exception("Could not parse", i)

            skip = True
            range_start = -1

    return p, values


if __name__ == "__main__":
    from sys import argv

    p, values = parse_input(argv)

    iterations = 10**6
    main(p, values, iterations)
