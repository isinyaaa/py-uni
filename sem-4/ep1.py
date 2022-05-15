#!/usr/bin/env python3

from random import random


# returns: door index
def get_random_door(chosen_doors: list[int], num_doors: int) -> int:
    if len(chosen_doors) >= num_doors:
        raise Exception("Cannot have more chosen doors than number of doors\
                        available")
    random_num = random() * num_doors // 1
    while random_num in chosen_doors:
        random_num = random() * num_doors // 1

    return random_num


def monty_hall(change_doors: bool, num_doors: int) -> bool:
    prize_door = get_random_door([], num_doors)
    user_chose = get_random_door([], num_doors)
    presenter_chose = get_random_door([prize_door, user_chose], num_doors)
    if change_doors:
        user_chose = get_random_door([presenter_chose, user_chose], num_doors)
    else:
        user_chose = get_random_door([presenter_chose], num_doors)
    return user_chose == prize_door


def main(change_doors: bool, num_doors: int, samples: int) -> None:
    print(sum([
        monty_hall(change_doors, num_doors)
        for i in range(samples)
    ])/samples)


if __name__ == "__main__":
    import sys

    change_doors = True
    num_doors = 4
    num_args = len(sys.argv)
    samples = 100000

    if num_args > 1:
        if sys.argv[1] in "keep k".split():
            change_doors = False

    if num_args > 2:
        try:
            num_doors = int(sys.argv[2])
        except Exception:
            pass

    main(change_doors, num_doors, samples)
