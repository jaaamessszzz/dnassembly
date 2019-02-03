#! /usr/bin/env python3
from itertools import tee

# Thank you itertools cookbook!
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def cycle_in_frames(iterable, frame=3):
    """
    Cycle through iterable and yield all possible frames of size frame. Treats iterable as a circular sequence.
    Janky but it works...
    :return:
    """
    assert frame < len(iterable)
    cycle_iterable = iterable * 2

    for cycle_count, derp in enumerate(cycle_iterable):
        if cycle_count < len(iterable):
            yield cycle_iterable[cycle_count: cycle_count + frame]
            cycle_count += 1
