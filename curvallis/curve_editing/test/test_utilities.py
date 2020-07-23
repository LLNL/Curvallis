from collections import namedtuple


def make_args(**kwargs):
    Args = namedtuple('Args', kwargs.keys())
    return Args(*kwargs.values())
