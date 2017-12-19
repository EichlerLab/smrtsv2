"""
Utility functions for use by any part the project.
"""


def as_bool(val):
    """
    Convert a value to a boolean. `val` is `True` if it is a boolean that contains True, is a string
    of 'true' or 't' (case-insensitive), or is a non-zero numeric (or string of non-zero numeric). False
    in all other cases including if `val` is `None`.

    :param val: Value to test.

    :return: `True` if `val` represents a True value, or `False` otherwise.
    """

    # None is False
    if val is None:
        return False

    # val if val is boolean
    if type(val) == bool:
        return val

    # True strings
    if val.lower in ('true', 't'):
        return True

    # Non-zero int or float is true
    try:
        return float(val) != 0.0

    except:
        pass

    # Everything else is False
    return False
