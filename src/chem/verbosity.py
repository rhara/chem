import functools
import inspect
import os
import sys

_QUIET_FALSE_VALUES = {"0", "n", "false"}


def is_quiet():
    """True if CHEM_QUIETNESS is set to anything other than '0', 'N', or 'FALSE' (case-insensitive)."""
    value = os.environ.get("CHEM_QUIETNESS")
    if value is None:
        return False
    return value.strip().lower() not in _QUIET_FALSE_VALUES


def logged(func):
    """Print the called function name and arguments to stderr, unless is_quiet()."""

    sig = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if not is_quiet():
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            params = ", ".join(f"{k}={v!r}" for k, v in bound.arguments.items())
            print(f"{func.__name__}({params})", file=sys.stderr)
        return func(*args, **kwargs)

    return wrapper
