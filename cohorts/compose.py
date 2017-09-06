

from functools import wraps
import re

class ComposableBool:
    def __init__(self, func):
        self.func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__

    def __and__(self, nextfunc):
        if isinstance(nextfunc, ComposableBool):
            return compose_and(self.func, nextfunc.func)
        else:
            return compose_and(self.func, nextfunc)

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

def compose_and(x, y):
    def comp(filterable_effect, **kwargs):
        return x(filterable_effect, **kwargs) and y(filterable_effect, **kwargs)
    comp.__name__ = re.sub(string="_".join([x.__name__, y.__name__]),
                           pattern="is_", repl="")
    comp.__doc__ = ' & '.join([x.__name__, y.__name__])
    return ComposableBool(comp)

def composable(func):
    return ComposableBool(func)

