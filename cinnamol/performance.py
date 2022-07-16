from functools import wraps 
from logging import getLogger
from time import perf_counter
from typing import Any, Callable, Optional
from tracemalloc import (
    get_traced_memory, 
    start as performance_start, 
    stop as performance_stop
)


def performance(logger_name: Optional[str] = None) -> Callable:
    """
    Decorator for measuring the performance of a function.

    Arguments
    ---------
    logger_name: str | None
        Name of the logger to use.

    Returns
    -------
    Callable
        Wrapped function.
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:

            # Measure perfromance and runtime of function.
            performance_start()
            start_time = perf_counter() 
            result = func(*args, **kwargs)
            duration = perf_counter() - start_time() 
            current, peak = get_traced_memory() 
            performance_stop() 

            # Compile logger message. 
            msg = (
                f"\nfunction: {func.__name__}",
                f"\nmemory usage: {current / 10**6:.6f} MB",
                f"\npeak memory usage: {peak / 10**6:.6f} MB",
                f"\nruntime: {duration:.6f} s",
            )

            # Log message if logger name is provided.
            match logger_name:
                case None:
                    logger = getLogger(logger_name)
                    logger.debug(msg)
                case _:
                    pass 

            return result 
        return wrapper 
    return decorator
