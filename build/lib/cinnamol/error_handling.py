from dataclasses import dataclass 
from functools import partial 
from typing import Any, Callable, Union


@dataclass 
class Fail:
    """
    Class for representing a failed execution of a function.

    Attributes
    ---------
    exception : Exception
        Exception raised during execution of function.
    name : str
        Name of function that failed.
    args : tuple
        Arguments passed to function.
    kwargs : dict
        Keyword arguments passed to function.
    """
    exception: Exception | None = None 
    name: str | None = None 
    args: Any = None 
    kwargs: Any = None 


@dataclass
class Success:
    """
    Class for representing a successful execution of a function.

    Attributes
    ---------
    value : Any
        Value returned by function.
    """
    value: Any = None


def failed(result) -> bool:
    """
    Check if result is a failed execution.
    
    Arguments
    ---------
    result : Any
        Result of function execution.

    Returns
    -------
    bool
        True if result is a failed execution.
    """
    return isinstance(result, Fail)


def success(result) -> bool:
    """
    Check if result is a successful execution.
    
    Arguments
    ---------
    result : Any
        Result of function execution.
    
    Returns
    -------
    bool
        True if result is a successful execution.
    """
    return isinstance(result, Success)


def try_catch(
    func: Callable, 
    *args, **kwargs
) -> Union[Success, Fail]:
    """
    Try to execute a function and catch any exceptions.

    Arguments
    ---------
    func : Callable
        Function to execute.
    args : tuple
        Arguments passed to function.
    kwargs : dict
        Keyword arguments passed to function.

    Returns
    -------
    Union[Success, Fail]
        Result of function execution.
    """
    try:
        return Success(func(*args, **kwargs))
    except Exception as exc:
        return Fail(exc, func.__name__, args, kwargs)


def wrapper(
    error_handler: Callable, 
    func: Callable, 
    *args, **kwargs
) -> Union[Any, Fail]:
    """
    Wrapper for executing a function and handling any exceptions.
    
    Arguments
    ---------
    error_handler : Callable
        Function to handle any exceptions.
    func : Callable
        Function to execute.
    args : tuple
        Arguments passed to function.
    kwargs : dict
        Keyword arguments passed to function.
    
    Returns
    -------
    Union[Any, Fail]
        Result of function execution.
    """
    match args:
        case [result, *_] if len(args) and success(result):
            return result.value 
        case _:
            return error_handler(func, *args, **kwargs)
            

def error_handler(func: Callable) -> Callable:
    """
    Decorator for handling any exceptions raised by a function.
    
    Arguments
    ---------
    func : Callable
        Function to execute.
    
    Returns
    -------
    Callable    
        Wrapped function.
    """
    return partial(wrapper, try_catch, func)
