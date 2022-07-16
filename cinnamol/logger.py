from logging import getLogger, Logger, StreamHandler, FileHandler 
from typing import Optional, Union


def set_logger(
    logger_name: str,
    logger_level: Union[str, int],
    logger_output: Optional[str] = None,
) -> Logger:
    """
    Set up a logger.

    Arguments
    ---------
    logger_name: str
        Name of the logger.
    logger_level: str | int
        Logging level.
    logger_output: str | None
        Path to logger output file.
    """
    logger = getLogger(logger_name)
    logger.setLevel(logger_level)

    # Add stream handler to logger.
    stream_handler = StreamHandler()
    stream_handler.setLevel(logger_level)
    logger.addHandler(stream_handler)

    # Add file handler to logger if output path is provided.
    if logger_output is not None:
        file_handler = FileHandler(logger_output)
        file_handler.setLevel(logger_level)
        logger.addHandler(file_handler)

    return logger
    