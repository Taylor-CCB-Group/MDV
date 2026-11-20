import logging
import sys


class _Formatter(logging.Formatter):
    """Appends ``trace_id=...`` when telemetry has given the record one.

    ``mdvtools.telemetry`` sets ``trace_id`` on every record while telemetry is
    on, and only lines written inside a traced request get a non-empty value.
    """

    def format(self, record: logging.LogRecord) -> str:
        text = super().format(record)
        trace_id = getattr(record, "trace_id", "")
        if not trace_id:
            return text
        # Docker splits a traceback into one log line per line, so the ID goes
        # on the message line, where Grafana looks for it.
        first, newline, rest = text.partition("\n")
        return f"{first} trace_id={trace_id}{newline}{rest}"


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Set up a logger with consistent formatting across all mdvtools modules.
    
    Args:
        name: The logger name (usually __name__)
        level: The logging level (default: INFO)
    
    Returns:
        A configured logger instance
    """
    logger = logging.getLogger(name)
    
    # Only configure if not already configured
    if not logger.handlers:
        # Create a formatter with the desired format
        formatter = _Formatter(
            fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        
        # Create a console handler and set the formatter
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        console_handler.setLevel(level)
        
        # Add the handler to the logger
        logger.addHandler(console_handler)
        logger.setLevel(level)
    
    return logger


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Get a logger with consistent formatting. This is a convenience function
    that calls setup_logger.
    
    Args:
        name: The logger name (usually __name__)
        level: The logging level (default: INFO)
    
    Returns:
        A configured logger instance
    """
    return setup_logger(name, level) 