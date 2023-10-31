import logging
from rich.logging import RichHandler

def set_logging(level="INFO"):
    logging.basicConfig(level=level, format="%(message)s", handlers=[RichHandler(level=level, log_time_format="[%X]", show_path=False)], force=True)
