import logging
import rich

def set_logging(level="INFO"):
    logging.basicConfig(level=level, format="%(message)s", handlers=[RichHandler(level=level)])
