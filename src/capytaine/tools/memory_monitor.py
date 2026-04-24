import logging
import time
from threading import Thread

from capytaine.tools.optional_imports import silently_import_optional_dependency

LOG = logging.getLogger(__name__)

class MemoryMonitor(Thread):
    """Monitor the memory usage in a separate thread.
    from : https://joblib.readthedocs.io/en/stable/auto_examples/parallel_generator.html#sphx-glr-auto-examples-parallel-generator-py
    """

    def __init__(self):
        super().__init__()
        self.stop = False
        self.memory_buffer = [0]
        self.psutil = silently_import_optional_dependency("psutil")

    def get_memory(self):
        "Get memory of a process and its children."
        p = self.psutil.Process()
        memory = p.memory_info().rss
        for c in p.children():
            try:
                memory += c.memory_info().rss
            except self.psutil.NoSuchProcess:
                pass
        return memory

    def run(self):
        if self.psutil is not None:
            memory_start = self.get_memory()
            while not self.stop:
                self.memory_buffer.append(self.get_memory() - memory_start)
                time.sleep(0.2)

    def get_memory_peak(self):
        self.stop = True
        super().join()
        if self.psutil is None:
            return None
        else:
            return round(max(self.memory_buffer) / 1e9, 2)

    def __enter__(self):
        self.start()

    def __exit__(self, exc_type, exc_value, traceback):
        memory_peak = self.get_memory_peak()
        if memory_peak is None:
            LOG.info("Actual peak RAM usage: Not measured since optional dependency `psutil` cannot be found.")
        else:
            LOG.info(f"Actual peak RAM usage: {memory_peak} GB.")
