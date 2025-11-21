import time
from threading import Thread

from psutil import Process

class MemoryMonitor(Thread):
    """Monitor the memory usage in a separate thread.
    from : https://joblib.readthedocs.io/en/stable/auto_examples/parallel_generator.html#sphx-glr-auto-examples-parallel-generator-py
    """

    def __init__(self, psutil):
        super().__init__()
        self.stop = False
        self.memory_buffer = [0]
        self.psutil = psutil 
        self.start()

    def get_memory(self):
        "Get memory of a process and its children."
        p = Process()
        memory = p.memory_info().rss
        for c in p.children():
            memory += c.memory_info().rss
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