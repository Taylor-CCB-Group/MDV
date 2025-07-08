import os

def get_memory_usage():
    """Get current memory usage in GB."""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024 / 1024
    except ImportError:
        return 0.0  # Return 0 if psutil is not available

def log_memory_usage(stage=""):
    """Log current memory usage."""
    memory_gb = get_memory_usage()
    if memory_gb > 0:
        print(f"Memory usage {stage}: {memory_gb:.2f} GB")
    else:
        print(f"Memory monitoring {stage}: psutil not available")
