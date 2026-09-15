# Subprocess entry: python -m mdvtools.jobs.run_worker "module:fn" <workspace>
import importlib
import sys


def main():
    target, workspace = sys.argv[1], sys.argv[2]
    mod, fn = target.split(":")
    getattr(importlib.import_module(mod), fn)(workspace)


if __name__ == "__main__":
    main()
