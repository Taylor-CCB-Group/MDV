test:
	uv run -- pytest mdvtools/ -m "not performance"

test-performance:
	uv run -- pytest mdvtools/tests/test_stress_performance.py -v

type:
	uv run -- pyright

format:
	uv run -- ruff format .

lint:
	uv run -- ruff check . --fix

all: format lint type test