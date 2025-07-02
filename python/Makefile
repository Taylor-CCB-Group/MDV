test:
	poetry run pytest mdvtools/ -m "not performance"

test-performance:
	poetry run pytest mdvtools/tests/test_stress_performance.py -v

type:
	poetry run pyright

format:
	poetry run ruff format .

lint:
	poetry run ruff check . --fix

all: format lint type test