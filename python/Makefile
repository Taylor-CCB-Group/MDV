test:
	poetry run pytest mdvtools/tests -m "not performance"

test-performance:
	poetry run pytest mdvtools/tests/test_stress_performance.py -v

test-backend:
	poetry run pytest mdvtools/dbutils/test -v

test-auth:
	poetry run pytest mdvtools/auth/tests -v

type:
	poetry run pyright

format:
	poetry run ruff format .

lint:
	poetry run ruff check . --fix

all: format lint type test
