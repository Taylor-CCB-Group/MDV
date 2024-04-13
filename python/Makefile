test:
	poetry run pytest

type:
	poetry run pyright

format:
	poetry run ruff format .

lint:
	poetry run ruff check . --fix

all: format lint type test