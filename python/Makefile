test:
	poetry run pytest mdvtools/

type:
	poetry run pyright

format:
	poetry run ruff format .

lint:
	poetry run ruff check . --fix

all: format lint type test