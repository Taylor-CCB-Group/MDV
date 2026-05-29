test:
	uv run -- pytest mdvtools/tests -m "not performance"

test-performance:
	uv run -- pytest mdvtools/tests/test_stress_performance.py -v

test-backend:
	uv run -- pytest mdvtools/dbutils/test -v

test-auth:
	uv run -- pytest mdvtools/auth/tests -v

# Cursor/VS Code integrated terminals may set NODE_OPTIONS for js-debug; pyright invokes Node.
type:
	env -u NODE_OPTIONS uv run -- pyright

format:
	uv run -- ruff format .

lint:
	uv run -- ruff check . --fix

all: format lint type test
