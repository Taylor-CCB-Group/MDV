check-docker:
	@echo "Checking Docker daemon..."
	@docker info >/dev/null 2>&1 && echo "Docker daemon is running" || (echo "Docker daemon is not running"; exit 1)

docker-frontend-only: check-docker
	docker-compose -f docker-compose.yml up -d --no-deps frontend

docker-backend: check-docker
	docker-compose -f docker-compose.yml up -d

docker-down: check-docker
	docker-compose -f docker-compose.yml down

docker-build: check-docker
	docker-compose -f docker-compose.yml build

docker-restart-frontend: docker-down docker-frontend-only

docker-restart-backend: docker-down docker-backend