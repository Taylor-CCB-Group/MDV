#!/bin/bash

# Script to manage test environment for MDV Playwright tests
# This script helps you manage a separate Docker container and database for testing

set -e  # Exit on error

case "$1" in
  start)
    echo "🚀 Starting test containers..."
    docker compose -f docker-test.yml up -d
    echo "⏳ Waiting for services to be ready..."
    sleep 10
    echo "✅ Test environment ready!"
    echo "   App: http://localhost:5056"
    echo "   Dev app still running on: http://localhost:5055"
    ;;
  stop)
    echo "🛑 Stopping test containers..."
    docker compose -f docker-test.yml down
    echo "✅ Test containers stopped"
    ;;
  clean)
    echo "🧹 Cleaning test containers and volumes..."
    docker compose -f docker-test.yml down -v
    echo "✅ Test environment cleaned (all data removed)"
    ;;
  restart)
    echo "🔄 Restarting test environment with fresh database..."
    docker compose -f docker-test.yml down -v
    docker compose -f docker-test.yml up -d
    echo "⏳ Waiting for services to be ready..."
    sleep 10
    echo "✅ Test environment ready with fresh database!"
    echo "   App: http://localhost:5056"
    ;;
  logs)
    echo "📋 Showing test container logs (Ctrl+C to exit)..."
    docker compose -f docker-test.yml logs -f
    ;;
  status)
    echo "📊 Test container status:"
    docker compose -f docker-test.yml ps
    ;;
  *)
    echo "MDV Test Environment Manager"
    echo ""
    echo "Usage: $0 {start|stop|clean|restart|logs|status}"
    echo ""
    echo "Commands:"
    echo "  start    - Start test containers (app on port 5056)"
    echo "  stop     - Stop test containers"
    echo "  clean    - Stop and remove containers and volumes (fresh database)"
    echo "  restart  - Clean restart with fresh database"
    echo "  logs     - Show container logs"
    echo "  status   - Show container status"
    echo ""
    echo "Examples:"
    echo "  $0 start          # Start test environment"
    echo "  $0 restart        # Fresh start with clean database"
    echo "  $0 stop           # Stop when done testing"
    exit 1
    ;;
esac

