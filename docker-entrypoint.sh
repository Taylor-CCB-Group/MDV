#!/bin/bash
set -e

# Create logs dir inside container and fix permissions
mkdir -p /app/logs
chown -R pn:pn /app/logs

# Run the original command
exec "$@"
