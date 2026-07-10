from __future__ import annotations

import logging
from typing import Any
from urllib.parse import urlencode

import requests
from flask import Flask, make_response, request

from mdvtools.server_extension import MDVProjectServerExtension

logger = logging.getLogger(__name__)

ALLOWED_UCSC_HOSTS = {
    "genome.ucsc.edu",
    "genome-euro.ucsc.edu",
    "genome-asia.ucsc.edu",
}

DEFAULT_MAX_UCSC_BYTES = 10 * 1024 * 1024  # 10MB max for UCSC image responses
UCSC_CHUNK_SIZE = 8192


def register_ucsc_proxy_routes(app: Flask) -> None:
    """
    Register an app-wide `/ucsc_proxy` route.

    This endpoint is not project-specific; it proxies to UCSC's `hgRenderTracks`.
    """
    # Avoid endpoint/rule duplication if called from multiple entry points.
    if any(rule.rule == "/ucsc_proxy" for rule in app.url_map.iter_rules()):
        return

    @app.route("/ucsc_proxy")
    def ucsc_image():
        try:
            params = request.args.to_dict()
            ucsc_host = params.pop("ucscHost", "genome.ucsc.edu")
            if ucsc_host not in ALLOWED_UCSC_HOSTS:
                logger.warning("Rejected invalid UCSC host: %s", ucsc_host)
                return "Invalid host", 400

            base_url = f"https://{ucsc_host}/cgi-bin/hgRenderTracks"
            ucsc_url = f"{base_url}?{urlencode(params)}"

            headers: dict[str, str] = {
                "User-Agent": "Mozilla/5.0 (compatible; MDV-Proxy/1.0)",
            }

            response = requests.get(ucsc_url, headers=headers, timeout=30, stream=True)
            response.raise_for_status()

            content_type = response.headers.get("Content-Type", "")
            if not content_type.startswith("image/"):
                logger.warning("Non-image content type from UCSC: %s", content_type)
                response.close()
                return "Unsupported media type", 415

            max_bytes = int(app.config.get("MAX_UCSC_BYTES", DEFAULT_MAX_UCSC_BYTES))
            content_chunks: list[bytes] = []
            total_bytes = 0
            for chunk in response.iter_content(chunk_size=UCSC_CHUNK_SIZE):
                total_bytes += len(chunk)
                if total_bytes > max_bytes:
                    logger.warning(
                        "UCSC response exceeded size limit: %s bytes", total_bytes
                    )
                    response.close()
                    return "Upstream response too large", 502
                content_chunks.append(chunk)

            content = b"".join(content_chunks)
            logger.info("Fetched UCSC image from: %s, size: %s bytes", ucsc_host, total_bytes)

            image_response = make_response(content)
            image_response.headers["Content-Type"] = content_type
            return image_response

        except requests.exceptions.RequestException as e:
            logger.error("Upstream fetch failed for UCSC: %s", e)
            return "Upstream fetch failed", 502
        except Exception as e:
            logger.exception("Internal error in ucsc_image: %s", e)
            return "Internal server error", 500


class UcscProxyServerExtension(MDVProjectServerExtension):
    """
    Server extension that registers the app-wide `/ucsc_proxy` route.
    """

    # Always registered for now. To make optional later:
    # - add "ucsc_proxy" to extension_classes in server_options.py
    # - gate on app.config / MDVServerOptions.extensions
    # - frontend could follow a similar pattern (e.g. get_frontend_extensions()
    #   returning chart modules to register before new ChartManager())

    def register_routes(self, project: Any, project_bp: Any):
        # App-wide endpoint; nothing to register on the per-project blueprint.
        return None

    def mutate_state_json(self, state_json: dict, project: Any, app: Flask):
        return None

    def register_global_routes(self, app: Flask, config: dict):
        register_ucsc_proxy_routes(app)

