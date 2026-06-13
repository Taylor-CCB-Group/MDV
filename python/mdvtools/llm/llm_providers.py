"""LLM provider discovery and factory helpers for ChatMDV."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Literal, Optional

import requests
from dotenv import load_dotenv
from langchain_core.callbacks import BaseCallbackHandler
from langchain_openai import ChatOpenAI, OpenAIEmbeddings

from mdvtools._optional import require_extra

require_extra("app", "langchain_community")

from langchain_community.embeddings import OllamaEmbeddings  # noqa: E402

ProviderKind = Literal["openai", "ollama"]
ModelKind = Literal["chat", "embedding"]

OPENAI_CHAT_MODELS = ("gpt-4.1", "gpt-4o", "gpt-4o-mini")
OPENAI_EMBEDDING_MODEL = "text-embedding-3-large"

DEFAULT_OLLAMA_BASE_URL = "http://localhost:11434"
DEFAULT_OLLAMA_EMBEDDING_MODEL = "nomic-embed-text"

# Common Ollama embedding model name substrings when capabilities are absent.
_EMBEDDING_NAME_HINTS = ("embed", "nomic", "mxbai", "bge")


@dataclass(frozen=True)
class ModelSpec:
    id: str
    label: str
    provider: ProviderKind
    model: str
    kind: ModelKind
    available: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "label": self.label,
            "provider": self.provider,
            "model": self.model,
            "kind": self.kind,
            "available": self.available,
        }


def _model_id(provider: ProviderKind, model: str, kind: ModelKind) -> str:
    return f"{provider}:{kind}:{model}"


def _default_ollama_base_url() -> str:
    """Default Ollama URL; inside Docker, localhost is the container, not the host."""
    if os.path.exists("/.dockerenv"):
        return "http://host.docker.internal:11434"
    return DEFAULT_OLLAMA_BASE_URL


def _ollama_base_url() -> str:
    load_dotenv()
    explicit = os.getenv("OLLAMA_BASE_URL")
    if explicit:
        return explicit.rstrip("/")
    return _default_ollama_base_url()


def _ollama_embedding_model_name() -> str:
    load_dotenv()
    return os.getenv("OLLAMA_EMBEDDING_MODEL", DEFAULT_OLLAMA_EMBEDDING_MODEL)


def _has_openai_api_key() -> bool:
    load_dotenv()
    return bool(os.getenv("OPENAI_API_KEY"))


def _fetch_ollama_tags(base_url: str) -> list[dict[str, Any]]:
    try:
        response = requests.get(f"{base_url}/api/tags", timeout=3)
        response.raise_for_status()
        payload = response.json()
        models = payload.get("models", [])
        if isinstance(models, list):
            return [m for m in models if isinstance(m, dict)]
    except (requests.RequestException, ValueError):
        pass
    return []


def _ollama_capabilities(entry: dict[str, Any]) -> set[str]:
    raw = entry.get("capabilities")
    if isinstance(raw, list):
        return {str(item).lower() for item in raw}
    return set()


def _is_ollama_embedding_model(name: str, entry: dict[str, Any]) -> bool:
    caps = _ollama_capabilities(entry)
    if caps:
        return "embedding" in caps and "completion" not in caps
    lowered = name.lower()
    return any(hint in lowered for hint in _EMBEDDING_NAME_HINTS)


def _is_ollama_chat_model(name: str, entry: dict[str, Any]) -> bool:
    caps = _ollama_capabilities(entry)
    if caps:
        if "embedding" in caps and "completion" not in caps:
            return False
        return "completion" in caps or "vision" in caps or not caps
    return not _is_ollama_embedding_model(name, entry)


def _openai_chat_models() -> list[ModelSpec]:
    if not _has_openai_api_key():
        return []
    return [
        ModelSpec(
            id=_model_id("openai", model, "chat"),
            label=f"OpenAI {model}",
            provider="openai",
            model=model,
            kind="chat",
        )
        for model in OPENAI_CHAT_MODELS
    ]


def _openai_embedding_model() -> ModelSpec:
    return ModelSpec(
        id=_model_id("openai", OPENAI_EMBEDDING_MODEL, "embedding"),
        label=f"OpenAI {OPENAI_EMBEDDING_MODEL}",
        provider="openai",
        model=OPENAI_EMBEDDING_MODEL,
        kind="embedding",
        available=_has_openai_api_key(),
    )


def _ollama_models_from_tags(tags: list[dict[str, Any]]) -> tuple[list[ModelSpec], list[ModelSpec]]:
    chat_models: list[ModelSpec] = []
    embedding_models: list[ModelSpec] = []
    for entry in tags:
        name = entry.get("name")
        if not isinstance(name, str) or not name:
            continue
        if _is_ollama_embedding_model(name, entry):
            embedding_models.append(
                ModelSpec(
                    id=_model_id("ollama", name, "embedding"),
                    label=f"Ollama {name}",
                    provider="ollama",
                    model=name,
                    kind="embedding",
                )
            )
        if _is_ollama_chat_model(name, entry):
            chat_models.append(
                ModelSpec(
                    id=_model_id("ollama", name, "chat"),
                    label=f"Ollama {name}",
                    provider="ollama",
                    model=name,
                    kind="chat",
                )
            )
    return chat_models, embedding_models


@dataclass
class DiscoveredModels:
    chat_models: list[ModelSpec]
    embedding_models: list[ModelSpec]
    default_model_id: Optional[str]

    def all_models(self) -> list[ModelSpec]:
        return [*self.chat_models, *self.embedding_models]

    def get_chat_model(self, model_id: Optional[str]) -> ModelSpec:
        if model_id:
            for spec in self.chat_models:
                if spec.id == model_id:
                    return spec
        if self.default_model_id:
            for spec in self.chat_models:
                if spec.id == self.default_model_id:
                    return spec
        if self.chat_models:
            return self.chat_models[0]
        raise ValueError(
            "No chat models are available. Set OPENAI_API_KEY or start Ollama "
            f"({DEFAULT_OLLAMA_BASE_URL}). In Docker, set OLLAMA_BASE_URL to reach the host."
        )


def discover_models() -> DiscoveredModels:
    """Probe OpenAI env and Ollama /api/tags for available models."""
    chat_models = _openai_chat_models()
    embedding_models: list[ModelSpec] = []
    if _has_openai_api_key():
        embedding_models.append(_openai_embedding_model())

    ollama_tags = _fetch_ollama_tags(_ollama_base_url())
    ollama_chat, ollama_embed = _ollama_models_from_tags(ollama_tags)
    chat_models.extend(ollama_chat)
    embedding_models.extend(ollama_embed)

    load_dotenv()
    default_env = os.getenv("CHATMDV_DEFAULT_MODEL")
    default_model_id: Optional[str] = None
    if default_env:
        for spec in chat_models:
            if spec.id == default_env or spec.model == default_env:
                default_model_id = spec.id
                break
    if default_model_id is None and chat_models:
        default_model_id = chat_models[0].id

    return DiscoveredModels(
        chat_models=chat_models,
        embedding_models=embedding_models,
        default_model_id=default_model_id,
    )


def ensure_any_provider_available() -> None:
    """Raise if neither OpenAI nor Ollama chat models are reachable."""
    discovered = discover_models()
    if not discovered.chat_models:
        raise ValueError(
            "No LLM providers are available. Set OPENAI_API_KEY or start Ollama "
            f"at {_ollama_base_url()}. In Docker, set OLLAMA_BASE_URL to reach the host."
        )


def resolve_embedding_model(
    chat_model: ModelSpec,
    discovered: Optional[DiscoveredModels] = None,
) -> ModelSpec:
    """Pair a chat model with its embedding backend for RAG."""
    if discovered is None:
        discovered = discover_models()

    if chat_model.provider == "openai":
        spec = _openai_embedding_model()
        if not spec.available:
            raise ValueError(
                "OpenAI embeddings require OPENAI_API_KEY. "
                "Select an Ollama chat model for a fully local setup."
            )
        return spec

    preferred_name = _ollama_embedding_model_name()
    for spec in discovered.embedding_models:
        if spec.provider == "ollama" and spec.model == preferred_name:
            return spec
    for spec in discovered.embedding_models:
        if spec.provider == "ollama":
            return spec
    raise ValueError(
        f"No Ollama embedding model is available for RAG. "
        f"Pull one with: ollama pull {preferred_name}"
    )


def create_chat_llm(
    model: ModelSpec,
    callbacks: Optional[list[BaseCallbackHandler]] = None,
) -> ChatOpenAI:
    """Create a LangChain chat model for the given spec."""
    cb = callbacks or []
    if model.provider == "openai":
        if not _has_openai_api_key():
            raise ValueError("OPENAI_API_KEY is not set.")
        return ChatOpenAI(
            temperature=0.1,
            model=model.model,
            callbacks=cb,
        )
    base_url = f"{_ollama_base_url()}/v1"
    return ChatOpenAI(
        temperature=0.1,
        model=model.model,
        base_url=base_url,
        api_key="ollama",
        callbacks=cb,
    )


def create_embeddings(model: ModelSpec):
    """Create a LangChain embeddings instance for the given spec."""
    if model.provider == "openai":
        if not _has_openai_api_key():
            raise ValueError("OPENAI_API_KEY is not set.")
        return OpenAIEmbeddings(model=model.model)
    return OllamaEmbeddings(
        model=model.model,
        base_url=_ollama_base_url(),
    )
