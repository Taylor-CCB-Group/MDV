"""Tests for ChatMDV LLM provider discovery and pairing."""

from unittest.mock import patch

import pytest

from mdvtools.llm import llm_providers


def test_ollama_base_url_defaults_to_host_docker_internal_in_container(monkeypatch):
    monkeypatch.delenv("OLLAMA_BASE_URL", raising=False)
    monkeypatch.setattr(llm_providers, "load_dotenv", lambda: None)
    with patch.object(llm_providers.os.path, "exists", lambda path: path == "/.dockerenv"):
        assert llm_providers._ollama_base_url() == "http://host.docker.internal:11434"


def test_ollama_base_url_explicit_env_overrides_docker_default(monkeypatch):
    monkeypatch.setenv("OLLAMA_BASE_URL", "http://custom:11434")
    monkeypatch.setattr(llm_providers, "load_dotenv", lambda: None)
    with patch.object(llm_providers.os.path, "exists", lambda path: path == "/.dockerenv"):
        assert llm_providers._ollama_base_url() == "http://custom:11434"


def test_discover_models_merges_openai_and_ollama(monkeypatch):
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")

    ollama_tags = {
        "models": [
            {"name": "gemma:26b", "capabilities": ["completion"]},
            {"name": "nomic-embed-text", "capabilities": ["embedding"]},
        ]
    }

    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=ollama_tags["models"]):
        discovered = llm_providers.discover_models()

    chat_ids = {m.id for m in discovered.chat_models}
    embed_ids = {m.id for m in discovered.embedding_models}

    assert "openai:chat:gpt-4.1" in chat_ids
    assert "ollama:chat:gemma:26b" in chat_ids
    assert "openai:embedding:text-embedding-3-large" in embed_ids
    assert "ollama:embedding:nomic-embed-text" in embed_ids
    assert discovered.default_model_id == "openai:chat:gpt-4.1"


def test_discover_models_ollama_only(monkeypatch):
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("CHATMDV_DEFAULT_MODEL", "ollama:chat:llama3.2")

    ollama_tags = [
        {"name": "llama3.2", "capabilities": ["completion"]},
        {"name": "nomic-embed-text", "capabilities": ["embedding"]},
    ]

    with patch.object(llm_providers, "_has_openai_api_key", return_value=False):
        with patch.object(llm_providers, "_fetch_ollama_tags", return_value=ollama_tags):
            discovered = llm_providers.discover_models()

    assert len(discovered.chat_models) == 1
    assert discovered.chat_models[0].model == "llama3.2"
    assert discovered.default_model_id == "ollama:chat:llama3.2"


def test_resolve_chat_model_id_by_full_id(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=[]):
        assert (
            llm_providers.resolve_chat_model_id("openai:chat:gpt-4.1")
            == "openai:chat:gpt-4.1"
        )


def test_resolve_chat_model_id_by_bare_name(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=[]):
        assert llm_providers.resolve_chat_model_id("gpt-4o-mini") == "openai:chat:gpt-4o-mini"


def test_resolve_chat_model_id_empty_uses_default(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=[]):
        assert llm_providers.resolve_chat_model_id(None) is None
        assert llm_providers.resolve_chat_model_id("") is None
        assert llm_providers.resolve_chat_model_id("   ") is None


def test_resolve_chat_model_id_unknown_raises(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=[]):
        with pytest.raises(ValueError, match="Unknown chat model"):
            llm_providers.resolve_chat_model_id("not-a-real-model")


def test_resolve_embedding_model_openai(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    discovered = llm_providers.discover_models()
    chat = discovered.get_chat_model("openai:chat:gpt-4.1")
    embedding = llm_providers.resolve_embedding_model(chat, discovered)
    assert embedding.provider == "openai"
    assert embedding.model == "text-embedding-3-large"


def test_resolve_embedding_model_ollama(monkeypatch):
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("OLLAMA_EMBEDDING_MODEL", "nomic-embed-text")

    ollama_tags = [
        {"name": "gemma:26b", "capabilities": ["completion"]},
        {"name": "nomic-embed-text", "capabilities": ["embedding"]},
    ]

    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=ollama_tags):
        discovered = llm_providers.discover_models()
        chat = discovered.get_chat_model("ollama:chat:gemma:26b")
        embedding = llm_providers.resolve_embedding_model(chat, discovered)

    assert embedding.provider == "ollama"
    assert embedding.model == "nomic-embed-text"


def test_resolve_embedding_model_ollama_missing_raises(monkeypatch):
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("OLLAMA_EMBEDDING_MODEL", "missing-embed-model")

    ollama_tags = [
        {"name": "gemma:26b", "capabilities": ["completion"]},
    ]

    with patch.object(llm_providers, "_fetch_ollama_tags", return_value=ollama_tags):
        discovered = llm_providers.discover_models()
        chat = discovered.get_chat_model("ollama:chat:gemma:26b")
        with pytest.raises(ValueError, match="No Ollama embedding model"):
            llm_providers.resolve_embedding_model(chat, discovered)


def test_rag_retriever_cache_keyed_by_embedding_model(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")

    from mdvtools.llm.langchain_mdv import _get_rag_retriever, _rag_retriever_cache

    _rag_retriever_cache.clear()

    openai_embed = llm_providers.ModelSpec(
        id="openai:embedding:text-embedding-3-large",
        label="OpenAI text-embedding-3-large",
        provider="openai",
        model="text-embedding-3-large",
        kind="embedding",
    )
    other_embed = llm_providers.ModelSpec(
        id="openai:embedding:text-embedding-3-small",
        label="OpenAI text-embedding-3-small",
        provider="openai",
        model="text-embedding-3-small",
        kind="embedding",
    )

    fake_retriever_a = object()
    fake_retriever_b = object()

    with patch("mdvtools.llm.langchain_mdv.crawl_local_repo", return_value=[]):
        with patch("mdvtools.llm.langchain_mdv.create_embeddings", return_value=object()):
            with patch("mdvtools.llm.langchain_mdv.FAISS") as mock_faiss:
                mock_faiss.from_documents.return_value.as_retriever.return_value = fake_retriever_a
                retriever_a = _get_rag_retriever(openai_embed)
                mock_faiss.from_documents.return_value.as_retriever.return_value = fake_retriever_b
                retriever_b = _get_rag_retriever(other_embed)

    assert retriever_a is fake_retriever_a
    assert retriever_b is fake_retriever_b
    assert len(_rag_retriever_cache) == 2
    assert openai_embed.id in _rag_retriever_cache
    assert other_embed.id in _rag_retriever_cache

    _rag_retriever_cache.clear()
