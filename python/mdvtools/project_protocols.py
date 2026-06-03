"""Structural typing helpers shared across MDVProject and test fakes."""

from __future__ import annotations

from typing import Any, Protocol


class CreateProjectPromptProject(Protocol):
    """Subset of MDVProject used when building create-project RAG prompts."""

    @property
    def dir(self) -> str: ...

    @property
    def datasources(self) -> list[Any]: ...

    def get_datasource_metadata(self, name: str) -> dict[str, Any]: ...

    def get_datasource_names(self) -> list[str]: ...


class ProjectViewsLike(Protocol):
    """Subset of MDVProject used when patching or preparing generated view code."""

    @property
    def views(self) -> dict: ...


class ProjectViewLookup(Protocol):
    """Subset of MDVProject used when checking whether a view was persisted."""

    def get_view(self, view: str) -> Any | None: ...


class RowsAsColumnsLinkSource(Protocol):
    def get_links(self, datasource: str, filter: str | None = None) -> list | None: ...
