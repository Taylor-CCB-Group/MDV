export type ColumnMetadataLike = Record<string, unknown> & {
    field?: string;
    deleted?: boolean;
    data?: unknown;
    buffer?: unknown;
    getValue?: unknown;
    originalData?: unknown;
};

export function isSoftDeletedColumn(
    column: ColumnMetadataLike | null | undefined,
): boolean {
    return Boolean(column?.deleted);
}

export function toColumnMetadata(
    column: ColumnMetadataLike,
): ColumnMetadataLike {
    const metadata = { ...column };
    delete metadata.data;
    delete metadata.buffer;
    delete metadata.getValue;
    delete metadata.originalData;
    if (isSoftDeletedColumn(metadata)) {
        metadata.deleted = true;
    } else {
        delete metadata.deleted;
    }
    return metadata;
}

export function normalizeColumnsMetadata(
    columns: ColumnMetadataLike[] | undefined,
): ColumnMetadataLike[] {
    return (columns ?? []).map((column) => toColumnMetadata(column));
}

export function findColumnMetadataIndex(
    columns: ColumnMetadataLike[] | undefined,
    field: string,
): number {
    return (columns ?? []).findIndex((item) => item.field === field);
}

export function updateColumnMetadata(
    columns: ColumnMetadataLike[] | undefined,
    column: ColumnMetadataLike,
): ColumnMetadataLike[] {
    const metadata = toColumnMetadata(column);
    const nextColumns = [...(columns ?? [])];
    const index = findColumnMetadataIndex(nextColumns, String(metadata.field));

    if (index === -1) {
        nextColumns.push(metadata);
    } else {
        nextColumns[index] = {
            ...metadata,
            deleted: false,
        };
    }

    return normalizeColumnsMetadata(nextColumns);
}

export function removeColumnMetadata(
    columns: ColumnMetadataLike[] | undefined,
    field: string,
): ColumnMetadataLike[] {
    return normalizeColumnsMetadata(
        (columns ?? []).filter((column) => column.field !== field),
    );
}
