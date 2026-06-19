export type TableAssociation =
    | { status: "none" }
    | { status: "ambiguous" }
    | { status: "resolved"; tableName?: string };

export const NO_TABLE_ASSOCIATION: TableAssociation = { status: "none" };
