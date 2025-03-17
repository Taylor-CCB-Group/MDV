import type { useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers";

export type RowsAsColsProps<T extends CTypes, M extends boolean> = ReturnType<typeof useRowsAsColumnsLinks>[0] & ColumnSelectionProps<T, M>;
