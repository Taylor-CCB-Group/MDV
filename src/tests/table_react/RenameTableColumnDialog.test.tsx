import RenameTableColumnDialog from "@/react/components/RenameTableColumnDialog";
import { fireEvent, render, screen } from "@testing-library/react";
import type { ReactNode } from "react";
import { describe, expect, test, vi } from "vitest";

vi.mock("@/catalog/ProjectRenameModal", () => ({
    DialogCloseIconButton: ({ onClose }: { onClose: () => void }) => (
        <button onClick={onClose} type="button">
            Close
        </button>
    ),
}));

vi.mock("@mui/material", () => ({
    Button: ({
        children,
        onClick,
        color,
    }: {
        children: ReactNode;
        onClick?: () => void;
        color?: string;
    }) => (
        <button type="button" data-color={color} onClick={onClick}>
            {children}
        </button>
    ),
    Dialog: ({
        children,
        open,
    }: {
        children: ReactNode;
        open: boolean;
    }) => (open ? <div>{children}</div> : null),
    DialogActions: ({ children }: { children: ReactNode }) => <div>{children}</div>,
    DialogContent: ({ children }: { children: ReactNode }) => <div>{children}</div>,
    DialogTitle: ({ children }: { children: ReactNode }) => <div>{children}</div>,
    Stack: ({ children }: { children: ReactNode }) => <div>{children}</div>,
    TextField: ({
        label,
        value,
        onChange,
        onKeyDown,
        helperText,
    }: {
        label: string;
        value: string;
        onChange: (event: { target: { value: string } }) => void;
        onKeyDown?: (event: { key: string; preventDefault: () => void }) => void;
        helperText?: string | null;
    }) => (
        <label>
            {label}
            <input
                aria-label={label}
                value={value}
                onChange={(event) => onChange({ target: { value: event.target.value } })}
                onKeyDown={(event) =>
                    onKeyDown?.({
                        key: event.key,
                        preventDefault: () => event.preventDefault(),
                    })}
            />
            {helperText ? <span>{helperText}</span> : null}
        </label>
    ),
    Typography: ({ children }: { children: ReactNode }) => <div>{children}</div>,
}));

describe("RenameTableColumnDialog", () => {
    test("prefills the current column name", () => {
        render(
            <RenameTableColumnDialog
                open
                columnField="age"
                initialName="Age"
                onClose={vi.fn()}
                onSubmit={vi.fn()}
            />,
        );

        expect(screen.getByDisplayValue("Age")).toBeDefined();
    });

    test("submits a valid rename", () => {
        const onSubmit = vi.fn();
        render(
            <RenameTableColumnDialog
                open
                columnField="age"
                initialName="Age"
                onClose={vi.fn()}
                onSubmit={onSubmit}
            />,
        );

        fireEvent.change(screen.getByLabelText("Column Name"), {
            target: { value: "Age label" },
        });
        fireEvent.click(screen.getByRole("button", { name: "Rename" }));

        expect(onSubmit).toHaveBeenCalledWith({
            columnField: "age",
            newName: "Age label",
        });
    });

    test("closes without submitting when the name is unchanged", () => {
        const onClose = vi.fn();
        const onSubmit = vi.fn();
        render(
            <RenameTableColumnDialog
                open
                columnField="age"
                initialName="Age"
                onClose={onClose}
                onSubmit={onSubmit}
            />,
        );

        fireEvent.click(screen.getByRole("button", { name: "Rename" }));

        expect(onClose).toHaveBeenCalledTimes(1);
        expect(onSubmit).not.toHaveBeenCalled();
    });
});
