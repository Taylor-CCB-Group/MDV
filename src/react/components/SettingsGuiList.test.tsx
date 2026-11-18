import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, test, vi } from "vitest";
import type { AnyGuiSpec } from "@/charts/charts";
import { g } from "@/lib/utils";
import { SettingsGuiList } from "./SettingsDialogComponent";

describe("SettingsGuiList", () => {
    test("controls show what they wrote, as in the settings dialog", () => {
        const func = vi.fn();
        const specs: AnyGuiSpec[] = [
            g({
                type: "radiobuttons",
                label: "Action on Filter",
                current_value: "hide",
                choices: [
                    ["Hide Points", "hide"],
                    ["Gray Out Points", "grey"],
                ],
                func,
            }),
        ];
        render(<SettingsGuiList specs={specs} />);
        const radio = (label: string) => screen.getByLabelText<HTMLInputElement>(label);
        const grey = radio("Gray Out Points");
        expect(grey.checked).toBe(false);
        fireEvent.click(grey);
        expect(func).toHaveBeenCalledWith("grey");
        // re-rendered from the written `current_value`, not left on the old choice
        expect(radio("Gray Out Points").checked).toBe(true);
        expect(radio("Hide Points").checked).toBe(false);
    });
});
