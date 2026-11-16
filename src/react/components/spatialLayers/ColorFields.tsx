import { Slider, Typography } from "@mui/material";

export type Rgb = [number, number, number];
export type Rgba = [number, number, number, number];

export const rgbToHex = (value: readonly number[]) =>
    `#${value
        .slice(0, 3)
        .map((channel) =>
            Math.round(Math.max(0, Math.min(255, channel)))
                .toString(16)
                .padStart(2, "0"),
        )
        .join("")}`;

export const hexToRgb = (hex: string): Rgb => [
    Number.parseInt(hex.slice(1, 3), 16),
    Number.parseInt(hex.slice(3, 5), 16),
    Number.parseInt(hex.slice(5, 7), 16),
];

/** A colour swatch, with an alpha slider (0–255) beside it when `alpha` is given. */
export function ColorFields({
    label,
    rgb,
    alpha,
    onChange,
    caption,
}: {
    label: string;
    rgb: Rgb;
    alpha?: number;
    onChange: (rgb: Rgb, alpha: number | undefined) => void;
    caption?: string;
}) {
    return (
        <div className="grid gap-1">
            <Typography variant="caption">{label}</Typography>
            <div className="flex items-center gap-3">
                <input
                    type="color"
                    aria-label={label}
                    value={rgbToHex(rgb)}
                    onChange={(event) => onChange(hexToRgb(event.target.value), alpha)}
                />
                {alpha !== undefined && (
                    <Slider
                        size="small"
                        min={0}
                        max={255}
                        step={1}
                        value={alpha}
                        aria-label={`${label} opacity`}
                        onChange={(_, next) => {
                            if (typeof next === "number") onChange(rgb, next);
                        }}
                    />
                )}
            </div>
            {caption && (
                <Typography variant="caption" color="text.secondary">
                    {caption}
                </Typography>
            )}
        </div>
    );
}
