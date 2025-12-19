import type { FieldName } from "@/charts/charts";
import { oklch2rgb } from "@/utilities/oklch2rgb";

/**
 * Simple hash function to convert a string to a number.
 * This provides a stable, deterministic mapping from field identifiers to numbers.
 */
function hashString(str: string): number {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
        const char = str.charCodeAt(i);
        hash = ((hash << 5) - hash) + char;
        hash = hash & hash; // Convert to 32-bit integer
    }
    return Math.abs(hash);
}

/**
 * Returns a stable RGB color for a given field identifier.
 * The same field will always get the same color, regardless of order or other fields.
 * 
 * Uses OKLCH color space with:
 * - Lightness: 200 (0-255 scale, converted to 0-1 for OKLCH)
 * - Chroma: 220 (0-255 scale, converted to 0-1 for OKLCH)
 * - Hue: Determined by hashing the field identifier, distributed across 0-360 degrees
 * 
 * @param field - The field identifier (FieldName) to get a color for
 * @returns RGB color as [r, g, b] array with values 0-255
 */
export function getFieldColor(field: FieldName): [number, number, number] {
    // Hash the field identifier to get a stable number
    const hash = hashString(field);
    
    // Map the hash to a hue value (0-360 degrees)
    // Using modulo to wrap around the hue circle
    const hue = hash % 360;
    
    // Use consistent lightness and chroma values matching the original implementation
    // Original used [200, 220, ...] - matching that pattern exactly
    // Based on the original code: oklch2rgb([200, 220, 360 * index/n])
    // Note: The scale for L and C in oklch2rgb appears non-standard (using 0-255 instead of 0-1)
    // but we match the original pattern to maintain consistency
    const lightness = 200;
    const chroma = 220;
    
    // Convert OKLCH to RGB - matching the original pattern exactly
    const rgb = oklch2rgb([lightness, chroma, hue]);
    
    // The original code uses the result directly: colorRange: [oklch2rgb(...)]
    // oklch2rgb may return values out of bounds, so we clamp them
    // Handle both 0-1 and 0-255 output scales
    const r = Math.max(0, Math.min(255, Math.round(rgb[0] > 1 ? rgb[0] : rgb[0] * 255)));
    const g = Math.max(0, Math.min(255, Math.round(rgb[1] > 1 ? rgb[1] : rgb[1] * 255)));
    const b = Math.max(0, Math.min(255, Math.round(rgb[2] > 1 ? rgb[2] : rgb[2] * 255)));
    
    return [r, g, b];
}

