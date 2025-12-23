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
 * - Lightness: 0.7 (70% - good visibility)
 * - Chroma: 0.2 (moderate saturation for vibrant but not overly intense colors)
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
    
    // Use standard OKLCH values for better color quality
    // Lightness: 0.7 (70%) provides good visibility
    // Chroma: 0.2 provides vibrant but not overly saturated colors
    const lightness = 0.7;
    const chroma = 0.2;
    
    // Convert OKLCH to RGB
    // oklch2rgb returns RGB values in 0-1 range
    const rgb = oklch2rgb([lightness, chroma, hue]);
    
    // Convert from 0-1 range to 0-255 range and clamp to valid RGB values
    // oklch2rgb may return values out of bounds, so we clamp them
    const r = Math.max(0, Math.min(255, Math.round(rgb[0] * 255)));
    const g = Math.max(0, Math.min(255, Math.round(rgb[1] * 255)));
    const b = Math.max(0, Math.min(255, Math.round(rgb[2] * 255)));
    
    return [r, g, b];
}

