/** An irregular polygon with a handle at each corner, drawn in lucide's stroke style to match the other tools. */
const PolygonIcon = ({ size = 24 }: { size?: number }) => (
    <svg
        width={size}
        height={size}
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        strokeWidth={2}
        strokeLinecap="round"
        strokeLinejoin="round"
        aria-hidden="true"
    >
        <path d="M5.92 6.43L12.08 4.57M15.03 5.71L18.97 12.29M18.24 14.96L10.76 19.04M8.28 18.13L4.72 8.87" />
        <circle cx="4" cy="7" r="2" />
        <circle cx="14" cy="4" r="2" />
        <circle cx="20" cy="14" r="2" />
        <circle cx="9" cy="20" r="2" />
    </svg>
);

export default PolygonIcon;
