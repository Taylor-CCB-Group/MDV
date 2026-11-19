import { Link, Paper, Typography } from "@mui/material";
import { useQuery } from "@tanstack/react-query";
import clsx from "clsx";
import { Fragment, type KeyboardEvent, type MouseEvent, useEffect, useRef, useState } from "react";
import {
    FEATURE_INFO_SOURCES,
    type FeatureInfoLink,
    type FeatureInfoSource,
    type Organism,
} from "./featureInfoSources";

interface FeatureInfoCardProps {
    id: string;
    source: FeatureInfoSource;
    organism: Organism;
    isHighlighted?: boolean;
    onCardClick?: (event: MouseEvent<HTMLDivElement>) => void;
    onCardKeyDown?: (event: KeyboardEvent<HTMLDivElement>) => void;
}

const clamp = (lines: number) => ({
    overflow: "hidden",
    textOverflow: "ellipsis",
    display: "-webkit-box",
    WebkitLineClamp: lines,
    WebkitBoxOrient: "vertical",
});

function Links({ links }: { links: FeatureInfoLink[] }) {
    return (
        <Typography variant="body2" sx={{ mt: "auto", lineHeight: 1.3 }}>
            {links.map((link, i) => (
                <Fragment key={link.href}>
                    {i > 0 && " · "}
                    <Link
                        href={link.href}
                        target="_blank"
                        rel="noopener noreferrer"
                        tabIndex={-1}
                        sx={{ "&:hover": { fontWeight: 600 } }}
                    >
                        {link.label}
                    </Link>
                </Fragment>
            ))}
        </Typography>
    );
}

/**
 * Card showing what `source` knows about an identifier (gene or protein).
 * Only fetches once scrolled near the viewport.
 */
export function FeatureInfoCard({
    id,
    source,
    organism,
    isHighlighted = false,
    onCardClick,
    onCardKeyDown,
}: FeatureInfoCardProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const [isVisible, setIsVisible] = useState(false);

    useEffect(() => {
        if (!containerRef.current || isVisible) return;
        const observer = new IntersectionObserver(
            (entries) => {
                if (entries[0].isIntersecting) setIsVisible(true);
            },
            { rootMargin: "100px" }, // start loading slightly before it's visible
        );
        observer.observe(containerRef.current);
        return () => observer.disconnect();
    }, [isVisible]);

    const sourceDef = FEATURE_INFO_SOURCES[source];
    const {
        data: info,
        isLoading,
        error,
    } = useQuery({
        queryKey: ["featureInfo", source, organism, id],
        queryFn: ({ signal }) => sourceDef.fetchInfo(id, { organism, signal }),
        enabled: isVisible && !!id,
        retry: false,
        staleTime: Number.POSITIVE_INFINITY,
    });

    const handleCardClick = (event: MouseEvent<HTMLDivElement>) => {
        // links inside the card shouldn't change the selection
        if (event.target instanceof Element && event.target.closest("a")) return;
        onCardClick?.(event);
    };

    const handleKeyDown = (event: KeyboardEvent<HTMLDivElement>) => {
        // Space selects; stop it scrolling the list
        if (event.key === " " || event.key === "Spacebar") event.preventDefault();
        onCardKeyDown?.(event);
    };

    return (
        <Paper
            ref={containerRef}
            variant="outlined"
            tabIndex={onCardClick || onCardKeyDown ? 0 : undefined}
            className={clsx(
                "h-[150px] w-full overflow-hidden flex flex-col rounded-md transition-shadow",
                onCardClick ? "cursor-pointer" : "cursor-default",
                "focus:outline-none focus:ring-2 focus:ring-offset-1",
            )}
            sx={{
                p: 1.5,
                ...(isHighlighted && { borderColor: "rgb(59, 130, 246)" }),
                "&:focus": { borderWidth: "2px" },
            }}
            onClick={handleCardClick}
            onKeyDown={handleKeyDown}
        >
            {!isVisible && (
                <Typography variant="body2" color="text.secondary">
                    Scroll to load info...
                </Typography>
            )}
            {isVisible && isLoading && <Typography variant="body2">Loading info for {id}...</Typography>}
            {isVisible && error && (
                <>
                    <Typography variant="subtitle1" fontWeight={600} sx={{ lineHeight: 1.2 }}>
                        {id}
                    </Typography>
                    <Typography variant="body2" color="error">
                        {error instanceof Error ? error.message : "Unknown error"}
                    </Typography>
                    <Links links={sourceDef.fallbackLinks(id, { organism })} />
                </>
            )}
            {isVisible && info && (
                <>
                    <Typography variant="subtitle1" fontWeight={600} sx={{ lineHeight: 1.2, ...clamp(1) }}>
                        {info.title}
                    </Typography>
                    {info.subtitle && (
                        <Typography variant="body2" color="text.secondary" sx={{ lineHeight: 1.3, ...clamp(1) }}>
                            {info.subtitle}
                            {info.note && (
                                <Typography component="span" variant="body2" color="warning.main">
                                    {` · ${info.note}`}
                                </Typography>
                            )}
                        </Typography>
                    )}
                    {info.facts?.map((fact) => (
                        <Typography key={fact.label} variant="body2" sx={{ lineHeight: 1.3, ...clamp(1) }}>
                            <strong>{fact.label}:</strong> {fact.value}
                        </Typography>
                    ))}
                    {info.description && (
                        <Typography variant="body2" sx={{ mt: 0.5, lineHeight: 1.3, ...clamp(info.facts ? 2 : 3) }}>
                            {info.description}
                        </Typography>
                    )}
                    <Links links={info.links} />
                </>
            )}
        </Paper>
    );
}
