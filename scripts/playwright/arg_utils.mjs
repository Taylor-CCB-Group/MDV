const FLAGS_WITH_VALUES = new Set([
    "--project",
    "-p",
    "--workers",
    "-j",
    "--grep",
    "-g",
    "--grep-invert",
    "--reporter",
    "--config",
    "--retries",
    "--shard",
    "--timeout",
    "--trace",
    "--output",
]);

export function hasPositionalArgs(args) {
    for (let index = 0; index < args.length; index += 1) {
        const arg = args[index];
        if (arg === "--") {
            return index < args.length - 1;
        }
        if (!arg.startsWith("-")) {
            return true;
        }
        if (FLAGS_WITH_VALUES.has(arg)) {
            index += 1;
        }
    }
    return false;
}
