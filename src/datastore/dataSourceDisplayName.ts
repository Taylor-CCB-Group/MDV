/**
 * A datasource has two names. `name` is its identifier: the h5 file, views, links, the API and
 * every lookup in the frontend use it, and it never changes. `display_name`, set by renaming
 * the datasource, is what users see. Anything that shows a datasource's name to a user goes
 * through `getDataSourceDisplayName`; anything that refers to the datasource keeps `name`.
 */

/** The parts of a datasource these helpers read: a `DataStore`, or a `datasources.json` entry. */
type NamedDataSource = { name: string; config?: { display_name?: unknown }; display_name?: unknown };

function usableDisplayName(value: unknown) {
    return typeof value === "string" && value.trim() ? value.trim() : null;
}

/**
 * The name a datasource is shown under: its `display_name` if it has a usable one, otherwise its
 * `name`. Matches `datasource_display_name` in `mdvproject.py`.
 */
export function getDataSourceDisplayName(dataSource: NamedDataSource): string {
    return (
        usableDisplayName(dataSource.config?.display_name) ??
        usableDisplayName(dataSource.display_name) ??
        dataSource.name
    );
}

/**
 * Why a datasource cannot be shown as `newName`, or null if it can. Mirrors the checks in
 * `MDVProject.rename_datasource`, which has the final say: the name must not be empty, and no
 * other datasource may be called or shown as it, ignoring case. Other datasources' identifiers
 * stay reserved after they are renamed.
 */
export function getDataSourceNameProblem(
    newName: string,
    dataSource: NamedDataSource,
    allDataSources: NamedDataSource[],
): string | null {
    const trimmed = newName.trim();
    if (!trimmed) return "Name is required";
    const normalized = trimmed.toLowerCase();
    const clash = allDataSources.some(
        (other) =>
            other.name !== dataSource.name &&
            [other.name, getDataSourceDisplayName(other)].some((taken) => taken.trim().toLowerCase() === normalized),
    );
    return clash ? `A datasource called '${trimmed}' already exists` : null;
}
