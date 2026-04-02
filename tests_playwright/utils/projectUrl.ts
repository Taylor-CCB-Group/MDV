import { PATH } from './routes';

export function getProjectUrl(projectId: string): string {
  if (PATH !== '/') {
    // When testing a non-root dashboard alias, use the smart root entrypoint for project mode.
    return `/?dir=/project/${projectId}`;
  } else {
    // For the default app route, use direct project paths.
    return `/project/${projectId}`;
  }
}
