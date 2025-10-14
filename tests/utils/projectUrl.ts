import { PATH } from './routes';

export function getProjectUrl(projectId: string): string {
  if (PATH !== '/') {
    // For /catalog_dev, use query param format at root
    return `/?dir=/project/${projectId}`;
  } else {
    // For root, use direct path
    return `/project/${projectId}`;
  }
}
