import { useEffect, useState } from 'react';
import './catalog_index.css';

/**
 * We should have better ways of documenting & specifying the types of the data we're working with.
 * We might use zod or similar to validate the data we're working with in front-end.
 */
type ProjectMetadata = {
    id: string;
    name: string;
}

const ProjectTile = ({ name, id }: ProjectMetadata) => {
    return (
        <a href={`/project/${id}`} className="project-tile">
            {name}
        </a>
    );
}


export default function App() {
    const [projects, setProjects] = useState<ProjectMetadata[]>([]);
    const [error, setError] = useState<string | null>(null);
    useEffect(() => {
        (async function () {
            const response = await fetch('/projects');
            const data = await response.json();
            setProjects(data);
        })();
    }, []);
    return (
        <div>
            <h1 className="bg-red-600">App TBD</h1>
            <div className='flex flex-col items-center'>
                <button className='bg-slate-500' onClick={async () => {
                    const response = await fetch('/create_project', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json'
                        },
                        body: JSON.stringify({
                            id: prompt('Enter project name')
                        })
                    });
                    if (response.ok) {
                        location.reload();
                    } else {
                        setError(response.statusText);
                    }
                } }>Create Project</button>
                {projects.map((p, i) => <ProjectTile key={i} name={p.name} id={p.id} />)}
                {error && <div className='bg-red-500'>{error}</div>}
            </div>
        </div>
    );
}