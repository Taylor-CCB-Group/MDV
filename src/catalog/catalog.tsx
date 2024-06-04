import { useEffect, useState } from 'react';
import './catalog_index.css';

const ProjectTile = ({ name }: { name: string }) => {
    return (
        <a href={`/project/${name}`} className="project-tile">
            {name}
        </a>
    );
}


export default function App() {
    const [projects, setProjects] = useState<string[]>([]);
    const [error, setError] = useState<string | null>(null);
    useEffect(() => {
        (async function () {
            const response = await fetch('/projects');
            const data = await response.json();
            setProjects(data);
        })();
    });
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
                {projects.map((p, i) => <ProjectTile key={i} name={p} />)}
                {error && <div className='bg-red-500'>{error}</div>}
            </div>
        </div>
    );
}