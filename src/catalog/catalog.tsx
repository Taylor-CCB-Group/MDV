import { useEffect, useState } from 'react';
import './catalog_index.css';
import ProjectTemplates from './CreateProject';
import { BugPlay } from 'lucide-react';

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
        <div className="container w-32 h-32 p-5 outline rounded-xl">
            <a href={`/project/${id}`}>
                {name}
            </a>
            {import.meta.env.DEV && (<a href={`http://localhost:5170?dir=/project/${id}`}>
            <BugPlay />
            </a>)}
            <button className='text-red-500' onClick={async () => {
                await fetch(`/delete_project/${id}`, { method: 'DELETE' });
                window.location.reload(); //probably should use react-router-dom or generally consider flow here.
            }}>Delete</button>
        </div>
    );
}

const Projects = () => {
    const [projects, setProjects] = useState<ProjectMetadata[]>([]);
    const [error, setError] = useState<string | null>(null);
    const [filter, setFilter] = useState<string | null>(null);
    useEffect(() => {
        (async function () {
            const response = await fetch('/projects');
            const data = await response.json();
            setProjects(data);
        })();
    }, []);
    const filteredProjects = filter ? projects.filter(p => p.name.includes(filter)) : projects;
    return (
        <div className='p-10 outline-dashed rounded-3xl'>
            Filter:
            <input type='text' placeholder='Search projects...' 
            className='p-2 m-8 bg-slate-100 rounded-xl'
            onChange={e => setFilter(e.target.value)} value={filter}></input>
            <div className='grid grid-flow-row grid-cols-8 w-full items-center gap-4'>
                {filteredProjects.map(p => <ProjectTile key={p.id} name={p.name} id={p.id} />)}
                {error && <div className='bg-red-500'>{error}</div>}
            </div>
        </div>
    )
}


export default function App() {
    return (
        <div className='p-10'>
            <h1 className="text-6xl text-center m-10">MDV</h1>
            <ProjectTemplates />
            <Projects />
        </div>
    );
}