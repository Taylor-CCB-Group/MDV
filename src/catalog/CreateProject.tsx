

export default function ProjectTemplates() {
    return (
        <div className="grid p-5 m-5 outline-dashed rounded-3xl">
            <h2 className="text-4xl text-center">Create Project...</h2>
            <button 
            className="p-2 m-2 bg-blue-500 text-white rounded-xl"
            onClick={async () => {
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
                    // go to the new project from response link...
                    location.reload();
                } else {
                    alert(response.statusText);
                }
            }}>Empty project</button>
        </div>
    )
}