export default function ProjectTemplates() {
    return (
        <div className="grid p-5 m-5 outline-dashed rounded-3xl">
            <h2 className="text-4xl text-center">Create Project...</h2>
            <button
                type="button"
                className="w-32 h-32 p-2 m-2 bg-blue-500 text-white rounded-xl"
                onClick={async () => {
                    const response = await fetch("/create_project", {
                        method: "POST",
                        headers: {
                            "Content-Type": "application/json",
                        },
                        body: JSON.stringify({
                            id: prompt("Enter project name"),
                        }),
                    });
                    if (response.ok) {
                        // go to the new project from response link...
                        // location.reload();
                        const data = await response.json();
                        window.location.href = `/project/${data.id}`;
                    } else {
                        alert(response.statusText);
                    }
                }}
            >
                <div className="text-6xl mb-2">+</div>
                <div className="text-sm">Empty Project</div>
            </button>
        </div>
    );
}
