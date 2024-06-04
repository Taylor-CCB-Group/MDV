//console.log('this is the place to build the catalog/project browse interface...')

(async function () {
    const response = await fetch('/projects');
    const data = await response.json();
    console.log(data);
    const container = document.createElement('div');
    container.className = 'projects';
    document.body.appendChild(container);
    for (const p of data) {
        const link = document.createElement('a');
        link.href = `/project/${p}`;
        link.innerText = p;
        container.appendChild(link);
    }
})();
const button = document.createElement('button');
button.innerText = 'Create Project';
button.onclick = async () => {
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
        alert('Failed to create project');
    }
};
document.body.appendChild(button);