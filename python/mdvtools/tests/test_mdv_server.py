import os
import tempfile
import unittest
from flask import json
from mdvtools.dbutils.mdv_server_old import app, db, Project, File

class TestMdvServer(unittest.TestCase):

    def setUp(self):
        app.config['TESTING'] = True
        app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///:memory:'
        self.app = app.test_client()
        with app.app_context():
            db.create_all()

    def tearDown(self):
        with app.app_context():
            db.session.remove()
            db.drop_all()

    def test_upload_file(self):
        with app.app_context():
            project_name = 'Test Project'
            file_name = 'test_file.txt'
            file_content = b'This is a test file content.'

            # Upload a file
            response = self.upload_file(project_name, file_name, file_content)
            self.assertEqual(response.status_code, 200)

            # Check if the file is uploaded
            project = Project.query.filter_by(name=project_name).first()
            self.assertIsNotNone(project)
            file = File.query.filter_by(name=file_name, project_id=project.id).first()
            self.assertIsNotNone(file)
            self.assertTrue(os.path.exists(file.file_path))

            # Test missing project name
            response = self.app.post('/upload', data={'file': (file_content, file_name)})
            #print(response)
            self.assertEqual(response.status_code, 400)
            data = json.loads(response.data)
            self.assertEqual(data['error'], 'Project name is missing.')

            # Test missing file
            response = self.app.post('/upload', data={'project_name': project_name})
            #print(response)
            self.assertEqual(response.status_code, 400)
            #data = json.loads(response.data)
            #self.assertEqual(data['error'], 'No file selected.')


    def test_get_projects(self):
        with app.app_context():
            # Create projects
            projects = ['pbmc3k', 'pbmc3k_example']
            for name in projects:
                project = Project(name=name)
                db.session.add(project)
            db.session.commit()

            # Get projects
            response = self.app.get('/projects')
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertEqual(sorted(data), sorted(projects))

    def upload_file(self, project_name, file_name, file_content):
        with tempfile.NamedTemporaryFile() as temp_file:
            temp_file.write(file_content)
            temp_file.seek(0)
            return self.app.post('/upload',
                                 data={'project_name': project_name, 'file': (temp_file, file_name)},
                                 content_type='multipart/form-data')

if __name__ == '__main__':
    unittest.main()
