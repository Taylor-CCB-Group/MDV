import pytest
from datetime import datetime
from flask_testing import TestCase
from mdvtools.dbutils.dbmodels import db, Project, File
from mdvtools.dbutils.dbservice import ProjectService, FileService
from mdvtools.dbutils.mdv_server_app import create_flask_app  # Updated import path

class ServiceTestCase(TestCase):
    def create_app(self):
        # Create app with test configuration
        app = create_flask_app(config_name='test')
        return app

    def setUp(self):
        """Setup the test environment"""
        with self.app.app_context():
            # Create all tables for the test database
            db.create_all()

    def tearDown(self):
        """Teardown the test environment"""
        with self.app.app_context():
            # Remove all sessions and drop all tables
            db.session.remove()
            db.drop_all()

    def test_get_active_projects(self):
        # Add sample projects
        project1 = Project(name='Project 1', path='/path1', is_deleted=False)
        project2 = Project(name='Project 2', path='/path2', is_deleted=True)
        db.session.add(project1)
        db.session.add(project2)
        db.session.commit()

        # Call service method
        active_projects = ProjectService.get_active_projects()

        self.assertEqual(len(active_projects), 1)
        self.assertEqual(active_projects[0].name, 'Project 1')

    def test_get_next_project_id(self):
        # Add a sample project
        project = Project(name='Sample Project', path='/sample/path')
        db.session.add(project)
        db.session.commit()

        # Call service method
        next_id = ProjectService.get_next_project_id()

        self.assertEqual(next_id, 2)  # Since 1 is the existing ID

    def test_add_new_project(self):
        # Call service method
        new_project = ProjectService.add_new_project('/new/path', 'New Project')

        self.assertIsNotNone(new_project)
        self.assertEqual(new_project.name, 'New Project')
        self.assertEqual(new_project.path, '/new/path')

    def test_get_project_by_id(self):
        # Add a sample project
        project = Project(name='Find Me', path='/find/me')
        db.session.add(project)
        db.session.commit()

        # Call service method
        found_project = ProjectService.get_project_by_id(project.id)

        self.assertIsNotNone(found_project)
        self.assertEqual(found_project.name, 'Find Me')

    def test_soft_delete_project(self):
        # Add a sample project
        project = Project(name='Delete Me', path='/delete/me')
        db.session.add(project)
        db.session.commit()

        # Call service method
        success = ProjectService.soft_delete_project(project.id)

        self.assertTrue(success)
        project = Project.query.get(project.id)
        self.assertTrue(project.is_deleted)
        self.assertIsNotNone(project.deleted_timestamp)

    def test_add_or_update_file_in_project(self):
        # Add a sample project
        project = Project(name='Project Files', path='/files/project')
        db.session.add(project)
        db.session.commit()

        # Call service method to add a file
        file = FileService.add_or_update_file_in_project('file1.txt', '/files/project/file1.txt', project.id)

        self.assertIsNotNone(file)
        self.assertEqual(file.name, 'file1.txt')

        # Call service method to update the file
        updated_file = FileService.add_or_update_file_in_project('file1_renamed.txt', '/files/project/file1.txt', project.id)

        self.assertIsNotNone(updated_file)
        self.assertEqual(updated_file.name, 'file1_renamed.txt')

    def test_get_file_by_path_and_project(self):
        # Add a sample project and file
        project = Project(name='Project Files', path='/files/project')
        db.session.add(project)
        db.session.commit()
        
        file = File(
            name='file1.txt',
            file_path='/files/project/file1.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        db.session.add(file)
        db.session.commit()

        # Call service method
        found_file = FileService.get_file_by_path_and_project('/files/project/file1.txt', project.id)

        self.assertIsNotNone(found_file)
        self.assertEqual(found_file.name, 'file1.txt')

    def test_delete_files_by_project(self):
        # Add a sample project and files
        project = Project(name='Project Files', path='/files/project')
        db.session.add(project)
        db.session.commit()

        file1 = File(
            name='file1.txt',
            file_path='/files/project/file1.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        file2 = File(
            name='file2.txt',
            file_path='/files/project/file2.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        db.session.add(file1)
        db.session.add(file2)
        db.session.commit()

        # Call service method
        success = FileService.delete_files_by_project(project.id)

        self.assertTrue(success)
        files = File.query.filter_by(project_id=project.id).all()
        self.assertEqual(len(files), 0)

    def test_update_file_timestamp(self):
        # Add a sample project and file
        project = Project(name='Project Files', path='/files/project')
        db.session.add(project)
        db.session.commit()

        file = File(
            name='file1.txt',
            file_path='/files/project/file1.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        db.session.add(file)
        db.session.commit()

        # Call service method
        success = FileService.update_file_timestamp(file.id)

        self.assertTrue(success)
        updated_file = File.query.get(file.id)
        self.assertGreater(updated_file.update_timestamp, file.update_timestamp)
    
    def test_file_exists_in_project(self):
        # Add a sample project and file
        project = Project(name='Project Files', path='/files/project')
        db.session.add(project)
        db.session.commit()

        file = File(
            name='file1.txt',
            file_path='/files/project/file1.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        db.session.add(file)
        db.session.commit()

        # Call service method to check if file exists
        file_exists = FileService.file_exists_in_project('/files/project/file1.txt', project.id)
        self.assertTrue(file_exists)

        # Check for a file that doesn't exist
        file_does_not_exist = FileService.file_exists_in_project('/files/project/nonexistent.txt', project.id)
        self.assertFalse(file_does_not_exist)

    def test_get_files_by_project(self):
        # Add a sample project and files
        project = Project(name='Project Files', path='/files/project')
        db.session.add(project)
        db.session.commit()

        file1 = File(
            name='file1.txt',
            file_path='/files/project/file1.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        file2 = File(
            name='file2.txt',
            file_path='/files/project/file2.txt',
            project_id=project.id,
            upload_timestamp=datetime.now(),
            update_timestamp=datetime.now()
        )
        db.session.add(file1)
        db.session.add(file2)
        db.session.commit()

        # Call service method to get files by project
        files = FileService.get_files_by_project(project.id)

        self.assertEqual(len(files), 2)
        self.assertEqual(files[0].name, 'file1.txt')
        self.assertEqual(files[1].name, 'file2.txt')

        # Test with a project that has no files
        empty_project = Project(name='Empty Project', path='/empty/project')
        db.session.add(empty_project)
        db.session.commit()

        no_files = FileService.get_files_by_project(empty_project.id)
        self.assertEqual(len(no_files), 0)

if __name__ == '__main__':
    pytest.main()
