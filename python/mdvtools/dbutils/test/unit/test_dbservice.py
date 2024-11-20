import unittest
from unittest.mock import patch, MagicMock
from mdvtools.dbutils.dbservice import ProjectService, FileService
#from mdvtools.dbutils.dbmodels import Project, File
from datetime import datetime

class TestProjectService(unittest.TestCase):
    
    @patch('mdvtools.dbutils.dbmodels.Project.query')
    def test_get_active_projects_success(self, mock_query):
        mock_query.filter_by.return_value.all.return_value = [MagicMock(), MagicMock()]
        result = ProjectService.get_active_projects()
        self.assertEqual(len(result), 2)
    
    @patch('mdvtools.dbutils.dbmodels.db.session.query')
    @patch('mdvtools.dbutils.dbmodels.db.func.max')
    def test_get_next_project_id_success(self, mock_max, mock_query):
        mock_max.return_value = 5
        mock_query.return_value.scalar.return_value = 5
        result = ProjectService.get_next_project_id()
        self.assertEqual(result, 6)
    
    @patch('mdvtools.dbutils.dbmodels.db.session')
    @patch('mdvtools.dbutils.dbmodels.Project')
    def test_add_new_project_success(self, mock_project, mock_session):
        mock_session.add = MagicMock()
        mock_session.commit = MagicMock()
        mock_project.return_value = MagicMock(id=1)
        result = ProjectService.add_new_project("/path/to/project", "Test Project")
        self.assertIsNotNone(result)
        self.assertEqual(result.name, "Test Project")
    
    @patch('project_service.Project.query')
    def test_get_project_by_id_success(self, mock_query):
        mock_query.get.return_value = MagicMock()
        result = ProjectService.get_project_by_id(1)
        self.assertIsNotNone(result)
    
    @patch('project_service.Project.query')
    @patch('project_service.db.session')
    def test_soft_delete_project_success(self, mock_session, mock_query):
        mock_query.get.return_value = MagicMock(is_deleted=False)
        mock_session.commit = MagicMock()
        result = ProjectService.soft_delete_project(1)
        self.assertTrue(result)

    # Failure scenarios

    @patch('project_service.Project.query')
    def test_get_active_projects_failure(self, mock_query):
        mock_query.filter_by.side_effect = Exception("Database error")
        result = ProjectService.get_active_projects()
        self.assertEqual(result, [])
    
    @patch('project_service.db.session.query')
    @patch('project_service.db.func.max')
    def test_get_next_project_id_failure(self, mock_max, mock_query):
        mock_max.side_effect = Exception("Database error")
        mock_query.return_value.scalar.return_value = None
        result = ProjectService.get_next_project_id()
        self.assertIsNone(result)
    
    @patch('project_service.db.session')
    @patch('project_service.Project')
    def test_add_new_project_failure(self, mock_project, mock_session):
        mock_session.add.side_effect = Exception("Database error")
        result = ProjectService.add_new_project("/path/to/project")
        self.assertIsNone(result)
    
    @patch('project_service.Project.query')
    def test_get_project_by_id_failure(self, mock_query):
        mock_query.get.side_effect = Exception("Database error")
        result = ProjectService.get_project_by_id(1)
        self.assertIsNone(result)
    
    @patch('project_service.Project.query')
    @patch('project_service.db.session')
    def test_soft_delete_project_failure(self, mock_session, mock_query):
        mock_query.get.return_value = MagicMock(is_deleted=False)
        mock_session.commit.side_effect = Exception("Database error")
        result = ProjectService.soft_delete_project(1)
        self.assertFalse(result)

class TestFileService(unittest.TestCase):
    
    @patch('project_service.FileService.get_file_by_path_and_project')
    @patch('project_service.db.session')
    def test_add_or_update_file_in_project_success(self, mock_session, mock_get_file):
        mock_get_file.return_value = None
        mock_session.add = MagicMock()
        mock_session.commit = MagicMock()
        result = FileService.add_or_update_file_in_project("file.txt", "/path/to/file.txt", 1)
        self.assertIsNotNone(result)
        self.assertEqual(result.name, "file.txt")
    
    @patch('project_service.FileService.get_file_by_path_and_project')
    @patch('project_service.db.session')
    def test_add_or_update_file_in_project_update_success(self, mock_session, mock_get_file):
        existing_file = MagicMock(name="file.txt", update_timestamp=datetime.now())
        mock_get_file.return_value = existing_file
        mock_session.commit = MagicMock()
        result = FileService.add_or_update_file_in_project("new_name.txt", "/path/to/file.txt", 1)
        self.assertIsNotNone(result)
        self.assertEqual(result.name, "new_name.txt")
    
    @patch('project_service.File.query')
    def test_get_file_by_path_and_project_success(self, mock_query):
        mock_query.filter_by.return_value.first.return_value = MagicMock()
        result = FileService.get_file_by_path_and_project("/path/to/file.txt", 1)
        self.assertIsNotNone(result)
    
    @patch('project_service.File.query')
    def test_file_exists_in_project_success(self, mock_query):
        mock_query.filter_by.return_value.first.return_value = MagicMock()
        result = FileService.file_exists_in_project("/path/to/file.txt", 1)
        self.assertTrue(result)
    
    @patch('project_service.File.query')
    def test_get_files_by_project_success(self, mock_query):
        mock_query.filter_by.return_value.all.return_value = [MagicMock(), MagicMock()]
        result = FileService.get_files_by_project(1)
        self.assertEqual(len(result), 2)
    
    @patch('project_service.File.query')
    @patch('project_service.db.session')
    def test_delete_files_by_project_success(self, mock_session, mock_query):
        mock_query.filter_by.return_value.all.return_value = [MagicMock(), MagicMock()]
        mock_session.delete = MagicMock()
        mock_session.commit = MagicMock()
        result = FileService.delete_files_by_project(1)
        self.assertTrue(result)
    
    @patch('project_service.File.query')
    @patch('project_service.db.session')
    def test_update_file_timestamp_success(self, mock_session, mock_query):
        mock_query.get.return_value = MagicMock()
        mock_session.commit = MagicMock()
        result = FileService.update_file_timestamp(1)
        self.assertTrue(result)

    # Failure scenarios

    @patch('project_service.FileService.get_file_by_path_and_project')
    @patch('project_service.db.session')
    def test_add_or_update_file_in_project_failure(self, mock_session, mock_get_file):
        mock_get_file.return_value = None
        mock_session.add.side_effect = Exception("Database error")
        result = FileService.add_or_update_file_in_project("file.txt", "/path/to/file.txt", 1)
        self.assertIsNone(result)

    @patch('project_service.FileService.get_file_by_path_and_project')
    @patch('project_service.db.session')
    def test_add_or_update_file_in_project_update_failure(self, mock_session, mock_get_file):
        existing_file = MagicMock(name="file.txt", update_timestamp=datetime.now())
        mock_get_file.return_value = existing_file
        mock_session.commit.side_effect = Exception("Database error")
        result = FileService.add_or_update_file_in_project("new_name.txt", "/path/to/file.txt", 1)
        self.assertIsNone(result)
    
    @patch('project_service.File.query')
    def test_get_file_by_path_and_project_failure(self, mock_query):
        mock_query.filter_by.side_effect = Exception("Database error")
        result = FileService.get_file_by_path_and_project("/path/to/file.txt", 1)
        self.assertIsNone(result)
    
    @patch('project_service.File.query')
    def test_file_exists_in_project_failure(self, mock_query):
        mock_query.filter_by.side_effect = Exception("Database error")
        result = FileService.file_exists_in_project("/path/to/file.txt", 1)
        self.assertFalse(result)
    
    @patch('project_service.File.query')
    def test_get_files_by_project_failure(self, mock_query):
        mock_query.filter_by.side_effect = Exception("Database error")
        result = FileService.get_files_by_project(1)
        self.assertEqual(result, [])
    
    @patch('project_service.File.query')
    @patch('project_service.db.session')
    def test_delete_files_by_project_failure(self, mock_session, mock_query):
        mock_query.filter_by.return_value.all.return_value = [MagicMock(), MagicMock()]
        mock_session.delete.side_effect = Exception("Database error")
        result = FileService.delete_files_by_project(1)
        self.assertFalse(result)
    
    @patch('project_service.File.query')
    @patch('project_service.db.session')
    def test_update_file_timestamp_failure(self, mock_session, mock_query):
        mock_query.get.return_value = MagicMock()
        mock_session.commit.side_effect = Exception("Database error")
        result = FileService.update_file_timestamp(1)
        self.assertFalse(result)

if __name__ == '__main__':
    unittest.main()
