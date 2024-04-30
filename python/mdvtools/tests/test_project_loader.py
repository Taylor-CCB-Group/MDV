import unittest
from unittest.mock import patch, Mock
from flask import Flask
from mdvtools.dbutils.project_loader import load_projects_from_config, create_projects_from_db, create_all_projects

class TestProjectLoader(unittest.TestCase):

    def setUp(self):
        # Create a Flask application context for testing
        self.app = Flask(__name__)
        self.app.testing = True
        self.app_context = self.app.app_context()
        self.app_context.push()

    def tearDown(self):
        # Pop the Flask application context after testing
        self.app_context.pop()

    @patch('mdvtools.dbutils.project_loader.os.path.exists')
    @patch('mdvtools.dbutils.project_loader.open')
    @patch('mdvtools.dbutils.project_loader.MDVProject')
    def test_load_projects_from_config(self, mock_mdvproject, mock_open, mock_exists):
        # Mock the content of the project_config.json file
        mock_exists.return_value = True
        mock_open.return_value.__enter__.return_value.read.return_value = '{"projects": ["project1", "project2"]}'
        mock_mdvproject.side_effect = [Mock(), Mock()]

        # Call the function
        projects = load_projects_from_config('/base/dir')
        print(f"loaded {len(projects)} projects from mock config.")

        # Assert that MDVProject constructor was called twice
        self.assertEqual(mock_mdvproject.call_count, 2)

    @patch('mdvtools.dbutils.project_loader.Project.query')
    @patch('mdvtools.dbutils.project_loader.MDVProject')
    def test_create_projects_from_db(self, mock_mdvproject, mock_query):
        # Mock Project.query result
        mock_query.all.return_value = [Mock(name='project1'), Mock(name='project2')]
        mock_mdvproject.side_effect = [Mock(), Mock()]

        # Call the function
        projects = create_projects_from_db('/base/dir')
        print(f"created {len(projects)} projects from mock db.")

        # Assert that MDVProject constructor was called twice
        self.assertEqual(mock_mdvproject.call_count, 2)


    @patch('mdvtools.dbutils.project_loader.load_projects_from_config')
    @patch('mdvtools.dbutils.project_loader.create_projects_from_db')
    def test_create_all_projects(self, mock_create_from_db, mock_load_from_config):
        # Mock the return values of dependent functions
        mock_create_from_db.return_value = [Mock()]
        mock_load_from_config.return_value = [Mock()]

        # Call the function
        projects = create_all_projects('/base/dir')

        # Assert that the returned list contains one element
        self.assertEqual(len(projects), 2)


if __name__ == '__main__':
    unittest.main()
