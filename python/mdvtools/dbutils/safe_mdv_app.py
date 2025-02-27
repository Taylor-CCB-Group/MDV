from flask import Flask
print("safe_mdv_app.py starting")

# Import the app at the module level (outside any conditionals)
try:
    from mdvtools.dbutils.mdv_server_app import app
    print("imported mdv_app at module level")
except Exception as e:
    print(f'Error importing mdv_app: {e}')
    app = Flask(__name__)

# Keep the existing __main__ block for direct execution
if __name__ == '__main__':
    print("running as __main__")
    app.run(host='0.0.0.0', debug=True, port=5055)