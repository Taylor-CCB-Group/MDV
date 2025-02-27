from flask import Flask
print("safe_mdv_app.py starting")
## why can't gunicorn import our app?
app = Flask(__name__)
if __name__ == '__main__':
    print("runnnig as __main__")
    try:
        from mdvtools.dbutils.mdv_server_app import app as mdv_app
        print("imported mdv_app")
        app = mdv_app
        app.run(host='0.0.0.0', debug=True, port=5055)
    except Exception as e:
        # the idea is that we can serve a fallback app if the main app fails to import
        # this could allow user to report errors, or possibly configure environment variables etc...
        # pending frontend.
        print('Error: %s' % e)
        app = Flask(__name__)
        app.run(port=5055)