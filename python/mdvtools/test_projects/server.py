from flask import Flask, request, jsonify
import os
import psycopg2
from datetime import datetime

app = Flask(__name__)

UPLOAD_FOLDER = 'uploads'
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)


def get_db_connection():
    return psycopg2.connect(
        dbname="mydatabase",
        user="postgres",
        host="localhost",
        port="5432" 
    )


@app.route('/')
def index():
    return open('index.html').read()


@app.route('/upload', methods=['POST'])
def upload():
    try:
        file = request.files['file']
        if not file:
            return jsonify({'error': 'No file selected.'}), 400
        
        
        conn = get_db_connection()
        cur = conn.cursor()

        try:
            # Save the file to the uploads directory
            file_path = os.path.join(UPLOAD_FOLDER, file.filename)
            file.save(file_path)

            # Insert upload log into database
            insert_upload_log(cur, file.filename, file_path)

            # Commit transaction
            conn.commit()

            return jsonify({'message': 'File uploaded successfully.'}), 200
        except Exception as e:
            # Rollback transaction on error
            conn.rollback()
            return jsonify({'error': str(e)}), 500
        finally:
            # Close database connection
            cur.close()
            conn.close()
    except Exception as e:
        return jsonify({'error': str(e)}), 500


def insert_upload_log(cur, file_name, file_path):
    cur.execute(
        "INSERT INTO upload_logs (file_name, file_path, upload_timestamp) VALUES (%s, %s, %s)",
        (file_name, file_path, datetime.now())
    )

if __name__ == '__main__':
    app.run(debug=True, port=5054)
