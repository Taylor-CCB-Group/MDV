from flask import Flask, request, jsonify
from flask import send_from_directory
import psycopg2
import os

import psycopg2._psycopg


app = Flask(__name__)

# Database connection details
DATABASE_URL = os.getenv(
    "DATABASE_URL", "postgres://postgres@localhost:5432/mydatabase"
)


def get_db_connection():
    conn = psycopg2.connect(DATABASE_URL)
    return conn


@app.route("/")
def index():
    return send_from_directory(".", "index.html")


@app.route("/projects", methods=["GET"])
def get_projects():
    # Example projects, replace with actual query to fetch projects
    projects = ["Project1", "Project2", "Project3"]
    return jsonify(projects)


@app.route("/execute_query", methods=["POST"])
def execute_query():
    if not request.json:
        return jsonify({"error": "No data provided"}), 500
    query_type = request.json.get("query_type")
    print(f"Received query type: {query_type}")
    sql_query = translate_query_type_to_sql(query_type)
    print(f"Translated SQL query: {sql_query}")
    if not sql_query:
        return jsonify({"error": "Invalid query type provided"}), 400

    conn: psycopg2._psycopg.connection | None = None
    cur: psycopg2._psycopg.cursor | None = None
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        cur.execute(sql_query)
        if sql_query.strip().lower().startswith("select"):
            results = cur.fetchall()
            if not cur.description:
                raise Exception("Query results not available")
            column_names = [desc[0] for desc in cur.description]
            data = [dict(zip(column_names, row)) for row in results]
            print(f"Query results: {data}")
            return jsonify({"data": data})
        else:
            conn.commit()
            return jsonify({"message": "Query executed successfully"})
    except Exception as e:
        print(f"Error executing query: {e}")
        return jsonify({"error": str(e)}), 500
    finally:
        if cur:
            cur.close()
        if conn:
            conn.close()


def translate_query_type_to_sql(query_type):
    if query_type == "list_projects":
        return "SELECT * FROM project"
    elif query_type == "list_files":
        return (
            "SELECT * FROM file"  # Adjust the SQL query based on your table structure
        )
    # Add more query types as needed
    return None


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5085, debug=True)
