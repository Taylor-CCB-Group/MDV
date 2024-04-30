import requests

# Specify the URL where the Flask server is running
UPLOAD_URL = 'http://127.0.0.1:5000/upload'

# Path to the file you want to upload
FILE_PATH = '/Users/jayesh/mdv/pbmc3k/datafile.h5'

# Create a dictionary containing the file data
files = {'file': open(FILE_PATH, 'rb')}

# Send the POST request to upload the file
response = requests.post(UPLOAD_URL, files=files)

# Print the response
print(response.json())