import io
import sys
from http.server import HTTPServer, BaseHTTPRequestHandler
import MolDisplay

# Define our own handler class that extends BaseHTTPRequestHandler
class MyHandler(BaseHTTPRequestHandler):
  
    # Handle GET requests
    def do_GET(self):
        # If the path is just "/", send the HTML form to upload SDF files
        if self.path == "/":
            # Send a 200 OK response and set the Content-type header to text/html
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            # Set the length of the response body
            self.send_header("Content-length", len(home_page))
            self.end_headers()
            # Write the HTML form to the response body
            self.wfile.write(bytes(home_page, "utf-8"))
        # If the path is anything else, send a 404 Not Found error
        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(bytes("404: not found", "utf-8"))
  
    # Handle POST requests
    def do_POST(self):
        # If the path is "/molecule", handle the uploaded SDF file and send back an SVG file
        if self.path == "/molecule":
            # Get the length of the uploaded file from the Content-length header
            content_length = int(self.headers.get("Content-length"))
            # Read the file data from the request body
            file_data = self.rfile.read(content_length)
            # Convert the file data into a TextIoWrapper object
            bytes_io = io.BytesIO(file_data)
            text_io = io.TextIOWrapper(bytes_io)
            # Skip the first four lines of header information
            text_io.readline()
            text_io.readline()
            text_io.readline()
            text_io.readline()
            # Parse the SDF file and generate an SVG string
            mol = MolDisplay.Molecule()
            mol.parse(text_io)
            mol.sort()
            svg_string = mol._svg_()
            # Encode the SVG string as bytes
            svg_bytes = svg_string.encode()
            # Send a 200 OK response and set the Content-type header to image/svg+xml
            self.send_response(200)
            self.send_header("Content-type", "image/svg+xml")
            # Set the length of the response body
            self.send_header("Content-length", len(svg_bytes))
            self.end_headers()
            # Write the SVG data to the response body
            self.wfile.write(svg_bytes)
        # If the path is anything else, send a 404 Not Found error
        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(bytes("404: not found", "utf-8"))

# Define the HTML form to be sent in the response to GET requests for "/"
home_page = """
<html>
<head>
<title> File Upload </title>
</head>
<body>
<h1> File Upload </h1>
<form action="molecule" enctype="multipart/form-data" method="post">
<p>
<input type="file" id="sdf_file" name="filename"/>
</p>
<p>
<input type="submit" value="Upload"/>
</p>
</form>
</body>
</html>
"""

# Define the HTTP server to listen on the specified port and use our handler
httpd = HTTPServer(('localhost', int(sys.argv[1])), MyHandler)

# Start the server and keep it running indefinitely
httpd.serve_forever()

