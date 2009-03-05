# Python test script for ObitMess
# Ends with test of user input

from xmlrpclib import ServerProxy

# Wait for server to start
import time
time.sleep(2)

# Get server proxy - default port for ObitMess
url = "http://localhost:8777/RPC2"
server = ServerProxy(url)

# Test ping 
answer = server.ping(42)
print "ping answer=",answer

# Test Create message window
task = "myTask"
answer = server.CreateWindow(task)
taskID = answer["taskID"]
print "create answer=",answer

# Test Add message
message = "Something profound"
arg = {"taskID":taskID, "message":message+"\n"}
answer = server.DisplayMessage(arg)
message = "Something else profound"
arg = {"taskID":taskID, "message":message+"\n"}
answer = server.DisplayMessage(arg)

# Status
message = "task gone wonky"
arg = {"taskID":taskID, "status":message}
answer = server.SetStatus(arg)

# Test user input
message = "Enter some input"
arg = {"taskID":taskID, "message":message+"\n"}
answer = server.DisplayMessage(arg)
answer = server.UserResponse(taskID)
print "UserResponse",answer
