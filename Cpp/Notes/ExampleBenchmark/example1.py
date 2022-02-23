# This is a python script which sends commands to the shell.
# You can run this with `python3 example1.py`


# 1) This is how you print formatted strings:
print("Part 1:")
name = "User"
print("Hello, {}! Welcome to python!\n\n".format(name)) # the {} is a place holder


# This is the command to import the OS module.
# You need this to be able to send commands to the shell.
# You don't need to install this module, it comes with Python 3.
# You can read the doc https://docs.python.org/3/library/os.html for details.
import os


# Here is how you send commands to the shell using the OS Module:
print("Part 2:")
os_output = os.system('ls -l') # run command
print("OS output type: {}".format(type(os_output)))
print("Return code: {}\n\n".format(os_output))


# It's actualy better to use the Subprocess module for running C programs.
# It let's you capture the output.
import subprocess

# Here's the syntax for running a command.
print("Part 3:")
sp_output = subprocess.run(["ls","-l"],capture_output=1,universal_newlines=1)
print("SP output type: {}\n\n".format(type(sp_output))) # output is an object

# This is the raw output.
print("Part 4:")
print("Raw output: \n{}".format(sp_output.stdout)) # output is captured as a string in .stdout