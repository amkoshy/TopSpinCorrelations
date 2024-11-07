import os

resubmitfile = open("Resubmit.sh", "w")

for filename in os.listdir("."):
    if os.path.splitext(filename)[-1] == ".out":
        file = open(filename, "r")
        lines = file.readlines()
        if len(lines) > 1:
            resubmitfile.write(lines[0])
        file.close()

resubmitfile.close()
