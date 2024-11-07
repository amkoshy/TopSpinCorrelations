import subprocess

problemlogs = subprocess.check_output("grep -l 'expired\|Error in <TBasket::ReadBasketBuffers>\|Error in <TNetXNGFile::Open>\|Error in <TNetXNGFile::ReadBuffer>' nohuplogs/*.out", shell=True).split('\n')
#print(logs)

jobs = []

for log in problemlogs:

    if log == '': continue

    items =  subprocess.check_output("grep -h '"+log+"' scripts/*", shell=True).split('\n')

    for item in items:

        if item == '': continue
        else: jobs.append(item)

#print(jobs)



file = open("Resubmit.sh", "w")

for job in jobs:

    file.write(job)
    file.write("\n")

file.close()
