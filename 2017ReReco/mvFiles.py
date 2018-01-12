#!/evn/pyton

from os import listdir
from os.path import isfile, join
import shutil

mypath="./"
newpath="./OldFiles"
onlyfiles = [f for f in listdir(mypath) if (isfile(join(mypath, f))) and f.endswith("root") ]

for my_file in onlyfiles:
    folder=my_file.split("_")[1]
    source=join(mypath,folder,my_file)
    destination=join(newpath,my_file)
    print "moving",source,"to",destination
    if isfile(source):
        shutil.move(source,destination)
    else:
        print source,"does not exist, moving on!"
