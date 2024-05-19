import os 

dirs = {
    "src":os.path.abspath(os.path.dirname(os.path.abspath(__file__))+"/")
}

dirs["root"] = os.path.abspath(dirs["src"]+"/../")
dirs["data"] = os.path.abspath(dirs["root"]+"/data/")

# Safely remove a file
def rmsafe(file:str):
    if file in ["","."]:
        print("WARNING: an attempt was made to remove the current working directory!")
        return
    if os.path.exists(file):
        os.remove(file)
