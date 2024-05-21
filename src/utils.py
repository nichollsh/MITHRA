import os 

h_pl = 6.62607015e-34 # planck's constant
c_vac = 299792458.0
k_B = 1.380649e-23


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
