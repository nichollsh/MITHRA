import os 

h_pl = 6.62607015e-34   # planck's constant
c_vac = 299792458.0     # speed of light
k_B = 1.380649e-23      # boltzmann's constant

# Safely remove a file
def rmsafe(file:str):
    if file in ["","."]:
        print("WARNING: an attempt was made to remove the current working directory!")
        return
    if os.path.exists(file):
        os.remove(file)
    return 

# Default path
dirs = {
    "pkg":os.path.abspath(os.path.dirname(os.path.abspath(__file__))+"/")
}
dirs["data"] = os.path.abspath(dirs["pkg"]+"/data/")
dirs["output"] = os.path.abspath(dirs["pkg"]+"/../../output/")

# Environment path
dirs["data"] = os.path.abspath(os.environ.get("MITHRA_DATA",dirs["data"]))

# Create data path
if not os.path.exists(dirs["data"]):
    os.mkdir(dirs["data"])
