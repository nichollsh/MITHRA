import matplotlib.pyplot as plt 
import numpy as np

def plot_spectrum(wl,fl,writepath="",title=""):
    plt.close("all")
    fig,ax = plt.subplots()

    ax.plot(wl,fl, lw=0.8, color='black')

    ax.set(xscale="log", yscale="log")
    ax.set(xlabel="Wavelength [nm]", ylabel="Spectral flux [erg s-1 cm-2 nm-1]")
    if len(title) > 0:
        ax.set_title(title)

    if len(writepath) > 0:
        fig.savefig(writepath, dpi=200, bbox_inches='tight')
    else:
        plt.show()
    
