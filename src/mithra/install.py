# Script to be run when `poetry install` is called

from mithra.tracks import download_bhac

def main():
    # Download tracks
    download_bhac(use_cache=False)
