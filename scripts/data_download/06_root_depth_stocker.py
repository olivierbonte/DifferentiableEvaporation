# %% imports
import shutil
import subprocess

from conf import stocker_dir

stocker_dir.mkdir(exist_ok=True)

# %% Download
zenodo_doi = "10.5281/zenodo.10885724"
file_names = [
    "cwdx80.nc",
    "cwdx80_forcing.nc",
    "zroot_cwd80.nc",
    "zroot_cwdx80_forcing.nc",
]
subprocess.run(["zenodo_get", zenodo_doi])
for file in file_names:
    shutil.move(file, stocker_dir / file)

# %%
