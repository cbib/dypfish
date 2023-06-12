#
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski
#

wget -P data/ https://zenodo.org/record/5155127/files/original.zip?download=1 https://zenodo.org/record/5155127/files/chx.zip?download=1 https://zenodo.org/record/5155127/files/cytod.zip?download=1 https://zenodo.org/record/5155127/files/muscle.zip?download=1 https://zenodo.org/record/5155127/files/nocodazole.zip?download=1 https://zenodo.org/record/5155127/files/prrc2c.zip?download=1

unzip data/original.zip -d data/savulescu/
unzip data/chx.zip -d data/savulescu/
unzip data/cytod.zip -d data/savulescu/
unzip data/muscle.zip -d data/savulescu/
unzip data/nocodazole.zip -d data/savulescu/
unzip data/prrc2c.zip -d data/savulescu/

rm data/*.zip
