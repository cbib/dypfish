#
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski
#

wget -P data/ https://zenodo.org/record/5155127/files/original.zip https://zenodo.org/record/5155127/files/chx.zip https://zenodo.org/record/5155127/files/cytod.zip https://zenodo.org/record/5155127/files/muscle.zip https://zenodo.org/record/5155127/files/nocodazole.zip https://zenodo.org/record/5155127/files/prrc2c.zip

unzip data/original.zip -d data/savulescu/
unzip data/chx.zip -d data/savulescu/
unzip data/cytod.zip -d data/savulescu/
unzip data/muscle.zip -d data/savulescu/
unzip data/nocodazole.zip -d data/savulescu/
unzip data/prrc2c.zip -d data/savulescu/

rm data/*.zip
