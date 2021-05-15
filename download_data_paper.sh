#
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski
#

wget -P data/ http://dypfish.org/material-kit-master/file/zip/original.zip http://dypfish.org/material-kit-master/file/zip/chx.zip http://dypfish.org/material-kit-master/file/zip/cytod.zip http://dypfish.org/material-kit-master/file/zip/nocodazole.zip http://dypfish.org/material-kit-master/file/zip/prrc2c.zip http://dypfish.org/material-kit-master/file/zip/muscle.zip

unzip data/original.zip -d data/savulescu/
unzip data/chx.zip -d data/savulescu/
unzip data/cytod.zip -d data/savulescu/
unzip data/muscle.zip -d data/savulescu/
unzip data/nocodazole.zip -d data/savulescu/
unzip data/prrc2c.zip -d data/savulescu/

rm data/*.zip