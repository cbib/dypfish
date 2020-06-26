

wget -P data/ http://dypfish.org/file/zip/micropatterned.zip http://dypfish.org/file/zip/chx.zip http://dypfish.org/file/zip/cytod.zip http://dypfish.org/file/zip/nocodazole.zip http://dypfish.org/file/zip/prrc2c.zip http://dypfish.org/file/zip/muscle.zip

unzip data/chx.zip -d data/
unzip data/cytod.zip -d data/
unzip data/muscle.zip -d data/
unzip data/nocodazole.zip -d data/
unzip data/prrc2c.zip -d data/

rm data/*.zip