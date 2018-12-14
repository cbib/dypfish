from setuptools import setup, find_packages

setup (
       name='DYPFISH',
       version='0.1',
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=['h5py'],

       # Fill in these to make your Egg ready for upload to
       # PyPI
       author='Hayssam Soueidan, Benjamin Dartigues and Macha Nikolski',
       author_email='',

       #summary = 'Just another Python package for the cheese shop',
       url='',
       license='',
       long_description='Long description of the package',

       # could also include long_description, download_url, classifiers, etc.
       )