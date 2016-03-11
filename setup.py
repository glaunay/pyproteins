from setuptools import setup, find_packages
setup(
  name = 'pyproteins',
  packages = find_packages(), # this must be the same as the name above
  version = '0.1',
  description = 'Toolbox to manipulate protein sequence data',
  author = 'Guillaume Launay',
  author_email = 'pitooon@gmail.com',
  url = 'https://github.com/glaunay/pyproteins', # use the URL to the github repo
  download_url = 'https://github.com/glaunay/pyproteins/tarball/0.1', # I'll explain this in a second
  keywords = ['protein', 'sequence'], # arbitrary keywords
  classifiers = [],
  install_requires=[
          'bs4', 'Bio', 'urllib2', 'numpy'
      ],
    dependency_links = [
        "http://dev.danse.us/trac/pathos"
    ]
)