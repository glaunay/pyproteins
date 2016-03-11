from distutils.core import setup
setup(
  name = 'pySeq',
  packages = ['pySeq'], # this must be the same as the name above
  version = '0.1',
  description = 'Toolbox to manipulate protein sequence data',
  author = 'Guillaume Launay',
  author_email = 'pitooon@gmail.com',
  url = 'https://github.com/glaunay/pySeq', # use the URL to the github repo
  download_url = 'https://github.com/glaunay/pySeq/tarball/0.1', # I'll explain this in a second
  keywords = ['protein', 'sequence'], # arbitrary keywords
  classifiers = [],
  install_requires=[
          'bs4', 'Bio', 'urllib2', 'numpy'
      ],
)