from setuptools import setup, find_packages
setup(
  name = 'pyproteins',
  packages = find_packages(), # this must be the same as the name above
  version = '0.2',
  license='BSD',
  description = 'Toolbox to manipulate protein sequence data',
  author = 'Guillaume Launay',
  author_email = 'pitooon@gmail.com',
  url = 'https://github.com/glaunay/pyproteins', # use the URL to the github repo
  packages=find_packages('src'),
  package_dir={'': 'src'},
  include_package_data=True,
  zip_safe=False,
  py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
  download_url = 'https://github.com/glaunay/pyproteins/tarball/0.2', # I'll explain this in a second
  keywords = ['protein', 'sequence'], # arbitrary keywords
  classifiers = [],
  install_requires=[
          'bs4', 'biopython', 'numpy', 'paramiko'
      ],
    dependency_links = [
        "http://dev.danse.us/trac/pathos"
    ]
)