import os
from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'TreeSAK', 'VERSION'))
    return version_file.readline().strip()

__long_description__ = ''' TreeSAK v%s ''' % version()

setup(name="TreeSAK",
      version=version(),
      long_description=__long_description__,
      license="GPL3+",
      author="Weizhi Song",
      author_email="songwz03@gmail.com",
      keywords="Bioinformatics",
      description="BioSAK",
      url="https://github.com/songweizhi/TreeSAK",
      packages=['TreeSAK'],
      package_data={'': ['*.r', '*.R', '*.py', '*.pl', '*.rb', '*.jar', 'VERSION', '*.hmm']},
      include_package_data=True,
      install_requires=['biopython', 'matplotlib', 'numpy', 'scipy', 'itolapi', 'networkx', 'seaborn', 'lxml', 'beautifulsoup4', 'ete3', 'arviz', 'plotly', 'kaleido', 'PyPDF3', 'dendropy'],  # reportlab
      scripts=['bin/TreeSAK'])
