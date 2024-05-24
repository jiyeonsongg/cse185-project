from setuptools import setup, find_packages
#from peakfind.__version__ import __version__
# peakfind version to import
setup(
  name = 'peak_a_view',
  version=0.1,#__version__,
  description = 'CSE185 DEMO PROJECT',
  author = 'Jiyeon, Jacob, Ivana',
  author_email = 'jis036@ucsd.edu, dketchum@ucsd.edu',
  packages = find_packages(),
  install_requires=[ # only if can be installed w pip!!
      'macs2',
    ],
  scripts=['scripts/script'],
  entry_points = {
    "console_scripts": [
      "peak_a_view = peak_a_view.peak_a_view:main"
    ],
  },
)
