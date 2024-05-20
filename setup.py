setup(
  name = 'mypileup', # change as ours
  version = VERSION,
  description = 'CSE185 DEMO PROJECT',
  author = 'Jiyeon', 'Jacob', 'Ivana',
  author_email = 'jis036@ucsd.edu',
  packages = find_packages(),
  entry_points = {
    "console_scripts": [
      "mypileup = mypileup.mypileup:main" # change as ours
    ]
  }
)
