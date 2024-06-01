from setuptools import setup, find_packages

setup(
    name='peak_a_view',
    version='0.1.0',
    description='CSE185 DEMO PROJECT',
    author='Jiyeon Song, Ivana Roque, Jacob Ketchum',
    author_email='jis036@ucsd.edu, iroque@ucsd.edu, dketchum@ucsd.edu',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'pysam',
        'macs2',
    ],
    scripts=['peak_a_view/peak_a_view_tool.py'],
    entry_points={
        'console_scripts': [
            'peak_a_view=peak_a_view.peak_a_view:main'
        ],
    },
)
