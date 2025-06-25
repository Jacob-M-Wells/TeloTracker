#!/usr/bin/env python3
"""
This is the TeloTracker installation script. Assuming you're in the same directory, it can be run
like this: `python3 setup.py install`, or (probably better) like this: `pip3 install .`

Copyright 2025 Jacob M. Wells (jacobwells1203@gmail.com)
https://github.com/Jacob-M-Wells/TeloTracker

This file is part of TeloTracker. TeloTracker is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. TeloTracker is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with TeloTracker.
If not, see <http://www.gnu.org/licenses/>.
"""

from setuptools import setup, find_packages

# Get the program version form the version.py file
with open('telotracker/version.py', 'r') as f:
    lines = f.readlines()
    __version__ = "0.0.0"
    for line in lines:
        if "__version__" in line:
            # Strip whitespace, then remove both single or double quotes
            __version__ = line.split("=")[1].strip().strip('"').strip("'")
            break

# Read long description from README
with open('README.md', 'r') as f:
    readme_description = f.read()

setup(
    name="TeloTracker",
    version=__version__,
    description="TeloTracker: A tool for telomere analysis in yeast",
    long_description=readme_description,
    long_description_content_type='text/markdown',
    url="https://github.com/Jacob-M-Wells/TeloTracker",
    author="Jacob M. Wells",
    author_email="jacobwells1203@gmail.com",
    license='GPLv3',
    packages=find_packages(),
    python_requires="3.11",
    install_requires=[
        "python-dateutil>=2.9.0",
        "python-edlib>=1.3.9",
        "numpy>=1.26.4",
        "pandas>=2.2.3",
        "pandas-flavor>=0.6.0",
        "scipy>=1.14.1",
        "matplotlib>=3.9.1",
        "seaborn>=0.13.2",
        "pingouin>=0.5.5",
        "biopython>=1.84",
        "pysam>=0.22.1",
        "parasail-python>=1.3.4",
    ],
    # External binaries not suitable for pip installation can be listed here
    # These will be installed via conda but not via pip
    extras_require={
        "bio_tools": [
            # This section can document the conda-only dependencies 
            # but doesn't actually install them via pip
        ]
    },
    entry_points={
        "console_scripts": [
            "telotracker=telotracker.__main__:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)