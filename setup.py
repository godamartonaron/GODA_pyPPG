from setuptools import setup, find_packages
import io
import os

with io.open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()


with open('docs'+os.sep+'requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='pyPPG',
    version='1.0.36',
    description='pyPPG: a python toolbox for PPG morphological analysis.',
    author='Marton A. Goda, PhD; Peter H. Charlton, PhD',
    author_email="marton.goda@campus.technion.ac.il",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/godamartonaron/GODA_pyPPG",
    project_urls={"Bug Tracker": "https://github.com/pypa/sampleproject/issues",},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
    ],
    packages={"pyPPG", "pyPPG"+os.sep+"ppg_bm", "pyPPG"+os.sep+"pack_ppg"},
    package_data={
        "pyPPG": ["*"],
        "pyPPG"+os.sep+"ppg_bm": ["*"],
        "pyPPG"+os.sep+"pack_ppg": ["*"],
    },

    install_requires=[required],
    python_requires=">=3.10",
    include_package_data=True,
)