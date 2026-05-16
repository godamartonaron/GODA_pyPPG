from setuptools import setup, find_packages
import io

with io.open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()


required = [
    "numpy>=1.26,<3",
    "scipy>=1.10",
    "pandas>=2.0",
    "dotmap>=1.3.30",
    "matplotlib>=3.6",
    "mne>=1.5",
]

setup(
    name='pyPPG',
    version='1.0.73.post1',
    description='pyPPG: a python toolbox for PPG morphological analysis.',
    author='Marton A. Goda, PhD; Peter H. Charlton, PhD',
    author_email="marton.goda@campus.technion.ac.il",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/godamartonaron/GODA_pyPPG",
    project_urls={
        "Documentation": "https://pyppg.readthedocs.io/en/latest/",
        "Source": "https://github.com/godamartonaron/GODA_pyPPG",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    package_data={
        "pyPPG": ["*"],
        "pyPPG.ppg_bm": ["*"],
        "pyPPG.pack_ppg": ["*"],
        "pyPPG.validation": ["pw_anal.py", "__init__.py"],
    },

    install_requires=required,
    python_requires=">=3.10",
    include_package_data=True,
)
