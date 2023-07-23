# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0,os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../pyPPG"))
sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))

sys.path.insert(0, os.path.abspath("../.."))

# sys.path.insert(0, os.path.abspath("../.."))
# sys.path.insert(0, os.path.abspath("./"))


project = 'pyPPG'
copyright = '2023, Marton A. GODA, PhD'
author = 'Marton A. GODA, PhD'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.todo","sphinx.ext.viewcode","sphinx.ext.autodoc"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_logo = "logo.png"
html_theme_options = {"logo_only": False}

exclude_patterns = ['test_annot_04', 'test_annot_BIDMC']