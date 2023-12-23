# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import pathlib
import sys
import os
#sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
sys.path.insert(0, "../src")

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Jelly'
copyright = '2023, Pietro Pustina'
author = 'Pietro Pustina'
release = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinxcontrib.matlab',
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.imgmath',
    'sphinxcontrib.bibtex',
    'sphinx.ext.githubpages',
]

templates_path = ['_templates']
exclude_patterns = []


# Configure the matlab extension
primary_domain = "mat"
this_dir = os.path.dirname(__file__)
matlab_src_dir = os.path.abspath(os.path.join(this_dir, "..", "..", "src"))
# matlab_src_dir = "../src"
matlab_short_links = True

# Configure the bibtex extension
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = 'unsrt'

#Configure latex
imgmath_latex_preamble="\\usepackage{xcolor} \\definecolor{myRed}{rgb}{0.7,0.125,0.309} \\definecolor{myBlue}{rgb}{0.396,0.560,0.917} \\color{myRed}"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_logo = 'logo.svg'
html_favicon = 'favicon.svg'
html_context = {
   "default_mode": "auto"
}
