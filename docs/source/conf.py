# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import pathlib
import sys
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Jelly'
copyright = '2023, Pietro Pustina'
author = 'Pietro Pustina'
release = '1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.imgmath',
    'sphinxcontrib.matlab',
    'sphinxcontrib.bibtex',
    'sphinx.ext.githubpages',
]

templates_path = ['_templates']
exclude_patterns = []


# Configure the matlab extension
primary_domain = "mat"
matlab_src_dir = "/media/Dati/Desktop/Sapienza/Dottorato/KaneToolboxV2/"
matlab_short_links = True

# Configure the bibtex extension
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = 'unsrt'

# Configure latex extension
# imgmath_embed = True
# imgmath_latex_preamble = '\\usepackage{siunitx}'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
