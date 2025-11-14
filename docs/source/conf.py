# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import os
import sys
sys.path.insert(0, os.path.abspath('../../src/'))


project = 'carmapy'
copyright = '2025, Wolf Cukier'
author = 'Wolf Cukier'
release = '0.5.13'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx_autodoc_typehints',
              'nbsphinx']

templates_path = ['_templates']
exclude_patterns = []


# autodoc_default_options = {
#     'members': True,
#     'undoc-members': True,
#     'show-inheritance': True,
#     'inherited-members': True,
# }


nbsphinx_allow_errors = False
autodoc_typehints = "description"
autosummary_generate = True
nbsphinx_execute = 'always'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_logo = "logo.png"
html_theme_options = {
    'logo_only': True,
}

nbsphinx_prolog = """
{% set docname = env.doc2path(env.docname, base=None) %}
.. note::  `Download full notebook here <https://github.com/wcukier/carmapy/tree/main/docs/source/{{ docname }}>`_
"""
