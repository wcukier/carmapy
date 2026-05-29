# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import os
import sys
from importlib.metadata import version as _get_version
sys.path.insert(0, os.path.abspath('../../src/'))

if os.environ.get("READTHEDOCS") == "True":
    _repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    _refdata = os.path.join(_repo_root, "picaso", "reference")
    os.environ.setdefault("picaso_refdata", _refdata)
    os.environ.setdefault("PYSYN_CDBS", os.path.join(_repo_root, "picaso", "reference", "stellar_grids"))

    # Patch Carma.run on RTD so tutorials can skip long runs without exposing
    # RTD-specific kwargs to readers. Convention: pass suppress_output=True in
    # a notebook cell to make that run a no-op on RTD (pre-committed output
    # files are read instead). Calls without suppress_output still execute but
    # always get suppress_output=True so Fortran output doesn't flood the build log.
    import carmapy as _carmapy
    _original_run = _carmapy.Carma.run
    def _rtd_run(self, *args, suppress_output=False, **kwargs):
        if suppress_output:
            return
        return _original_run(self, *args, suppress_output=True, **kwargs)
    _carmapy.Carma.run = _rtd_run
else:
    os.environ["PYSYN_CDBS"]     = "/Users/wcukier/Dropbox/Research/Projects/24-Brown Dwarfs/picaso/PYSYN_CDBS"
    os.environ["picaso_refdata"] = "/Users/wcukier/Dropbox/Research/Projects/24-Brown Dwarfs/picaso/reference"


project = 'carmapy'
copyright = '2025, Wolf Cukier'
author = 'Wolf Cukier'
release = _get_version('carmapy')

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


nbsphinx_allow_errors = True
autodoc_typehints = "description"
autosummary_generate = True
nbsphinx_execute = 'auto'
nbsphinx_timeout = 1800

nbsphinx_custom_formats = {
    ".py": ["jupytext.reads", {"fmt": "py:percent"}],
}

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
