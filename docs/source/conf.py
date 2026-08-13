# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import os
import sys
import glob
import gzip
import shutil
from importlib.metadata import version as _get_version
sys.path.insert(0, os.path.abspath('../../src/'))


def _decompress_notebook_data():
    """Decompress ``*.gz`` data files shipped alongside the notebooks.

    Large CARMA result files (e.g. ``my_first_carma.txt``,
    ``rates_my_first_carma.txt``) are committed in gzip-compressed form to keep
    them well under GitHub's file-size limits without needing git LFS. The
    notebooks read the raw ``.txt`` files, so decompress any ``*.gz`` whose
    decompressed counterpart is missing or stale before Sphinx executes them.
    """
    notebooks_dir = os.path.join(os.path.dirname(__file__), "notebooks")
    for gz_path in glob.glob(os.path.join(notebooks_dir, "**", "*.gz"),
                             recursive=True):
        out_path = gz_path[:-len(".gz")]
        if (os.path.exists(out_path)
                and os.path.getmtime(out_path) >= os.path.getmtime(gz_path)):
            continue
        with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)


_decompress_notebook_data()

if os.environ.get("READTHEDOCS") == "True":
    _repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    _refdata = os.path.join(_repo_root, "picaso", "reference")
    os.environ.setdefault("picaso_refdata", _refdata)

    # The STScI synphot tarballs extract under a ``grp/redcat/trds/`` prefix, so
    # the stellar grids do not sit directly under ``stellar_grids/``. Locate the
    # phoenix catalog and point PYSYN_CDBS at the dir that contains ``grid/`` so
    # pysynphot can find ``grid/phoenix/catalog.fits``.
    _stellar_grids = os.path.join(_refdata, "stellar_grids")
    _catalogs = glob.glob(
        os.path.join(_stellar_grids, "**", "grid", "phoenix", "catalog.fits"),
        recursive=True,
    )
    if _catalogs:
        # .../grid/phoenix/catalog.fits -> .../  (parent of grid/)
        _pysyn_cdbs = os.path.dirname(
            os.path.dirname(os.path.dirname(_catalogs[0]))
        )
    else:
        _pysyn_cdbs = _stellar_grids
    os.environ.setdefault("PYSYN_CDBS", _pysyn_cdbs)

else:
    os.environ["PYSYN_CDBS"]     = "/Users/wcukier/Library/CloudStorage/Dropbox/Research/Projects/2026/2026-Brown Dwarfs/picaso/PYSYN_CDBS"
    os.environ["picaso_refdata"] = "/Users/wcukier/Library/CloudStorage/Dropbox/Research/Projects/2026/2026-Brown Dwarfs/picaso/reference"


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
