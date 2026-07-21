# Sphinx configuration for the PyMC user documentation (pymcpp).
#
# Local build:
#   pip install pymcpp -r docs/requirements.txt && sphinx-build docs docs/_build
# pymcpp must be importable for the API reference (Read the Docs builds it
# from source; see .readthedocs.yaml).

import os

# Make autodoc import the compiled pymcpp module (real docstrings) instead of
# executing the installed .pyi stub, which is not runnable Python.
os.environ["SPHINX_AUTODOC_IGNORE_NATIVE_MODULE_TYPE_STUBS"] = "1"

project = "PyMC"
author = "OMEGA Research Group, Imperial College London"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "myst_nb",
    "sphinx_copybutton",
]

# The tutorial notebooks live in notebook/ at the repository root; they are
# reached through the docs/notebooks symlink (Sphinx requires all sources
# under its root directory).

# Render the outputs stored in the notebooks; never execute them (some, like
# ffmlp, need a torch-enabled build).
nb_execution_mode = "off"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
}

# The notebooks symlink exposes everything in notebook/, including local
# virtualenvs; keep Sphinx away from anything that is not a tutorial.
exclude_patterns = ["notebooks/.venv", "notebooks/**/.*"]

html_theme = "pydata_sphinx_theme"
html_title = "PyMC"
