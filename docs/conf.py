# Sphinx configuration for the PyMC++ user documentation (pymcpp).
#
# Local build:
#   pip install pymcpp -r docs/requirements.txt && sphinx-build docs docs/_build
# pymcpp must be importable for the API reference (Read the Docs builds it
# from source; see .readthedocs.yaml).

import os

# Make autodoc import the compiled pymcpp module (real docstrings) instead of
# executing the installed .pyi stub, which is not runnable Python.
os.environ["SPHINX_AUTODOC_IGNORE_NATIVE_MODULE_TYPE_STUBS"] = "1"

project = "PyMC++"
author = "OMEGA Research Group, Imperial College London"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
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

# Render $...$ / $$...$$ LaTeX in the notebooks' markdown cells, and bare
# \begin{align} environments.
myst_enable_extensions = ["dollarmath", "amsmath"]

# MathJax: load the mathtools extension for macros like \coloneqq used in the
# notebooks.
mathjax3_config = {
    "loader": {"load": ["[tex]/mathtools"]},
    "tex": {"packages": {"[+]": ["mathtools"]}},
}

autosummary_generate = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
}

# The notebooks symlink exposes everything in notebook/, including local
# virtualenvs; keep Sphinx away from anything that is not a tutorial.
exclude_patterns = ["notebooks/.venv", "notebooks/**/.*"]

html_theme = "pydata_sphinx_theme"
html_title = "PyMC++"

# Logo of MC++ (doc/ at the repository root): the theme shows the light or the
# dark variant with its colour mode; the icon follows the browser's colour
# scheme by itself.
html_theme_options = {
    "logo": {
        "image_light": "../doc/mcpp-logo.svg",
        "image_dark": "../doc/mcpp-logo-dark.svg",
        "text": "PyMC++",
        "alt_text": "PyMC++ documentation - Home",
    },
}
html_favicon = "../doc/mcpp-logo-55.svg"
