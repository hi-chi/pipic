# Configuration file for Sphinx documentation builder.
# For full options, see the Sphinx documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import subprocess

# Add parent directory to path for autodoc
sys.path.insert(0, os.path.abspath('..'))

# -- Project information --
project = 'π-PIC'
copyright = '2026, LWFA Team'
author = 'LWFA Team'

# Full version
release = '1.0'
version = '1.0.0'

# -- General configuration --
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.duration',
    'sphinx.ext.autosummary',
]

# Source file suffix
source_suffix = '.rts'

# Master doc (in docs/ root)
master_doc = 'index'

# List of patterns to exclude
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Language
language = 'en'

# -- Options for HTML output --
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for LaTeX output --
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '11pt',
    'preamble': r'''
\usepackage{amsmath}
\usepackage{amssymb}
''',
}

latex_documents = [
    ('index', 'pipic.tex', 'π-PIC Documentation',
     'LWFA Team', 'manual'),
]

# -- Options for reStructuredText --
rst_prolog = r"""
.. |pi-pic| replace:: π-PIC
.. |cgs| replace:: CGS
"""

# -- Math configuration (MathJax 3) --
mathjax3_config = {
    'tex': {
        'inlinemath': [('$', '$')],
        'displaymath': [('$$', '$$')],
    },
}

# -- Autodoc configuration --
autodoc_typehints = 'description'
autodoc_member_order = 'bysource'

# -- Intersphinx mapping --
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}

# -- Options for EPUB output
epub_show_urls = 'footnote'