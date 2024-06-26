# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import sys

#-----------------------------------------------------------------------------

#Adding this to generate API reference on readthedocs website
def run_apidoc(_):
	from sphinx.ext.apidoc import main
	sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
	cur_dir = os.path.abspath(os.path.dirname(__file__))
	module = os.path.join(cur_dir,"..","meps")
	main(['-e', '-o', cur_dir, module, '--force'])

def setup(app):
	app.connect('builder-inited', run_apidoc)

# ----------------------------------------------------------------------------

sys.path.insert(0,os.path.abspath('../src/epspy/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'epspy'
copyright = '2024, Shikhar Mittal'
author = 'Shikhar Mittal'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc','sphinx.ext.mathjax']

source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


