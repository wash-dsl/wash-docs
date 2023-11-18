from cgitb import html
import subprocess, os

extensions = ["breathe", 'sphinx.ext.todo', 'exhale']

html_theme = "sphinx_rtd_theme"


breathe_projects = {
    "my_project": "../../Doxygen/gen_docs/xml"
}
# Breathe configuration
breathe_default_project = "my_project"
# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
}
# TODO: CHANGE TO /src

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'
