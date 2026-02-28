# Configuration file for the Sphinx documentation builder.

project = "xgalois"
copyright = "2026, amirmulla"
author = "amirmulla"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "breathe",
    "exhale",
]

# -- Breathe (Doxygen bridge) ------------------------------------------------
breathe_projects = {"xgalois": "./doxyoutput/xml"}
breathe_default_project = "xgalois"

# -- Exhale (auto-generates RST from Doxygen XML) ----------------------------
exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "doxygenStripFromPath": "..",
    "rootFileTitle": "C++ API Reference",
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "INPUT = ../xgalois\nRECURSIVE = YES\n",
}

primary_domain = "cpp"
highlight_language = "cpp"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "doxyoutput", ".venv"]

# -- HTML Theme ---------------------------------------------------------------
html_theme = "furo"

html_title = "xgalois Documentation"
html_short_title = "xgalois"
html_logo = "_static/logo.png"
html_favicon = "_static/logo.png"

html_theme_options = {
    "navigation_with_keys": True,
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#1F425F",
        "color-brand-content": "#0d4b8e",
    },
    "dark_css_variables": {
        "color-brand-primary": "#7eb8da",
        "color-brand-content": "#7eb8da",
    },
}

html_static_path = ["_static"]
html_css_files = [
    "css/custom.css",
]
