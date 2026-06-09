# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from pathlib import Path
import datetime
import importlib
import sys
import tomllib

from sphinx.ext.autodoc import AttributeDocumenter
from stpipe import Step

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

with open(Path(__file__).parent.parent / "pyproject.toml", "rb") as metadata_file:
    metadata = tomllib.load(metadata_file)["project"]

project = metadata["name"]
author = "Space Telescope Science Institute (`STScI <stsci.edu>`_), Association of Universities for Research in Astronomy (`AURA <https://www.aura-astronomy.org>`_)"
copyright = f"{datetime.datetime.today().year}, {author}"

package = importlib.import_module(metadata["name"])
try:
    version = package.__version__.split("-", 1)[0]
    # The full version, including alpha/beta/rc tags.
    release = package.__version__
except AttributeError:
    version = "dev"
    release = "dev"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "numfig",
    "numpydoc",
    "sphinxcontrib.jquery",
    "pytest_doctestplus.sphinx.doctestplus",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.automodsumm",
    "sphinx_automodapi.autodoc_enhancements",
    "sphinx_automodapi.smart_resolver",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# reST default role used for single backticks (`text`)
default_role = "obj"

# we're overriding directives in sphinx.ext.autodoc_enhancements;
# TODO remove this when https://github.com/sphinx-doc/sphinx/pull/1843 is released
suppress_warnings = [
    "app.add_directive",
]

rst_epilog = """.. _jwst: high-level_API.html"""

# -- HTML output configuration ----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": True,
    "sticky_navigation": False,
    "style_external_links": True,
}
html_logo = "_static/jwst_logo.png"
html_last_updated_fmt = "%b %d, %Y"
html_sidebars = {"**": ["globaltoc.html", "relations.html", "searchbox.html"]}
html_domain_indices = True
html_use_index = True

html_css_files = ["custom.css"]
htmlhelp_basename = "jwstdoc"

# -- EPUB output configuration -----------------------------------------------

epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright
epub_show_urls = "footnote"
epub_exclude_files = ["search.html"]

# -- LaTeX output configuration ----------------------------------------------

latex_elements = {
    "papersize": "letterpaper",  #'letterpaper' or 'a4paper'
    "pointsize": "14pt",  #'10pt', '11pt' or '12pt'
    "preamble": r"""\usepackage{enumitem} \setlistdepth{99}""",
}
latex_documents = [
    (
        "index",  # source start file
        f"{project}.tex",  # target name
        "JWST Calibration Pipeline",  # title
        author,  # author
        "manual",  # documentclass [howto, manual, or own class]
    ),
]
latex_show_urls = "True"
latex_domain_indices = True

# -- Texinfo output configuration -------------------------------------------

texinfo_documents = [
    (
        "index",  # source start file
        project,  # target name
        "JWST Calibration Pipeline",  # title
        author,  # author
        project,  # dir menu entry
        "James Webb Space Telescope (JWST) Calibration Pipeline",  # description
        "Miscellaneous",  # category
    ),
]
texinfo_domain_indices = True
texinfo_show_urls = "inline"  # 'footnote', 'no', or 'inline'

# If true, do not generate a @detailmenu in the "Top" node's menu.
# texinfo_no_detailmenu = False

# -- manpage output configuration ---------------------------------------

man_pages = [
    (
        "index",  # source start file
        project,  # name
        "James Webb Space Telescope (JWST) Calibration Pipeline",  # description
        [author],  # authors
        1,  # manual section
    )
]
man_show_urls = True

# -- linkcheck configuration -------------------------------------------------

linkcheck_retry = 5
linkcheck_ignore = [
    "http://stsci.edu/schemas/fits-schema/",  # Old schema from CHANGES.rst
    "https://stsci.edu",  # CI blocked by service provider
    "https://jwst-docs.stsci.edu",  # CI blocked by service provider
    "https://outerspace.stsci.edu",  # CI blocked by service provider
    "https://jira.stsci.edu/",  # Internal access only
    r"https://.*\.readthedocs\.io",  # 429 Client Error: Too Many Requests
    "https://doi.org",  # CI blocked by service provider (timeout)
    r"https://github\.com/spacetelescope/jwst/(?:issues|pull|blob)",
]
linkcheck_timeout = 180
linkcheck_anchors = False
linkcheck_report_timeouts_as_broken = True
linkcheck_allow_unauthorized = False

# Enable nitpicky mode - which ensures that all references in the docs resolve.
nitpicky = True

# -- numpydoc configuration --------------------------------------------------

# Don't show summaries of the members in each class along with the class' docstring
numpydoc_show_class_members = False

# -- sphinx-autodoc configuration --------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#configuration

autoclass_content = "both"  # combine class and __init__ docstrings

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath("jwst/"))
sys.path.insert(0, os.path.abspath("exts/"))


class StepSpecDocumenter(AttributeDocumenter):
    def should_suppress_value_header(self):
        if self.name == "spec" and issubclass(self.parent, Step):
            # if this attribute is named "spec" and belongs to a "Step"
            # don't show the value, it will be formatted in add_context below
            return True
        return super().should_suppress_value_header()

    def add_content(self, more_content):
        super().add_content(more_content)
        if self.name != "spec" or not issubclass(self.parent, Step):
            return
        if not self.object.strip():
            return

        # format the long "Step.spec" string to improve readability
        source_name = self.get_sourcename()
        self.add_line("::", source_name, 0)
        self.add_line("  ", source_name, 1)
        txt = "\n".join((l.strip() for l in self.object.strip().splitlines()))
        self.add_line(f"  {txt}", source_name, 2)


def setup(app):
    # add a custom AttributeDocumenter subclass to handle Step.spec formatting
    def register_documenter(app, config):
        app.add_autodocumenter(StepSpecDocumenter, True)

    # register it with a high priority so it behaves with the built-in autodoc
    app.connect("config-inited", register_documenter, priority=9000)


# -- sphinx-automodapi configuration --------------------------------------------
# https://sphinx-automodapi.readthedocs.io/en/latest/automodapi.html#automatically-generating-module-documentation-with-automodapi

automodapi_toctreedirnm = "api"

# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]

# If true, '()' will be appended to :func: etc. cross-reference text.
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# syntax highlighting
pygments_style = "default"

# -- sphinx.ext.autosummary configuration ----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html

autosummary_generate = True

# -- sphinx.ext.intersphinx configuration ------------------------------------

# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration

intersphinx_mapping = {
    "asdf": ("https://asdf.readthedocs.io/en/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "drizzle": ("https://spacetelescope-drizzle.readthedocs.io/en/latest/", None),
    "gwcs": ("https://gwcs.readthedocs.io/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "numpy": ("https://numpy.org/devdocs", None),
    "photutils": ("https://photutils.readthedocs.io/en/stable/", None),
    "python": ("https://docs.python.org/3/", None),
    "requests": ("https://requests.readthedocs.io/en/latest/", None),
    "scipy": ("https://scipy.github.io/devdocs", None),
    "stcal": ("https://stcal.readthedocs.io/en/latest/", None),
    "stdatamodels": ("https://stdatamodels.readthedocs.io/en/latest/", None),
    "stpipe": ("https://stpipe.readthedocs.io/en/latest/", None),
    "synphot": ("https://synphot.readthedocs.io/en/latest/", None),
    "tweakwcs": ("https://tweakwcs.readthedocs.io/en/latest/", None),
}
intersphinx_disabled_domains = ["std"]
