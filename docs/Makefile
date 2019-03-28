# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -m sphinx
SPHINXPROJ    = capytaine
SOURCEDIR     = .
BUILDDIR      = _build

all: html

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile all show test autodoc

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

show: html
	@firefox _build/html/index.html &

test: show

autodoc:
	rm -rf developer_manual/api/capytaine*.rst
	sphinx-apidoc -o developer_manual/api/ ../capytaine --separate --doc-project "API documentation"


REMOTE_PATH=/var/www/html/capytaine
REMOTE_NAME=vps
VERSION=$(shell python -c 'import capytaine; print(capytaine.__version__)')
upload: html
	ssh $(REMOTE_NAME) "cd $(REMOTE_PATH) && mkdir $(VERSION) && rm latest && ln -s $(VERSION) latest"
	rsync -a _build/html/* $(REMOTE_NAME):$(REMOTE_PATH)/$(VERSION)/

