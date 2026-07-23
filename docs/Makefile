SPHINXOPTS ?= -W --keep-going
SPHINXBUILD ?= sphinx-build
PYTHON ?= python
SOURCEDIR = source
BUILDDIR = _build

.PHONY: help clean assets html linkcheck

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

clean:
	rm -rf "$(BUILDDIR)"

assets:
	@$(PYTHON) tools/generate_elasticity_assets.py

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

linkcheck:
	@$(SPHINXBUILD) -M linkcheck "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)
