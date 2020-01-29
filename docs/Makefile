# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXAPIDOC  = sphinx-apidoc
PYTHONAPISRC  = ../pdb_prot_align
SPHINXBUILD   = python3 -msphinx
SPHINXPROJ    = pdb_prot_align
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	$(SPHINXAPIDOC) -E -T -M -e -f -o . $(PYTHONAPISRC)
	cp -r ../notebooks/input_files .  # copy files for notebook examples
	cp -r ../notebooks/output_files .  # copy files for notebook examples
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
