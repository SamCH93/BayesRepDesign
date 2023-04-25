## ----------------------------------------------------------------------------
## Document/Build/Check/Install R package (depends on devtools and roxygen2)
## Author: Samuel Pawel
## adapted from Manuela Ott, Sebastian Meyer, Florian Gerber
## ----------------------------------------------------------------------------

PACKAGE = BayesRepDesign
VERSION = 0.42
TAR = $(PACKAGE)_$(VERSION).tar.gz

all: build

description:
	sed -i -r -- 's/^Version:.*/Version: '$(VERSION)'/g' DESCRIPTION ;      
	sed -i -r -- 's/^Date:.*/Date: '`date +'%F'`'/g' DESCRIPTION ;
	R -e 'roxygen2::roxygenize()'

document: description
	R -e 'devtools::document()'

manual: document
	R -e 'devtools::build_manual(path = "out/")'

$(TAR): manual
	R -e 'devtools::build(path = "out/")'

build: $(TAR)

install: $(TAR)
	R -e 'devtools::install_local(path = "out/$(TAR)")'

check: $(TAR)
	R -e 'devtools::check_built(path = "out/$(TAR)", cran = FALSE)'

cran: $(TAR)
	R -e 'devtools::check_built(path = "out/$(TAR)", cran = TRUE, remote = TRUE, force_suggests = TRUE, run_dont_test = TRUE)'

test:
	R -e 'devtools::test()'

.PHONY: all document manual build install check cran description test
