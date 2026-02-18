PKGNAME := WGCNAplus

.PHONY: doc build install test check clean

doc:
	Rscript -e "roxygen2::roxygenise('.')"

build: doc
	cd .. && R CMD build $(PKGNAME)

install: doc
	R CMD INSTALL .

test:
ifdef filter
	Rscript -e "testthat::test_file(list.files('tests/testthat', pattern='$(filter)', full.names=TRUE))"
else
	Rscript -e "testthat::test_local('.')"
endif

check: doc
	cd .. && R CMD check $(PKGNAME) --no-manual

clean:
	rm -rf man/*.Rd NAMESPACE
	rm -rf ../$(PKGNAME)_*.tar.gz
	rm -rf ../$(PKGNAME).Rcheck
