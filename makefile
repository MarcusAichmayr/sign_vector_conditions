.PHONY: install

install:
	sage -pip install --upgrade .

test:
	sage -t sign_crn/ examples/

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
