.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage tests/test_functions.py -v
	sage tests/test_conditions_injectivity.py -v
	sage tests/test_conditions_bijectivity.py -v
	sage tests/test_conditions_bijectivity_robust.py -v
	# Todo: tests for utility

