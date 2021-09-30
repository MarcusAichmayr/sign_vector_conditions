.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t bijectivity_exponential_maps/functions.py
	sage tests/test_conditions_injectivity.py -v
	sage -t bijectivity_exponential_maps/conditions_bijectivity.py
	sage tests/test_conditions_bijectivity_robust.py -v
	# Todo: tests for utility

