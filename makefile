.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t bijectivity_exponential_maps/functions.py
	sage -t bijectivity_exponential_maps/conditions_injectivity.py
	sage tests/test_conditions_injectivity.py -v # remove this later
	sage -t bijectivity_exponential_maps/conditions_bijectivity.py
	sage tests/test_conditions_bijectivity_robust.py -v
	# TODO: tests for utility

