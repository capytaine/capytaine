# How to contribute

Thank you for considering to contribute to Capytaine.
Here are some things you should know before you go further.

## Contact

The best place to report a bug or to propose an enhancement is the
[issue tracker](https://github.com/mancellin/capytaine/issues)
of Capytaine's Github repository.

The [discussion page](https://github.com/mancellin/capytaine/discussions)
can also be used to ask questions and submit ideas.

## Submitting changes

Pull requests are welcome. Before submitting one, please take a look at the
following points:

* The code follows the [PEP8](https://pep8.org/) guidelines for code style.
* Indentation is done with four spaces.
* The Github repository includes tests in the `pytest` directory.
  Test can be run with the `pytest` module (`python -m pytest`).
  Before submitting a change, make sure that all tests are passing.
  If you'd like to add a feature, please add the relevant tests to 
  the `pytest` directory.
* All modifications of the code are listed in the `docs/changelog.rst`.
  If you submit a bug fix or a new feature, please add a line in
  the changelog.
* New features should also be documented in the user manual (`docs/user_manual/`).
  Please add a short documentation of the feature in the relevant page.

