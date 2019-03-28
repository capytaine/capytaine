# How to contribute

Thank you for considering to contribute to Capytaine.
Here are some things you should know before you go further.

## Contact

The best place to report a bug or to propose an enhancement is the
[issue tracker](https://github.com/mancellin/capytaine/issues)
of Capytaine's Github repository.

For questions or comments that are not directly related to the present Python
implementation (theory, interpretation of the results, ...), you can also leave
a message on [Nemoh's forum](http://lheea.ec-nantes.fr/redmine/projects/nemoh)
([Registration](http://lheea.ec-nantes.fr/redmine/account/register))

## Submitting changes

Pull requests are welcome. Before submitting one, please take a look at the
following points:

* The `master` branch in the git repository is the latest stable version of the code.
  Developments are done on the `develop` branch. Please prefer this branch for
  pull requests.
* The code follows the [PEP8](https://pep8.org/) guidelines for code style.
* Indentation is done with four spaces.
* The Github repository includes tests in the `pytest` directory.
  Test can be run with the `pytest` module (`python -m pytest`).
  Before submitting a change, make sure that all tests are passing.
  If you'd like to add a feature, please add the relevant tests to 
  the `pytest` directory.

