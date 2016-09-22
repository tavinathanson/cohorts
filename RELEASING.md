# Releasing Cohorts

This document explains what do once your [Pull Request](https://www.atlassian.com/git/tutorials/making-a-pull-request/) has been reviewed and all final changes applied. Now you're ready merge your branch into master and release it to the world:

1. Assign a version to the release you are preparing. `cohorts` uses [versioneer](https://github.com/warner/python-versioneer), so rather
than explicitly editing a `__version__` variable in the code you must instead tag the branch with your new version number (e.g. `git tag 1.2.3`).

Warning: whether you should do `git tag v1.2.3` vs. `git tag 1.2.3` depends on versioneer settings, which can be found in `setup.cfg` and `_version.py` (because `_version.py` was generated from `setup.cfg` configuration, the two should match.) In `cohorts`, `git tag 1.2.3` (without the `v`) is what's expected.

2. Once your candidate release branch is tagged you must then run `git push --tags` and merge the branch.

3. After the `cohorts` unit tests complete successfully on Travis then the latest version
of the code (with the version specified above) will be pushed to [PyPI](https://pypi.python.org/pypi) automatically. If you're curious about how automatic deployment is achieved, see our [Travis configuration](https://github.com/hammerlab/cohorts/blob/master/.travis.yml#L51).

4. Finally, run [github-changelog-generator](https://github.com/skywinder/github-changelog-generator) from the `cohorts` directory (first, install it if needed).
