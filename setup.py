# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import os
from setuptools import setup
import versioneer

current_directory = os.path.dirname(__file__)
readme_filename = "README.md"
readme_path = os.path.join(current_directory, readme_filename)

readme = ""
try:
    with open(readme_path, "r") as f:
        readme = f.read()
except IOError as e:
    print(e)
    print("Failed to open %s" % readme_path)

try:
    import pypandoc
    readme = pypandoc.convert(readme, to="rst", format="md")
except ImportError as e:
    print(e)
    print("Failed to convert %s to reStructuredText", readme_filename)
    pass


if __name__ == "__main__":
    setup(
        name="cohorts",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        description="Utilities for analyzing mutations and neoepitopes in patient cohorts",
        author="Tavi Nathanson",
        author_email="tavi {dot} nathanson {at} gmail {dot} com",
        url="https://github.com/tavinathanson/cohorts",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Environment :: Console",
            "Operating System :: OS Independent",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: Apache Software License",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
        install_requires=[
            "pandas>=0.15",
            "seaborn>=0.7.0",
            "scipy>=0.17.0",
            "topiary>=0.0.15",
            "six>=1.10.0",
            "lifelines>=0.9.1.0",
            "isovar>=0.0.2",
        ],
        dependency_links=[
            "git+git://github.com/hammerlab/isovar",
        ],
        long_description=readme,
        packages=["cohorts"],
    )
