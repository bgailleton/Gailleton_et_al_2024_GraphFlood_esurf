========
scabbard
========


.. image:: https://img.shields.io/pypi/v/scabbard.svg
        :target: https://pypi.python.org/pypi/scabbard

.. image:: https://img.shields.io/travis/bgailleton/scabbard.svg
        :target: https://travis-ci.com/bgailleton/scabbard

.. image:: https://readthedocs.org/projects/scabbard/badge/?version=latest
        :target: https://scabbard.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




high-level python package for the DAGGER suite


* Free software: MIT license
* Documentation: https://scabbard.readthedocs.io.


Features
--------

* TODO


Experimental
------------

Among the experimental features is a tentative link to Blender to make nice 3D plots. The latter will only get activated if called from blender and will happily be ignored by most.

Only tested on Ubuntu, it requires the following steps:

1. identify the `python` version of your blender version
2. create a new env with this python version
3. install all the packages you need either with conda or pip.
4. AFTER EVERY NEW ADDITIONS (of one or multiple packages at once) YOU'LL NEED TO LINK THE PACKAGES TO BLENDER:

`ln -s /path/to/your/environment/site-packages/* ~/.config/blender/3.X/scripts/addons/modules/`

Where you obviously need to adapts the paths and version number (yes I had multiple times the question "urgh it cannot find /path/to/your/environment")

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
