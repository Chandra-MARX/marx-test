marxtest.database module
========================

In addition to pretty pictures, many |marx| tests summary the test result
in a single number or string. After each test run, these numbers are added
to a permanent database of test results for later inspection.
This database can be used to e.g. compare test results between different |marx|
versions, different `CIAO`_ versions or for different computers to compare
performance.

The database is stored in the `SQLite <http://www.sqlite.org>` format through
the `sqlite3` package in the Python standard library.
The path to the database file is set in the cfg file and the methods in this
module provide the interface to that database.

The following tables summarize the database scheme. There are three
tables in the database with the following fields:

Table ``config``
----------------
This table holds information about the setup of the tests, like the
marx version number. The fields are:

:marx: |marx| version
:ciao: `CIAO`_ version
:caldb: CalDB version
:saotrace: `SAOTrace`_ version
:compiler: Set this field to compare results from different compilers.
    Since this cannot be determined automatically from the |marx| binaries,
    it has to be set correctly be hand.
:compilerversion: As before, this has to be set by hand.
:compilerflags: As before, this has to be set by hand.
:host: For comparing performance on different computers

Table ``expresult``
-------------------
This table contains the expect results of a test. For example, when fitting
a spectrum, the expected result would be a certain spectral slope:

:testclass: name of the test
:name: A single test might output more than one number. So, in addition to the
    name of the test class, a "name" for the result can be set.
:version: If the definition of a test changes so that results cannot be
    compared to previous runs any longer, the version needs to be incremented.
:title: A title used when results are displayed on websites etc.
    Limited to 20 characters.
:description: Describe the test if more details are needed.
:value: expected value
:unit: Unit of value. Used for display only.
:acceptable: Values are rarely exact. This field defined the maximum deviation
   from ``value`` that still makes a test pass.

Table ``result``
----------------
This table stores test results:

:configid: ``ROWID`` in table ``config`` that describes the configuration of
    the test system.
:expid: ``ROWID`` in table ``expresults`` matching the test that was run
:value: results of the test
:sigma: uncertainty on ``value`` e.g. for fit results. Can be ``NONE``.
:sigma2: In some cases a single value for ``sigma`` might not be sufficient,
    e.g. if error up and down a different. In this case, ``sigma2`` hold the
    second number.


Accessing the database
----------------------
The definition of the database is kept simple on purpose to allow for easy
treatment of special cases. As one example, it is possible to write test
results that do not contain a ``value`` with `~marxtest.database.TBD`. In the
current implementation that is not a meaningful entry and one could define the
database with an ``NOT NULL`` constraint to prevent empty value, but we could
imagine a future use where we want to record that a test was run, even if it
produced no output because the |marx| version in question does not support the
commands on anything similar. Similarly, it is possible to manually add rows to
the ``result`` that have a ``configid`` value with no counterpart in the
``config`` table. (A ``FOREIGN KEY`` constraint would prevent that, but this
was added only recently to `SQLite <http://www.sqlite.org>`.)

On the downside, this means that caution is required when writing to the
database outside of the functions defined in this module.


.. warning:: Certain settings are cached!

   The database is initialized and the version information for |marx|,
   `CIAO`_, etc. is cached when they are first used.
   To reset this cache, exit Python and start again or call
   `~marxtest.database.reset_cache()` explicitly.



Reference/API
-------------

.. automodule:: marxtest.database
    :members:
    :undoc-members:
    :show-inheritance:
