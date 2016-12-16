---------------------------
|marx| accuracy and testing
---------------------------

No simulation is ever perfect. In this section, we show tests for the internal
consistency of |marx| and compare simulations to real observations.
This demonstrates the limits to which |marx| simulations can be
trusted. Anybody using |marx| for scientific analysis needs to consult the
following pages to determine if the simulations are realistic enough for
their purpose.

For example, the |marx| HRMA module uses a simplified geometry to speed up
simulations so that the user can simulate millions of rays in just
a few minutes on a typical desktop computer. Thus, the simulated PSF will
certainly differ from the observed PSF. `SAOTrace`_ has more details, but does 
not the
`"hook" feature of the PSF
<http://cxc.harvard.edu/ciao/caveats/psf_artifact.html>`_ seen on very small
scales either. Only a comparison between simulations and observed data can tell
how well the simulations are doing.


.. toctree::
   :maxdepth: 2

   {% for module in modulelist %}{{ module }}
   {% endfor %}

The code to run all these tests is available and linked below. Advanced users
may wish to inspect the code of the tests for some more ideas on how to use
|marx|. The code is lightly commented, but not in as much details as :ref:`examples`.
It will not execute as-is on your computer because it depends on the
local ``$PATH`` and other environment variables. For example, we use Python
scripts to set up a directory structure, download Chandra
data, initialize `CIAO`_ etc. See `<https://github.com/Chandra-MARX/marx-test>`_
for the full code to run all tests.


.. toctree::
   :maxdepth: 2

   listofcode

