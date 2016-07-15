----------------
Tests for |marx|
----------------

No simulation is ever perfect. In this section, we shows tests for the internal
consistency of |marx| and compare simulations to real observations.
This demonstrates the limits to which |marx| simulations can be
trusted. Anybody using |marx| for scientific analysis needs to consult the
following pages to determine if the simulations are realistic enough for
their purpose.

|marx| was originally written to help with the calibration of the HETGS
instrument and later expanded for other purposes. Today, it includes modules
for all gratings and detectors on Chandra. Different modules have a different
level of detail (see :ref:`sect-modelsforthespacecraft`) and more details are
added with new |marx| versions. To explain how this impacts the accuracy of
|marx| simulations, we discuss a few of the design decisions behind the |marx|
code. 

For example, the |marx| HRMA module has fewer details than `SAOTRace`_,
but it runs much faster, allowing the user to simulate millions of rays in just
a few minutes on a typical desktop computer. Thus, a PSF simulated by |marx| should
differ from the observed PSF more than a PSF simulated with
`SAOTrace`_. However, we know that neither |marx| nor `SAOTrace`_ reproduce the
`"hook" feature of the PSF <http://cxc.harvard.edu/ciao/caveats/psf_artifact.html>`_ seen on very small scales, so both simulators are not
perfect. (Note that it is possible to combine |marx| and `SAOTrace`_ to use the
strength of both. See :ref:`sect-ex-simobs`.) In the end, only a comparison
between simulations and observed data can tell how well the simulations are
really doing.

A second example is the ACIS CCD. |marx| traces each ray to the exact
intersection point. Physically, the photon then interacts with the silicon of
the detector and this charge cloud might actually spread over several
pixels. This process `can be modeled <http://space.mit.edu/ACIS/ps_files/gyp_model_spie.pdf>`_ 
and it is possible to extract the flight grade of the detected event form the
model results, but this model is relatively slow and still carries some
uncertainty. Instead, |marx| draws flight grades from a fixed table included in
the code.



.. toctree::
   :maxdepth: 2

   {% for module in modulelist %}{{ module }}
   {% endfor %}

   listofcode

