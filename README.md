# marx-test
The marx software itself is written in C and distributed using the traditional configure/make/make install
tool chain based on GNU make (and automake to make the makefiles).

However, tests for marx often require interaction with CIAO or with real Chandra data. Those tests run longer
and have requirements that marx itself does not have. Thus, we have a separate test suite for marx.

This splits into several parts:

- `CI`: Relatively quick test that could (though that's not set up right now) be run on github actions for every PR to marx itself. These test will check marx help messages, run a few simple simulations, etc. but do not include long simulations or grids of simulations wth millions of photons.
- `speed`: A framework to compare marx runtimes with different compilers and compiler optimizations.
- `notebooks`: A set of Jupyter notebooks that demonstrate how to use marx and how to compare marx with real data.
- `tests`: Tests written using a separate test harness (which is defined in "marxtest"). Plan: Convert those to notebooks
- `marxtest`: Testing harness. In the past, it was not possible to run CIAO, Python, and Astropy in the same environment, nor was there a good way to control CIAO form Python or Python from CIAO. This testing framework is a workaround for that. Given that CIAO now easily runs in a conda environment, this may be obsolete and can be spimfied/removed.

## Notes on the CI tests
The tests in CI are meant to be run as continuous integration. They are in no way complete, but tests get added
when work on the C code happens, so they are testing any changes compared to marx circa marx 5.0.
The are to be run with the pytest framework, and they need marx in the path as well as an activated CIAO in 
the environment. 
Several of those tests go directly into the C code and use [cffi](https://cffi.readthedocs.io/en/stable/) to test individual C functions in marx or jdfits. This is a bit brittle, in the sense that it depends on the OS and compiler options. A few incomplete notes are below:

- On Linux with `gcc`, marx itself needs to be compiled with `-fPIC` to be able to load the shared library into the test. Even after doing that, I still run into `ImportError: /tmp/tmpp2lqsq9a/_acis_cffi.cpython-311-x86_64-linux-gnu.so: undefined symbol: PF_Errno` despite the `pfile` library being linked into the shared library. I've not debugged this further. I suspect that the age of the marx makefiles leads to settings that are non-standard today (and marx was never meant to be used as shared library), but since the tests pass on MacOS and the tests that call marx as a binary (instead of linking to it as a shared library) work, I have not pursued this further.
- On MacOS with `clang`, the tests work with default compiler options.