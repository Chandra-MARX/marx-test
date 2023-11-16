# marx-test
The marx software itself is written in C and distributed using the traditional configure/make/make install
tool chain based on GNU make (and automake to make the makefiles).

However, tests for marx often require interaction with CIAO or with real Chandra data. Those tests run longer
and have requirements that marx itself does not have. Thus, we have a separate test suite for marx.

This splits into several parts:

- `CI`: Relatively quick test that could (though that's not set up right now) be run on github actions for every PR to marx itself.
- `speed`: A framework to compare marx runtimes with different compilers and compiler optimizations.
- `notebooks`: A set of Jupyter notebooks that demonstrate how to use marx and how to compare marx with real data.
- `tests`: Tests written using a separate test harness (which is defined in "marxtest"). Plan: Convert those to notebooks
- `marxtest`: Testing harness. In the past, it was not possible to run CIAO, Python, and Astropy in the same environment, nor was there a good way to control CIAO form Python or Python from CIAO. This testing framework is a workaround for that. Given that CIAO now easily runs in a conda environment, this may be obsolete and can be spimfied/removed.