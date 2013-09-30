hitran
======

hitran - a set of Python utilities for reading and manipulating HITRAN 2012 files.

A brief word of warning: although the ``hitran.py`` module file is quite small, the repository contains all of the HITRAN 2012 "parameter" files, so that the complete repository is ~500MB in size!

If run as an executable, then ``hitran.py`` will read in one of the ``.par`` files from the folder of HITRAN 2012 data, parse the file, convert the result to absorption cross-section, and display as a plot. (If one of the zipped par-files is selected, then it will automatically read directly from the zipped version.)
