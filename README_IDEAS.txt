NOBLEFIT

NOBLEFIT is an flexible toolbox for quantitative interpretation of aquieous (noble) gas concentrations in terms of environmental processes.

NOBLEFIT is distributed as free software under the GNU General Public License (see LICENSE.txt).

Copyright (C) 2013 Matthias S. Brennwald.
Contact: matthias.brennwald@eawag.ch





Concepts and ideas

The idea behind NOBLEFIT is very similar to the noblegas fitter by Frank Peeters, but with some important differences:

1. The models are more modular, and users can write and use their own modules. This is thought to provide a more flexible system where special cases can be treated by specific modules.

2. There is no graphical user interface (GUI). This allows calling the fitter program from the shell, e.g., to run the fits with different data sets or to use the fitter results for further analyses without going through the hassle of reading data files in complicated formats.

Further ideas:
- NOBLEFIT should run on Matlab, but also on open source alternatives (e.g., Octave).
- NOBLEFIT may also include gases that are not noble gases, such as N2, SF6 or CFCs, etc.
- ...