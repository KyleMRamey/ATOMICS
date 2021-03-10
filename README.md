AHTR Thermal Behavior Modeling and Iterative Criticality Suite

Code is written in C++. See C++ file "ATOMICS.cpp" for code.

Triangular superimposed mesh for detector tallies based on the superimposed hexagonal mesh within Serpent has its own directory. The changes were implemented within Serpent 2.31. This is no longer the most recent Serpent release. The capability is not tested with newer releases. See PhD dissertation documentation (Appendix B) for explicit changes to files. For Serpent source files which have not changed since the Serpent 2.31 release, the modified files can just be carried over to the new version. For files which have been changed since the Serpent 2.31 release, they will need to be manually updated with the documented code changes to achieve the same functionality. Note that line changes are likely for these files, so might need to refer to the Serpent 2.31 file for the relative location of the code changes.
