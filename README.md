# CausalExplanationSimplicity

## causal_structure_R_data_to_pdf_graph.py

Python script that generates the possible causal diagramms from the R output file r_output.txt, using the latex template Latex_Template.tex. Both files must be in the same folder as the script file. 

**This is our main file (to be executed).**

This requires besides a Python3 environment, a latex distribution installed on the system running this script, as well as the following Python3 libraries:
codecs,
re,
jinja2.


## Latex_Template.tex

This template sets the environment for plotting the obtained mechanisms. It can be modified at will but must not be renamed.


## r_output.txt

The script requires input data in form of the CNA output standard. This file has to be generated or modified prior to running the Python script.


## R scripts

A folder that contains some sample R-scripts that make use of the CNA R-package. Some of the examples model multi level mechanisms. (not required for running the script)


## cna documentation

For information purposes the CNA documentation file and application examples are added in the respective folder.
