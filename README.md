# CausalExplanationSimplicity

## causal_structure_R_data_to_pdf_graph.py

Python script that generates the possible causal diagrams from the R output file r_output.txt, using the latex template Latex_Template.tex. Both files must be in the same folder as the script file. 

**This is the main file (to be executed).**

This requires besides a Python3 environment, a latex distribution installed on the system running this script, as well as the jinja2 Python3 library and the TikZ library for latex.


## Latex_Template.tex

This template sets the environment for plotting the obtained mechanisms. It can be modified at will but must not be renamed.


## r_output.txt

The script requires input data in form of the CNA output standard. This file has to be generated or modified prior to running the Python script.


## R_test_frames

A folder that contains some sample R-scripts that make use of the CNA-package. Some of the examples model multi level mechanisms. The content of this folder is not required for running the script.

- **test_1.r**  
causal structure from (Harbecke, 2018) with an added upper level

- **test_2.r**  
a causal chain with two levels, a very simple multi-level structure 

- **test_3.r**  
a variation of test_2.r which involves disjunctions, conjunctions and negations

- **test_4.r**  
a simple single level structure

- **test_5.r**  
the example data set "women" from CNA

- **test_6.r**  
the example data set "volatile" from CNA (circular causal structure)

- **test_7.r**  
the example data set "highdim" from CNA (it is not intended to be executed without further reduction of the data set)

- **test_8.r**  
an adjustable test frame, specify a formula and generate the respective truth table and Coincidence Analysis

- **test_9.r**  
causal structure from the introductory part of "An Algorithm for Boolean Constitutive Inference for Mechanistic Explanation and Potential Measures of Simplicity"


## cna documentation

For information purposes the CNA documentation file and application examples are added in the respective folder.
