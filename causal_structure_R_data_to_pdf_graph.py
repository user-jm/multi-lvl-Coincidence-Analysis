#!/usr/bin/env python3

# file: causal_structure_R_data_to_pdf_graph.py

# This script reads the output of a Coincidence Analysis (CNA) of some given truth table of a tiered list
# of causal factors and attempts to generate causal graphs for each unique solution. These are finally exported
# into a graph in Latex TikZ-code. Therefore latex and the tikz library have to be installed on the system running this script. 
# It assumes that the causal factors pertain to different levels, which are separated by the cna causal ordering relation "<".
# Relations between factors of different levels are not considered as causal but constitution relations. Accordingly the 
# script separates causal and constitution relations in the complex solution formulae of cna.

# It proceeds in eight steps:
# step 1: generates a list of the causal factors that appear in the cna-generated solutions (part of the function read_R_file)
# step 2: identifies the equivalence formulae in the cna-output (part of the function read_R_file)
# step 3: categorises the formulae into constitution relations and causal relations of different levels,
# already discards all formulae that do not fit in any of these categories (part of the function read_R_file)
# step 4: sorts out those cna-generated solutions that make no sense given the imposed level hierarchy (function cull_complex_solutions)
# step 5: identifies the incoming factors of each level - those that are not equivalent to complex expressions of other factors
# (part of the function determine_factor_order)
# step 6: defines the causal order on each level iteratively, starting with incoming factors a total causal order is imposed
# (part of the function determine_factor_order)
# step 7: minimises the constitution relations, while the causal relations are obtained from cna's complex solution formulae,
# the constitution relations are generated from its atomic solution formulae, there still many to be sorted out, since all
# possible equivalences are listed (function find_structure)
# step 8: translates the obtained structures into a graphical output via Latex


import os                       # operating system interfaces is required to find the files of the own path
import codecs                   # for en- and decoding of strings (esp. to get tex-files in utf-8)
import re                       # regex for complex search patterns in strings
import jinja2                   # Latex interface
from jinja2 import Template

# syntax definitions for Latex expressions
latex_jinja_env = jinja2.Environment(
	block_start_string = '\BLOCK{',
	block_end_string = '}',
	variable_start_string = '\VAR{',
	variable_end_string = '}',
	comment_start_string = '\#{',
	comment_end_string = '}',
	line_statement_prefix = '%%',
	line_comment_prefix = '%#',
	trim_blocks = True,
	autoescape = False,
	loader = jinja2.FileSystemLoader(os.path.abspath('.')))


# auxiliary function:
def find_causal_factors(st) :
    # looks in string st for denominators of causal factors, which are separated by ", " or " < "
    # returns a (possibly empty) list of causal factors
    
    # deletes "Factors: " from line (if it occurs)
    st = st.replace("Factors: ","")
    
    # deletes end-of-line-symbol ("\n") and spaces at the end of line if necessary
    st = re.sub("\r?\n","",st).rstrip()
    
    # returns the list of components of st that were separated by ", " or " < "
    return re.split(",\s*|\s*<\s*", st)
    
def convert_negation_syntax(st) :
    # conversion of the negation syntax (in cna by minuscle) such that "a" -> "~A"
    # 1. step: add "~" before each minuscle, which is either
    # a) at the beginning of a formula
    # b) follows a conjunctor
    # c) follows a disjunctor
    # d) follows an opening bracket
    st = re.sub(r'^([a-z])',  r'~\1', st)
    # explanation:  "sub" replaces each instance of a minuscle (expressed by "[a-z]")
    # by itself plus the prefix "~",
    # if it has been found at the first position of the string (implicated by "^")

    # b) if following a "*", the letter will be placed behind "*~"
    st = re.sub(r'\*([a-z])',  r'*~\1', st)
    # The regex expression "\*" picks the star symbol "*" from the string.

    # c) if following " + ", the letter will be placed behind "+ ~"
    st = re.sub(r'\s\+\s([a-z]+)',  r' + ~\1', st)
    # in regex "\s" corresponds to spaces, "\+" to "+"
    
    # d) if following a "(", the letter will be placed behind "(~"
    st = re.sub(r'\(([a-z])',  r'(~\1', st)
    # The regex expression "\(" picks the star symbol "(" from the string.

    # 2. step replacement of the minuscle that follow to "~" by majuscle
    st = re.sub(r'(~[a-z]+)', lambda pat: pat.group(1).upper(), st)
    
    return st
    
    
def get_equiv_formula(st) :
    # returns the leftside and rightside partial formulae of the "alpha <-> beta" from the input string st
    # as pairs of strings (a,b)
    a = re.split(" <-> ",st)[0].strip()          # strip() removes leading spaces
    b = re.split(" <-> ",st)[1].strip()
    
    a = convert_negation_syntax(a)
    
    # The lines of the cna output contain further stuff, we can get rid off it:
    b = re.split("[ \t]",b)[0]
    return (a,b)

def get_components_from_formula(st, factor_list):
    # returns a list of the elements of factor_list that appear in the input string st
    # the returned list is empty if no factor from factor_list appears in st or factor_list is empty, no list of strings or no list at all
    
    component_list = []                          # declaration of the list
    
    # since we use several nested lists of causal factors, we have to treat all cases separately
    # - the factors are elements of the list as in factor_list from main -> first case
    # - the factors are elements of the elements of the list as in level_factor_list from main -> second case
    # - the factors are elements of elements of the elements of the list as in level_factor_list_order from main -> third case
    
    if isinstance(factor_list[0], str) :
        # first case: check the elements from factor_list
        for element in factor_list :
            if st.find(element) > -1 :
                component_list.append(element)
    
    elif isinstance(factor_list[0], list) :
        if isinstance(factor_list[0][0], str) :
            # second case: traverse the sublists of factor_list and check for occurrences in st 
            for m in range(len(factor_list)) :
                for element in factor_list[m] :
                    if st.find(element) > -1 :
                        component_list.append(element)
                        
        elif isinstance(factor_list[0][0], list) :
            if isinstance(factor_list[0][0][0], str) :
                # third case: traverse the subsublists of factor_list and check for occurrences in st
                for m in range(len(factor_list)) :
                    for o in range(len(factor_list[m])) :
                        for element in factor_list[m][o] :
                            if st.find(element) > -1 :
                                component_list.append(element) 
    
    return component_list

def get_formula_level(st, level_factor_list):
    # if all factors in st are of the same level, get_formula_level returns this level,
    # otherwise it returns -1
    
    inequal = False
    factors = get_components_from_formula(st, level_factor_list) # list of factors that occur in st
    level = -1
    
    # Since we use different forms of nested lists, we have to distinguish the different cases:
    # case 1: multi_order = True -> level_factor_list is subdivided into levels and causal orders
    # case 2: multi_order = False -> level_factor_list is only subdivided into levels
    multi_order = isinstance(level_factor_list[0], list)
    
    
    # first case
    if multi_order :
        # get the level of factors[0]
        for m in range(len(level_factor_list)) :
            for i in range(len(level_factor_list[m])) :
                if factors[0] in level_factor_list[m][i] :
                    level = m
            
            # check if all remaining factors have the same level as factors[0]
            for k in range(1,len(factors)):
                inequal = True
                # for-loop over the orders
                for i in range(len(level_factor_list[level])) :
                    if factors[k] in level_factor_list[level][i] :
                        inequal = False
                        break  # after we have found the factor in the right level list, we can stop the loop over the orders
        
            if inequal:
                level = -1
    
    # second case
    else :
        # get the level of factors[0]
        for i in range(len(level_factor_list)) :
            if factors[0] in level_factor_list[i] :
                level = i
            
        # check if all remaining factors have the same level as factors[0]
        for k in range(1,len(factors)):
            if not(factors[k] in level_factor_list[level]) :
                inequal = True
                break  # we can stop after we have found the first factor that is not of the same level as factors[0]
        
        if inequal:
            level = -1
        
    return level

def get_factor_order(factor, factor_list):
    # determines the order of a causal factor relative to a factor list of the form
    # list[LEVEL][ORDER][FACTORS], if factor is contained in FACTORS, return the corresponding value of ORDER in the respective level
    # otherwise return -1
    order = -1
    for m in range(len(factor_list)):
        for o in range(len(factor_list[m])):
            if factor in factor_list[m][o]:
                order = o
                break       # we can stop after we have 
        if o != -1 : break  # determined the order
        
    return order
    
def get_formula_order(formula, factor_list):
    # determines the order of a formula = the highest order of the factors occurring in it
    # The order of the factors is determined by factor_list using the function get_factor_order.
    # if all factors have a determinable order, the maximum value is returned,
    # otherwise -1
    
    order = -1
       
    
    for fac in get_components_from_formula(formula, factor_list):
        if order < get_factor_order(fac, factor_list) :
            order = get_factor_order(fac, factor_list)
                
    return order


def create_term_list(solution_list) :
    # for a given list of formulae, this function splits conjunctions into a nested list of conjuncts for each formula
    # expected form of the elements of the input list is: "(...)*(...)* .... *(...)"
    
    solution_term_list = []
    
    for i in range(len(solution_list)) :        # for each solution from cna
        solution_term_list.append([])           # add a new sublist to solution_term_list
        
        # CNA's csf output comes in form of a table -- the formula has to be separated
        if solution_list[i].find("(") > -1 :
            string = "(" + solution_list[i].split("(",1)[1]  # cut everything before the first bracket
            
        else :
            string = solution_list[i]
        
        # remove evertything after the last closing bracket
        if solution_list[i].find(")") > -1 :
            string2 = ""
            sub_string = ""
            for st in string.split(")") :
                string2 = string2 + sub_string
                sub_string = st + ")"
            string = string2
         
        # the entries of solution_list have the form "(A*B + C <-> E)*(D*E + B + ~C <-> F)* ... "
        # split these formulae into their terms (split by occurrences of the string ")*(")
        # and add the parts to solution_term_list[i]
        solution_term_list[i] = re.split("\)\*\(", string)
        
        # remove first opening bracket of the first element
        if solution_term_list[i][0][0] == "(" :
            solution_term_list[i][0] = solution_term_list[i][0][1:] # "[1:]" gives the whole string up from the second character
        
        # remove last closing bracket of the last element
        if solution_term_list[i][len(solution_term_list[i])-1][len(solution_term_list[i][len(solution_term_list[i])-1])-1] == ")" :
            # this reads: if the last character of the last element of solution_term_list[i] is ")"
            
            solution_term_list[i][len(solution_term_list[i])-1] = solution_term_list[i][len(solution_term_list[i])-1][:-1]
            # "[:-1]" returns the input string minus its last character    
            
    return solution_term_list


def get_causal_prefactors(factor, formula_list, factor_list) :
    # returns a list of direct and indirect causal prefactors of a given factor relative to a list of causal relations
    
    return_list = []
    
    for formula in formula_list :
        if formula[1] == factor :
            # if the factor we are interested in is the target factor of the relation formula,            
            # then add all factors that appear on the left side of formula to the return list
            return_list.extend(get_components_from_formula(formula[0], factor_list))
            
            for pfac in get_components_from_formula(formula[0], factor_list) :
                # get the indirect prefactors recursively
                return_list.extend(get_causal_prefactors(pfac, formula_list, factor_list))
    
    return_list = list(set(return_list))  # get rid of duplicates
    return return_list

def convert_causal_relation(formula, level_factor_list_order, tex_code) :
    # translates a formula of causal relations into TikZ-Latex code
    # returns the code as string
    
    st = ""  # output code
    
    ###########################################
    # determine the syntactic type of formula #
    ###########################################
    # assumption: the cna output contains only the following types of formula:
    # 1) equivalence between causal factors (A <-> B) 
    # 2) non-equivalence between causal factors (~A <-> B)
    # 3) equivalence of one causal factor with conjunctions of factors ( A*B* ... <-> E)
    # 4) equivalence of one causal factor with conjunctions of factors and negated factors ( A*~B* ... <-> E)
    # 5) equivalence of one causal factor with disjunctions of factors ( A + B + ... <-> E)
    # 6) equivalence of one causal factor with disjunctions of factors and negated factors ( A + ~B + ... <-> E)
    # 7) equivalence of one causal factor with disjunctions of at least one conjunct of factors and possibly negated factors
    #    ( A + ~B*C + ... <-> E)
    
    
    level = get_formula_level(formula[0], level_factor_list_order)
    
    #######################################
    # equivalences with (negated) factors #
    #######################################
    
    for o in range(len(level_factor_list_order[level])) :
        for fac in level_factor_list_order[level][o] :
            if fac == formula[0] :  # the left side of formula equals one factor from the factor list
                st = "\draw[->] (" + fac + ".east) -- (" + formula[1] + ".west);"
                break    # stop after the factor has been found
            elif formula[0] == "~" + fac : # the left side of formula equals the negation of one factor from the factor list
                st = "\\node[neg] (" + fac + "neg) at ([xshift=\LNeg]" + fac + ".south east) {};\n\draw[->] (" + fac + "neg) -- (" + formula[1] + ");"
                break # stop after the factor has been found
                
                
    if st == "" :
        if formula[0].find("+") > -1 :
        
            ################
            # disjunctions #
            ################
            
            # create a list of the possibly complex disjuncts of formula
            disjunctor_list = re.split("\s*\+\s*", formula[0])  
            
            
            for disj in disjunctor_list : 
                # running through this list
                # each disjunct is tested whether it is
                # A) a simple causal factor
                # B) a conjunction
                # C) a negated causal factor
                # each case is treated separately

                
                if disj in get_components_from_formula(formula[0], level_factor_list_order) :
                    # case A: the discunct is one causal factor
                    
                    # a straight arrow is drawn from source factor to target factor
                    st = st + "% simple disjunction\n"
                    st = st + "\draw[->] (" + disj + ".east) to (" + formula[1] + ".west);\n"
                
                elif formula[0].find("*") > -1 :
                    # case B: the disjunct is a conjunction
                    st = st + "% complex disjunction\n"
                    
                    # first, the conjuncts are connected by curved lines meeting in one junction point
                    # placed one right (with a slight upward shift) to factor of the highest causal order
                    # second, this junction point is connected with the target factor by a straight arrow like
                    # in case A)
                    conjunctor_list = re.split("\*", disj)
                    
                    # set the junction of the conjuncts
                    # place it beside the (first) conjunct of the highest causal order (the factor that is most to the right in the graph)
                    # find the corresponding node of that conjunct -- f_fac
                    cross_point = ""                    # name of node of the junction of the conjunctions
                    for fac in get_components_from_formula(disj, level_factor_list_order) :
                        cross_point = cross_point + fac
                        if fac == get_components_from_formula(disj, level_factor_list_order)[0] :
                            f_fac = fac
                        else :
                            if get_factor_order(fac, level_factor_list_order) > get_factor_order(f_fac, level_factor_list_order) :
                                f_fac = fac
                    
                    cross_point = cross_point + formula[1]
                    
                    if get_factor_order(f_fac, level_factor_list_order) < get_factor_order(formula[1], level_factor_list_order) :
                        # this is the normal non-circular case
                        position = "at ([xshift=\hDisjConj, yshift=\\vDisjConj]" + f_fac + ".east)"
                        circular = False
                    else :
                        # circular case, f_fac has the same horizontal coordinate as the target factor,
                        # so the junction should not be placed to the right of f_fac but to the left
                        
                        position = "at ([xshift=\hcDisjConj, yshift=\\vDisjConj]" + f_fac + ".west)"
                        circular = True
                    
                    # Attention: It might happen that several disjuncts of conjuncts meet at the same factor f_fac,
                    # therefore we have to check whether the position of the junction node has to be shifted.
                    q = 1
                    while (tex_code.find(position) > -1) or (st.find(position) > -1) :  
                        # this position has already been specified in earlier vertices (tex_code) or this one
                        # -> shift it above by \tDisjConj
                        if circular :
                            position = "at ([xshift=\hcDisjConj, yshift={\\vDisjConj + " + str(q) + "*\\tDisjConj}]" + f_fac + ".west)"
                        else :
                            position = "at ([xshift=\hDisjConj, yshift={\\vDisjConj + " + str(q) + "*\\tDisjConj}]" + f_fac + ".east)"
                        
                        q = q + 1
                        
                    st = st  + "% junction of the conjuncts\n\\node[aux] (" + cross_point + "aux) " + position + " {};\n% partial arrows from the conjuncts to the junction\n"
                        
                    
                    
                        
                    

                    
                    for conj in conjunctor_list :
                    
                        # now connect the conjuncts with the junction
                        if conj[0] == "~" and conj[1:] in get_components_from_formula(disj, level_factor_list_order) :
                            # case B i) the conjunct is a negated factor
                            st = st  + "\\node[neg] (" + conj[1:] + "neg) at ([xshift=\LNeg]" + conj[1:] + ".south east) {};\n"
                            st = st + "\draw[conjunctonsegment] (" + conj[1:] + "neg) to (" + cross_point + "aux);\n"
                            
                        elif conj in get_components_from_formula(disj, level_factor_list_order) :
                            # case B ii) the conjunct is a mere factor
                            st = st + "\draw[conjunctonsegment] (" + conj + ".east) to (" + cross_point + "aux);\n"
                            
                        else :
                            # Is there anything else that might happen??
                            print("The disjunction " + formula[0] + "  could not be plotted because its substructure was not recognized(1).")

                    # connect the junction with the target factor
                    
                    # experimental label above the connection for very convoluted graphs
                    st_conj = "$"
                    for conj in conjunctor_list :
                        if conj[0] == "~" :
                            st_conj = st_conj + "\\neg " + conj[1:] + "\cdot "                    
                        else :
                            st_conj = st_conj + conj + "\cdot "
                        
                    st_conj = st_conj[:-6] + "$"
                    st = st + "% arrow from junction to target factor\n\draw[->] (" + cross_point + "aux) -- (" + formula[1] + ".west) node[draw=none,fill=none,font=\\tiny,pos=0,sloped,above=\LabelDist] {\scalebox{.3}{" + st_conj + "}};\n"
                
                elif (disj[0] == "~") and (disj[1:] in get_components_from_formula(formula[0], level_factor_list_order)) :
                    # case C: the disjunction is a negated factor
                    # = the first character is "~" and the further characters correspond to one element from level_factor_list_order
                    
                    # assumption: one disjunctive chain can contain either a mere or its negation (otherwise above "elif" has to been
                    # changed to "if")
                    
                    st = st + "% negated disjunct\n"
                    st = st  + "\\node[neg] (" + disj[1:] + "neg) at ([xshift=\LNeg]" + disj[1:] + ".south east) {};\n"
                    st = st + "\draw[->] (" + disj[1:] + "neg) to (" + formula[1] + ".west);\n"
                
                
                else :
                    # disjunction is of a different form (is there any further possible??)
                    
                    print("The disjunction " + formula[0] + "  could not be plotted because its substructure was not recognized.")

            
            
        elif formula[0].find("*") > -1 :
            
            ################
            # conjunctions #
            ################
            
            # procedure is the same as conjuncts in disjunctions
            # the only difference is that we here do not deal with a subformula of formula[0], but the whole
            
            # place junction of the conjuncts
            st = "% junction of the conjuncts\n\\node[aux] (" + formula[1] + "aux) at ([xshift=\LConj]" + formula[1] + ".west) {};\n% partial arrows from the conjuncts to the junction\n"
            
            # plot the arrows from the conjuncts to the junction
            for conj in get_components_from_formula(formula[0], level_factor_list_order) :

                if formula[0].find("~" + conj) > -1 :
                    # case A: the factor conj appears negated n formula
                    
                    # assumption: a conjunction chain can only contain a factor or its negation
                    # (otherwise the subsequent "else" has to be changed into a new "if")
                    
                    st = st  + "\\node[neg] (" + conj + "neg) at ([xshift=\LNeg]" + conj + ".south east) {};\n"
                    st = st + "\draw[conjunctonsegment] (" + conj + "neg) to (" + formula[1] + "aux);\n"
                else :
                    # case B: factor conj occurs non-negated
                    
                    st = st + "\draw[conjunctonsegment] (" + conj + ".east) to (" + formula[1] + "aux);\n"
            
            # draw arrow from junction to target factor
            # experimental with tiny label above vertex
            st = st + "% arrow from junction to target factor\n\draw[->] (" + formula[1] + "aux) -- (" + formula[1] + ") node[draw=none,fill=none,font=\\tiny,above=\LabelDist,pos=0,sloped] {\scalebox{.3}{$" + formula[0].replace("*", "\cdot ").replace("~", "\\neg ") + "$}};\n"
            
        else :
            # this should never happen
                        
            ##############################
            # structure not recognisable #
            ##############################
            
            print(formula[0] + " -> " + formula[1] + "  could not be drawn in because its structure was not recognized.")
            
    return st
    
def convert_constitution_relation(formula, level_factor_list_order, constitution_relation_list) :
    # converts a formula of constitution relations into TikZ-Latex code
    # returns the code as string
    
    st = ""
    
    # constitution relations are drawn differently depending on whether they are to the left or to the right of the upper level factor
    c_left = False
    c_right = False
    
    # check whether it is a left- or rightside relation
    for f in constitution_relation_list :
        if (formula[1] == f[1]) and (formula[0] != f[0]) :
            # is there a further constitution relation to the same causal factor which includes factors of higher causal order
            # than those from formula? -> if true it is a leftside relation
            # if there is no further constitution relation it is neighter left- nor rightside
            # if there further relations but of lower order -> rightside relation
            if get_formula_order(formula[0], level_factor_list_order) < get_formula_order(f[0], level_factor_list_order) :
                c_left = True
            elif get_formula_order(formula[0], level_factor_list_order) > get_formula_order(f[0], level_factor_list_order) :
                c_right = True
    
    # draw one connecting line toward formula[1] for each causal factor in formula[0]
    # (usually there should only be one factor in formula[0])
    for fac in get_components_from_formula(formula[0], level_factor_list_order) :           
        if c_left and not(c_right) :
            # case 1: leftside relation
            st = st + "\draw[crelationleft] (" + fac + ".north west) to (" + formula[1] + ".south);\n"
    
        elif not(c_left) and c_right :
            # case 2: rightside relation
            st = st + "\draw[crelationright] (" + fac + ".north east) to (" + formula[1] + ".south);\n"
    
        else: 
            # case 3: otherwise
            st = st + "\draw[crelationstraight] (" + fac + ".north) to (" + formula[1] + ".south);\n"
    
    return st

# main functions start here
def read_R_file(file_name):
    # gets the relevant data of the cna output and returns in some lists
    # level_factor_list -- list of causal factors separated into one sublist for each constitution level,
    # level_equiv_list -- list of causal relations separated into one sublist for each constitution level,
    # constitution_relation_list -- list of constitution relations
    
    
    file_lines = []                              # declaration of a list to include the lines of text of the R output file
    with open (file_name, 'rt') as text_file:    # open file_name
        for next_line in text_file:              # for each line in that file do
            file_lines.append(next_line)         # append it to the list file_lines
            
    abort = False                                # since corrupted/unexpected input cannot be used, we check at several stages
                                                 # whether we may continue or not - this what the variable abort is for
           
    ######################################## 
    # step 1: determine the causal factors #
    ########################################
    factor_list = []                             # declaration of factor list
    
    
    # Attention the following might change if the formatting of the cna output changes
    for i in range(len(file_lines)) : # search for the list of causal factors in the R output
        if file_lines[i].find("Causal ordering:") > -1:
            # case 1: if the factors are divided into several levels, cna prints "Causal ordering:"
            # then the factors are listed in the subsequent line
            factor_list = find_causal_factors(file_lines[i+1]) 
            file_line_factors = i + 1            # Later it will be helpful to know the line where the factors are listed. 
            break                                # leave for-loop after the line has been found
        
        elif file_lines[i].find("Factors:") > -1:
            # case 2: the R input does not include a separation of causal factors into different levels,
            # then the output contains "Factors:" followed by the causal factors in the same line
            factor_list = find_causal_factors(file_lines[i])
            file_line_factors = i
            break
    
    level_count = file_lines[file_line_factors].count("<") + 1 # level_count = number of constitution levels
    
    if not(factor_list) :  # factor_list is empty
        print("Abort no causal factors have been found in " + file_name)
        abort = True
        return abort, "", "", "", ""
    
    # continue if factor_list is non-empty
    else :
        #########################
        # step 2: find formulae #
        #########################
        
        equiv_list = []                          # declaration of the list for atomic solution formulae
        solution_list = []                       # declaration of the list for the complex solutions
            
        for line in file_lines :
            if line.count("<->") == 1 :          # cna's atomic solution formulae 
                # exactly one "<->" has been found in the line
                # read the partial formulae on its left and right side and add them to equiv_list
                equiv_list.append(get_equiv_formula(line))
                
            elif line.count("<->") > 1 :         # cna's complex solution formulae
                # more than one "<->" has been found in the line
                
                # delete end-of-line-symbol ("\n") and spaces at the end of line (with rstrip()) if necessary, as well as
                # removes leading spaces (with strip())
                # then add this line to solution_list                
                solution_list.append(convert_negation_syntax(re.sub("\r?\n","",line).rstrip().strip()))
                
                
                
        same_as_asf = False
        if not(solution_list) :  # if solution_list is empty (no line with more than one occurrence of "<->" has been found)
            for line in file_lines :
                if line.find("Same as asf") > -1 : # If the complex solutions are equal to atomic solutions, cna returns "Same as asf".
                    same_as_asf = True
                    
                    st = ""
                    for i in range(len(equiv_list)) :
                        if st == "" :  # first formula
                            st = equiv_list[i][0] + " <-> " + equiv_list[i][1] # recompose the logical formula from equiv_list
                        elif st[0] != "(" : # second formula
                            # In case that there are more than one atomic solution formulae, we connect them via conjuntion 
                            # and have to put them into brackets and add some brackets for the first formula
                            st = "(" + st + ")*(" + equiv_list[i][0] + " <-> " + equiv_list[i][1] + ")"
                        else : # after the second formula
                            st = st + "*(" + equiv_list[i][0] + " <-> " + equiv_list[i][1] + ")"
                        
                    
                    solution_list.append(st)
                    break  # break for-loop running over the text lines
                    
        if not(equiv_list) or (not(solution_list) and not(same_as_asf)):  # if equiv_list is empty or
            # solution_list is empty while not flagged as same_as_asf, then abort
            print("Abort no formula has been found in " + file_name)
            abort = True
            return abort, "", "", "", ""
        else:                                    # otherwise continue
            # check whether each causal factor appears in some atomic formula, otherwise remove it from factor_list

            for k in range(len(factor_list)-1,-1,-1):
                i = 0
                found = False
                while not(found) and i < len(equiv_list):
                    found = (equiv_list[i][0].find(factor_list[k]) > -1 or equiv_list[i][1] == factor_list[k])
                    i = i + 1
                    
                if not(found):  # since the for-loop is regressive, it should be no problem to remove the elements from the list
                    # within the loop
                    print("Factor " + factor_list[k] + " has been discarded, since it does not occur in any formula.")
                    del factor_list[k]
                    
            
            
            #####################################
            # step 3: categorising the formulae #
            #####################################
            
            # separating into the constitutive levels
            
            
            
            level_factor_list = []               # declaration of new lists
            level_equiv_list = []
            constitution_relation_list = []
            
            for i in range(level_count):
                if level_count > 1: # multi-level case
                    st = re.split(" < ",file_lines[file_line_factors])[i].strip()
                else :   # single-level case
                    st = file_lines[file_line_factors]
                    
                level_factor_list.append(find_causal_factors(st))
                level_equiv_list.append([])
                
            for formula in equiv_list : 
                # add all formulae that contain only factors from level i to level_equiv_list[i]
                if get_formula_level(formula[1], level_factor_list) == get_formula_level(formula[0], level_factor_list) :
                    level_equiv_list[get_formula_level(formula[1], level_factor_list)].append(formula)
                    
                elif get_formula_level(formula[0], level_factor_list) > -1 :
                    # all formulae, which relate the element on the right side with factors of one different level on the left side
                    # are added to constitution_relation_list
                    # assumption: only constitution relations with level difference of one are maintained
                    
                    if get_formula_level(formula[0], level_factor_list) == get_formula_level(formula[1], level_factor_list) - 1 :
                        constitution_relation_list.append(formula)
                        
                           
                # all further formulae will not be considered any longer
    return abort, level_factor_list, level_equiv_list, constitution_relation_list, solution_list
    

def cull_complex_solutions(solution_list, level_factor_list, level_equiv_list, constitution_relation_list) :
    # This functions rejects cna-solutions that do not meet the specifities of multi-level structures.
    # It returns a reduced list of elements from solution_list that 
    # A) contain formulae which exihibit a mixture of causal factors from different levels within one disjunction or conjunction
    # B) are incomplete after constitution relations have been separated -- some causal factors are then absent in some solutions
    # C) are duplicates of other solutions (originally they might have differed in terms that are constitution relations)
    # It constitutes the fourth step of the overall procedure.
    
    ##########################################################################
    # step 4: single out valid and unique solutions for the causal structure #
    ##########################################################################
    
    solution_term_list = [] # elements of solution_list are conjunctions of equivalences, these will be separated and added to this list
    
    # fill this list
    solution_term_list = create_term_list(solution_list)        
            
    #####################################
    # step 4A: discard invalid solutions #
    #####################################
    
    # sorting out solutions that mix causal factors from different levels on the left side of an equivalence
    print("Number of solutions before filtering " + str(len(solution_list)))
    for i in range(len(solution_list)-1,-1,-1) : # regressive for-loop over all elements of solution_list
        # (i goes from [len(solution_list) - 1] to 0)
        
        valid = True
        for term in solution_term_list[i] :
            # does term not match any valid atomic solution formula?
            if not(any(get_equiv_formula(term) in j for j in level_equiv_list)) and not(get_equiv_formula(term) in constitution_relation_list) :
                # explanation: any(a in b for b in c) checks whether a is an element of a sublist of c

                valid = False  # then flag the corresponding solution as invalid
                break          # one invalid term is sufficient, leave the for-loop over the terms if one has been found
        
        
        if not(valid) :
            del solution_list[i]      # since the for-loop is regressive list elements can safely be deleted
            del solution_term_list[i]
           
            

    print("Number of solutions after sorting out invalid terms " + str(len(solution_list)))
    
    ##############################################
    # step 4B: remove incomplete causal relations #
    ##############################################
    
    # Since cna does not distinguish between causal and constitution relations, we might end up with solutions in which
    # some factors do not appear anymore. These solutions will be discarded next.
    
    # first remove constitution relations from solution_term_list
    
    for j in range(len(solution_list)) :
        for i in range(len(solution_term_list[j]) - 1, -1, -1) :
            if get_equiv_formula(solution_term_list[j][i]) in constitution_relation_list :
                del solution_term_list[j][i]
                
    # now check that every factor still appears for every solution
    for i in range(len(solution_list)-1,-1,-1) : # regressive for-loop over all elements of solution_list
        complete = True
        for lvl in level_factor_list : # go through all levels
            
            if len(lvl) > 1 : # if that level contains only one factor, there is no causal relation on that level
                # possible problem: a level that contains incoming factors only
                # at this point these cannot be detected
                
                
                for fac in lvl :       # go through all factors of that level
                    found = False
                    for term in solution_term_list[i] :
                        if fac in get_components_from_formula(term, level_factor_list) : 
                            # if fac occurs in term, we can stop looking for that factor
                            found = True
                            break      # hence, break the term-loop
                    
                    if not(found) :    # if one factor hasn't been found in any term, this solution is to be discarded
                        complete = False
                        print(fac)
                        break         # break from fac-loop
        
            if not(complete) : break # break from lvl-loop if one absent factor has been found
        
        if not(complete) :
            del solution_list[i]      # since the for-loop is regressive list elements can safely be deleted
            del solution_term_list[i]  
            
        
    print("Number of solutions after sorting out causally incomplete structures " + str(len(solution_list)))
    
    # reconstruct solution_list without constitution relations
    new_solution_list = []
    
    for sol in solution_term_list :
        st = "("
        for term in sol :
            st = st + term + ")*("
            
        # removing the surplus "*(" at the end of st
        st = st[:-2]
        
        #################################
        # step 4C: discarding duplicates #
        #################################
        
        if not(st in new_solution_list) : # checking for duplicates
            new_solution_list.append(st)
     
    
    print("Number of solutions after discarding duplicates " + str(len(new_solution_list)))
    
    return new_solution_list
    
def determine_factor_order(level_factor_list, level_equiv_list) :    
    # uses the list of causal relations in level_equiv_list to determine a total causal ordering of the causal factors
    # in level_factor_list for each level separately
    # returns a Boolean value whether the process has been aborted due to some error, another Boolean whether this
    # ordering permits a unique order and the new factor list level_factor_list_order that is nested twice by level and causal order
    
    level_count = len(level_factor_list)
    abort = False
    unique = True 
    
            
    level_factor_list_order = []  # declaration of a new list for causal factors with one sublist for each level that contains one
    # sublist for each causal order in this level, wherein we find the causal factors
    # order = 0 -> incoming factors, 
    # order = i -> this is a target factor, all of its source factors are of order < i and at least one is of order i - 1
            
    # starting from here, all steps will be executed separately for each level
    for m in range(level_count):
                
        ##########################################
        # step 5: determine the incoming factors #
        ##########################################
                
                
        level_factor_list_order.append([]) 
        level_factor_list_order[m].append([])  # add the empty list for factors of order 0
                

        for element in level_factor_list[m]:
            # incoming factors do never appear on the right side of equivalence formulae
            i = 0
            found = False
            while not(found) and i < len(level_equiv_list[m]):
                found = (level_equiv_list[m][i][1] == element)
                i = i + 1
                    
            if not(found):
                # add element to incoming_factor_list
                level_factor_list_order[m][0].append(element)
        
        if not(level_factor_list_order[m][0]):
            # no incoming factors have been found
            # Exist circular relations? (A<->B; B<->A)
            # yes -> structure is non-unique
            # no -> cannot continue with input data
            circular = False
            
            for element in level_factor_list[m] :
                e_circular = True
                for formula in level_equiv_list[m] :
                    if formula[1] == element :
                        if not((formula[1],formula[0]) in level_equiv_list[m]) :
                            e_circular = False
                  
                if e_circular :
                    circular = True
            
            unique = not(circular)
            if unique :
                # no incoming factors found but no circularity determined -> abort
                print("Abort, no incoming causal factors have been found in the data of " + file_name + ".")
                abort = True
                    
           
                
        else:
            # level_factor_list_order[m][0] is not empty
                        
                
            ###########################################
            # step 6: define causal order iteratively #
            ###########################################
                    
            # create a total ordering of all causal factors of one level into their causal order
            # this is done iteratively starting with the yet obtained order zero by categorising all further factors
            
            # create a list of all yet non-categorised factors (= those that are causally down_stream to the considered factors)
            downstream_factor_list = []
                
            # initially downstream_factor_list consists of those elements of level_factor_list[m], that are not of order 0
            for element in level_factor_list[m]:
                if not (element in level_factor_list_order[m][0]):
                    downstream_factor_list.append(element)
                        
           
                    
            # successively add to the list level_factor_list_order[m] those factors which appear on the right side of causal relation
            # whose left side factors are all already contained in level_factor_list_order[m]
            # this is done by four nested lists:
            # 1) a while loop that runs over the indexes of the elements of downstream_factor_list
            # it may have to pass the same element multiple times since it might be necessary that other factors are to
            # be classified first
            # 2) a for loop over all causal equivalence formulae
            # checks whether in all formulae where the considered element stands on the right side (as a singular term)
            # all terms on the left side have an order assigned, therefore a third loop has to been run over
            # 3) a for loop over the causal factors that appear on the left side of the current formula
            # 4) a for loop over all orders in level_factor_list_order to check whether the elements from 3) are already contained
            # in one sublist

            # traverse downstream_factor_list regressively (loop 1)
            j = len(downstream_factor_list) - 1
            while j > -1 :
                classifiable = False    # Boolean value whether a causal order can be assigned to downstream_factor_list[j] given
                # the current level_factor_list_order
                order = 0               # the order that will be given to downstream_factor_list[j]
                
                # loop 2 - for loop over all causal formulae
                for formula in level_equiv_list[m] :
                    if formula[1] == downstream_factor_list[j] :
                        # go through all formulae where element j is the right side term
                        # then check whether in all of these formulae every factor on the left side has already an order assigned
                        # if so, element j gets the max order + 1
                        # if not, continue with the next element of downstream_factor_list
                        
                        order = 0 # order has to be reset to zero
                        
                        # loop 3 - for loop over all factors on the left side of formula
                        for fac in get_components_from_formula(formula[0], level_factor_list) :
                            fac_is_listed = False      # is there already an order assigned to the currently considered factor fac?
                            
                            # loop 4 - for loop over all orders
                            for o in range(len(level_factor_list_order[m])) :
                                if fac in level_factor_list_order[m][o] :
                                    fac_is_listed = True
                                    if order < o + 1 :
                                        order = o + 1  # the order of downstream_factor_list[j] is at least one higher than that of fac
                                        
                                    break              # break from loop over orders, if fac has already been found
                            # end of loop 4 over orders        
                            
                            if not(fac_is_listed) :
                                classifiable = False   # if one source factor has no assigned order, the target factor is not (yet)
                                # classifiable
                                break                  # break from loop over factors, since one non-categorised factor suffices
                                
                            else :
                                classifiable = True    # an order can be assigned to target factor j (given the current information)
                                
                        # end of loop 3 over factors of formula[0]
                        
                        if not(classifiable) :
                            break                      # break from loop over formulae after one has been found that turns out that
                            # factor j is unclassifiable by now
                            
                # end of loop 2 over all causal relations
                
                if classifiable :
                    if len(level_factor_list_order[m]) < order :
                        # the determined order of j is by at least 2 higher than that of the highest elements of the current list
                        # this should never happen - and possibly can't
                        print("Error in determining the order of the causal factors")
                        abort = True
                        break                          # break from while loop, something went wrong with level_factor_list_order
                    else :
                        if len(level_factor_list_order[m]) == order :
                            # factor j is the first one of its order
                            # a new sublist is to be created
                            level_factor_list_order[m].append([])
                        
                        # j to the ordered factor list    
                        level_factor_list_order[m][order].append(downstream_factor_list[j])     
                        
                        # and delete it from the unordered one
                        del downstream_factor_list[j]
                        
                        # proceed with the last element of downstream_factor_list
                        j = len(downstream_factor_list) - 1
                            
                else :
                    j = j - 1                          # proceed with the next factor in downstream_factor_list
                        
            # end of loop 1 over the elements of downstream_factor_list
            
            # output as console text     
            #print("causal factors of level " + str(m) + " with their causal order:")
            #for p in range(len(level_factor_list_order[m])):
            #    print("order " + str(p))
            #     print(level_factor_list_order[m][p])
            
            
                       
            if downstream_factor_list:
                # If downstream_factor_list is not empty, there are some causal factors that cannot be categorised
                # into level_factor_list_order.
                
                # test-wise set the causal order of all remaining factors to max order + 1
                circular = True
                
                level_factor_list_order[m].append([])
                
                for fac in downstream_factor_list :
                    level_factor_list_order[m][len(level_factor_list_order[m]) - 1].append(fac)
                    
            else :
                circular = False 
                
               
    return abort, unique, circular, level_factor_list_order
    
    
    
    
    
    
def find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list) :
    # this functions is intended to minimise constitution_relation_list and to
    # return the minimal list of constitution relations
    
    # In most cases constitution_relation_list contains many relations that should not be plotted. 
    # These are either wrong, since they include the causal preconditions of a complex that constitutes a higher level factor,
    # or dispensable, since the are middle terms of the complex.
    # e.g. if a lower level causal structure (A*B -> C) * (D*E -> F) * (C*F -> G) constitutes F1 on an upper level,
    # we only want the connectives A <-> F1, B <-> F1, D <-> F1, E <-> F1 and G <-> F1,
    # also relations between F1 and further preconditions of A,B,D,E should be discarded
    
    ################################################################
    # step 7: discard wrong and dispensable constitution realtions #
    ################################################################
    
    
    level_count = len(level_factor_list_order)
    new_constitution_list = []
    return_list = []  
    
    for m in range(1,level_count) :
        # start at level one, since the first level (not the zeroth) is the first one that might have constitution relations 
        
        for o in range(len(level_factor_list_order[m])) :
            # loop over all causal orders of level m

            for fac in level_factor_list_order[m][o] :
                # loop over all factors of order o and level m
                
                # create an auxiliary list of constitution relations with respect to fac
                auxiliary_list = []
                for c in constitution_relation_list :
                    # loop over all constitution relations
                    if c[1] == fac :
                        for l_fac in get_components_from_formula(c[0], level_factor_list_order) :
                            if not(l_fac in auxiliary_list) :
                                # add the factors on the left side of c to the auxiliary list if fac stands on its right side
                                auxiliary_list.append(l_fac) 
                    
            
                if len(auxiliary_list) > 2 :
                    # further steps are only necessary if auxiliary_list contains more than two entries
                    # 1) find the highest order -> factors of this order a kept for sure
                    # 2) determine the lowest order that should be kept
                    # a) if fac is of order 0 (on level m), factors of downto order 0 (on level m-1) should be kept
                    # b) if fac is of a higher order, its lowest factors must not be in constitution relation
                    # with fac's causal pre-factors
                    
                    max_order = 0
                    
                    # determine the value of max_order
                    for l_fac in auxiliary_list :
                        if get_factor_order(l_fac, level_factor_list_order) > max_order :
                            max_order = get_factor_order(l_fac, level_factor_list_order)

                            
                    if o > 0 :
                        # if the considered factor fac on the higher level is not an incomming factor, check whether its alleged
                        # constitution factors already constitute a causally upstream factor on fac's level 
                        # -> if it has been found remove it from the prospective list for fac
                        
                        # a list of all direct and indirect causal prefactors of fac
                        prefactor_list = get_causal_prefactors(fac, level_equiv_list[m], level_factor_list_order[m])
                        
                        for l in range(len(auxiliary_list) - 1, -1, -1) :
                            for pfac in prefactor_list :
                                pair = (auxiliary_list[l], pfac)
                                if pair in new_constitution_list :
                                    del auxiliary_list[l]
                                    break   # break from loop over prefactors, since we are done with this element fron auxiliary list
                                                   
                    # otherwise if the considered factor on the higher level is an incomming factor, do not restrict the lowest order
                    # of factors from auxiliary_list                         

                # carry over the entries from auxiliary_list into new_constitution_list    
                for lfac in auxiliary_list :
                    entry = (lfac, fac)
                    if not(entry in new_constitution_list) :
                        new_constitution_list.append(entry) 
                
                # clear the auxiliary list    
                auxiliary_list.clear()                    
                    
                    
        # now discard constitution relations to terms that are middle terms of causal chains whose
        # upstream and downstream factors are also in a constitution relation with the considered higher level factor
        
        # completely new loops are necessary in order to make the first one run with correct results
        # (some relations that will be discarded now, are necessary before in order to test whether some constitution relations
        # overlap with those of causal pre-factors
        
        for o in range(len(level_factor_list_order[m])) :
            # loop over all causal orders of level m

            for fac in level_factor_list_order[m][o] :
                # loop over all factors of order o and level m
                
                # create a new auxiliary list of constitution relations with respect to fac
                auxiliary_list = []
                                    
                for entry in new_constitution_list :
                    # obtain the results from the loop above
                    if entry[1] == fac :
                        auxiliary_list.append(entry[0])
                    
                if len(auxiliary_list) > 2 :
                    for l in range(len(auxiliary_list) - 1, -1, -1) :
                        
                        # get the level of auxiliary_list[l]
                        level = 0
                        for lvl in range(len(level_factor_list_order)) :
                            if any(auxiliary_list[l] in o for o in level_factor_list_order[lvl]) :
                                # in this if-clause we check whether auxiliary_list[l] is element in any sublist of
                                # level_factor_list_order[lvl] (the sublists are the order lists, which contain the causal factors)
                                
                                level = lvl
                                break          # break from the lvl loop, after it has been obtained

                        if get_causal_prefactors(auxiliary_list[l], level_equiv_list[lvl], auxiliary_list) :
                            # if factor l has causal prefactors within auxiliary_list
                            
                            # and if it is also a causal prefactor of another element in auxiliary list,
                            # it should be removed
                            for lfac in auxiliary_list :
                                if auxiliary_list[l] in get_causal_prefactors(lfac, level_equiv_list[lvl], auxiliary_list) :
                                    del auxiliary_list[l]
                                    break      # break from loop over the elements of auxiliary_list
                            

                # carry over the entries from auxiliary_list into return_list    
                for lfac in auxiliary_list :
                    entry = (lfac, fac)
                    if not(entry in return_list) :
                        return_list.append(entry) 
                
                # clear the auxiliary list    
                auxiliary_list.clear()
                                
                                
                
    #print("minimal list of constitution relations:")
    #for c in constitution_relation_list :
    #    print(c[0] + " -- " + c[1]) 
      
    return return_list
# end of find_structure                   


def print_structure_in_tikz_plot(level_factor_list_order, level_equiv_list, constitution_relation_list) :
    # Prepares the TikZ code for plotting one solution
    # a) places the causal factors as nodes separated by causal order (horizontally) and constitution level (vertically)
    # b) adds the causal and constitution relations as vertices between nodes
    
    ######################################
    # step 8: preparing the output files #
    ######################################
    
    tex_code = "% placement of the nodes\n"
    
    ######################################################
    # step 8a) - placement of the nodes = causal factors #
    ######################################################
    
    # [possible improvement 1]
    # - start with the level thas has the highest number of causal orders (highest horizontal length)
    # here: start with level 1
    
    placement = ""

    for m in range(len(level_factor_list_order)) :
        # add some tex-comments in order to increase the readability of the tex-code
        tex_code = tex_code  + "% factors of level " + str(m) + ":\n"
        
        max_num_factors_order = len(level_factor_list_order[m][0])
        for o in range(len(level_factor_list_order[m])) :
            
            # placement of the factors in their causal order
            # factors of the same order are placed on top of each other
            
            # [possible improvement 2]
            # - group factors that appear in the same constitution relations
            # here: procede according to the sequence in factor_list
            
            tex_code = tex_code  + "% causal order " + str(o) + ":\n"
            for e in level_factor_list_order[m][o] :
                
                # add a line to tex_code in which the node is placed, its name is the same as the one of the causal factor
                # and it is displayed on a label
                tex_code = tex_code + "\\node" + placement + " (" + e + ") {$" + e +"$};\n"
                
                if o == 0 :
                    # highlight incoming factors
                    tex_code = tex_code + "\hilightsource{" + e + "};\n"
                
                # maybe outgoing factors = those that have no outgoing arrows instead of max order ??
                elif o == len(level_factor_list_order[m]) - 1 :
                    # highlight outgoing factors
                    tex_code = tex_code + "\hilighttarget{" + e + "};\n"
                
                # prepare the variable placement for the next factor
                placement = "[above= \LvDist of " + e + "]"
                # the next factor of the same level and causal order will be positioned above by \LvDist
                
            # factors of the subsequent order -> next factor will be placed to the right of the bottom factor of the current order
            
            placement = "[right= \LhDist of " + level_factor_list_order[m][o][0] + "]"
            
            if len(level_factor_list_order[m][o]) > max_num_factors_order :
                max_num_factors_order = len(level_factor_list_order[m][o])
                
        # factors of the subsequent level -> next factor will be shifted upwards by \iLvDist
        # more precisely: \iLvDist  + height of the current level (= max_num_factors_order * \LvDist) above the first factor of this level
        placement = "[above= {" + str(max_num_factors_order) + "*\LvDist + " + str(max_num_factors_order - 2) + "* \HeightNode  + \iLvDist} of " + level_factor_list_order[m][0][0] + "]"
        
    ##########################################################
    # step 8 b) - plot the causal and constitution relations #
    ##########################################################
    tex_code = tex_code  + "\n% causal relations\n"
    for m in range(len(level_equiv_list)) : 
        tex_code = tex_code  + "% of level "  + str(m) + "\n"
        for formula in level_equiv_list[m] :
            tex_code = tex_code  + "% formula: "  + formula[0] + " <-> " + formula[1] + "\n"
            tex_code = tex_code + convert_causal_relation(formula, level_factor_list_order, tex_code) + "\n\n"
    
    tex_code = tex_code  + "\n% constitution relations\n"        
    for formula in constitution_relation_list :
        tex_code = tex_code  + "% formula: "  + formula[0] + " <-> " + formula[1] + "\n"
        tex_code = tex_code + convert_constitution_relation(formula, level_factor_list_order, constitution_relation_list) + "\n"  
    
    
    return tex_code

                   
def create_pdf(tex_code_table, latex_template_file) :
   # creates a pdf of tex_code_table
   # Latex_Template assumes that tex_code_table is a table of pairs of a number (used to reference the index of the element) and
   # a string that contains valid TikZ code.
   
    output_file = "output_graph.tex"
    # defining the template (already prepared file)
    template = latex_jinja_env.get_template(latex_template_file)
    render = template.render(data = tex_code_table, maxnumber = len(tex_code_table))
    
    # save the generarted string as tex file
    f = open(output_file, 'wb')
    f.write(render.encode('utf-8'))
    f.close()
    
    
    # compile this tex file with pdflatex -- requires the Latex compiler to be installed on the executing system
    with codecs.open(str(output_file), "w","utf-8") as letter:
        letter.write(render);
        letter.close();
        os.system("pdflatex -interaction=batchmode " + str(output_file))
        print('File ' + str(output_file[:-4]) + '.pdf created.')
    
    # remove automatically generated log files
    os.remove(str(output_file[:-4]+".log"))
    os.remove(str(output_file[:-4]+".aux"))
                   
def main() :
    # main function
    
    input_file = "r_output.txt"
    latex_template_file = "Latex_Template.tex"
    
    level_factor_list = []               # declaration of the lists
    level_factor_order_list = []
    level_equiv_list = []
    constitution_relation_list = []
    solution_list = []
    
    # steps 1 - 3 start function read_R_file -> converts cna output into lists that are sorted by constitution level
    # of causal factors (level_factor_list), causal relations (level_equiv_list) and one list of constitution relations
    # (constitution_relation_list), if the cna output is not as expeceted stop the procedure with abort = True
    if os.path.exists(input_file) :
        abort, level_factor_list, level_equiv_list, constitution_relation_list, solution_list = read_R_file(input_file)
    else :
        abort = True
        print("Error: Expected input data file " + input_file + " has not been found in path folder.")
        
    if not(abort) :
        # step 4 will sort out the correct solutions among those generated by cna
        solution_term_list = []
        solution_list = cull_complex_solutions(solution_list, level_factor_list, level_equiv_list, constitution_relation_list)
        
        if solution_list : # solution_list isn't empty after culling
            tex_table = []
            solution_term_list = [] # recreate a local solution_term_list for the reduced solution_list
    
            solution_term_list = create_term_list(solution_list)
            
            counter = 0
            
            
            if len(solution_term_list) > 100 :   # limiting the output to 100 solutions
                # remove this if-clause if more solutions are required
                print("More than 100 solutions obtained, plotting only the first 100.")
                for i in range(len(solution_term_list) - 1, 99, -1) :
                    del solution_term_list[i]
                        
            for sol in solution_term_list :

                # level_equiv_list has to be adapted to sol
                for i in range(len(level_factor_list)) :
                    level_equiv_list[i].clear() # clear its current elements
                    
                for term in sol :
                    formula = get_equiv_formula(term)
                    level_equiv_list[get_formula_level(formula[0], level_factor_list)].append(formula)

                
                # continue with steps 5 and 6 that determine the causal order of the factors
                abort, unique, circular, level_factor_list_order = determine_factor_order(level_factor_list, level_equiv_list)
                if not(abort) :
                    if unique :
                        # if the causal order of the factors is uniquely determinable:
                        # step 7 minimalisation of the constitution relations
                        new_constitution_relation_list = find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list)
                        
                        # step 8: graphical output as a graph in pdf
                        
                        # generating the tex-code
                        st = print_structure_in_tikz_plot(level_factor_list_order, level_equiv_list, new_constitution_relation_list)
                        
                        # subtitle of the graph will be the formula in tex-math syntax
                        subtitle = "\\tiny $"
                        for term in sol :
                            subtitle = subtitle + "(" + term.replace("*", " \cdot ").replace("~", "\\neg ").replace("<->", "\leftrightarrow ") + ")\cdot"
                        
                        subtitle = subtitle[:-5] + "$"
                        
                        if circular :
                            subtitle = "circular causal structure\\\\[4mm]" + subtitle
                        
                        # counter enumerates the solutions
                        counter = counter + 1
                        entry = (counter, st, subtitle)
                        tex_table.append(entry)
                        
                        
                    else :
                        # unique == False
                        # causal structure is non-unique -> no valid solution
                        
                        # print that this solution does not correspond with a single causal structure
                        st = ""
                        for term in sol :
                            st = st + "*(" + term + ")"
                        
                        st = st[1:] # remove first character (formula should not start with leading "*")
                        print("Solution " + st + " has a non-unique causal structure. It is discarded.")
             
            
            # after one entry for each solution has been generated in tex_table compile the pdf
            if tex_table :
                # if tex_table is non-empty
                if os.path.exists(latex_template_file) :
                    create_pdf(tex_table, latex_template_file)
                else :
                    abort = True
                    print("Error: Cannot plot the graphs since the expected template file " + latex_template_file + " does not exist in path folder.")
            
        else :
            # solution_list is empty
            # no solution survived selection of valid solutions
            print("No valid complex solution formula has been found in " + input_file + ".")   
            
                     
if __name__ == '__main__':
    # start main() when executing this script file
    main()
