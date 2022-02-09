#!/usr/bin/env python3

# file: causal_structure_R_data_to_pdf_graph.py

# proceeds in eight steps


import os                       # operating system interfaces is required to find the files of the own path
import codecs                   # for en- and decoding of strings (esp. to get tex-files in utf-8)
import re                       # regex for complex search patterns in strings
import jinja2                   # Latex interface
from jinja2 import Template

from search_formula import search_formula # own function that finds equivalent logical expressions relative to some formulae list

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
    
    
def get_equiv_formula(st):
    # returns the leftside and rightside partial formulae of the "alpha <-> beta" from the input string st
    # as pairs of strings (a,b)
    a = re.split(" <-> ",st)[0].strip()          # strip() removes leading spaces
    b = re.split(" <-> ",st)[1].strip()
    
    # conversion of the negation syntax (in cna by minuscle) such that "a" -> "~A"
    # 1. step: add "~" before each minuscle, which is either
    # a) at the beginning of a formula
    # b) follows a conjunctor
    # c) follows a disjunctor
    a = re.sub(r'^([a-z])',  r'~\1', a)
    # explanation:  "sub" replaces each instance of a minuscle (expressed by "[a-z]")
    # by itself plus the prefix "~",
    # if it has been found at the first position of the string (implicated by "^")

    # b) if following a "*", the letter will be placed behind "*~"
    a = re.sub(r'\*([a-z])',  r'*~\1', a)
    # The regex expression "\*" picks the star symbol "*" from the string.

    # c) if following " + ", the letter will be placed behind "+ ~"
    a = re.sub(r'\s\+\s([a-z]+)',  r' + ~\1', a)
    # in regex "\s" corresponds to spaces, "\+" to "+"

    # 2. step replacement of the minuscle that follow to "~" by majuscle
    a = re.sub(r'(~[a-z]+)', lambda pat: pat.group(1).upper(), a)
    
    
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
    # list[LEVEL][ORDER][FACTORS], if factor is contained in FACTORS, return the corresponding value of ORDER
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


def convert_causal_relation(formula, level_factor_list_order, tex_code):
    # translates a formula of causal relations into TikZ-Latex code
    # returns the code as string
    
    st = ""  # output code
    
    ###########################################
    # determine the syntactic type of formula #
    ###########################################
    # assumption: the cna output contains only the following types of formula:
    # 1) equivalence between causal factors (A <-> B) or non-equivalence between causal factors (~A <-> B)
    # 2) equivalence of one causal factor with conjunctions of factors ( A*B* ... <-> E)
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
            
            disjunctor_list = re.split("\s*\+\s*", formula[0])  # list of the possibly complex disjuncts of formula
            
            for disj in disjunctor_list : 

                
                if disj in get_components_from_formula(formula[0], level_factor_list_order) :
                    # case A: the discunct is one causal factor
                    st = st + "% simple disjunction\n"
                    st = st + "\draw[->] (" + disj + ".east) to (" + formula[1] + ".west);\n"
                
                elif formula[0].find("*") > -1 :
                    # case B: the disjunct is a conjunction
                    st = st + "% complex disjunction\n"
                    
                    conjunctor_list = re.split("\*", disj)
                    
                    # first set the junction of the conjuncts
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
                    position = "at ([xshift=\hDisjConj, yshift=\\vDisjConj]" + f_fac + ".east)"
                    # Attention: It might happen that several disjuncts of conjuncts meet at the same factor f_fac,
                    # therefore we have to check whether the position of the junction node has to be shifted.
                    q = 1
                    while (tex_code.find(position) > -1) or (st.find(position) > -1) :
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
                            st = st + "\draw[conjunctonsegment] (" + conj + ") to (" + cross_point + "aux);\n"
                            
                        else :
                            # Is there anything else that might happen??
                            print("The disjunction " + formula[0] + "  could not be plotted because its substructure was not recognized(1).")

                    # connect the junction with the target factor
                    st = st + "% arrow from junction to target factor\n\draw[->] (" + cross_point + "aux) -- (" + formula[1] + ".west);\n"
                
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
                    
                    st = st + "\draw[conjunctonsegment] (" + conj + ") to (" + formula[1] + "aux);\n"
            
            # draw arrow from junction to target factor
            st = st + "% arrow from junction to target factor\n\draw[->] (" + formula[1] + "aux) -- (" + formula[1] + ");"
        else :
            
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
                solution_list.append(re.sub("\r?\n","",line).rstrip().strip())
                
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
                    factor_list.remove(factor_list[k])  
                    
            
            
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
                    # alle Formeln, bei denen rechts ein Element aus einer anderen Ebene mit Elementen links (alle aus gleicher Ebene)
                    # verbunden wird, werde zu constitution_relation_list hinzugefuegt
                    # assumption: only constitution relations with level difference of one are maintained
                    
                    if get_formula_level(formula[0], level_factor_list) == get_formula_level(formula[1], level_factor_list) - 1 :
                        constitution_relation_list.append(formula)
                        
                           
                # all further formulae will not be considered any longer
    return abort, level_factor_list, level_equiv_list, constitution_relation_list, solution_list
    
    
    
def determine_factor_order(level_factor_list,level_equiv_list) :    
    # description
    
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
        # step 4: determine the incoming factors #
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
                        
                
            ###########################################
            # step 5: define causal order iteratively #
            ###########################################
                    
            # dazu zunaechst Zerlegung von level_factor_list in Faktoren, die nur rechts in den verbleibenden
            # Formeln stehen -> downstream_factor_list
            downstream_factor_list = []
                
            # downstream_factor_list besteht zunaechst aus den Elementen von level_factor_list[m], die nicht
            # in von Ordnung 0 sind
            for element in level_factor_list[m]:
                if not (element in level_factor_list_order[m][0]):
                    downstream_factor_list.append(element)
                        
           
                    
            # jetzt wird level_factor_list_order[m] sukzessiv ergaenzt um Faktoren, die 
            # nur von Faktoren abhaengen, die bereits in level_factor_list_order[m] stehen
            # dazu wird downstream_factor_list von hinten durchgegangen
            for j in range(len(downstream_factor_list)-1,-1,-1):
                stop = False
                k = 0
                candidate = False
                # k verweist auf den aktuellen Kandidaten aus downstream_factor_list
                while not(candidate) and k < len(downstream_factor_list):
                    i = 0
                    candidate = True
                    order = 0           # Ordnung des betrachteten Faktors

                    # es werden alle Formeln dieser Ebene durchsucht, wobei
                    # i auf aktuelle Formel verweist
                    while candidate and i < len(level_equiv_list[m]):
                        # steht der aktuelle Kandidat rechts in der Formel?
                        if downstream_factor_list[k] == level_equiv_list[m][i][1]:
                            # dann wird geprueft ob alle Faktoren links Elemente von 
                            # level_factor_list_order[m] sind, dabei:
                            # level_equiv_list[m][i][0] -> m = Level der Formel
                            # i = Nummer des Eintrags, 0 = linke Seite der Formel
                            for element in get_components_from_formula(level_equiv_list[m][i][0], level_factor_list):
                                candidate = False
                                for o in range(len(level_factor_list_order[m])):     # o laeuft ueber die Anzahl der Ordnungen
                                    # in level_factor_list_order
                                    if ( element in level_factor_list_order[m][o]):
                                        candidate = True
                                        if order < o + 1:
                                            order = o + 1
                                
                        i = i + 1
                                
                    k = k + 1
                if len(level_factor_list_order[m]) == order :
                    level_factor_list_order[m].append([])  # neue Ordnung zur Liste hinzufuegen
                            
                # falls in allen Formeln des Levels, alle Faktoren auf der linken Seite
                # level_factor_list_order[m] sind, wird downstream_factor_list[k-1] in einen solchen
                # umgewandelt mit der entsprechenden Ordnung (k-1, weil zwischenzeitlich k um eins erhoeht wurde)    
                if candidate :
                    level_factor_list_order[m][order].append(downstream_factor_list[k-1])                    
                            
                    downstream_factor_list.remove(downstream_factor_list[k-1])
                else:
                    stop = True
            
            # output as console text     
            print("causal factors of level " + str(m) + " with their causal order:")
            for p in range(len(level_factor_list_order[m])):
                print("order " + str(p))
                print(level_factor_list_order[m][p])
            
                       
            if downstream_factor_list:
                # falls downstream_factor_list nicht leer ist, ist Kausalstruktur nicht eindeutig
                print("No unique structure")
                unique = False
               
    return abort, unique, level_factor_list_order
    
    
    
    
    
    
def find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list) :
    # bestimmt Kausal- und Konstitutionsstruktur aus den eingegebenen Listen von Kausalfaktoren, Kausalrelationen und 
    # Konstitutionsbeziehungen
    # Ausgabe erfolgt in Form von (moeglichst) minimalen Kausal- und Konstitutionsrelationen
    
    level_count = len(level_factor_list_order)
    
    for m in range(level_count) :
        #########################################################
        # Schritt 6: falsche Konstitutionsbeziehungen verwerfen #
        #########################################################
        # vor den Kausalrelationen durchfuehren um volle Zahl (auch redundante) an Aequivalenzformeln
        # fuer Umformungen nutzen zu koennen
                    
        if m > 0 :
            # erst ab 1. Level gibt es Konstitutionsbeziehungen
            for o in range(len(level_factor_list_order[m])) :
                for fac in level_factor_list_order[m][o] :
                    # zaehle Konstitutionsformeln, die fac mit unterer Ebene verbinden
                    # bei mehr als zwei (eine fuer incomming, andere outgoing), muessen falsche Formeln
                    # geloescht werden
                    auxiliary_list = [] # Hilfsliste fuer Konstitutionsbeziehungen bzgl. fac
                    for c in constitution_relation_list :
                        if c[1] == fac :
                            auxiliary_list.append(c)
                            
                    if len(auxiliary_list) < 1 :
                        # Meldung das weniger Konstitutionsbeziehungen als erwartet gefunden worden sind
                        print("Für " + fac + " wurden nur " + str(len(auxiliary_list)) + " Konstitutionsbeziehungen gefunden.")
                                
                    elif len(auxiliary_list) > 2 :
                        # loesche ueberzaehlige Konstitutionsbeziehungen
                        discard_list = [] # Liste der zu verwerfenden Formeln
                                
                        if o == 0 :
                            # wenn Ordnung des betrachteten Faktors = 0, dann koennen falsche Formeln nur aus Zwischenfaktoren
                            # bestehen, d.h. incomming Beziehung muss aus Faktoren der Ordnung 0 (bzgl. der unteren Ebene) bestehen,
                            # outgoing Beziehung besteht aus Faktor hoechster Ordnung (bzgl. der unteren Ebene)
                            max_order = 0
                            for ac in auxiliary_list :
                                if get_formula_order(ac[0], level_factor_list_order) > max_order :
                                    max_order = get_formula_order(ac[0], level_factor_list_order)

                            # alle Formeln zwischen Ordnung = 0 und max_order koennen verworfen werden
                            for ac in auxiliary_list :
                                if (get_formula_order(ac[0], level_factor_list_order) < max_order) and (get_formula_order(ac[0], level_factor_list_order) > 0) :   
                                    discard_list.append(ac)
                                            
                            # pruefen, ob es nun nur noch zwei uebrige Formeln sind
                            if len(auxiliary_list) - len(discard_list) == 2 :
                                # falls ja, loesche Elemente von discard_list aus constitution_relation_list
                                for e in discard_list :
                                    constitution_relation_list.remove(e)
                                            
                            elif len(auxiliary_list) - len(discard_list) < 2 :
                                print("Für den Faktor " + fac + " gibt es zu wenig Konstitutionsbeziehungen.")
                                # !!! Was tun ? !!!
                                        
                            elif len(auxiliary_list) - len(discard_list) > 2 :
                                # loesche zunaechst klar falsche Formeln
                                for e in discard_list :
                                    constitution_relation_list.remove(e)
                                # Meldung dass immer noch zu viele Formeln vorliegen
                                print("Für den Faktor " + fac + " gibt es zu viele Konstitutionsbeziehungen.")
                                # !!! Was tun ? !!!
                                        
                        else : # Term der hoeheren Ebene ist von groesserer Kausalordnung als 0
                        
                            sec_auxiliary_list = [] # zweite Hilfsliste, hier werden alle bekannten "Uebersetzungen" von
                            # Faktoren der Ebene m durch ihre Konstituenten der Ebene m-1 eingetragen
                            
                                        
                            for j in level_equiv_list[m] :
                                if j[1] == fac :                                                
                                    # suche unter den Aequivalenzen der Ebene m nach solchen, die fac kausal aufloesen
                                    for k in get_components_from_formula(j[0],level_factor_list_order[m]) :
                                        # fuer jeden kausalen Faktor bzgl. fac pruefe, sammle dessen Konstitutionsbeziehungen
                                        # in sec_auxiliary_list
                                        for c in constitution_relation_list :
                                            if c[1] == k :
                                                sec_auxiliary_list.append(c)
                                                            
                                                            
                                                            
                                                           

                            sec_auxiliary_list.extend(auxiliary_list) # alle Listen werden fuer search_formula benoetigt
                            sec_auxiliary_list.extend(level_equiv_list[m]) # dort werden sie zum Teil wieder separiert
                            
                            for ac in auxiliary_list :
                                if search_formula(ac, sec_auxiliary_list) :
                                    print("redundant:")
                                    print(ac)
                                    discard_list.append(ac)
                                else : #hierhier
                                    print("NICHT redundant")
                                    print(ac)
                             
                                    
                            # loesche Duplikate    
                            discard_list = list(dict.fromkeys(discard_list))
                                        
                            # pruefen, ob es nun nur noch zwei uebrige Formeln sind
                            if (len(auxiliary_list) - len(discard_list)) == 2 :
                                # falls ja, loesche Elemente von discard_list aus constitution_relation_list
                                    for e in discard_list :
                                        constitution_relation_list.remove(e)
                            elif (len(auxiliary_list) - len(discard_list)) == 1 :
                                for e in discard_list :
                                    constitution_relation_list.remove(e)
                                    # !!! Ist das so erwuenscht ? !!!
                                                
                            elif (len(auxiliary_list) - len(discard_list)) < 1 :
                                print("Für den Faktor " + fac + " gibt es zu wenig Konstitutionsbeziehungen.")
                                # !!! Was tun ? !!!
                                        
                            elif (len(auxiliary_list) - len(discard_list)) > 2 :
                                # loesche zunaechst klar falsche Formeln
                                for e in discard_list :
                                    constitution_relation_list.remove(e)

                                # Meldung dass immer noch zu viele Formeln vorliegen
                                print("Für den Faktor " + fac + " gibt es zu viele Konstitutionsbeziehungen.")
                                    
                            sec_auxiliary_list.clear()
                                        
                                    
                                        
                                    
                        discard_list.clear()
                                
                          
                                
                                
                auxiliary_list.clear()
                                
                                
                
                    
        #################################################
        # Schritt 7: redundante Kausalformeln streichen #
        #################################################
        # Kausalformeln, die durch Einsetzen von anderen Formeln in einander gebildet werden,
        # sollen verworfen werden  
                    
        #discard_list = [] # Liste der zuloeschenden Formeln
                    
        for j in range(len(level_equiv_list[m])-1,-1,-1):
            # fuer jede Formel j wird separat geprueft, ob sie redundant ist, d.h. es gibt eine Formel mit weniger
            # Kausalfaktoren (dann notwendigerweise hoeherer Ordnung), die dieselbe Relation ausdrueckt
            
            if get_factor_order(level_equiv_list[m][j][1],level_factor_list_order) > 1:
                # nur Formeln, mit Ordnung > 1 koennen redundant sein, da nur solche durch
                # Einsetzungen entstehen koennen
                
                
                print(level_equiv_list[m][j])            
                
                new_equiv_list = []
                
                for eq in level_equiv_list[m] : # nur Aequivalenzformeln gleicher oder niedrigerer Ordnung sind brauchbar fuer
                    # Suche nach Redundanzen
                    if get_factor_order(level_equiv_list[m][j][1],level_factor_list_order) >= get_factor_order(eq[1],level_factor_list_order) :  
                        new_equiv_list.append(eq)
                    
                if search_formula(level_equiv_list[m][j], new_equiv_list) : # verwendet Funktion aus search_formula.py
                    del level_equiv_list[m][j]  # loescht j-ten Eintrag von level_equiv_list[m]
                
                
        
                            
        print("minimale Kausalrelationen der Konsitutionsebene " + str(m) + ":")
        for formula in level_equiv_list[m] : 
            print(formula[0] + " -> " + formula[1])
                        
    
    # Schleife ueber Konstitutionsebenen endet hier
                
    print("Liste der minimalen Konstitutionsbeziehungen:")
    for c in constitution_relation_list :
        print(c[0] + " -- " + c[1]) 
      
    return level_equiv_list, constitution_relation_list
# Ende find_structure                   


def print_structure_in_tikz_plot(level_factor_list_order, level_equiv_list, constitution_relation_list) :
    # description to follow
    
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
        
        max_num_factors_order = 0
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
                tex_code = tex_code + "\\node" + placement + " (" + e + ") {" + e +"};\n"
                
                if o == 0 :
                    # highlight incoming factors
                    tex_code = tex_code + "\hilightsource{" + e + "};\n"
                elif o == len(level_factor_list_order[m]) - 1 :
                    # highlight outgoing factors
                    tex_code = tex_code + "\hilighttarget{" + e + "};\n"
                
                # prepare the variable placement for the next factor
                placement = "[above= \LvDist of " + e + "]"
                # the next factor of the same level and causal order will be positioned above by \LvDist
                
            # factors of the subsequent order -> next factor will be placed to the right of the bottom factor of the current order
            
            placement = "[right= \LhDist of " + level_factor_list_order[m][o][0] + "]"
            
            if len(level_factor_list_order[m]) > max_num_factors_order :
                max_num_factors_order = len(level_factor_list_order[m])
                
        # factors of the subsequent level -> next factor will be shifted upwards by \iLvDist
        # more precisely: \iLvDist  + height of the current level (= max_num_factors_order * \LvDist) above the first factor of this level
        placement = "[above= {" + str(max_num_factors_order) + "*\LvDist + " + str(max_num_factors_order) + "* \HeightNode  + \iLvDist} of " + level_factor_list_order[m][0][0] + "]"


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
    
    
    output_file = "output_graph.tex"
    # defining the template (already prepared file)
    template = latex_jinja_env.get_template('Latex_Template.tex')
    render = template.render(tex_formula = tex_code)
    
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
    
    
# end of print_structure_in_tikz_plot    
                   
def main() :
    # main function
    
    
    level_factor_list = []               # declaration of the lists
    level_factor_order_list = []
    level_equiv_list = []
    constitution_relation_list = []
    solution_list = []
    
    # steps 1 - 3 start function read_R_file -> converts cna output into lists that are sorted by constitution level
    # of causal factors (level_factor_list), causal relations (level_equiv_list) and one list of constitution relations
    # (constitution_relation_list), if the cna output is not as expeceted stop the procedure with abort = True
    abort, level_factor_list, level_equiv_list, constitution_relation_list, solution_list = read_R_file("r_output.txt")
    
    if not(abort) :
        # continue with steps 4 and 5 that determine the causal order of the factors
        abort, unique, level_factor_list_order = determine_factor_order(level_factor_list,level_equiv_list)
        if not(abort) :
            if unique :
                # if the causal order of the factors is uniquely determinable:
                # steps 6 and 7 minisation of the causal and constitution relations
                level_equiv_list, constitution_relation_list = find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list)
                # step 8: graphical output as a graph in pdf
                print_structure_in_tikz_plot(level_factor_list_order, level_equiv_list, constitution_relation_list)
                
            else :
                print("non-unique result")                     #[hierhier]
                level_equiv_list, constitution_relation_list = find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list)
          
                    
if __name__ == '__main__':
    main()
