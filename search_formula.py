#!/usr/bin/env python3

# provides a function to identify a redundant logical formula relative to a list of formulae
# Assumption: input formulae are of type of disjunctive normal forms (dnf)
# search_formula(f,list_f) transforms formula f until all logically legitimite transformations are 
# obtained or a formula from list_f is matches one of the transformations
#
# The performed transformations are:
# 1.) distributive law: A*B*C + A*B*D*E + B ... -> A*B*(C + D*E) + B ... / but also -> A*(B*C + B*D*E) + B ...
# 2.) commutations
# - a) of disjunctions A*B + B*C + E + ... -> yields the set of all permutations of the discuncts
# - b) of conjunctions A*B + B*C + E + ... -> yields the set of all permutations of the conjuncts within on discunct
# 3.) idempotences
# - a) of disjuncts A*B*C + A*B*C + ... -> A*B*C + ... - looks only at first and second disjunct
# - b) of conjuncts A*A*B + ... -> A*B + ... - looks only at first and second conjunct of the first disjunct
# 4.) Substitutions of formulae from list_f 
#     -- each formula is of the form alpha <-> beta, with beta being a causal factor and alpha a discunctive normal form
#        if the leftside expression is found in formula f or one of its transformations, then the rightside term beta is substituted
#        in its place
# 
# The six transformation operations are performed parallelised.
# Valid transformations are filled into a tree structure below the current node, if they are not already nodes of the tree.
# Dazu wird eine Abgleichsmenge gefuehrt, an der jede neue
# Umformung abgeglichen wird. For this purpose, a list of all tree nodes is kept, against which each new transformation is adjusted.
# The tree structure is traversed in a deep search, which is terminated by finding one of the formulas from list_f  is terminated,
# yielding the output: True.
# Each node is traversed only once, so that the search process is also ended after all nodes have been traversed, yielding: False.


from tree import Tree     # own file tree.py - offers a simple tree structure with a deep search function
import itertools          # itertools provides functions to obtain all permutations of a string and Cartesian products of lists
import re                 # regex for complex search patterns in strings
import multiprocessing    # multiprocessing and functools for multicore usage
import functools          


# rules of formula transformation
def de_morgan_disj(formula) :
    # De Morgan rule for disjunctions  ~alpha*~beta -> ~(alpha + beta)
    # Depending on whether alpha and beta are complex expressions or individual causal factors, three cases must be distinguished
    # 1) both are complex expressions
    # 2) both are individual causal factors
    # 3) alpha is a complex expression, beta a causal factor (We do not have to care about the reversed case, since we have
    # the commutation operation.)
    # For each type we have to discern the cases of ~alpha*~beta; ~~alpha*~beta; ~alpha*~~beta; ~~alpha*~~beta
    
    return_list = [] # the list of obtained equivalences will be returned
    
    # first case: complex replacements, the de Morgan terms are in brackets
    if (formula.count("~(") > 1) and (formula.count(")*~(") > 0) :
        # is only executed if formula contains at least two negated brackets,
        # at least one of which is followed by a conjunction and a negated bracket,
        
        # now check that the first opening bracket "~(" belongs to the closing one in ")*~(" 
        if formula.find(")*~(") == formula.replace("(","r",1).find("~(") - 2 :
            if (formula.find("~(") == 0) or (formula[formula.find("~(")-1] != "~") :
                #1.A there is only one negator infront of the first bracket
    
                # If this is the case, replace "~(alpha)*~(beta)" by  "~(alpha + beta)".
                string = re.sub(r'~\(([a-zA-Z\+~\*\s]+)\)\*~\(([a-zA-Z\+~\*\s]+)\)', r'~(\1 + \2)', formula)
                # explanation: looks for a partial expression in formula, which starts with "~(",
                # is followed by an arbitrary quantity (but at least one) of letters, "*", "+" and " ", 
                # this sequence will be saved as \1, which is indicated by "( ... )"
                # thereupon must follow ")*~(", which is followed by another sequence of letters, "*", "+" and " ",
                # which will be stored as \2,
                # finally the expression has to end with ")"
                # If such an expression has been found in formula, replace it by "~(" + \1 + " + " + \2 + ")".
                # ATTENTION: The brackets must not contain subordinate brackets, otherwise the groups \1 and \2 do not capture.
                # This shouldn't be a problem(??)
        
            else :
                #1.B We have the form ~~(alpha)*~(beta) -> ~(~(alpha) + beta)
                string = re.sub(r'~~\(([a-zA-Z\+~\*\s]+)\)\*~\(([a-zA-Z\+~\*\s]+)\)', r'~(~(\1) + \2)', formula)

            return_list.append(string)     # add the string to the list of results
            
    if (formula.count("~(") > 1) and (formula.count(")*~~(") > 0) :
        # 1.C or 1.D we have ~(alpha)*~~(beta) or ~~(alpha)*~~(beta)
        
        # now check that the first opening bracket "~(" belongs to the closing one in ")*~~(" 
        if formula.find(")*~~(") == formula.replace("(","r",1).find("~~(") - 2 :
            if (formula.find("~(") == 0) or (formula[formula.find("~(")-1] != "~") :
                #1.C there is only one negator infront of the first bracket
    
                # If this is the case, replace "~(alpha)*~~(beta)" by  "~(alpha + ~(beta))".
                string = re.sub(r'~\(([a-zA-Z\+~\*\s]+)\)\*~~\(([a-zA-Z\+~\*\s]+)\)', r'~(\1 + ~(\2))', formula)
                
        
            else :
                #1.D We have the form ~~(alpha)*~~(beta) -> ~(~(alpha) + ~(beta))
                string = re.sub(r'~~\(([a-zA-Z\+~\*\s]+)\)\*~~\(([a-zA-Z\+~\*\s]+)\)', r'~(~(\1) + ~(\2))', formula)

            return_list.append(string)     # add the string to the list of results    
    
    
    # second case: the de Morgan terms are negated causal factors
    if (formula.count("~") > 1) and (formula.count("*") > 0) :
        # There are at least two negators and one conjunctor.

        match = re.search(r'~([a-zA-Z]+)\*~([a-zA-Z]+)', formula)
        match_2  = re.search(r'~([a-zA-Z]+)\*~~([a-zA-Z]+)', formula)
        
        if match and ((match.start() == 0) or (formula[match.start() -1] != "~")) :
            # 2.A the simplest case ~A*~B -> ~(A + B)
            string = re.sub(r'~([a-zA-Z]+)\*~([a-zA-Z]+)', r'~(\1 + \2)', formula)
            # explanation: looks for a partial expression in formula, which starts with "~",
            # is followed by an arbitrary quantity (but at least one) of letters, 
            # this sequence will be saved as \1, which is indicated by "( ... )"
            # thereupon must follow "*~", which is followed by another sequence of letters, that will be saved as \2
            # If such an expression has been found in formula, replace it by "~(" + \1 + " + " + \2 + ")".
            
            return_list.append(string)  # add the string to the list of results
        
        elif match and ((match.start() > 0) and (formula[match.start() -1] == "~")) :
            # 2.B the first conjunct has two negators: ~~A*~B -> ~(~A + B)
            string = re.sub(r'~~([a-zA-Z]+)\*~([a-zA-Z]+)', r'~(~\1 + \2)', formula)
            return_list.append(string)  # add the string to the list of results
        
        elif match and ((match_2.start() == 0) or (formula[match_2.start() -1] != "~")) :
            # 2.C only the second conjunct has to negators ~A*~~B -> ~(A + ~B)
            string = re.sub(r'~([a-zA-Z]+)\*~~([a-zA-Z]+)', r'~(\1 + ~\2)', formula)
            return_list.append(string)  # add the string to the list of results
            
        elif match and ((match_2.start() > 0) and (formula[match_2.start() -1] == "~")) :
            # 2.D both conjuncts have two negators: ~~A*~~B -> ~(~A + ~B)
            string = re.sub(r'~~([a-zA-Z]+)\*~~([a-zA-Z]+)', r'~(~\1 + ~\2)', formula)
            return_list.append(string)  # add the string to the list of results

    # third case: the first de Morgan term is a complex expression, the second is a negated causal factor
    if (formula.count("~(") > 0) and (formula.count(")*~") > 0) :
        # formula contains at least one negation of a bracket and the combination of a closing a bracket with
        # a conjunction of a negated expression
        
        if (formula.count("(") == 1) or ((formula.find(")*~")) < formula.replace("(","r",1).find("(") - 3) :
            # if ")*~" isn't followed by another bracket (which should have been handled by case one), ...
            
            # again starting with the most simple case - both conjuncts have only one negator
            if (formula.find("~(") == 0) or (formula[formula.find("~(")-1] != "~") :
                # ... then replace the expression "~(alpha)*~A" by  "~(alpha + A)" (case 3.A)
                string = re.sub(r'~\(([a-zA-Z\+~\*\s]+)\)\*~([a-zA-Z]+)', r'~(\1 + \2)', formula)
                # explanation: looks for a partial expression in formula, which starts with "~(",
                # is followed by an arbitrary quantity (but at least one) of letters, "*", "+" and " ", 
                # this sequence will be saved as \1, which is indicated by "( ... )"
                # thereupon must follow ")*~", which is followed by a sequence of letters, which will be saved as \2
                # If such an expression has been found in formula, replace it by "~(" + \1 + " + " + \2 + ")".
                # ATTENTION: The brackets must not contain subordinate brackets, otherwise the groups \1 and \2 do not capture.
                # This shouldn't be a problem(??)
                return_list.append(string)     
                
                
                # case 3.C ~alpha*~~B -> ~(alpha + ~B)
                string_2 = re.sub(r'~\(([a-zA-Z\+~\*\s]+)\)\*~~([a-zA-Z]+)', r'~(\1 + ~\2)', formula)
                return_list.append(string_2)     
           
            elif formula[formula.find("~(") - 1] == "~" :
                # case 3.B ~~alpha*~B -> ~(~(alpha) + B)
                string = re.sub(r'~~\(([a-zA-Z\+~\*\s]+)\)\*~~([a-zA-Z]+)', r'~(~(\1) + \2)', formula)
                return_list.append(string)     
               
                # case 3.D ~~alpha*~~B -> ~(~(alpha) + ~B)
                string_2 = re.sub(r'~~\(([a-zA-Z\+~\*\s]+)\)\*~~([a-zA-Z]+)', r'~(~(\1) + ~\2)', formula)
                return_list.append(string_2)

    if return_list : return return_list # If return_list isn't empty, return it,
    else : return ""                    # otherwise return "".

def commutation_conj(formula) :
    # commutation of conjuncts
    
    q = 0
    start_formula = formula
    
    if (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
        # If the first conjunction is in brackets, only the bracket term is transformed.
        
        while (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
            
            # determining the string inside the brackets
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            bracket_term = ""
            for i in range(left_bracket_pos + 1, right_bracket_pos) :
                bracket_term = bracket_term + formula[i]
        
            formula = bracket_term
    
    elif formula.count("(") > 0 :
        # other bracket expressions are substituted by pseudo-literals during the commutation operation
        
        while formula.find("(") > -1 :
            # put the expression in parenthesis into a dictionary -> dictionary identifier functions as replacement in the formula
            # 1. step: get "( ... )" expression and save it in the dictionary
        
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            if (left_bracket_pos == 0) or (formula[left_bracket_pos - 1] != "~") :
                # no negation before the bracket
            
                # .join() is faster than a for-loop over all characters
                slist = [formula[i]  for i in range(left_bracket_pos, right_bracket_pos + 1)]
       
            else :
                # There is a negation right before the bracket.
                # -> The replacement has to capture the negation symbol also.
                if (left_bracket_pos >= 2) and (formula.find("~~") == left_bracket_pos -2) :
                    # double negation
                    # .join() is faster than a for-loop over all characters
                    slist = [formula[i]  for i in range(left_bracket_pos - 2, right_bracket_pos + 1)]
                
                else :
                    # single negation
                    # .join() is faster than a for-loop over all characters
                    slist = [formula[i]  for i in range(left_bracket_pos - 1, right_bracket_pos + 1)]
            
            bracket_term = "".join(slist)            
        
            formula = formula.replace(bracket_term, "bracket_term" + str(q), 1)
            
            if q == 0 : # create dictionary
                dictionary = {str(q) : bracket_term}
            else :     
                dictionary[str(q)] = bracket_term

            q = q + 1
    
    if formula.find("*") > -1 :
        # procedure: 1) split the formula into disjuncts
        # 2) form the possible permutations of the conjuncts for each disjunct separately
        # 3) combine each combination of disjuncts with each other
        
        return_list = []  # output - yields all *-commutative formulae of the input formula
        
        disj_list = re.split("\s*\+\s*", formula) # list of disjuncts in formula - step 1) is done
        
        
        # start with step 2):    
        partial_permutation_list = [] # partial_permutation_list will contain for each disjunct, one list of all permutations of its
        # conjuncts
        
        for disj in disj_list :
            if disj.find("*") > -1 :
                # creating that list of permutations
                # list(itertools.permutations(...) -> lists all permuations of the conjuncts of the disjunct disj
                partial_permutation_list.append(list(itertools.permutations(re.split("\*", disj))))
                
            else :    
                # in case that the discunct disj does not contain any conjunction,
                # record disj itself into the list of permutations
                partial_permutation_list.append(list(disj))
            
            
        auxiliary_list = [] # merges related conjuncts back into one formula (eg. ("A","~B") -> "A*~B")
        for i in range(len(partial_permutation_list)) :                
            auxiliary_list.append([])
            for term in partial_permutation_list[i] :
                
                st = ""
                
                for konj in term :
                    if st == "" :
                        st = st + konj
                    else :
                        st = st  + "*" + konj

                if i < len(partial_permutation_list) - 1 :
                    st = st + " + "
                auxiliary_list[i].append(st)


        # step 2) is done
        
        # start with step 3): form the Cartesian products from the variants of the disjuncts
        for item in itertools.product(*auxiliary_list):
            # .join() is faster than a for-loop over all characters
            #st = ""
            slist = [t for t in item]
            st = "".join(slist)
            
            
            if q > 0 :     # if we have substituted terms in brackets from the original formula
                # its placeholders "bracket_term1", ... have to be rereplaced pursuant to the dictionary entries
                for i in range(q) :
                    bracket_term = dictionary[str(i)]
                    st = st.replace("bracket_term" + str(i), bracket_term, 1) 
            
            elif formula != start_formula :
                # if only a partial formula was evaluated, st must be transferred to start_formula
                st = start_formula.replace(formula,st,1)
            
            return_list.append(st)

        auxiliary_list.clear()
        partial_permutation_list.clear()
            
        return return_list
    else : return ""

def commutation_disj(formula) :
    # commutation of disjuncts  [hierhier]
    q = 0
    start_formula = formula
    
    if (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
        # falls sich die erste Disjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        while (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos + 1, right_bracket_pos)]
            formula = "".join(slist)
        
            
            
    elif formula.count("(") > 0 :
        # andere Klammerausdruecke werden zur Bearbeitung durch Pseudo-Literale substituiert
        
        while formula.find("(") > -1 :
            # trage Klammerausdruck in ein dictionary ein -> dictionary Bezeichner als Ersatz in der Formel
            # 1. Schritt "( ... )" auslesen und in dictionary speichern 
        
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            if (left_bracket_pos == 0) or (formula[left_bracket_pos - 1] != "~") :
                # keine Negation vor der Klammer
            
                # .join() ist schneller als Schleife ueber Buchstaben
                slist = [formula[i]  for i in range(left_bracket_pos, right_bracket_pos + 1)]
       
            else :
                # There is a negation right before the bracket.
                # -> The replacement has to capture the negator also.
                if (left_bracket_pos >= 2) and (formula.find("~~") == left_bracket_pos -2) :
                    # double negation
                    # .join() is faster than a for-loop over all characters
                    slist = [formula[i]  for i in range(left_bracket_pos - 2, right_bracket_pos + 1)]
                
                else :
                    # single negation
                    # .join() is faster than a for-loop over all characters
                    slist = [formula[i]  for i in range(left_bracket_pos - 1, right_bracket_pos + 1)]
            
            bracket_term = "".join(slist)
            
            formula = formula.replace(bracket_term, "bracket_term" + str(q), 1)
            
            if q == 0 : # lege dictionary an
                dictionary = {str(q) : bracket_term}
            else :     
                dictionary[str(q)] = bracket_term

            q = q + 1
    
    
    
    if formula.find(" + ") > -1 :
        # Vorgehen: bestimme alle Disjunkten, bilde alle mgl. Permutationen dieser Disjunkten
        
        # Liste aller mgl. Permutationen von Disjunkten aus Formeln von eq_list
        permutation_list = list(itertools.permutations(re.split("\s*\+\s*", formula)))
        
        return_list = []
        
        # zusammensetzen der Permutationen
        for permut in permutation_list :
            st = ""
            for disj in permut :
                if st == "" :
                    st = st + disj
                else :
                    st = st  + " + " + disj
            
            
            if q > 0 :     # aus Ursprungsformel wurden Klammern substituiert,
                # also Platzhalter "bracket_term1", ... wieder gemaess dictionary ersetzen
                for i in range(q) :
                    bracket_term = dictionary[str(i)]
                    st = st.replace("bracket_term" + str(i), bracket_term, 1) 
            
            elif formula != start_formula :
                # falls nur eine Teilformel ausgewertet wurde, muss st auf start_formula uebertragen werden
                st = start_formula.replace(formula,st,1)    
                
            
            # eintragen in die Ausgabeliste
            return_list.append(st)
        
        permutation_list.clear()
        
        return return_list
    else : return ""

def distribution(formula) :
    # Distributivitaet
     
    q = 0
    while formula.find("(") > -1 :
        # innerhalb von Klammern koennen keine neuen Formeln durch Distribution gefunden werden
        # daher koennen Klammerausdruecke einstweilen durch Pseudo-Literale ersetzt werden
        # trage Klammerausdruck in ein dictionary ein -> dictionary Bezeichner als Ersatz in der Formel
        # 1. Schritt "( ... )" auslesen und in dictionary speichern 
        
        left_bracket_pos = formula.find("(")
        right_bracket_pos = formula.find(")")
        
        if (left_bracket_pos == 0) or (formula[left_bracket_pos - 1] != "~") :
            # keine Negation vor der Klammer
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos, right_bracket_pos + 1)]
            
       
        else :
            # There is a negation right before the bracket.
            # -> The replacement has to capture the negation symbol also.
            if (left_bracket_pos > 2) and (formula.find("~~") == left_bracket_pos -2) :
                # double negation
                # .join() is faster than a for-loop over all characters
                slist = [formula[i]  for i in range(left_bracket_pos - 2, right_bracket_pos + 1)]
                
            else :
                # single negation
                # .join() is faster than a for-loop over all characters
                slist = [formula[i]  for i in range(left_bracket_pos - 1, right_bracket_pos + 1)]
            
        
        bracket_term = "".join(slist)
        
        formula = formula.replace(bracket_term, "bracket_term" + str(q), 1)
        
        if q == 0 :  # lege dictionary an
            dictionary = {str(q) : bracket_term}
        else :     
            dictionary[str(q)] = bracket_term

        q = q + 1
    
    # Formel muss mindestens einen Disjunktor, zwei Konjunktoren enthalten
    if (formula.find(" + ") > -1) and (formula.count("*") > 1) :
                
        # Vorgehen: prueft ob erste n Konjunkte der ersten m Disjunkten identisch sind, wobei n < Anzahl der Konjunkte jeder dieser
        # Disjunkten gelten muss
        # falls ja: Klammer die n Konjunkten aus (A*B*C + A*B*D*E + B ... -> A*B*(C + D*E) + B ...)
        # andere Distributivitaeten koennen ignoriert werden, da sie durch Kommutationen gefunden werden
        
        return_list = []
        
        disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
        counter = 0
        max_disj = len(disj_list)
        conj_list = []
        while counter < len(disj_list) and (max_disj == len(disj_list)) :
            if disj_list[counter].find("*") > -1 :
                conj_list.append(re.split("\*", disj_list[counter]))
            else :
                max_disj = counter
                conj_list.append([])
            counter = counter + 1
        
        # conj_list ist nun eine Liste der Form [[c11, c12, ...], [c21, c22, ... ],  ...] wobei cij die jte Konjunkte der iten 
        # Disjunkten ist
        
        
        # nun: ruecklaufend was ist das groesste n (=Anzahl an gemeinsamen, ersten Konjunkten) der groesstmoeglichen Anzahl m
        # an ersten Disjunkten
        conj_counter = len(conj_list[0]) - 1
        while conj_counter > 0 :
            disj_counter = max_disj
            while disj_counter > 1 :
                condition = True
                for i in range(disj_counter) :
                    if conj_counter < len(conj_list[i]) :  # besitzt die Disjunktion ueberhaupt genuegend viele Konjunkte?
                        for j in range(conj_counter) :
                            if conj_list[0][j] != conj_list[i][j] :
                                condition = False
                    else : condition = False            
                        
                if condition :
                    # stimmen die ersten n Konjunkten der ersten m Disjunkten ueberein -> fuege Distributionsaequivalente zur 
                    # Ausgabeliste hinzu
                    
                    # konjunktiver Vorfaktor
                    
                    # .join() ist schneller als Schleife ueber Buchstaben
                    slist = [conj_list[0][i] + "*(" if i == conj_counter -1 else conj_list[0][i] + "*" for i in range(conj_counter)]
                    st = "".join(slist)
                    
                              
                    # disjunktiver Klammerausdruck
                    for i in range(disj_counter) :
                        for j in range(conj_counter,len(conj_list[i])) :
                            st = st + conj_list[i][j]
                            if j < len(conj_list[i]) - 1:
                                st = st + "*"
                            
                        if i < disj_counter - 1 :
                            st = st + " + "    
                        else :
                            st = st + ")"
                            
                    # weiteren Disjunkte
                    for i in range(disj_counter,len(disj_list)) :
                        st = st + " + " + disj_list[i]
                    
                                
                    if q > 0 : # Ursprungsformel enthielt keine Klammern, daher 
                        for i in range(q) : # Platzhalter "bracket_term1", ... wieder gemaess dictionary ersetzen
                            bracket_term = dictionary[str(i)]
                            st = st.replace("bracket_term" + str(i), bracket_term, 1) 
                        
                    
                    return_list.append(st)
                
                # versuche es mit einer Disjunkten weniger              
                disj_counter = disj_counter - 1
            # versuche es mit einer Konjunkten weniger    
            conj_counter = conj_counter - 1    
        if return_list :
            return return_list
        else : return ""
    else : return ""
    
def substitution_reduce(formula, equivalences) :
    # reducing substitutions replace complex expressions by a single causal factor if one of the formulae in equivalences
    # contains a formula alpha <-> beta, with alpha appearing as term in formula
    
    return_list = []
    
    for equiv in equivalences :
        
        if formula.find(equiv[0]) > -1 :   # falls die linke Seite einer Aequivalenz in der Formel auftritt
            
            aux_formula = ""
            
            if formula.find(equiv[0]) == 0 : # falls der gefundene Ausdruck ganz links in der Formel steht,
                # kann die Ersetzung ohne weiteres durchgefuehrt werden
                
                if (equiv[0].find("+") == -1) or (len(formula) == len(equiv[0])) or (formula[formula.find(equiv[0]) + len(equiv[0])] != "*") :
                    # further condition: either equiv[0] is no disjunction,
                    # or it is equal to formula,
                    # or the match of equiv[0] in formula does not end with "*", which implies that either ")" or " + " is following
                    
                    aux_formula = formula.replace(equiv[0], equiv[1], 1) # erstmaliges Vorkommen durch die rechte Seite ersetzen
                    # weitere Vorkommen koennen durch weitere Aufrufe ersetzt werden
            elif formula[formula.find(equiv[0]) - 1] != "~" : 
                # sonst muss gefordert werden, dass kein Negationssymbol vor dem Ausdruck steht
                # A*B <-> C darf nicht in ~A*B zu ~C eingesetzt werden
                if (equiv[0].find("+") == -1) or (len(formula) == len(equiv[0])) or (formula[formula.find(equiv[0]) + len(equiv[0])] != "*") :
                    # moreover: if the match is part of a conjunction, equiv[0] must not be a disjunction
                    # since we must not substitute in cases like this one:
                    # A + B <-> C
                    # D*A + B <->  E 
                    
                    aux_formula = formula.replace(equiv[0], equiv[1], 1) # erstmaliges Vorkommen durch die rechte Seite ersetzen
            
            if aux_formula.find("(" + equiv[1] + ")") > -1 :
                # falls nun die Form "(A)" gewonnen wurde, loesche die Klammern um "A"
                aux_formula = aux_formula.replace("(" + equiv[1] + ")", equiv[1])
            
            return_list.append(aux_formula)  
    
    if return_list : return return_list
    else : return ""


def idempotence_disj(formula) :
    # formt den ersten Term der Form A*B*C + A*B*C + ... -> A*B*C + ... um
    
    start_formula = formula
    
    # Umgang mit Klammern
    if (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
        # falls sich die erste Disjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        while (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos + 1, right_bracket_pos)]
            formula = "".join(slist)
        
    
    if formula.find("+") > -1 :
        # Idempotenz ist nur moeglich, wenn es mindestens einen Disjunktor gibt
        disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
        
        if disj_list[0] == disj_list[1] :
            # erste und zweite Disjunkte sind identisch -> entferne erste Disjunkte und setze Formel wieder aus disj_list zusammen
            
            slist = [disj_list[i] if i == len(disj_list) - 1 else disj_list[i] + " + " for i in range(1, len(disj_list))]
            st = "".join(slist)
            
            # falls nur erster Klammerausdruck umgeformt wurde, setzen ihn wieder in Gesamtformel ein
            if formula != start_formula :
                st = start_formula.replace(formula,st,1)  
            
            return_list = []
            return_list.append(st)
            return return_list # Rueckgabe in Liste
            
        else : return "" # keine Idempotenz zwischen ersten beiden Disjunkten
        
    else :  return "" # keine Disjunktoren in formula, also gib "" zurueck
    

def idempotence_conj(formula) :
    # formt den ersten Term der Form A*A*B + ... -> A*B + ... um
    
    start_formula = formula
    
    # Umgang mit Klammern
    if (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
        # falls sich die erste Konjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        while (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos + 1, right_bracket_pos)]
            formula = "".join(slist)
        
    
    if formula.find("*") > -1 :
        # Idempotenz ist nur moeglich, wenn es mindestens einen Konjunktor gibt
        
        valid = True
        if formula.find("+") > -1 :
            # falls die Formel eine Disjunktion ist, spalte in Disjunkte auf
            disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
            
            if disj_list[0].find("*") > -1 : # erste Disjunkte muss einen Konjuktor enthalten
                conj_list_one = re.split("\*", disj_list[0]) # nur erste Disjunkte ist von Relevanz
            else : valid = False     # sonst kann die Formel nicht in der richtigen Weise idempotent sein
            
        else : # Fall 2: formula ist eine Konjunktion, keine Disjunktion
            conj_list_one = re.split("\*", formula) 
        
        
        if valid :
            # formula enthaelt an der relevanten Stele einen Konjunktor
            
            if conj_list_one[0] == conj_list_one[1] :
                # erste und zweite Konjunkte sind identisch -> entferne erste Konjunkte und setze Formel wieder aus conj_list_one
                # und disj_list zusammen
                
                # Zusammensetzung des Restes der ersten Disjunkten
                slist = [conj_list_one[i] if i == len(conj_list_one) - 1 else conj_list_one[i] + "*" for i in range(1, len(conj_list_one))]
                st = "".join(slist)
                
                if formula.find("+") > -1 :
                    # Hinzufuegen der anderen Disjunkten, falls Ausgangsformel eine Disjunktion ist
                    slist = [" + " + disj_list[i] for i in range(1, len(disj_list))]
                    st = st + "".join(slist)
                
                # falls nur erster Klammerausdruck umgeformt wurde, setzen ihn wieder in Gesamtformel ein
                if formula != start_formula :
                    st = start_formula.replace(formula,st,1)  
            
                return_list = []
                return_list.append(st)
                return return_list # Rueckgabe in Liste
            
            
            else : return "" # keine Idempotenz zwischen ersten beiden Konjunkten

        else : return "" # Konjunktoren an falschen Stellen
        
    else :  return "" # keine Konjunktoren in formula, also gib "" zurueck

def double_negation_elimination(formula) :
    # eliminates double negations ("~~"), returns the list of possible formulae with successively eliminated double negations
    # if formula does not contain double negations, return ""
    
    return_list = []
    while formula.find("~~") > -1 :
        formula = formula.replace("~~","",1)
        return_list.append(formula)
    
    if return_list : return return_list 
    else : return ""    

def double_negation_introduction(formula) :
    # introduces a double negation ("~~"), this might be necessary for performing de Morgan transformations
    # fairly limited application restriction in order not to inflate the space of possible formulae to much
    # only applicable at the first disjunct inside the first bracket, if this bracket is a disjunction and the first term is not already
    # negated
    
    return_list = []
    
    if (len(formula) > 4) and (formula[0] == "~") and ((formula[0] == "(" and formula[1] != "~") or (re.match(r'^[A-Z]',formula))) :
        # checks whether formula starts with a bracket, which is not followed by a negation

        if (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
            formula = formula[0] + "~~" + formula[1:]
            return_list.append(formula)
            
    if return_list : return return_list
    else : return "" 

# auxiliary function for parallel processing
def smap(f):
    return f()

def search_formula(formula,formula_list) :
    # determines all equivalent expressions that can be generated from formula by means of six transformation operations and 
    # checks whether a formula from formula_list can be generated from them
    # returns the corresponding truth value
    #
    # Valid transformation operations are:
    # 1) substitution by means of the equivalence formulae from formula_list
    # 2) distributive transformation: A*B*C + A*B*~C*D + A*B*G + B ... -> A*B*(C + ~C*D + G) + B ...
    # 3) commutativity of disjunction
    # 4) commutativity of conjunction
    # 5) idempotence of disjunction
    # 6) idempotence of conjunction
    
    # create tree with input formula
    root = Tree(formula[0]) 
    
    # filtering of formula_list for useful formulae and separation into two types
    formula_list_eq = []
    formula_trans_list = []
    for f in formula_list :
        # a) formulae that have the same term on the right side as formula
        # These will be used to check the redundancy of formula.
        if f[1] == formula[1] :
            formula_list_eq.append(f)
        
        # b) formulae with rightside terms of lower causal order than formula
        # these will be used for substitutions
        else :  # (the conditions are already checked in the main function, so we know that all formulae that do not fullfill
            # the first condition comply to this one)
            formula_trans_list.append(f)
    
    formula_list.remove(formula)
    
    # set of all nodes of our tree, such that we do not add the same formula twice
    # In python, sets are more efficient in checking "a in list/set" than lists, this is why we do not use a
    # proper list here.
    full_set = {formula[0]}

    
    print("Checking tree " + str(root))
    
    # terminating the loops as soon as we have a result
    stop = not(formula_list_eq)   # If formula_list_eq is empty, we will never find an equivalent of formula therein, so we may stop
    # immediately.


    for tree in root.deep_search() :
        # start deep search of tree, while iteratively deepening it
        
        for formula_t in formula_list :
            # comparison of our obtained formula with those from formula_list
            if formula_t[0] == str(tree) : # If they are equal, we have finished.
                 print(formula_t[0] + " found.")
                 stop = True
                 break  # Break from inner for-loop
            
        if stop : break # Break from outer for-loop
       
      
        # If we haven't finished yet, we have to generate more tree nodes.
        # Perform all transformation operations on the current node at once:
        func1 = functools.partial(substitution_reduce, str(tree),formula_trans_list)
        func2 = functools.partial(de_morgan_disj, str(tree))
        func4 = functools.partial(distribution, str(tree))
        func6 = functools.partial(commutation_disj, str(tree))
        func7 = functools.partial(commutation_conj, str(tree))
        func8 = functools.partial(idempotence_disj, str(tree))
        func9 = functools.partial(idempotence_conj, str(tree))
        func10 = functools.partial(double_negation_elimination, str(tree))
        func11 = functools.partial(double_negation_introduction, str(tree))
        
        pool = multiprocessing.Pool(processes=9)
        res = pool.map(smap, [func1, func2, func4, func6, func7, func8, func9, func10, func11])
        # end process pooling and join results
        pool.close()
        pool.join()
                      
        
        
        
        for sub_list in res :
            if not(sub_list == "") : # All Operations are defined to yield "" if they aren't applicable.
                for j in sub_list :  
                    if not(j in full_set) : # If they are, check (via the check set) whether transformed formulae are already in the tree.
                        full_set.add(j)     # If they are new, add them to tree and check set.
                        tree.add_child(Tree(j))
        
                                
    return stop # True, if an equivalent formula has been found, False, if all nodes have been traversed without any match

def main() :
    # just a dummy function to test search_formula
    # We define one test formula three further formulae for substitutions and comparisons.

    
    #start_f = "E*G*~D*G*H + E*G*H + E*A" # our test formula
    #start_f = "~A*~D + ~B*~D + ~C*~D + ~A*E + ~B*~~E + ~C*E"
    start_f = "~(~~A*B*C)*~(D*~E)"
    start_eq = "J"
    formula = (start_f, start_eq)
    
    # a couple of further formulae for the formula_list
    st = "H*I"
    st1 = "G"
    st2= "E*A"
    st3 = "D"
    
    # Each formula is an equivalence alpha <-> beta, so form pairs (alpha, beta).
    formula_1 =(st,st1)
    formula_2 = (st2,st3)
    
    st4 = "D + E*G*~D*H + E*G*H" # our test formula should be transformed to this one
    st5 = "J"
    formula_3 = (st4,st5) 
    
    st6 = "A*B*C + D*~E"
    st7 = "C"
    formula_4 = (st6,st7) 
    
    st8 = "~C"
    st9 = "J"
    formula_5 = (st8,st9) 

    formula_list = []
    formula_list.append(formula_1)
    formula_list.append(formula_2)
    formula_list.append(formula_3)
    formula_list.append(formula_4)
    formula_list.append(formula_5)
    formula_list.append(formula)
    print(search_formula(formula, formula_list)) # start search_formula and print the result
    
if __name__ == '__main__':
    main()
