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
    
    # deletes end-of-line-symbol ("\n") and spaces at the end of line of necessary
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
    # gibt eine Liste der im String st vorkommenden Elemente aus factor_list zurueck
    component_list = []                          # Deklaration der Ausgabeliste
    
    # Schritt 1: auf welcher Ebene von factor_list befinden sich die Faktoren?
    # - factor_list -> 1
    # - level_factor_list -> 2
    # - level_factor_list_order -> 3
    
    if isinstance(factor_list[0], str) :
        # Fall 1: gehe Elemente von factor_list durch
        for element in factor_list :
            if st.find(element) > -1 :
                component_list.append(element)
    
    elif isinstance(factor_list[0], list) :
        if isinstance(factor_list[0][0], str) :
            # Fall 2: gehe Elemente der Elemente von factor_list durch
            for m in range(len(factor_list)) :
                for element in factor_list[m] :
                    if st.find(element) > -1 :
                        component_list.append(element)
                        
        elif isinstance(factor_list[0][0], list) :
            if isinstance(factor_list[0][0][0], str) :
                # Fall 3: gehe Elemente der Elemente der Elemente von factor_list durch
                for m in range(len(factor_list)) :
                    for o in range(len(factor_list[m])) :
                        for element in factor_list[m][o] :
                            if st.find(element) > -1 :
                                component_list.append(element) 
    
    return component_list

def get_formula_level(st, level_factor_list):
    # falls alle Faktoren in st vom gleichen Level sind, wird dieses Level zurueck gegeben
    # anderenfalls ist der Rueckgabewert -1
    
    inequal = False
    factors = get_components_from_formula(st, level_factor_list)
    level = 0
    
    # je nachdem ob level_factor_list eine Liste von Elementen oder eine Liste von Listen von Elementen (nach Ordnungen separiert) ist.
    # muss verschieden vorgegangen werden
    # daher Unterscheidung nach multi_level
    multi_order = isinstance(level_factor_list[0], list)
    
    if multi_order :
        for m in range(len(level_factor_list)) :
            for i in range(len(level_factor_list[m])) :
                if factors[0] in level_factor_list[m][i] :
                    level = m
            
            # pruefe ob alle anderen factors vom gleichen Level sind
            for k in range(1,len(factors)):
                inequal = True
                # for-Schleife ueber die Ordnungen (fuer diese Untersuchung irrelevant)
                for i in range(len(level_factor_list[level])) :
                    if factors[k] in level_factor_list[level][i] :
                        inequal = False
        
            if inequal:
                level = -1
    
    else :
        # suche das Level von factors[0]
        for i in range(len(level_factor_list)) :
            if factors[0] in level_factor_list[i] :
                level = i
            
        # pruefe ob alle anderen factors vom gleichen Level sind
        for k in range(1,len(factors)):
            if not(factors[k] in level_factor_list[level]) :
                inequal = True
        
        if inequal:
            level = -1
        
    return level

def get_factor_order(factor, factor_list):
    # Bestimmt die Ordnung eines Kausalfaktors bezueglich einer Faktorenliste der Form 
    # liste[LEVEL][ORDNUNG][FAKTOREN], sobald factor in FAKTOREN enthalten ist, gebe ORDNUNG zurueck
    # ansonsten -1
    order = -1
    for m in range(len(factor_list)):
        for o in range(len(factor_list[m])):
            if factor in factor_list[m][o]:
                order = o
    return order
    
def get_formula_order(formula, factor_list):
    # Bestimmt die Ordnung einer Formel = die hoechste Ordnung eines darin auftretenden Kausalfaktors bezueglich
    # einer Faktorenliste der Form, gemaess get_factor_order
    # falls Ordnung wohl definiert ist, wird diese zurueck gegeben
    # ansonsten -1
    order = -1
       
    
    for fac in get_components_from_formula(formula, factor_list):
        if order < get_factor_order(fac, factor_list) :
            order = get_factor_order(fac, factor_list)
                
    return order


def convert_causal_relation(formula, level_factor_list_order, tex_code):
    # ueberfuehrt formula in TikZ-Latex-Code um und gibt diesen zurueck als String zurueck
    st = ""
    
    #########################################
    # Ermittle syntaktischen Typ der Formel #
    #########################################
    # Annahme: cna gibt nur Formeltypen aus:
    # 1) Aequivalenz (A <-> B)
    # 2) Nicht-Aequivalenz (~A <-> B)
    # 3) Konjunktionen von Faktoren ( A*B* ... <-> E)
    # 4) Konjunktionen von Faktoren oder Negaten ( A*~B* ... <-> E)
    # 5) Disjunktionen von Faktoren ( A + B + ... <-> E)
    # 6) Disjunktionen von Faktoren oder Negaten ( A + ~B + ... <-> E)
    # 7) Disjunktionen von Faktoren oder Negaten mit Konjunktionen ( A + ~B*C + ... <-> E)
    
    
    level = get_formula_level(formula[0], level_factor_list_order)
    
    #################
    # Aequivalenzen #
    #################
    
    for o in range(len(level_factor_list_order[level])) :
        for fac in level_factor_list_order[level][o] :
            if fac == formula[0] :
                st = "\draw[->] (" + fac + ".east) -- (" + formula[1] + ".west);"
            elif formula[0] == "~" + fac :
                st = "\\node[neg] (" + fac + "neg) at ([xshift=\LNeg]" + fac + ".south east) {};\n\draw[->] (" + fac + "neg) -- (" + formula[1] + ");"
                
                
    if st == "" :
        # weiter mit Verknuepfungen zwischen Kausalfaktoren (und, oder)
        if formula[0].find("+") > -1 :
        
            #################
            # Disjunktionen #
            #################
            
            disjunctor_list = re.split("\s*\+\s*", formula[0])
            
            for disj in disjunctor_list : #get_components_from_formula(formula[0], level_factor_list_order) :

                
                if disj in get_components_from_formula(formula[0], level_factor_list_order) :
                    # Fall 1: Disjunktor ist ein einfacher Kausalfaktor
                    st = st + "% einfache Disjunktion\n"
                    st = st + "\draw[->] (" + disj + ".east) to (" + formula[1] + ".west);\n"
                
                elif formula[0].find("*") > -1 :
                    # Fall 2: disj ist eine Konjunktion
                    st = st + "% komplexe Disjunktion\n"
                    
                    conjunctor_list = re.split("\*", disj)
                    
                    # zuerst setze gemeinsamen Schnittpunkt
                    # Platziere diesen neben (erster) Konjunkten der hoechsten Ordnung in der Liste
                    # ermittle Node der dieser Konjunkten -- f_fac
                    cross_point = ""                    # Name der Node des Schnittpunkts der Konjunkten
                    for fac in get_components_from_formula(disj, level_factor_list_order) :
                        cross_point = cross_point + fac
                        if fac == get_components_from_formula(disj, level_factor_list_order)[0] :
                            f_fac = fac
                        else :
                            if get_factor_order(fac, level_factor_list_order) > get_factor_order(f_fac, level_factor_list_order) :
                                f_fac = fac
                    
                    cross_point = cross_point + formula[1]
                    position = "at ([xshift=\hDisjConj, yshift=\\vDisjConj]" + f_fac + ".east)"
                    # Vorsicht: Es kann passieren, dass mehrere Disjunkten aus Konjunktionen beim gleichen Faktor f_fac
                    # zusammentreffen, daher muss die Position des Konjunktionsknotenpunktes ggf. verschoben werden
                    q = 1
                    while (tex_code.find(position) > -1) or (st.find(position) > -1) :
                        position = "at ([xshift=\hDisjConj, yshift={\\vDisjConj + " + str(q) + "*\\tDisjConj}]" + f_fac + ".east)"
                        q = q + 1
                        
                    st = st  + "% Treffpunkt der Konjunkten\n\\node[aux] (" + cross_point + "aux) " + position + " {};\n% Teilpfeile von den Konjunkten zum Schnittpunkt\n"
                        
                    
                    
                        
                    

                    
                    for conj in conjunctor_list :
                    
                        # verbinde Konjunkten mit dem Schnittpunkt
                        if conj[0] == "~" and conj[1:] in get_components_from_formula(disj, level_factor_list_order) :
                            # Konjunkte ist negierter Faktor
                            st = st  + "\\node[neg] (" + conj[1:] + "neg) at ([xshift=\LNeg]" + conj[1:] + ".south east) {};\n"
                            st = st + "\draw[conjunctonsegment] (" + conj[1:] + "neg) to (" + cross_point + "aux);\n"
                            
                        elif conj in get_components_from_formula(disj, level_factor_list_order) :
                            # Konjunkte ist einfacher Faktor
                            st = st + "\draw[conjunctonsegment] (" + conj + ") to (" + cross_point + "aux);\n"
                            
                        else :
                            # was koennte sonst passieren(??)
                            print("Die Disjunktion " + formula[0] + " -> " + formula[1] + "  konnte nicht eingezeichnet werden, da ihre Substruktur nicht erkannt wurde(1).")

                    # verbinde Schnittpunkt mit Zielfaktor
                    st = st + "% Pfeil vom Schnittpunkt zum Zielfaktor\n\draw[->] (" + cross_point + "aux) -- (" + formula[1] + ".west);\n"
                
                elif (disj[0] == "~") and (disj[1:] in get_components_from_formula(formula[0], level_factor_list_order)) :
                    # Fall 3: disj ist ein negierter Faktor (erstes Zeichen = "~" und weitere Zeichen entsprechen
                    # einem bekannten Kausalfaktor
                    
                    # ANNAHME: in einer Disjunktionskette kann derselbe Faktor nur entweder negiert oder 
                    # ohne Negation auftreten (sonst muss aus "elif" ein neues "if" gemacht werden)
                    
                    st = st + "% negierte Disjunkte\n"
                    st = st  + "\\node[neg] (" + disj[1:] + "neg) at ([xshift=\LNeg]" + disj[1:] + ".south east) {};\n"
                    st = st + "\draw[->] (" + disj[1:] + "neg) to (" + formula[1] + ".west);\n"
                
                
                else :
                    # disj besitzt eine andere Struktur
                    
                    print("Die Disjunktion " + formula[0] + " -> " + formula[1] + "  konnte nicht eingezeichnet werden, da ihre Substruktur nicht erkannt wurde.")

            
            
        elif formula[0].find("*") > -1 :
            
            #################
            # Konjunktionen #
            #################
            
            # setze Treffpunkt der Konjunkten
            st = "% Treffpunkt der Konjunkten\n\\node[aux] (" + formula[1] + "aux) at ([xshift=\LConj]" + formula[1] + ".west) {};\n% Teilpfeile von den Konjunkten zum Schnittpunkt\n"
            
            # zeichne Pfeile von Konjunkten bis Schnittpunkt
            for conj in get_components_from_formula(formula[0], level_factor_list_order) :

                if formula[0].find("~" + conj) > -1 :
                    # falls Faktor conj negiert in der Formel auftaucht
                    
                    # ANNAHME: in einer Konjunktionskette kann derselbe Faktor nur entweder negiert oder 
                    # ohne Negation auftreten (sonst muss aus dem folgenden "else" ein neues "if" gemacht werden)
                    
                    st = st  + "\\node[neg] (" + conj + "neg) at ([xshift=\LNeg]" + conj + ".south east) {};\n"
                    st = st + "\draw[conjunctonsegment] (" + conj + "neg) to (" + formula[1] + "aux);\n"
                else :
                    # Faktor conj taucht unnegiert auf
                    
                    st = st + "\draw[conjunctonsegment] (" + conj + ") to (" + formula[1] + "aux);\n"
            
            # zeichne Pfeil von Schnittpunkt bis Zielfaktor
            st = st + "% Pfeil vom Schnittpunkt zum Zielfaktor\n\draw[->] (" + formula[1] + "aux) -- (" + formula[1] + ");"
        else :
            
            ##########################
            # Struktur nicht erkannt #
            ##########################
            
            print(formula[0] + " -> " + formula[1] + "  konnte nicht eingezeichnet werden, da ihre Struktur nicht erkannt wurde.")
            
    return st
    
def convert_constitution_relation(formula, level_factor_list_order, constitution_relation_list) :
    # wandelt Formelsyntax in TikZ-Latex-Code um und gibt diesen zurueck als String zurueck
    st = ""
    
    # Konstitutionsverbindungen sind je nach dem, ob der Faktor die Substruktur nach rechts oder links begrenzt, verschieden gestaltet
    c_left = False
    c_right = False
    
    # Pruefung, ob untersuchte Formel nach links oder rechts begrenzend ist    
    for f in constitution_relation_list :
        if (formula[1] == f[1]) and (formula[0] != f[0]) :
            # gibt es eine andere Konstitutionsformel bezgl. den gleichen Faktor, die aus Faktoren hoeherer
            # (Kausal-)Ordnung besteht?
            if get_formula_order(formula[0], level_factor_list_order) < get_formula_order(f[0], level_factor_list_order) :
                c_left = True
            elif get_formula_order(formula[0], level_factor_list_order) > get_formula_order(f[0], level_factor_list_order) :
                c_right = True
    
    # fuer jeden Faktor in der Konstitutionsrelation wird eine Verbindung zum hoeherstufigen Faktor gezogen
    for fac in get_components_from_formula(formula[0], level_factor_list_order) :           
        if c_left and not(c_right) :
            # nach links begrenzend
            st = st + "\draw[crelationleft] (" + fac + ".north west) to (" + formula[1] + ".south);\n"
    
        elif not(c_left) and c_right :
            # nach rechts begrenzend
            st = st + "\draw[crelationright] (" + fac + ".north east) to (" + formula[1] + ".south);\n"
    
        else: 
            # weder noch (beispielsweise weil ein Kausalfaktor allein konstituierend fuer einen hoeherstufigen Faktor ist)
            st = st + "\draw[crelationstraight] (" + fac + ".north) to (" + formula[1] + ".south);\n"
    
    return st

# ab hier gehts zur Sache
def read_R_file(file_name):
    # liest cna-Ausgabe aus und gibt relevante Daten in Form von listen zurueck
    # level_factor_list -- nach Konstitutionslevel separierte Listen von Kausalfaktoren,
    # level_equiv_list -- nach Konstitutionslevel separierte Listen von Kausalrelationen, sowie
    # constitution_relation_list -- Liste von Konstitutionsbeziehungen
    
    
    file_lines = []                              # Deklaration der Liste fuer die Textzeilen
    with open (file_name, 'rt') as text_file:    # Oeffne file_name um Textdatei auszulesen
        for next_line in text_file:              # fuer jede Zeile in file_name
            file_lines.append(next_line)         # fuege sie zur Liste file_lines hinzu
    abort = False
           
    ################################################ 
    # Schritt 1: bestimme Menge der Kausalfaktoren #
    ################################################
    factor_list = []                             # Deklaration der Liste der Kausalfaktoren
    
    # offenbar steht in Zeile 3 der R-Ausgabe immer entweder "Causal ordering:" oder "Factors:"
    # das sollte spaeter robuster gestaltet werden
    if file_lines[2].find("Causal ordering:") > -1:
        # im Falle von "Causal ordering:" sind Faktoren in Zeile 4 gelistet
        factor_list = find_causal_factors(file_lines[3]) 
        multi_level = True                       # multi_level, wenn dies in R angegeben worden ist
        
    elif file_lines[2].find("Factors:") > -1:
        # im Falle von "Factors:" sind Faktoren in Zeile 3 gelistet
        factor_list = find_causal_factors(file_lines[2])
        multi_level = False
    else:
        print("Programmabbruch keine Kausalfaktoren gefunden")
        abort = True
    
    # weiter falls Kausalfaktoren gefunden worden sind (Liste factor_list ist nicht leer)
    if factor_list:
        #################################
        # Schritt 2: suche nach Formeln #
        #################################
        
        equiv_list = []                          # Deklaration der Liste fuer Aequivalenzformeln
            
        for line in file_lines:
            if line.find("<->") > 0:
                # falls "<->" gefunden wurde, lies Teilformeln links und rechts aus
                # und fuege sie zu equiv_list hinzu
                equiv_list.append(get_equiv_formula(line))
                    
        if not(equiv_list):                      # falls equiv_list leer Programmabbruch
            print("Programmabbruch keine Formeln gefunden")
            abort = True
            return abort, "", "", ""
        else:                                    # sonst weiter mit
            # pruefe ob jeder Kausalfaktor in mindestens einer Formel vorkommt, ansonsten loesche ihn
            # aus factor_list
            for k in range(len(factor_list)-1,-1,-1):
                i = 0
                found = False
                while not(found) and i < len(equiv_list):
                    found = (equiv_list[i][0].find(factor_list[k]) > -1 or equiv_list[i][1] == factor_list[k])
                    i = i + 1
                    
                if not(found):  # das Loeschen in der laufenden for-Schleife sollte aufgrund des Rueckwaertslaufens
                    # keine Probleme machen
                    print("Faktor " + factor_list[k] + " gelöscht, da in keiner Formel auftretend")
                    factor_list.remove(factor_list[k])  
                    
            
            # Separierung der Ebenen in den Formeln
            # Vorsicht die folgende Zeile haengt kritisch von der Struktur der Eingabedatei ab!!!
            # Kausalfaktoren muessen in 4. Zeile stehen
            level_count = file_lines[3].count("<")
            # level_count = Anzahl der Levels -1
            
            level_factor_list = []               # Deklaration der separater Faktorlisten
            level_equiv_list = []
            constitution_relation_list = []
            for i in range(level_count + 1):
                if multi_level:
                    st = re.split(" < ",file_lines[3])[i].strip()
                    level_factor_list.append(find_causal_factors(st))
                    level_equiv_list.append([])
                else :
                    level_factor_list.append(find_causal_factors(file_lines[2]))
                    level_equiv_list.append([])

                
            #######################################
            # Schritt 3: Entflechtung der Formeln #
            #######################################
            for formula in equiv_list : 
                # fuege alle Formeln, die allein aus Elementen der Ebene i bestehen, zu level_equiv_list[i] hinzu
                if get_formula_level(formula[1], level_factor_list) == get_formula_level(formula[0], level_factor_list) :
                    level_equiv_list[get_formula_level(formula[1], level_factor_list)].append(formula)
                    
                elif get_formula_level(formula[0], level_factor_list) > -1 :
                    # alle Formeln, bei denen rechts ein Element aus einer anderen Ebene mit Elementen links (alle aus gleicher Ebene)
                    # verbunden wird, werde zu constitution_relation_list hinzugefuegt
                    # !!! ANNAHME: Es werden nur Ebenen mit Leveldifferenz von 1 miteinander verbunden !!!!!!
                    if get_formula_level(formula[0], level_factor_list) == get_formula_level(formula[1], level_factor_list) - 1 :
                        constitution_relation_list.append(formula)
                        
                           
                # alle Formeln mit Ebenenvermischung auf der linken Seite sind damit verworfen
    return abort, level_factor_list, level_equiv_list, constitution_relation_list
    
    
    
def determine_factor_order(level_factor_list,level_equiv_list) :    
    level_count = len(level_factor_list)
    abort = False
    unique = True 
    
    # rekonstruiere factor_list aus level_factor_list
    factor_list = []
    zaehler = -1
    for level in level_factor_list :
        zaehler = zaehler + 1
        if level_factor_list[zaehler] :
            for element in level_factor_list[zaehler] :
                factor_list.append(element)
            
    level_factor_list_order = []  # Deklaration einer Liste (je Level) von Listen (je Ordnung) deren Eintraege die Faktoren sind
    # Ordnung = 0 -> Eingangsfaktor, Ordnung = i -> alle abhaengigen
    # Faktoren sind Ordnung < i und mindestens einer der Ordnung i-1
            
    # ab hier folgend alle Schritte nach Ebenen separiert
    for m in range(level_count):
                
        ########################################
        # Schritt 4: bestimme Eingangsfaktoren #
        ########################################
                
                
        level_factor_list_order.append([]) 
        level_factor_list_order[m].append([])  # ein Element hinzufuegen fuer Liste der Faktoren mit Ordnung 0
                
        # zunaechst Bestimmung der Eingangsfaktoren (Ordnung = 0)
        for element in level_factor_list[m]:
            # Eingangsfaktoren stehen nie rechts in Formeln
            i = 0
            found = False
            while not(found) and i < len(level_equiv_list[m]):
                found = (level_equiv_list[m][i][1] == element)
                i = i + 1
                    
            if not(found):
                # fuege element zu incoming_factor_list
                level_factor_list_order[m][0].append(element)
        
        if not(level_factor_list_order[m][0]):
            # keine Eingangsfaktoren identifiziert
            # zu Testen gibt es zirkulaere Beziehungen (A<->B; B<->A) -- Struktur als nicht eindeutig klassifizieren
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
                # keine Eingangsfaktoren gefunden und keine uneindeutige Struktur -> Programmabbruch
                print("Programmabbruch keine eingehenden Kausalfaktoren gefunden")
                abort = True
                    
           
                
        else:
                        
                
            ############################################################
            # Schritt 5: lege Ordnung der Kausalfaktoren iterativ fest #
            ############################################################
                    
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
                            for element in get_components_from_formula(level_equiv_list[m][i][0], factor_list):
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
            
            # vorerst Ausgabe als Text in Konsole -> spaeter graphisch als Diagramm        
            print("Kausalfaktoren vom Level " + str(m) + " mit ihrer Ordnung:")
            for p in range(len(level_factor_list_order[m])):
                print("Ordnung " + str(p))
                print(level_factor_list_order[m][p])
            
                       
            if downstream_factor_list:
                # falls downstream_factor_list nicht leer ist, ist Konstitutionsstruktur nicht eindeutig
                print("Keine eindeutige Struktur")
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
    # Beschreibung folgt
    ######################################
    # step 8: preparing the output files #
    ######################################
    
    tex_code = "% Platzierung der Nodes\n"
    
    ######################################################
    # step 8a) - placement of the nodes = causal factors #
    ######################################################
    
    # [Verbesserungsmoeglichkeit 1]
    # - mit jener Ebene starten, die die hoechste Zahl an Ordnungen besitzt (horizontal am laengsten ist)
    # hier: Start auf unterster Ebene
    
    placement = ""
    
    for m in range(len(level_factor_list_order)) :
        # zur besseren Arbeit mit der tex-Ausgabe werden tex-Kommentare gesetzt
        tex_code = tex_code  + "% Faktoren der Ebene " + str(m) + ":\n"
        
        max_num_factors_order = 0
        for o in range(len(level_factor_list_order[m])) :
            
            # Platzierung der Faktoren in kausaler Reihenfolge
            # Faktoren gleicher Ordnung stehen uebereinander
            
            # [Verbesserungsmoeglichkeit 2]
            # - Faktoren, die in Konstitutionsrelationen gemeinsam auftreten zusammmen gruppieren
            # hier: Vorgehen nach Reihenfolge in List
            
            tex_code = tex_code  + "% Kausalordnung " + str(o) + ":\n"
            for e in level_factor_list_order[m][o] :
                
                # fuege Zeile zu tex_code hinzu, der Node platziert, mit Namen des Faktors bezeichnet 
                # und mit gleichem Label versieht
                tex_code = tex_code + "\\node" + placement + " (" + e + ") {" + e +"};\n"
                
                if o == 0 :
                    # Eingangsfaktoren als solche kennzeichnen
                    tex_code = tex_code + "\hilightsource{" + e + "};\n"
                elif o == len(level_factor_list_order[m]) - 1 :
                    # Ausgangsfaktoren kennzeichnen
                    tex_code = tex_code + "\hilighttarget{" + e + "};\n"
                
                # Anpassung von placement fuer naechsten Faktor
                placement = "[above= \LvDist of " + e + "]"
                # naechster Kausalfaktor gleicher Ebene und Ordnung befinden sich um \LvDist oberhalb des aktuellen
                
            # naechste Ordnung beginnt -> naechster Faktor wird rechts auf Hoehe des untersten Faktors dieser Ordnung gesetzt
            
            placement = "[right= \LhDist of " + level_factor_list_order[m][o][0] + "]"
            
            if len(level_factor_list_order[m]) > max_num_factors_order :
                max_num_factors_order = len(level_factor_list_order[m])
                
        # naechste Ebene beginnt -> naechster Faktor wird gegenueber der aktuellen Ebene um \iLvDist nach oben verschoben
        # genauer: \iLvDist  + Hoehe der aktuellen Ebene (= max_num_factors_order * \LvDist) oberhalb des ersten Elements
        # der ersten Ordnung dieser Ebene 
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
    
    # pdf fertig generiert und aufgeraeumt
    
# Ende print_structure_in_tikz_plot    
                   
def main() :
    # main function
    
    
    level_factor_list = []               # declaration of the lists
    level_factor_order_list = []
    level_equiv_list = []
    constitution_relation_list = []
    
    # steps 1 - 3 start function read_R_file -> converts cna output into lists that are sorted by constitution level
    # of causal factors (level_factor_list), causal relations (level_equiv_list) and one list of constitution relations
    # (constitution_relation_list), if the cna output is not as expeceted stop the procedure with abort = True
    abort, level_factor_list, level_equiv_list, constitution_relation_list = read_R_file("r_output.txt")
    
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
