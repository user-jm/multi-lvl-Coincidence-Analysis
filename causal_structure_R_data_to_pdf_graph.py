#!/usr/bin/env python3

# Datei: causal_structure_R_data_to_pdf_graph.py

import re                                        # Regex fuer Suche in Strings

# erst einmal nur Hilfsfunktionen:
def find_causal_factors(st):
    # sucht in String st nach durch ", ", " < " oder "*" separierten Bezeichnern fuer Kausalfaktoren
    # gibt die Liste an Kausalfaktoren zurueck
    
    # loesche "Factors: " aus Textzeile (sofern auftretend)
    st = st.replace("Factors: ","")
    # loesche Zeilenendensymbol und ggf. Leerzeichen am Zeilenende
    st = re.sub("\r?\n","",st).rstrip()
    # gib Liste der durch ", ", "< " oder "*" separierten Eintraege zurueck
    return re.split(",\s*|\s*<\s*|\*", st)

def get_equiv_formula(st):
    # gibt die Teilformeln links und rechts von "<->" in Eingabestring st als Paar von Strings (a,b) zurueck
    a = re.split(" <-> ",st)[0].strip()          # strip() loescht fuehrende Leerzeichen
    b = re.split(" <-> ",st)[1].strip()
    # in cna-Ausgabe haengen rechts noch weitere Informationen, diese von b abtrennen:
    b = re.split("[ \t]",b)[0]
    return (a,b)

def get_components_from_formula(st, factor_list):
    # gibt eine Liste der im String st vorkommenden Elemente aus factor_list zurueck
    component_list = []                          # Deklaration der Ausgabeliste
    for element in factor_list:
        if st.find(element) > -1:
            component_list.append(element)
    
    return component_list

def get_formula_level(st, level_factor_list):
    # falls alle Faktoren in st vom gleichen Level sind, wird dieses Level zurueck gegeben
    # anderenfalls ist der Rueckgabewert -1
    
    inequal = False
    factors = find_causal_factors(st)
    level = 0
    
    # suche das Level von factors[0]
    for i in range(len(level_factor_list)):
        if factors[0] in level_factor_list[i]:
            level = i
            
    # pruefe ob alle anderen factors vom gleichen Level sind
    for k in range(1,len(factors)):
        if not(factors[k] in level_factor_list[level]):
            inequal = True
        
    if inequal:
        level = -1
        
    return level

def get_factor_order(faktor, liste):
    # Bestimmt die Ordnung eines Kausalfaktors bezueglich einer Faktorenliste der Form 
    # liste[LEVEL][ORDNUNG][FAKTOREN], sobald faktor in FAKTOREN enthalten ist, gebe ORDNUNG zurueck
    # ansonsten -1
    ordnung = -1
    for m in range(len(liste)):
        for o in range(len(liste[m])):
            if faktor in liste[m][o]:
                ordnung = o
    return ordnung
    
def get_formula_order(formula, liste):
    # Bestimmt die Ordnung einer Formel = die hoechste Ordnung eines darin auftretenden Kausalfaktors bezueglich
    # einer Faktorenliste der Form, gemaess get_factor_order
    # falls Ordnung wohl definiert ist, wird diese zurueck gegeben
    # ansonsten -1
    ordnung = -1
    
    auxiliary_list = [] # unstrukturierte Hilfsliste aller Elemente von liste, da die Funktion get_components_from_formula
    for m in range(len(liste)):  # nur mit unstrukturierten Listen umgehen kann
        for o in range(len(liste[m])):
            for e in liste[m][o]:
                auxiliary_list.append(e)    
    
    for fac in get_components_from_formula(formula, auxiliary_list):
        if ordnung < get_factor_order(fac, liste):
            ordnung = get_factor_order(fac, liste)
                
                
    auxiliary_list.clear()
    return ordnung

def substitute_formula(formel_1,formel_2):
    # substituitert falls moeglich im linken Term von formel_1 alle Vorkommnisse des Ausdrucks der rechts in formel_2 
    # durch den linken in formel_2
    # gibt entsprechend die linke Seite der eingesetzten Formel zurueck
    # falls keine Einsetzung moeglich, gib "" zurueck
    if formel_1[0].find(formel_2[1]) > -1:
        st = formel_1[0].replace(formel_2[1],formel_2[0])
    else:
        st = ""
    
    return st

def convert_formula(st):
    # wandelt Formelsyntax von String st in Latex-Code um und
    # gibt diesen zurueck
    
    
    return st

# ab hier gehts zur Sache
def read_R_file(file_name):
    file_lines = []                              # Deklaration der Liste fuer die Textzeilen
    with open (file_name, 'rt') as text_file:    # Oeffne file_name um Textdatei auszulesen
        for next_line in text_file:              # fuer jede Zeile in file_name
            file_lines.append(next_line)         # fuege sie zur Liste file_lines hinzu
    abbruch = False
           
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
        abbruch = True
    
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
            abbruch = True
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
            
            
            level_factor_list_order = []  # Deklaration einer Liste (je Level) von Listen (je Ordnung) deren Eintraege die Faktoren sind
            # Ordnung = 0 -> Eingangsfaktor, Ordnung = i -> alle abhaengigen
            # Faktoren sind Ordnung < i und mindestens einer der Ordnung i-1
            
            # ab hier folgend alle Schritte nach Ebenen separiert
            for m in range(level_count + 1):
                
                ########################################
                # Schritt 4: bestimme Eingangsfaktoren #
                ########################################
                
                
                level_factor_list_order.append([]) 
                
                # zunaechst Bestimmung der Eingangsfaktoren (Ordnung = 0)
                for element in level_factor_list[m]:
                    level_factor_list_order[m].append([]) 
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
                    print("Programmabbruch keine eingehenden Kausalfaktoren gefunden")
                    abbruch = True
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
                            ordnung = 0           # Ordnung des betrachteten Faktors

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
                                                if ordnung < o + 1:
                                                    ordnung = o + 1
                                
                                i = i + 1
                                
                            k = k + 1
                        
                        # falls in allen Formeln des Levels, alle Faktoren auf der linken Seite
                        # level_factor_list_order[m] sind, wird downstream_factor_list[k-1] in einen solchen
                        # umgewandelt mit der entsprechenden Ordnung (k-1, weil zwischenzeitlich k um eins erhoeht wurde)    
                        if candidate :
                            level_factor_list_order[m][ordnung].append(downstream_factor_list[k-1])                    
                            
                            downstream_factor_list.remove(downstream_factor_list[k-1])
                        else:
                            stop = True
                    
                    print("Kausalfaktoren vom Level " + str(m) + " mit ihrer Ordnung:")
                    for p in range(len(level_factor_list_order[m])):
                        if level_factor_list_order[m][p]:
                            print("Ordnung " + str(p))
                            print(level_factor_list_order[m][p])
                        
                    if downstream_factor_list:
                    # zu Testzwecken, downstream_factor_list sollte eigentlich leer sein, wenn nicht gib dessen Elemente aus
                        print("downstream_factor_list auf Level " + str(m) + " enthält noch die folgenden Elemente (sollte eigentlich leer sein):")
                        print(downstream_factor_list)
               
                    
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
                                        
                                    else :
                                        sec_auxiliary_list = [] # zweite Hilfsliste, hier werden alle bekannten "Uebersetzungen" von
                                        # Faktoren der Ebene m durch ihre Konstituenten der Ebene m-1 eingetragen
                                        
                                        for j in level_equiv_list[m] :
                                            if j[1] == fac :                                                
                                                # suche unter den Aequivalenzen der Ebene m nach solchen, die fac kausal aufloesen
                                                for k in get_components_from_formula(j[0],level_factor_list[m]) :
                                                    # fuer jeden kausalen Faktor bzgl. fac pruefe, sammle dessen Konstitutionsbeziehungen
                                                    # in sec_auxiliary_list
                                                    for c in constitution_relation_list :
                                                        if c[1] == k :
                                                            sec_auxiliary_list.append(c)
                                                            
                                                            
                                                            
                                                           

              
                                        for ac in auxiliary_list :
                                            
                                            ter_auxiliary_list = []
                                            for j in level_equiv_list[m] :
                                                for k in get_components_from_formula(j[0], level_factor_list[m]) :
                                                    for sc in sec_auxiliary_list :
                                                        
                                                        sc_transpose = (sc[1],sc[0])
                                                        if substitute_formula(ac,sc_transpose) != "" :
                                                            ele = (substitute_formula(ac,sc_transpose),fac)
                                                            ter_auxiliary_list.append(ele)
                           
                                                for tc in ter_auxiliary_list :
                                                    for k in get_components_from_formula(j[0], level_factor_list[m]) :
                                                        for sc in sec_auxiliary_list :
                                                            if substitute_formula(tc,sc) != "" :
                                                                if substitute_formula(tc,sc) == ac[0] :
                                                                    discard_list.append(ac)
                                 
                                                            
                                            
                                            ter_auxiliary_list.clear()
                                    
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
                    
                    discard_list = [] # Liste der zuloeschenden Formeln
                    
                    
                    for j in range(len(level_equiv_list[m])-1,-1,-1):

                        if get_factor_order(level_equiv_list[m][j][1],level_factor_list_order) > 1:
                            # nur Formeln, mit Ordnung > 1 koennen durch Einsetzungen entstehen
                            
                            # Hilfsliste in die zu pruefende Formeln aufgenommen werden
                            auxiliary_list = []
                            for fac in level_equiv_list[m]:
                                auxiliary_list.append(fac)
                            
                            
                            auxiliary_list.remove(level_equiv_list[m][j])
                            
                            
                            # Vorgehen muss mehrstufig erfolgen um mehrfache Einsetzungen zu ermoeglichen
                            # hoechste Ordnung einer Einsetzung entspricht Ordnung(Formel j)-1
                            # minimale Ordnung ist 1 (jede Formel hat rechts einen Term der Ordnung >0)
                            
                            q = get_factor_order(level_equiv_list[m][j][1],level_factor_list_order) - 1 # Zaehlvariable = Ordnung(Formel j)-1
                            formula_redundant = False
                            while not(formula_redundant) and q > 0:
                                                    
                                for g in auxiliary_list:
                                    # nur Formeln deren rechte Seite gleich jener der zu pruefenden Formel sind, sind von belang
                                    if level_equiv_list[m][j][1] == g[1] :                            
                                        for h in auxiliary_list:
                                            # ist Ordnung von h = q?
                                            if get_factor_order(h[1],level_factor_list_order) == q:
                                                # nur positiver Fall ist sinnvoll
                                                if not((substitute_formula(g,h),g[1]) in auxiliary_list):
                                                    # ist die durch Substitution generierte Formel noch nicht in auxliary_list enthalten?
                                                    if substitute_formula(g,h) == level_equiv_list[m][j][0]:
                                                        # wenn sie gleich der zu pruefenden Formel j ist, dann ist diese redundant
                                                        formula_redundant = True
                                                                                                           
                                                    elif substitute_formula(g,h) != "":
                                                        # ist diese Formel verschieden von Formel j und nicht leer,
                                                        # dann nimm sie in auxiliary_list auf
                                                        ele = (substitute_formula(g,h),g[1])
                                                        auxiliary_list.append(ele)
                                   
                                q = q - 1   # reduziere Zaehlvariable um eins
                            
                            # loesche Inhalte der Hilfsliste
                            auxiliary_list.clear()
                                    
                            if formula_redundant:
                                # setze Formel auf Loeschliste
                                discard_list.append(level_equiv_list[m][j])       
                            
                    for f in discard_list:
                        # loesche als redundant markierte Formeln 
                        level_equiv_list[m].remove(f)       
                        
                    discard_list.clear() 
                            
                    print("minimale Kausalrelationen der Konsitutionsebene " + str(m) + ":")
                    for formula in level_equiv_list[m] : 
                        print(formula[0] + " -> " + formula[1])
                        
    
            # Schleife ueber Konstitutionsebenen endet hier
                
            print("Liste der minimalen Konstitutionsbeziehungen:")
            for c in constitution_relation_list :
                print(c[0] + " -- " + c[1]) 
            
                    
                   
                    
                    
if __name__ == '__main__':
    read_R_file("r_output.txt")
