import re
######################################################################
#               TP de Python - Alignement de Séquence                #
#                                                                    #
#                                 ~                                  #
#                                                                    #
#                Implémentation de l'algorithme de                   #
#                      Needleman et Wunsh                            #
#                                                                    #
#                  Avec et sans pénalité de gap                      #
#                                                                    #
# Version: 2016.04.03                                                #
# Version Python: 3.4.3                                              #
#                                          Auteur: Quentin BONENFANT #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Description:                                                       #
#                                                                    #
# Ce programme contient les fonctions nécessaires pour aligner deux  #
# séquences d'ADN et afficher le résultat.                           #
# Elle est basée sur l'algorithme de Needleman et Wunsch.            #
# Le calcule des scores est fait à partr d'une matrice de score      #
# déclarée en variable globale, et indexé par les nucléotides.       #
# A coté de cette  matrice est déclaré le score de gap par défaut    #
# ainsi que le score de gap ajusté pour les prolongation de gap.     #
# Le but est d'avoir un alignement optimal plus proche de la réalité.#
#                                                                    #
######################################################################

#       Sommaire
#
# 1/ Déclaration et initialisation
# 2/ Implémentation de l'algorithme simple
# 3/ Implémentation de l'algorithme avec gap Affine
# 4/ Outils d'affichage et de manipulation
# 5/ Procédure de test
# 6/ Menu de navigation


                              
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #     #
 ##    #            Déclaration & 
  #   #                       Initialisation 
 ### #    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
######################################################################
# Création d'objets facilement manipulables pour le calcule de score #
# et la procédure de callback (chemin inverse). Cela nous permettra  #
# de savoir d'où vient un node, et si l'alignement précédent était   #
# un gap ou non.                                                     #
######################################################################

class node:
    ''' Noeud d'un graphe de score d'alignement selon l'algorithme
        de Needleman Wunsch. Un noeud est défini par son score
        et une liste de noeuds parents qui peuvent êtres nulles'''
    def __init__(self, s,p):
         
        self.score  = None
        self.parent = None
        if (s!=None):
            self.score = s
        if (p!=None):
            self.parent = p
         
    # Fonction de modificiation du score d'un node
    def setScore(self,n):
        self.score=n
    
    # Fonction de modificiation des parents d'un node
    def setParent(self,p):
        self.parent=p
         
         
class parent:
    ''' Parent d'un noeud. Il s'agit d'objets 'node' associés à
        des alignements (ex: ('A','-') ).'''
        
    def __init__(self, n , al):
         self.nodes = []
         if(n!=None):
            self.nodes.append(n)
            
         self.alignements = []
         if(al!=None):
            self.alignements.append(al)
    
    def addNode(self,n):
        self.nodes.append(n)
    
    def addAlign(self,al):
        self.alignements.append(al)
         
######################################################################
#                                                                    #
#                   Creation des variables globales                  #
#                                                                    #
######################################################################

# Valeurs de gap par défaut
d= -2  # Coût de gap / d'ouverture de gap
k= -1  # Coût de prolongation de gap

# Matrice pour le score de match / mismatch par défaut
costmat= { "A":{"A": 3 ,"T":-1 ,"C":-1 ,"G":-1 },
           "T":{"A":-1 ,"T": 3 ,"C":-1 ,"G":-1 },
           "C":{"A":-1 ,"T":-1 ,"C": 3 ,"G":-1 },
           "G":{"A":-1 ,"T":-1 ,"C":-1 ,"G":3  }}

######################################################################
#                                                                    #
#                   Creation d'une regex de controle                 #
#                                                                    #
######################################################################

# Afin de vérifier que les séquences soient bien de l'ADN et ne fassent
# pas planter le programme, les séquences seront testées avant d'être
# alignée.
DNA= re.compile(r"^[ATCG]+$",re.IGNORECASE)

def isDNA(string):
    if DNA.search(string):
        return(True)
    return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##      #
 #  #    #            Implémentation
   #    #                     de
 ####  #                     l'algorithme simple
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


######################################################################
# Implementation de l'algorithme de Needleman-Wunsch simple          #
######################################################################

def alignementSimple(seq1,seq2,d,cost=costmat):
    # Implementation de l'aglgorythme de Needleman Wunch utilisant une valeur fixe de d et une matrice de score.
    
    # Pré traitement des séquences ( passage en majuscule )
    seq1=seq1.upper()
    seq2=seq2.upper()
    
    # Récupération de la taille de la séquence
    l1= len(seq1)
    l2= len(seq2)
    
    # Création d'une matrice contenant des noeuds "vides".
    # Cela permet au passage de "zeroter" le noeud (0,0).
    matriceAlignement = [[node(0,parent(None,None)) for i in range(l1+1)] for j in range(l2+1)] 
    
    # Initialisation de la matrice (et remplissage des lignes de gap)
    
    # Remplissage de la première ligne
    for i in range(1,l1+1):
        matriceAlignement[0][i]= node( matriceAlignement[0][i-1].score + d, parent( matriceAlignement[0][i-1] , (seq1[i-1],"-") ) )
    
    # Remplissage de la première colone
    for i in range(1,l2+1):
        matriceAlignement[i][0]= node( matriceAlignement[i-1][0].score + d, parent( matriceAlignement[i-1][0] , ("-",seq2[i-1]) ) )
    
    # Remplissage de la matrice
    for j in range(1,l2+1):
        for i in range(1, l1+1):
            
            # Récupération des objets node précédent.
            ant1= matriceAlignement[j-1][i-1] # Diagonal     
            ant2= matriceAlignement[j-1][i]   # Haut
            ant3= matriceAlignement[j][i-1]   # Gauche
        
            # Calcule du meilleur score
            maximum= max(ant1.score + cost[seq1[i-1]][seq2[j-1]], ant2.score + d , ant3.score + d)
            
            # Determination des noeuds parents
            # Ces noeuds sont liés à l'alignement correspondant.
            
            # Création d'un parent vide
            parents= parent(None,None)
            
            if (ant1.score+cost[seq1[i-1]][seq2[j-1]] == maximum):
                parents.addNode(ant1)
                parents.addAlign( (seq1[i-1],seq2[j-1]) )
            
            if(ant2.score + d == maximum):
                parents.addNode(ant2)
                parents.addAlign( ("-",seq2[j-1]) )
                
            if(ant3.score + d == maximum):
                parents.addNode(ant3)
                parents.addAlign( (seq1[i-1],"-") )

            # création du nouveau noeud
            matriceAlignement[j][i]= node(maximum, parents)
    
    # Revoi de la matrice complétée
    return(matriceAlignement)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ###      #
    #    #            Implémentation
  ##    #                     de 
    #  #                   l'algorithme affine
 ###  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################################
# Implementation de l'algorithme de Needleman-Wunsch + gap affine    #
######################################################################

def alignementAffine(seq1,seq2,d,k,cost=costmat):
    
    # Implementation de l'aglgorithme de Needleman Wunsch utilisant
    # une valeur ajustée du score de gap d avec un nombre k tel que d<k<0.
    # Cette méthode utilise la même matrice de score.
    
    
    # INITIALISATION
    
    # Pré traitement des séquences ( passage en majuscule )
    seq1=seq1.upper()
    seq2=seq2.upper()
    
    # Récupération de la taille de la séquence
    l1= len(seq1)
    l2= len(seq2)
    
    # On défini un "infini" suffisamment grand
    infini= 42*(-10**6)
    
    # Création de trois matrices vide contenant des noeuds "vides".
    # Le zerotage des matrices en  0,0 est alors déjà effectué.
    # Cela peut se faire en une seule ligne.
    
    matriceA, matriceB, matriceC = [[[node(0,parent(None,None)) for i in range(l1+1)] for j in range(l2+1)] for nombre_matrice in range(3)]
    
    # On ajoute le score d'ouverture de gap au en position [1][0] de la matrice  B
    matriceB[1][0].setScore(d+k)
    matriceB[1][0].setParent(parent(matriceB[0][0],("-" , seq2[0])))
    
    
    # On ajoute le score d'ouverture de gap au en position [0][1] de la matrice C
    matriceC[0][1].setScore(d+k)
    matriceC[0][1].setParent(parent(matriceC[0][0] , (seq1[0],"-")))
    


    # REMPLISSAGE DES LIGNES

    # remplissage de la première ligne des matrices A et B
    for i in range(1,l1+1):
        matriceA[0][i].setScore(infini)
        matriceB[0][i].setScore(infini)
    
    # remplissage de la première ligne de C
    for i in range(2,l1+1):
        matriceC[0][i].setScore(matriceC[0][i-1].score+k)
        matriceC[0][i].setParent(parent(matriceC[0][i-1],(seq1[i-1],"-")))
    
    # REMPLISSAGE DES COLONES 
    
    # remplissage de la première colone de A et C
    for i in range(1,l2+1):
        matriceA[i][0].setScore(infini)
        matriceC[i][0].setScore(infini)
        
    # Remplissage de la première colonne de B
    for i in range(2,l2+1):
        matriceB[i][0].setScore(matriceB[i-1][0].score+k)
        matriceB[i][0].setParent(parent(matriceB[i-1][0],("-",seq2[i-1])))
    

    # Affichage des matrices initialisée
    # Décommenter la ligne suivante pour afficher.

    #printMatrices(matriceAlignement)                                           



    # REMPLISSAGE DES MATRICES                                 

    for j in range(1,l2+1):
        for i in range(1, l1+1):
            
            # MATRICE A
            # Pour remplire la case "A" du node, il faut
            # tester quel est le(s) meilleur(s) score(s) venant de
            # (i-1, j-1)
            
            # Récupération des objets node d'intérêt.
            
            antA= matriceA[j-1][i-1]
            antB= matriceB[j-1][i-1]
            antC= matriceC[j-1][i-1]
            
            # Recherche du meilleurs score
            maximum= max(antA.score, antB.score, antC.score)
            
            # calcul de la valeur de match / mismatch
            sub=cost[seq1[i-1]][seq2[j-1]]
            
            # Determination des noeuds parents
            # Ces noeuds sont liés à l'alignement correspondant.
            if(j!=1 or i!=1):
                parents= parent(None,None)            
                if (antA.score == maximum):
                    parents.addNode(antA)
                parents.addAlign( (seq1[i-1],seq2[j-1]) )
        
                if(antB.score == maximum):
                    parents.addNode(antB)
                    parents.addAlign( (seq1[i-1],seq2[j-1]) )
                    
                if(antC.score == maximum):
                    parents.addNode(antC)
                    parents.addAlign( (seq1[i-1],seq2[j-1]) )
                    
                # création du nouveau noeud)
                matriceA[j][i].setScore(maximum + sub)
                matriceA[j][i].setParent(parents) 
            
            
            # Cas spéciale de la "première" (1,1) case de A, où les origines 
            # des trois matrices  serais considérée comme parents.
            # Cela triplerais les résultats lors de l'affichage des alignements.
            elif(i==1 and j==1):
                matriceA[1][1]= node(sub, parent(matriceA[0][0],(seq1[0],seq2[0])))
            
            
            # MATRICE B
            
            # Récupération des objets node d'intérêt.
            antA= matriceA[j-1][i]       
            antB= matriceB[j-1][i]       
            antC= matriceC[j-1][i]
            
            # Recherche du meilleurs score
            maximum= max(antA.score + d + k , antB.score + k , antC.score + d + k)
            
            # Determination des noeuds parents
            # Ces noeuds sont liés à l'alignement correspondant.
            parents= parent(None,None)
            
            if (antA.score + d + k == maximum):
                parents.addNode(antA)
                parents.addAlign( ("-",seq2[j-1]) )
            
            if(antB.score + k == maximum):
                parents.addNode(antB)
                parents.addAlign( ("-",seq2[j-1]) )
                
            if(antC.score + d + k== maximum):
                parents.addNode(antC)
                parents.addAlign( ("-",seq2[j-1]) )

            # création du nouveau noeud
            matriceB[j][i].setScore(maximum)
            matriceB[j][i].setParent(parents) 
            
            
            # MATRICE C
            
            # Récupération des objets node d'intérêt.
            antA= matriceA[j][i-1]      
            antB= matriceB[j][i-1]
            antC= matriceC[j][i-1]
            
            # Recherche du meilleurs score
            maximum= max(antA.score + d + k , antB.score + d + k , antC.score + k)
            
            # Determination des noeuds parents
            # Ces noeuds sont liés à l'alignement correspondant.
            parents= parent(None,None)
            
            if (antA.score + d + k == maximum):
                parents.addNode(antA)
                parents.addAlign( (seq1[i-1],"-") )
            
            if(antB.score + d + k == maximum):
                parents.addNode(antB)
                parents.addAlign( (seq1[i-1],"-") )

            if(antC.score + k== maximum):
                parents.addNode(antC)
                parents.addAlign( (seq1[i-1],"-") )

            # création du nouveau noeud
            matriceC[j][i].setScore(maximum)
            matriceC[j][i].setParent(parents) 
            
    # On retourne les trois matrices remplies
    return(matriceA,matriceB,matriceC)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #        #
  #        #       Outils d'affichage 
 #  #     #                     et 
 #####   #                   de manipulation
    #   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################################
#  Méthode récursive de callback, permettant de récupérer            #
#  tout les meileurs alignements de deux séquences. Cette méthode    #
#  s'appuie sur les algorithmes de récupération de la liste des mots #
#  d'un arbre binaire.                                               #
######################################################################

def recallback(node,seq):
    # Condition d'arret: si on a plus de nodes (arrive en  0,0 )
    if(node.parent.nodes==[]):
        return([seq])
    
    # On determine le nombre de parents
    l=len(node.parent.nodes)
    
    #Si on a qu'un seul parent, on rappelle la fonction sur le node parent
    if (l==1):
        return( recallback(node.parent.nodes[0],seq+[node.parent.alignements[0]]) )
      
    # Sinon, on le fait pour chaque parent
    elif (l==2):
        return( recallback(node.parent.nodes[0],seq+[node.parent.alignements[0]])+
                recallback(node.parent.nodes[1],seq+[node.parent.alignements[1]]) )
    
    # Pour les parents multiples, on retourne une liste pour chaques alignements
    # possible sous la forme [ [(X,X),(Y,Y)] , [(X,X),(Y,Y)] ]
    elif (l==3):    
        return( recallback(node.parent.nodes[0],seq+[node.parent.alignements[0]])+
                recallback(node.parent.nodes[1],seq+[node.parent.alignements[1]])+
                recallback(node.parent.nodes[2],seq+[node.parent.alignements[2]]) )


######################################################################
#                                                                    #
#             Affichage des matrices et des alignements              #
#                                                                    #
######################################################################

# Fonction d'affichage des matrices
# Prend en paramètre un tableau (longueur max=3)
# de matrices contenant des objets de type "node"
# L'affichage se fait sur des colonnes de 10 caractères de largeur.

def printMatrices(mat):
    for n,m in enumerate (mat):
        print("Matrice "+ "ABC"[n])
        for line in m:
            l=""
            for objects in line:
                l+="% 10d" % objects.score
                
            print(l)
        print("\n")
        
# Fonction d'affichage "propre" des différents alignements de séquences.

def printSequence (couples):
    for al in couples:
        al.reverse()
        sequences=list(zip(*al))
        print("".join(sequences[0]))
        print("".join(sequences[1]))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        
        
# Procédure de recherche et d'affichage de l'alignement optimal

def lancerAlignementSimple(s1="42",s2="42"):
    global costmat,d
    print("Alignement simple")
    while (not isDNA(s1) or not isDNA(s2)):
        s1=input("Entrez la première séquence: \n")
        s2=input("Entrez la seconde séquence: \n")

    print("\nLes séquences qui seront alignées sont:\n")
    print("Séquence 1: "+s1+"\nSéquence 2: "+s2+"\n")
    print(" Score de gap d = " +str(d))
    print("\n Matrice d'alignement\n")
    resultat=alignementSimple(s1,s2,d,costmat)
    for line in resultat:
        l=""
        for objects in line:
            l+="% 5d" % objects.score
        print(l)
    print("\n")
    print("Le meilleur score d'alignement est: "+str(resultat[-1][-1].score))
    # on lance le retour sur trace depuis le coin en bas à droite
    printSequence(recallback(resultat[-1][-1],[]))


# Procédure de recherche et d'affichage de l'alignement optimal
# pour la méthode avec gap affine
# Par défaut, les valeurs s1 et s2 sont invalide.
# Si on ne passe pas de paramètre, un couple de séquence
# valide sera demandé.

def lancerAlignementAffine(s1="42",s2="42"):
    # les couts sont des valeurs globales
    global costmat,d,k 
    
    print(" Alignement Affine")
    while (not isDNA(s1) or not isDNA(s2)):
        s1=input("Entrez la première séquence: \n")
        s2=input("Entrez la seconde séquence: \n")
    
    print("\nLes séquences qui seront alignées sont:\n")
    print("Séquence 1: "+s1+"\nSéquence 2: "+s2+"\n")
    print(" Score de gap d = " +str(d))
    print(" et k = " +str(k))
    print("\n Matrices d'alignement\n")
    
    resultat=alignementAffine(s1,s2,d,k,costmat)
    printMatrices(resultat)
    maximum=max( [ resultat[0][-1][-1].score,
                resultat[1][-1][-1].score,
                resultat[2][-1][-1].score ] )
    print("Le meilleur score est: "+str(maximum))
    # on lance le retour sur trace à partir du meilleur
    # node terminal (test entre chaques matrice)
    depart=[]
    for i in range(3):
        if(resultat[i][-1][-1].score==maximum):
            depart.append(resultat[i][-1][-1])
    # parcours depuis les meilleurs noeuds
    for node in depart:
        printSequence(recallback(node,[]))
        

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ####      #
  #        #           Procédure
  ###     #                      de
     #   #                          Test
  ###   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def lancerTest():
    print("##################################################################")
    print("#                                                                #")
    print("#                       -=~ TESTS ~=-                            #")
    print("#                                                                #")
    print("##################################################################")
    
    # Les première séquences testée viennent de l'exercice 6 du TD
    # Ici, ce sont les séquences 1 et 2.
    # Le résultat est supposé être un alignement simple avec un score finale de 3
    # et l'alignement suivant:
    #   TACGATGA
    #   TCCGAT-A
    print(" TEST N° 1")
    print("Test des séquences S1 et S2 de l'exercice 6 du TD, NW simple")
    print("Alignement attendu: ")
    print("TACGATGA")
    print("TCCGAT-A")
    print("Test: \n")
    
    s1="TACGATGA"
    s2="TCCGATA"
    lancerAlignementSimple(s1,s2)

    print(" \nTEST N° 2")
    print("Test des séquences S2 et S3 de l'exercice 6 du TD, NW simple")
    print("De multiples alignements doivent être trouvés")
    print("\n")
    # Sequences 3 et 4 du TD, comprenant des alignements optimaux multiples
    # Le programme se chargera de retrouver ces alignements et de les afficher
    s3="ACGACGA"
    lancerAlignementSimple(s2,s3)
    
    print(" \nTEST N° 3")
    print("Test sur l'alignement gap affine")
    print("Séquence TTATT vs TT, doit donner soit:")
    print("TT---, soit ---TT, soit T---T ")
    print("C'est à dire toujours le gap le plus long possible.")
    print("L'algorithme classique donnerai plus de réponses (mais moins probables).")
    print("\n")
    s1="TTATT"
    s2="TT"
    lancerAlignementAffine(s1,s2)
    
    print(" \nTEST N° 4")
    print("Test des deux algorithme sur une séquence comportant la possibilité")
    print("d'aligne soit un grand gap, soit de le couper en deux.")
    print("Ce test doit donner un seul alignement possible pour la version affine")
    print("et deux alignements pour la version simple (sans gestion du gap) ")
    
    print("\n")
    s1="ATGTGACGA"
    s2="ATACGA"
    lancerAlignementAffine(s1,s2)
    
    print("\n")
    lancerAlignementSimple(s1,s2)
    
    print("\n")
    print(r"/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!\ ")
    print("ATTENTION, le test génère beacoup de texte, remontez bien jusqu'en haut")
    print(r"/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!\ ")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #       #
   #       #           
  ###     #                   Menu
 #   #   #       
  ###   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("##################################################################")
print("#                                                                #")
print("#                      Needleman Wunsch                          #")
print("#                                                                #")
print("##################################################################")
print(""" Ce programme vous permettra d'aligner deux séquences d'ADN 
          en utilisant l'algorithme de Needleman et Wunsch. Deux versions 
          vous seront proposés ici: avec gap fixe, ou gap affine """)

# On stock les objets fonctions dans un dictionnaire
# Il n'y a pas de switch/case en python.
routine={'1':lancerAlignementSimple, '2':lancerAlignementAffine, '3':lancerTest}

# on continue tant qu'on a pas de demande pour quitter
continuer = True
while(continuer):
    print("\n Veuillez selecionner l'algorithme que vous souhaitez utiliser")
    print(" ou lancer des alignements de test")
    choix="-1"
    # véfification des inputs
    while(choix not in ['1','2','3','0']):
        print("1- Classique   2- Gap Affine  3- Tests ou 0- Quitter")
        choix=input()
        
    # on quitte si l'utilisateur le demande
    if (choix=='0'):
        continuer=False
        break
    # on execute la fonction appropriée
    routine[choix]()

# Fin
###################################################################################B
#                           TGUgcHl0aG9uIGMnZXN0IHN5bXBh                           #
#6#################################################################################4