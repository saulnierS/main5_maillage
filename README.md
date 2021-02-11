# main5_maillage

Ce projet est implémenté avec python3, il utilise principalement gmsh, matplotlib, scipy, numpy

#Fichiers

omega:

	appartement.msh : généré par gmsh au moment de faire le maillage avec la fonction dans appartement.py  
	resolution.py : permet d'appeler toutes les fonctions nécéssaires à la résolution du problème avec la méthode des élement finis.

	omega.msh : généré par gmsh au moment de faire le maillage avec la fonction dans omega.py  
	resolution.py : permet d'appeler toutes les fonctions nécéssaires à la résolution du problème avec la méthode des élement finis.

	appartement.py : crée et maille un appartement (qui correspond à l'espace omega)

	omega.py : crée et maille un carré simple (qui correspond à l'espace omega) (permet de tester)

maillage:

	maillage.py : permet de retranscrire les informations gmsh dans une classe plus simple à utiliser par la suite

solveur:

	common.py : classe triplets pour le format des matrices creuses

	fem_p1.py : permet de calculer la matrice de masse, de rigidité, les points de quatratures et d'appliquer la condition de dirichlet

	resolution.py: appelle l'ensemble des fonctions pour résoudre le problème avec les éléments finis

resultat:
	
	U.png: solution après la résolution

#Execution:

	dans un terminal, entrez la commande : ``python3 resolution.py``
