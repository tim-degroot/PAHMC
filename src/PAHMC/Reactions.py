import numpy as np
import logging

logger = logging.getLogger(__name__)
  

def Do_scramble(reactionkey, molecule):
    # Get the index of the aliphatic site
    if len(molecule.index[0]) == 1:
        i = molecule.index[0][0]
        j = 0

    elif len(molecule.index[0]) == 2:
        i,j = molecule.index[0]
    
    # Convert the reaction key into something usable
    
    # Save the atom to move
    atom = reactionkey[0] 
    # Copy key and remove the atom to move (first character of string)
    key = reactionkey.replace(reactionkey[0],'')
    # Split the remaining key into the current site number, and the next site number
    current, next = key.split('to')
    
    # Remove the atom from the current site:
    if molecule.al_place == 'e':
        molecule.edges[i][j] = molecule.edges[i][j].replace(atom, '', 1)
    if molecule.al_place == 'l':
        molecule.links[i][j] = '0'
    
    
    # If next site is on the edge add the atom to that site
    for k in range(len(molecule.edges)):
        for l in range(len(molecule.edges[k])):
            if molecule.edge_numbers[k][l] == next:
                molecule.edges[k][l] += atom
                molecule.index[0] = [k,l]
                x,y = k,l
                molecule.al_place = 'e'

        
    # If next site is on the links add the atom to that site
    for m in range(len(molecule.links)):
        for n in range(len(molecule.links[m])):
            if molecule.link_numbers[m][n] == next:
                molecule.links[m][n] = atom
                molecule.index[0] = [m,n]
                x,y = k,l
                molecule.al_place = 'l'
                

    # Check if the new aliphatic site has a deuterium
    # if moving D, of course yes
    if atom == 'D':
        molecule.al_deuterium = True
    # if not moving D, depends on the next site
    elif atom == 'H':
        if molecule.al_place == 'l':
            molecule.al_deuterium = False
        elif molecule.al_place == 'e':
            if 'D' not in molecule.edges[x][y]:
                molecule.al_deuterium = False
            elif 'D' in molecule.edges[x][y]:
                molecule.al_deuterium = True
    

def Do_dissociation(atom, molecule):

    if len(molecule.index[0]) == 1:
        i = molecule.index[0][0]
        j = 0

    elif len(molecule.index[0]) == 2:
        i,j = molecule.index[0]

    # Remove the dissociating atom
    if molecule.al_place == 'e':
        molecule.edges[i][j] = molecule.edges[i][j].replace(atom, '', 1)
            
    if molecule.al_place == 'l':
        molecule.links[i][j] = '0'

    molecule.index = None
    molecule.al_place = None
