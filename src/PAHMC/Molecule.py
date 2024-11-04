# Import relevant libraries
import numpy as np
import copy as c

# Import files from PAH-MC
import PAHMC.Functions as Functions

class Molecule():

    def __init__(self, input):
        # Make copies 
        self.edges = c.deepcopy(input.mol_edge)
        self.edge_numbers = c.deepcopy(input.mol_edge_numbers)
        self.links = c.deepcopy(input.mol_links)
        self.link_numbers = c.deepcopy(input.mol_links_numbers)
        self.Deuterium = False
    
        self.index, self.al_deuterium = Functions.Find_aliphatic_site(self.edges)
        
        self.al_place = 'e' # for now we can only start with aliphatic site on edge, not on link
        
        
        if self.al_deuterium:
            self.Deuterium = True
        
        else:
            for i in range(len(self.edges)):
                for j in range(len(self.edges[i])):
                    if 'D' in self.edges[i][j]:
                        self.Deuterium = True
        
        self.possible_reactions = [None] 
        self.al_type = [None] 
        self.HH_time = 0
        self.HD_time = 0
        self.DD_time = 0
        self.positions = []
        
#        for m in range(len(self.edge_numbers)):
#            for edgenumber in self.edge_numbers[m]:
#                self.position_times[edgenumber] = 0
#            for linknumber in self.link_numbers[m]:
#                self.position_times[linknumber] = 0
#
