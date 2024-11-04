# Import relevant libraries
import os, re, sys
import numpy as np
import logging




logger = logging.getLogger(__name__)

# Import other relevant files
# (none currently)

# Create a class object to contain the parameters of the simulation
class Input_reader():
    """Takes the input file and initializes a parameter class for the simulation to be run"""
    
    def __init__(self, filename):
    

    
        # Try to open file, raise an error if file was not found
        try:
            file = open(filename, 'r')
        except FileNotFoundError:
            logger.error('The provided input file \''+filename+'\'  was not found.')
            sys.exit(2)
        
        # Define some counters to keep track of reading the input
        linecounter = 0 # keep track of all lines of file
        readingcounter = 0 # keep track of lines actually containing input
        ratefilecounter = 0
        # Declare some variables that are needed
        self.reactions = [] # link reaction keys with rate filenames

        
        for line in file:
            # Remove unnecessary whitespace from the end of each line
            line = line.rstrip()
            
            linecounter += 1
            
            # Ignore lines that are empty or start with # (using regular expression)
            if re.match(r'^#|^\s*$', line):
                continue
            
            if readingcounter == 0:
                # Read in the molecule/simulation name
                self.mol_name = line
                
            elif readingcounter == 1:
                # Read in and save the starting/initial edge structure
                self.mol_edge = line.split()
                for i, edge in enumerate(self.mol_edge):
                    edge = edge[1:-1].split(',') 
                    if len(edge) not in (1,2,3,4):
                        logger.error('Individual edges can only have up to 4 components. Line: '+ str(linecounter)+ '.')
                        sys.exit(2)
                
                    for j, atoms in enumerate(edge):
                        
                        # need to check for correct atom weights in edges
                        #print(atoms)
                        if atoms not in ('H','D','HH', 'HD', 'DD'):
                            logger.error('Incorrect atoms defined, only no atom (0), hydrogen (H), deuterium (D), or aliphatic groups (HH, HD, DD) are currently supported. Line: '+ str(linecounter)+ '.')
                            sys.exit(2)
                            
                        
                        # Assemble the seperate parts back together into a filled edge
                        edge[j] = atoms
                    # Edges to use for the simulation
                    self.mol_edge[i] = edge
                
                # keep a record of the initial edge structure just in case    
                self.init_edge = [0] * len(self.mol_edge)
               
                for n in range(0,len(self.init_edge)):
                    self.init_edge[n] = tuple(self.mol_edge[n]) 
                
                self.init_edge = tuple(self.init_edge)

            elif readingcounter == 2:
                # Read in and save the starting/initial edge structure
                self.mol_edge_numbers = line.split()
                for i, edge in enumerate(self.mol_edge_numbers):
                    edge = edge[1:-1].split(',') 
                    if len(edge) is not len(self.mol_edge[i]):
                        logger.error('Edges and edge numbering need to have the same shape. Line: '+ str(linecounter)+ '.')
                        sys.exit(2)
                
                    for j, number in enumerate(edge):
                        
                        # Assemble the seperate parts back together into integer filled edge
                        edge[j] = number
                    # Edges to use for the simulation
                    self.mol_edge_numbers[i] = edge
                    
                if len(self.mol_edge_numbers) is not len(self.mol_edge):
                    logger.error('Edges and edge numbering need to have the same size. Line: '+ str(linecounter)+ '.')
                    sys.exit(2)
                
                # keep a record of the initial edge structure just in case    
                self.init_edge_numbers = [0] * len(self.mol_edge_numbers)
               
                for n in range(0,len(self.init_edge_numbers)):
                    self.init_edge_numbers[n] = tuple(self.mol_edge_numbers[n])
                
                self.init_edge_numbers = tuple(self.init_edge_numbers)

            elif readingcounter == 3:
                # Read in and save the initial connections between the edges
                self.mol_links = line.split()
                for i, link in enumerate(self.mol_links):
                    link = link[1:-1].split(',')
                    
                    if len(link) > 2:
                        logger.error('Individual links can only have up to two components currently. Line: '+ str(linecounter)+ '.')
                        sys.exit(2)
                    for j, link_ in enumerate(link):
                        if link_ != '0':
                            logger.error('Incorrect link specified, only empty tertiary carbons (0) are accepted. Line: '+ str(linecounter)+ '.')
                            sys.exit(2)
                        link[j] = link_
                    self.mol_links[i] = link
                if len(self.mol_links) != len(self.mol_edge):
                    logger.error('Edges and links need to have the same size. Line: '+ str(linecounter)+ '.')
                    sys.exit(2)
                
                self.init_links = tuple(self.mol_links)

            elif readingcounter == 4:
                # Read in and save the initial connections between the edges
                self.mol_links_numbers = line.split()
                for i, link in enumerate(self.mol_links_numbers):
                    link = link[1:-1].split(',')
                    
                    if len(link) is not len(self.mol_links[i]):
                        logger.error('Links and link numbering need to have the same shape. Line: '+ str(linecounter)+ '.')
                        sys.exit(2)
                    
                    
                    for j, link_number in enumerate(link):
                        link[j] = link_number
                    self.mol_links_numbers[i] = link
                if len(self.mol_links_numbers) is not len(self.mol_links):
                    logger.error('Links and link numbering need to have the same size. Line: '+ str(linecounter)+ '.')
                    sys.exit(2)
                
                self.init_links_numbers = tuple(self.mol_links_numbers)
                    
            elif readingcounter == 5:
                # Read in the range of the simulation, and convert to float/integer
                self.energy_range = line.split(',')
                if len(self.energy_range) == 1:
                    try:
                        self.energy_range = [float(self.energy_range[0])]
                    except ValueError:
                        logger.error('The single energy should be a floating point value. Line: '+str(linecounter)+'.')
                        sys.exit(2)
                elif len(self.energy_range) == 3:
                    try:
                        self.energy_range[0] = float(self.energy_range[0])
                        self.energy_range[1] = float(self.energy_range[1])
                        self.energy_range[2] = int(self.energy_range[2])
                    except ValueError:
                        logger.error('Cannot understand the simulation range given. Line: '+str(linecounter)+'.')
                        sys.exit(2)
                else:
                    logger.error('Cannot understand the simulation range given, it should be either a single energy or an energy range (min, max, nsteps). Line: '+str(linecounter)+'.')
                
            elif readingcounter == 6:
                # Read in the number of simulations for each energy
                try:
                    self.iterations = int(line)
                except ValueError:
                    logger.error('Invalid number of iterations given. Line: '+ str(linecounter)+ '.')
                    sys.exit(2)
                            
            elif readingcounter == 7:
                try:
                    self.t_max = float(line)
                except ValueError:
                    logger.error('Invalid maximum time given. Line: '+ str(linecounter)+ '.')
                    sys.exit(2)
            
            elif readingcounter == 8:
                self.handling = line.lstrip()
                if self.handling not in ('w','q'):
                    logger.error('Invalid error handling mode given. Supported are: do nothing (n), display warnings (w), display warnings and abort (q). Line: '+ str(linecounter)+ '.')
                    sys.exit(2)
            
#            elif readingcounter >=7:
#                # read in the reaction keys with corresponding rate filenames 
#                ratefilecounter += 1
#                reac = line.split()
#                self.reactions.append(reac)
#                
#                if len(reac) ==1:
#                    raise AttributeError('No rate file specified for reaction '+reac+'. Line: '+ str(linecounter)+ '.')
#                elif len(reac) > 2:
#                    raise AttributeError('Can only take one rate file per reaction. Line: '+ str(linecounter)+ '.')
#                if reac[0] not in Reactions.implemented_reactions:
#                    raise AttributeError('Invalid reaction key given. Line: '+ str(linecounter)+ '.')
                    
            readingcounter += 1
#        print('Number of reactions with rate files specified: '+str(ratefilecounter))
        
        file.close()    

