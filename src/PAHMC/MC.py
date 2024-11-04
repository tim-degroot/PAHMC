# Contains all MC functions

import numpy as np
import random
import sys
import logging

from PAHMC.Functions import Possible_reactions
import PAHMC.Reactions as React
from PAHMC.Functions import Choose_reaction
from PAHMC.Input import Input_reader
from PAHMC.RateDef import Rates
from PAHMC.Molecule import Molecule
from PAHMC.Ratecheck import Check_available_rates
from PAHMC.Output import STD_Output

logger = logging.getLogger(__name__)

def Single_MC(E, max_time, molecule, rates):
    """Run a single MC"""
    
    # Make a list of all reaction keys that have rates specified
    specified_rates = list(rates.reactionrates.keys())
    
    energy = E
    time = 0
    total_hops = 0
    diss_atom = None
    diss_position = None
    
    while time < max_time:
        
        # Determine the possible reactions (and aliphatic site surroundings)
        molecule.possible_reactions = Possible_reactions(molecule, specified_rates)
        
        # Choose reaction from the possibilities
        reactionkey, dt = Choose_reaction(E, molecule.possible_reactions, rates)
        
        # Error if there's no dissociation before the rates go to 0
        if reactionkey == None:
            logger.error('Ran out of reaction rates before dissociation.')
            break
        
        # Update time
        time += dt
        
        # Update energy
        E -= rates.dE[reactionkey]
        

#        print(reactionkey)
#        if total_hops%1000 == 0:
#            print(time)

        # Carry out the reaction
        if 'to' in reactionkey:
            # First some bookkeeping (saving the position of the aliphatic site and such)
            # Copy key and remove the atom to move (first character of string)
            key = reactionkey.replace(reactionkey[0],'')
            # Split the remaining key into the current site number, and the next site number
            current = key.split('to')[0]
            molecule.positions.append(current)
            molecule.position_times[current] += dt
        
            # Do the scrambling
            React.Do_scramble(reactionkey, molecule)
        elif 'diss' in reactionkey:
            # Save the atom that dissociates
            diss_atom = reactionkey[0]
            # Save the site the dissociation happened
            diss_position = reactionkey.replace(reactionkey[0], '')
            diss_position = diss_position.replace('diss', '')
            
            # Do the dissociation
            React.Do_dissociation(diss_atom, molecule)
            #print(total_hops, diss_atom, diss_position)
            break
       

        # Keep track of how many 'hops' are done by the aliphatic site
        total_hops +=1
        
    
    return diss_atom, diss_position, time, total_hops



def Do_MC(inputfile, ratedef_file, rate_dir, outputfile):
    """Function to perform multiple Monte Carlo simulations"""
    
    input = Input_reader(inputfile)
    warn_setting = input.handling
    
    rates = Rates(ratedef_file, rate_dir)
    
    molecule = Molecule(input)
    
    Check_available_rates(rates, molecule, warn_setting)
    
    # Initialize some variables, arrays, and dictionaries
    
    if len(input.energy_range) == 1:
        Energy = input.energy_range
    else:
        # Make a linear spaced array of energies in range specified (nsteps +1 because linspace includes the stop value specified.)
        Energy = np.linspace(input.range[0], input.range[1], num=input.range[2]+1)
    
    dissociation_atoms = {}
    dissociation_times = {}
    dissociation_positions = {}
    N_scramble_hops = {}
    
            
    # Run multiple MC events
    for k in range(len(Energy)):

        dissociation_atoms[Energy[k]] = {'H': 0, 'D': 0, 'None': 0}
        N_scramble_hops[Energy[k]] = []
        dissociation_times[Energy[k]] = []
        dissociation_positions[Energy[k]] = {}
    
        # Initialize the dissociation positions dictionary to zero
        for i in range(len(molecule.edge_numbers)):
            for e_num in molecule.edge_numbers[i]:
                dissociation_positions[Energy[k]][e_num] = 0
            for l_num in molecule.link_numbers[i]:
                dissociation_positions[Energy[k]][l_num] = 0
    
    
    
        for iter in range(input.iterations):
            logger.info('Running iteration '+str(iter+1)+' out of '+str(input.iterations)+'...')
            
            diss_atom, diss_position, time, hops = Single_MC(Energy[k], input.t_max, molecule, rates)

            if diss_atom == None:
                dissociation_atoms[Energy[k]]['None'] +=1
            else:
                dissociation_atoms[Energy[k]][diss_atom] += 1
                dissociation_positions[Energy[k]][diss_position] += 1
                dissociation_times[Energy[k]].append(time)
            N_scramble_hops[Energy[k]].append(hops)
            
            
#            print(diss_atom, diss_position, time, hops)
            # Reset the molecule for the next event
            molecule = Molecule(input)
            
            
    
    STD_Output(outputfile, dissociation_atoms, dissociation_positions, dissociation_times, N_scramble_hops)
    #output        
