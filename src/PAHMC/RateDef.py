import numpy as np
import os, re, sys
import logging

logger = logging.getLogger(__name__)

class Rates():

    def __init__(self, rate_def_file, rates_directory):
    
        self.reactionrates = {}
        self.dE = {}
        linecounter = 0

        # Try to open file, raise an error if file was not found
        try:
            file = open(rate_def_file, 'r')
        except FileNotFoundError:
            logger.error('The provided rate definition file \''+rate_def_file+'\'  was not found.')
            sys.exit(2)
            
        if not os.path.isdir(rates_directory):
            logger.error('The provided rates directory \''+rates_directory+'\'  was not found.')
            sys.exit(2)

        for line in file:
            # Remove unnecessary whitespace from the end of each line
            line = line.rstrip()
            
            linecounter += 1

            # Ignore lines that are empty or start with # (using regular expression)
            if re.match(r'^#|^\s*$', line):
                continue
                
            # Split line into the reactions and ratefile names
            reactions, ratefile = line.split('\t')
            
            # Make the path of the file
            filepath = os.path.join(rates_directory, ratefile)
            
            # Check if file exists
            if not os.path.isfile(filepath):
                logger.error('The provided rate file \''+filepath+'\'  was not found.')
                sys.exit(2)
            
            # Read in rates
            rates = np.loadtxt(filepath,unpack=True,skiprows=2)
            
            reactions = reactions.split(',')
            
            
            
            for key in reactions:
                
                
                # Remove spaces if they are present
                reac = key.strip()
                
                
                self.reactionrates[reac] = rates
            
                with open(filepath) as f:
                    firstline = next(f) # Read in first line of the file
                    # Then select the value of Delta out of the rate file by regular expression
                    # (this regex finds all numbers in the first line, then saves only the last one)
                    delta = float(re.findall(r"[-+]?\d+[\.]?\d*[eE]?[-+]?\d*\b",firstline)[-1])
            
                    self.dE[reac]= delta
        
        N_rates = len(self.reactionrates.keys())
        
        logger.info(str(N_rates)+' rate files sucessfully read.')
