import numpy
import random
from string import ascii_uppercase
from itertools import cycle

population = ['A']

# the rich get richer model
# takes:
#       rho, the probability of innovation at any one time step
#       len, how long to run the model for
def simons_model(rho, len):
    # naming of groups
    name = 'B'
    repeats = 0
    # time steps
    for t in range(len):
        # innovation with prob rho
        if random.random() <= rho:
            population.append(name)
            name = name_next_group(name)
        else:
            population.append(random.choice(population))
            
# get the next string in the sequence ...ABX, ABY, ABZ, ACA, ACB...
def name_next_group(curr_name):
    all_zs = True
    new_name = list(curr_name)
    # loop backwards over the current string
    for i, letter in enumerate(new_name[::-1]):
        # increment first non-Z letter by 1
        if letter != "Z":
            new_name[len(new_name)-(i+1)] = chr(ord(letter)+1)
            all_zs = False
            break
        # if we find a Z, turn it over to A
        else:
            new_name[len(new_name)-(i+1)] = "A"
    # if the whole string was Zs, the new string is all As with length one greater
    if all_zs:
        new_name = "A"*(len(curr_name)+1)
    return "".join(new_name)