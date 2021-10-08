import random
import numpy as np
from collections import Counter
from scipy.stats import linregress
import matplotlib.pyplot as plt

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

# the rich get richer model
# takes:
#       rho, the probability of innovation at any one time step
#       steps, how long to run the model for
def simons_model(rho, steps):
    population = [1]
    groups = 1
    # time steps
    for t in range(steps):
        # innovation with prob rho
        if random.random() <= rho:
            groups += 1
            population.append(groups)
        else:
            population.append(random.choice(population))
    return population
            
rhos = [0.1, 0.01, 0.001]
for rho in rhos:
    sims = 100
    steps = 100000
    results = []
    # run the simulation sims times and aggregate results
    for j in range(sims):
        results.extend(simons_model(rho, steps))
    
    # get counts for each group
    avg_counts = ({k: v/10 for k, v in dict(Counter(results)).items()})
    rank = list(avg_counts.keys())
    freq = list(avg_counts.values())

    mod = linregress(np.log10(rank), np.log10(freq))
    beta = -mod.slope
    alpha = 1 - rho

    plt.figure()
    plt.plot(np.log10(rank), np.log10(freq))
    plt.plot(np.log10(rank), mod.intercept+np.log10(rank)*mod.slope)
    plt.title("rho = "+ str(rho) + " beta = " + str(round(beta, 3)))
    plt.xlabel("log group number")
    plt.ylabel("log group size")
    plt.savefig("zipfian_simon_model_rho_"+str(rho).replace('.','_')+".png")