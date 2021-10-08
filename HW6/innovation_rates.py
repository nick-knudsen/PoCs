import pandas as pd
from collections import Counter
import numpy as np
import spacy
import en_core_web_sm

def ulysses():
    f = open("ulysses.txt",'r')
    words=[]
    counts=[]
    for line in f:
        line=line.split(':')
        words.append(line[0])
        counts.append(int(line[1]))          
    f.close()

    df = pd.DataFrame({'word': words,'count':counts})

    # 5b
    totalNumWords = df['count'].sum()
    uniqueWords = len(counts)
    rho_est = uniqueWords/totalNumWords
    print("Frac: ", rho_est)
    freq = Counter(counts)
    #print("Freq: ", freq)

    # 5c
    n1=freq[1]/sum(freq.values())
    n2=freq[2]/sum(freq.values())
    n3=freq[3]/sum(freq.values())
    print("n_1(g) emp.:{0}   n_2(g) emp.:{1}  n_3(g) emp.:{2}".format(n1,n2,n3))

    n1g_est = 1/(2-rho_est)
    n2g_est = (1-rho_est)/((2-rho_est)*(3-rho_est))
    n3g_est = 2*(1-rho_est)**2/((2-rho_est)*(3-2*rho_est)*(4-3*rho_est))
    print("n_1(g) est: {0} n_2(g) est.: {1} n_3(g) est.: {2}".format(n1g_est, n2g_est, n3g_est))

def pride_and_prej():

    try:
        nlp = spacy.load("en_core_web_sm")
    except:
        nlp = en_core_web_sm.load()

    def clean_text(doc):
        doc = [str(token).lower() for token in doc if token.is_punct !=
            True and token.is_space != True]
        return doc
    nlp = spacy.load('en_core_web_sm')
    nlp.add_pipe(clean_text, name="clean", last=True)
    filename = "pride_prejudice.txt"
    with open(filename) as text:
        words = text.read()
    clean_words = nlp(words)
    clean_words = np.array(clean_words[1:])
    (w_unique, w_counts) = np.unique(clean_words, return_counts=True)
    frac = len(w_counts)/w_counts.sum()
    print(frac)
    freq = Counter(w_counts)
    n1 = freq[1]/sum(freq.values())
    n2 = freq[2]/sum(freq.values())
    n3 = freq[3]/sum(freq.values())
    print(f'{n1:.3f}  {n2:.3f}  {n3:.3f}')

def monte_cristo():
    filename = "monte_cristo.txt"
    nlp = spacy.load('en_core_web_sm')
    with open(filename) as text:
        words = text.read()
    clean_words = nlp(words)
    clean_words = np.array(clean_words[1:])
    (f_unique, f_counts) = np.unique(clean_words, return_counts=True)
    frac = len(f_counts)/f_counts.sum()
    print(frac)
    freq = Counter(f_counts)
    n1 = freq[1]/sum(freq.values())
    n2 = freq[2]/sum(freq.values())
    n3 = freq[3]/sum(freq.values())
    print(f'{n1:.3f}  {n2:.3f}  {n3:.3f}')

pride_and_prej()
monte_cristo()