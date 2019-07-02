from __future__ import division
import numpy as np
from os.path import expanduser
import sys
from scipy import spatial, stats
from biom import load_table

#file_name = 'GPC_OTU_Table_97sim.biom'
file_name = 'emp_deblur_100bp.release1.biom'
#file_name = 'Global.Global2000Subset.Bacteria.EMP.biom'


mydir = expanduser("~/GitHub/DeepLife_16S")
path = mydir + '/data/'+file_name

import metrics


dat = load_table(path) 
ids1 = dat.ids().tolist()


sads = []
ids = []
for i in ids1:
    try:
        sad = dat.data(i)
        sads.append(sad)
        ids.append(i)
    except:
        continue
    
    
OUT = open(mydir + '/diversity_data/GPC_OTU_Table_97sim_Sample_Div_Measures.txt','w+')

print>>OUT, 'name N S Var Evar ESimp EQ O ENee EPielou EHeip BP SimpDom Nmax McN skew logmodskew chao1 ace jknife1 jknife2 margalef menhinick preston_a preston_S'

ct = 0
numRADs = len(sads)
for ind, sad in enumerate(sads):

    name = ids[ind]
    RAD = list([x for x in sad if x > 0])
    RAD = [int(i) for i in RAD]
    RAD.sort(reverse=True)

    N = sum(RAD)
    S = len(RAD)

    if S < 5: continue
    if max(RAD) == min(RAD): continue

    # Evenness/Unevenness
    Var = np.var(RAD, ddof = 1)
    Evar = metrics.e_var(RAD)
    ESimp = metrics.e_simpson(RAD)
    EQ = metrics.EQ(RAD)
    O = metrics.OE(RAD)
    
    #Camargo = 0.0 # metrics.camargo(RAD)   # Takes too long
    ENee = metrics.NHC(RAD)
    EPielou = metrics.e_pielou(RAD)
    EHeip = metrics.e_heip(RAD)


    # Dominance
    BP = metrics.Berger_Parker(RAD)
    SimpDom = metrics.simpsons_dom(RAD)
    Nmax = max(RAD)
    McN = metrics.McNaughton(RAD)

    # Rarity
    
    logmodskew = metrics.logmodskew(RAD)
    skew = stats.skew(sad)
    #p_ones = metrics.r_singletons(RAD)
    #p_zpt1 = metrics.p_ZPtOne(RAD)

    # Preston's alpha and richness, from Curtis and Sloan (2002).
    # Estimating prokaryotic diversity and its limits. PNAS.
    preston_a, preston_S = metrics.Preston(RAD)

    # Richness estimators
    chao1, ace, jknife1, jknife2 = metrics.EstimateS1(RAD)
    margalef = metrics.Margalef(RAD)
    menhinick = metrics.Menhinick(RAD)

    ct+=1

    print>>OUT, name, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logmodskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S
    

    print numRADs - ct
    
OUT.close()
