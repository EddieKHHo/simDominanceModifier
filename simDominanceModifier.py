##########----------Import packages
import os, sys, time, argparse
import numpy as np
import pandas as pd
import random
from multiprocessing import Pool

##########----------calc fractional occupany for haploid
def calcFO_Haploid(conc1, m1, k):
    numer = conc1 / (k ** m1)
    denom = 1 + conc1 / (k ** m1)
    return (numer / denom)

#########----------calc fractional occupancy
def calcFO(conc1, conc2, m1, m2, k):
    numer = conc1 / (k ** m1)
    denom = 1 + conc1 / (k ** m1) + conc2 / (k ** m2)
    return (numer / denom)

##########----------calc expression of B based on allele at A
##########  if A is 0 it reduce expression, if A is 1 it enhances expression
def calcExpB(A, conc1, conc2, m1, m2, k):
    FO = calcFO(conc1, conc2, m1, m2, k)
    if A==0:
        return (-1 * FO)
    else:
        return (FO)

#########----------return gt given two alleles
#########   assume values are 0 or 1 only
def calcGt(allele1, allele2):
        if allele1 == 0 and allele2 == 0:
            return 0
        elif allele1 == 1 and allele2 == 1:
            return 2
        else:
            return 1

##########----------return binding affinity for alpha based on sex
def calcMAlpha(SEX, ALPHA):
    if SEX=='m':
        return ALPHA
    else:
        return round(1-ALPHA, 1)

##########----------return binding affinity for beta based allele of A
def calcMBeta(A, BETA):
    if A==0:
        return BETA
    else:
        return round(1-BETA, 1)

##########----------calc absolute fitness based on sex, standardize concentration of B and selection regime
##########  SR = 0 for sexual antagonism
##########  SR = 1 for sexually concordant selection using 'female' fitness eq
##########  SR = 2 for sexually concordant selection using 'male' fitness eq

def calcAbsFitness(SEX, STCONC, SR, S_F, S_M, GAMMA_F, GAMMA_M):
    if SR == 0:
        if SEX == 'm':
            return (1 - (S_M * STCONC) ** GAMMA_M)
        else:
            return (1 - (S_F * (1 - STCONC) ** GAMMA_F))
    elif SR == 1:
        return (1 - (S_F * (1 - STCONC) ** GAMMA_F))
    elif SR == 2:
        return (1 - (S_M * STCONC) ** GAMMA_M)
'''
def calcAbsFitness(SEX, STCONC, SR, S_F, S_M, GAMMA_F, GAMMA_M):
    if SR == 0:
        if SEX == 'm':
            return 0.9
        else:
            return 0.1
    elif SR == 1:
        return 0.25
    elif SR == 2:
        return 0.75
'''
##########----------mutate binding values (alpha, beta)
##########   increments of 0.1 and value bound between 0 and 1
def mutateBindValue(mu, bValue):
    rBinom = np.random.binomial(1, mu)
    if rBinom == 1:
        newBValue = -1
        while (newBValue < 0 or newBValue > 1):
            epsilon = random.choices([-0.1, 0.1], k=1)
            newBValue = round(bValue + epsilon[0], 1)
    else:
        newBValue = bValue
    return (newBValue)

###########----------generate gamete through selection on parents
def selectionGamete(listIndiv, maxFit, popSize):
    #####-----while loop stop when get enough gametes
    listID, nGam = [], 0
    lGamAlpha, lGamA, lGamBeta = [], [], []
    while nGam < popSize:
        #####-----choose random parent and get relative fitness
        pID = int(random.randint(0, len(listIndiv) - 1))
        pFit = listIndiv[pID].getAbsFitness()
        relFit = pFit/maxFit
        #####-----selection by testing relFit against random uniform variate
        rVar = np.random.uniform(0, 1, 1)
        if relFit >= rVar:
            nGam += 1
            pGam = listIndiv[pID].getGamete()
            lGamAlpha.append(pGam[0])
            lGamA.append(pGam[1])
            lGamBeta.append(pGam[2])
    #####-----return list of gametes
    return ([lGamAlpha, lGamA, lGamBeta])

##########----------get highest fitness
def getMaxAbsFit(listIndiv):
    listAbsFit = [x.getAbsFitness() for x in listIndiv]
    maxAbsFit = np.max(listAbsFit)
    return(maxAbsFit)

##########----------get average fitness
def getAvgAbsFit(listIndiv):
    listStat = [x.getAbsFitness() for x in listIndiv]
    avgStat = np.mean(listStat)
    return (avgStat)

##########----------get average fitness by genotype of protein A
def getAvgAbsFit_GtA(listIndiv):
    listFit = [x.getAbsFitness() for x in listIndiv]
    listGt = [x.getGtA() for x in listIndiv]
    n0, n1, n2, fitA0, fitA1, fitA2 = 0, 0, 0, [], [], []
    for i in range(len(listIndiv)):
        gt, fit = listGt[i], listFit[i]
        if gt == 0:
            n0 += 1
            fitA0 = fitA0 + [fit]
        elif gt == 1:
            n1 += 1
            fitA1 = fitA1 + [fit]
        elif gt == 2:
            n2 += 1
            fitA2 = fitA2 + [fit]
    avgFitA0 = np.mean(fitA0) if n0 > 0 else 0
    avgFitA1 = np.mean(fitA1) if n1 > 0 else 0
    avgFitA2 = np.mean(fitA2) if n2 > 0 else 0
    return [avgFitA0, avgFitA1, avgFitA2]

##########----------get mean for total expression of A
def getAvgExpA(listIndiv):
    listStat = [sum(x.getExpA()) for x in listIndiv]
    avgStat = np.mean(listStat)
    return (avgStat)

##########----------get variance for total expression of A
def getVarExpA(listIndiv):
    listStat = [sum(x.getExpA()) for x in listIndiv]
    varStat = np.var(listStat)
    return (varStat)

##########----------get expression of allele A0 and A1 from heterozygotes averaged across individuals
def getAvgAlleleExpA_Het(listIndiv):
    hetConcA0, hetConcA1 = [], []
    for INDIV in listIndiv:
        gt = INDIV.getGtA()
        a1, a2 = INDIV.getA()
        concA1, concA2 = INDIV.getExpA()
        #####-----get conc of allele A0 and A1 in heterozygote
        if gt == 1:
            if a1 == 0:
                hetConcA0 += [concA1]
                hetConcA1 += [concA2]
            else:
                hetConcA0 += [concA2]
                hetConcA1 += [concA1]
    if len(hetConcA0) > 0:
        avgConcC0, avgConcA1 = np.mean(hetConcA0), np.mean(hetConcA1)
        return ([avgConcC0, avgConcA1])
    else:
        return ([-1, -1])

##########----------get average standardized experssion of B
def getAvgStExpB(listIndiv):
    listStat = [x.getStExpB() for x in listIndiv]
    avgStat = np.mean(listStat)
    return (avgStat)

##########----------get variance for standardized experssion of B
def getVarStExpB(listIndiv):
    listStat = [x.getStExpB() for x in listIndiv]
    varStat = np.var(listStat)
    return (varStat)


##########----------get frequency of A=1 allele
def getFreqA(listIndiv):
    listStat = [x.getGtA() for x in listIndiv]
    freqA = np.sum(listStat)/(2*len(listIndiv))
    return (freqA)

##########----------get genotype frequency at coding region of gene A
def getGtFreqA(listIndiv):
    listStat = [x.getGtA() for x in listIndiv]
    gtFreqA = [listStat.count(x)/len(listIndiv) for x in [0,1,2]]
    return (gtFreqA)

##########----------get average value of ALPHA
def getAvgAlpha(listIndiv):
    listStat = []
    [listStat.extend(x.getAlpha()) for x in listIndiv]
    avgStat = np.mean(listStat)
    return (avgStat)

##########----------get average value of BETA
def getAvgBeta(listIndiv):
    listStat = []
    [listStat.extend(x.getBeta()) for x in listIndiv]
    avgStat = np.mean(listStat)
    return (avgStat)

#########----------count occurence of each binding value
def countBindValue(listStat):
    lValues = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    lCounts = [listStat.count(x) for x in lValues]
    return (lCounts)

#########----------count occurence of values for ALPHA
def getDistAlpha(listIndiv):
    listStat = []
    [listStat.extend(x.getAlpha()) for x in listIndiv]
    listCounts = countBindValue(listStat)
    return (listCounts)

#########----------count occurence of values for BETA
def getDistBeta(listIndiv):
    listStat = []
    [listStat.extend(x.getBeta()) for x in listIndiv]
    listCounts = countBindValue(listStat)
    return (listCounts)


##########----------count num individuals with exp B in discrete bins
##########  bins are MIN <= x < MAX
##########      last bin uses 0.9 <= x < 1.1 because maximum expression of B is 1
def countBinnedExpB(listStat):
    lMin = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    lMax = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1]
    lCounts = []
    for MIN, MAX in zip(lMin, lMax):
        lCounts += [sum([1 for x in listStat if MIN <= x < MAX])]
    return (lCounts)

def getDistStExpB(listIndiv):
    listStat = [x.getStExpB() for x in listIndiv]
    listCounts = countBinnedExpB(listStat)
    return (listCounts)

##########----------calculate correlation between loci
def getCorrel(listIndiv):
    #####get lists of value from each indiv
    N = len(listIndiv)
    alpha1, alpha2, a1, a2, beta1, beta2 = [-1]*N, [-1]*N, [-1]*N, [-1]*N, [-1]*N, [-1]*N
    for i, x in enumerate(listIndiv):
        alpha1[i], alpha2[i] = x.getAlpha()
        a1[i], a2[i] = x.getA()
        beta1[i], beta2[i] = x.getBeta()
    #####combine list between homologous Cm for haplotype correlations
    alphaAll, aAll, betaAll = alpha1 + alpha2, a1 + a2, beta1 + beta2
    #####haplotype correlations (across loci)
    cAlphaA, cAlphaBeta, cABeta = float("nan"), float("nan"), float("nan")
    if np.var(alphaAll) > 0 and np.var(aAll) > 0:
        cAlphaA = np.corrcoef(alphaAll, aAll)[0, 1]
    if np.var(alphaAll) > 0 and np.var(betaAll) > 0:
        cAlphaBeta = np.corrcoef(alphaAll, betaAll)[0, 1]
    if np.var(aAll) > 0 and np.var(betaAll) > 0:
        cABeta = np.corrcoef(aAll, betaAll)[0, 1]
    #####return
    return [cAlphaA, cAlphaBeta, cABeta]



####################--------------------Object to represent individual
class indiv():
    #####-----Initizlize individual
    def __init__(self, SEX, CONC_D, GA, ALPHA1, ALPHA2, A1, A2, BETA1, BETA2, KA, KB, ASAT, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB):
        
        #####-----assign sex and allele
        #####   'alpha' is binding site allele for gene A (real number between 0 and 1)
        #####   'a' is protein coding allele for gene A (binary 0 or 1)
        #####   'beta' is binding site allele for gene B (real number between 0 and 1)
        self.sex = SEX
        self.uAlpha, self.uBeta = UA, UB
        self.alpha1, self.alpha2 = ALPHA1, ALPHA2
        self.a1, self.a2 = A1, A2
        self.beta1, self.beta2 = BETA1, BETA2
        self.gtA = calcGt(self.a1, self.a2)
        self.absFitness = 1

        #####-----assign binding affinity between sex TF and gene A binding site
        self.mSexAlpha1 = calcMAlpha(self.sex, self.alpha1)
        self.mSexAlpha2 = calcMAlpha(self.sex, self.alpha2)

        #####-----assign binding affinity between gene A product and gene B binding site
        #####   depends on value of 'a' (0 or 1) and 'beta'
        self.mA1Beta1 = calcMBeta(self.a1, self.beta1)
        self.mA1Beta2 = calcMBeta(self.a1, self.beta2)
        self.mA2Beta1 = calcMBeta(self.a2, self.beta1)
        self.mA2Beta2 = calcMBeta(self.a2, self.beta2)

        #####-----calc concentration of product from gene A (TF for gene B)
        self.concA1 = calcFO_Haploid(CONC_D, self.mSexAlpha1, KA)*GA
        self.concA2 = calcFO_Haploid(CONC_D, self.mSexAlpha2, KA)*GA

        #####-----calc expression of gene B
        #####   depends on concentration of gene A products (concA1, concA2) ,their identity (a1, a2) and binding sites affinity (mABeta)
        self.concB1 = 0.5 * (calcExpB(self.a1, self.concA1, self.concA2, self.mA1Beta1, self.mA2Beta1, KB) + calcExpB(self.a2, self.concA2, self.concA1, self.mA2Beta1, self.mA1Beta1, KB))
        self.concB2 = 0.5 * (calcExpB(self.a1, self.concA1, self.concA2, self.mA1Beta2, self.mA2Beta2, KB) + calcExpB(self.a2, self.concA2, self.concA1, self.mA2Beta2, self.mA1Beta2, KB))
        self.concB = self.concB1 + self.concB2

        #####-----calc max and min expression of B using saturating conc of A
        self.concMinB = -1 * calcFO_Haploid(ASAT, 0, KB)
        self.concMaxB = calcFO_Haploid(ASAT, 0, KB)

        #####-----standardize expression of B
        self.stConcB = (self.concB-self.concMinB)/( self.concMaxB-self.concMinB)

        #####-----calc absolute fitness using stConcB and selection parameters
        #self.absFitness = calcAbsFitness(SEX, self.stConcB, S_F, S_M, GAMMA_F, GAMMA_M)

    #####-----Propoerties of individual
    def getSex(self):
        return self.sex
    def getGtA(self):
        return self.gtA
    def getA(self):
        return [self.a1, self.a2]
    def getAlpha(self):
        return [self.alpha1, self.alpha2]
    def getBeta(self):
        return [self.beta1, self.beta2]
    def getMSAlpha(self):
        return [self.mSexAlpha1, self.mSexAlpha2]
    def getMABeta(self):
        return [self.mA1Beta1, self.mA1Beta2, self.mA2Beta1, self.mA2Beta2]
    def getExpA(self):
        return [self.concA1, self.concA2]
    def getExpB(self):
        return [self.concB1, self.concB2, self.concB]
    def getSatExpB(self):
        return [self.concMinB, self.concMaxB]
    def getStExpB(self):
        return self.stConcB
    def getAbsFitness(self):
        return self.absFitness
    def getSmy(self):
        SMY = [self.sex, self.gtA, self.a1, self.a2, self.alpha1, self.alpha2, self.beta1, self.beta2, self.concA1, self.concA2, self.stConcB, self.absFitness]
        return SMY

    #####-----randomly pass on one allele as gamete
    #####   assume complete linkage for gene A and free recomb between A and B
    def getGamete(self):
        #####-----get random 0 or 1
        rBinom1 = np.random.binomial(1, 0.5)
        rBinom2 = np.random.binomial(1, 0.5)
        #####-----rBinom1, rBinom2 controls inheritance of gene A and B, respectively
        if rBinom1 == 0:
            if rBinom2 == 0:
                return [mutateBindValue(self.uAlpha, self.alpha1), self.a1, mutateBindValue(self.uBeta, self.beta1)]
            else:
                return [mutateBindValue(self.uAlpha, self.alpha1), self.a1, mutateBindValue(self.uBeta, self.beta2)]
        else:
            if rBinom2 == 0:
                return [mutateBindValue(self.uAlpha, self.alpha2), self.a2, mutateBindValue(self.uBeta, self.beta1)]
            else:
                return [mutateBindValue(self.uAlpha, self.alpha2), self.a2, mutateBindValue(self.uBeta, self.beta2)]


####################--------------------Simulation type 1
####################  random assigment of ALPHA, A and BETA values
def runSim01(CONC_D, GA, KA, KB, ASAT, SR, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB, N, PROPF, NUMGEN):
    ##########----------initialize dfs to collect data
    SIMTYPE = 1
    COLNAME = [
        'SIMTYPE','CONC_D','GA','KA','KB','ASAT','SR','S_F','S_M','GAMMA_F','GAMMA_M','UA','UB','N','PROPF','NUMGEN','GEN',
        'maxFit_F','maxFit_M','maxFit','avgFit_F','avgFit_M','avgFit',
        'avgFit_F00','avgFit_F01','avgFit_F11','avgFit_M00','avgFit_M01','avgFit_M11','avgFit_00','avgFit_01','avgFit_11',
        'avgExpA_F', 'avgExpA_M', 'avgExpA', 'varExpA_F', 'varExpA_M', 'varExpA',
        'avgStExpB_F', 'avgStExpB_M', 'avgStExpB', 'varStExpB_F', 'varStExpB_M', 'varStExpB',
        'avgHetExpA0_F', 'avgHetExpA1_F', 'avgHetExpA0_M', 'avgHetExpA1_M',
        'freqA_F','freqA_M','freqA','gtA00_F','gtA01_F','gtA11_F','gtA00_M','gtA01_M','gtA11_M','gtA00','gtA01','gtA11',
        'avgAlpha_F','avgAlpha_M','avgAlpha','avgBeta_F','avgBeta_M','avgBeta',
        'cAlphaA_F', 'cAlphaBeta_F', 'cABeta_F', 'cAlphaA_M', 'cAlphaBeta_M', 'cABeta_M'
    ]
    dfGeneral = pd.DataFrame(columns=COLNAME)

    COLNAME = [
        'SIMTYPE', 'CONC_D', 'GA', 'KA', 'KB', 'ASAT', 'SR', 'S_F', 'S_M', 'GAMMA_F', 'GAMMA_M', 'UA', 'UB', 'N', 'PROPF', 'NUMGEN', 'GEN',
        'alF0', 'alF1', 'alF2', 'alF3', 'alF4', 'alF5', 'alF6', 'alF7', 'alF8', 'alF9', 'alF10',
        'alM0', 'alM1', 'alM2', 'alM3', 'alM4', 'alM5', 'alM6', 'alM7', 'alM8', 'alM9', 'alM10',
        'beF0', 'beF1', 'beF2', 'beF3', 'beF4', 'beF5', 'beF6', 'beF7', 'beF8', 'beF9', 'beF10',
        'beM0', 'beM1', 'beM2', 'beM3', 'beM4', 'beM5', 'beM6', 'beM7', 'beM8', 'beM9', 'beM10',
        'al0', 'al1', 'al2', 'al3', 'al4', 'al5', 'al6', 'al7', 'al8', 'al9', 'al10',
        'be0', 'be1', 'be2', 'be3', 'be4', 'be5', 'be6', 'be7', 'be8', 'be9', 'be10',
        'BF0', 'BF1', 'BF2', 'BF3', 'BF4', 'BF5', 'BF6', 'BF7', 'BF8', 'BF9',
        'BM0', 'BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7', 'BM8', 'BM9',
    ]
    dfDist = pd.DataFrame(columns=COLNAME)

    ##########----------calculate num of male and females
    N_F= int(N*(1-PROPF))
    N_M = N - N_F

    ##########----------Initialize alleles
    ##########  equal prob for each allele of each locus
    listBind = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #listBind = [1.0]
    lAlpha1 = random.choices(listBind, k=N)
    lAlpha2 = random.choices(listBind, k=N)
    lA1 = np.random.binomial(1, 0.5, N)
    lA2 = np.random.binomial(1, 0.5, N)
    lBeta1 = random.choices(listBind, k=N)
    lBeta2 = random.choices(listBind, k=N)
    random.shuffle(lAlpha1)
    random.shuffle(lAlpha2)
    random.shuffle(lA1)
    random.shuffle(lA2)
    random.shuffle(lBeta1)
    random.shuffle(lBeta2)

    ##########----------Initialize sexes based on N_F, N_M
    popSex = ['f'] * int(N_F) + ['m'] * int(N_M)
    random.shuffle(popSex)

    ##########----------Initilize N individuals in population
    POPUL_F, POPUL_M = [], []
    for i in range(N):
        SEX = popSex[i]
        if SEX == 'f':
            POPUL_F.append(
                indiv('f', CONC_D, GA, lAlpha1[i], lAlpha2[i], lA1[i], lA2[i], lBeta1[i], lBeta2[i], KA, KB, ASAT, S_F,
                      S_M, GAMMA_F, GAMMA_M, UA, UB))
        else:
            POPUL_M.append(
                indiv('m', CONC_D, GA, lAlpha1[i], lAlpha2[i], lA1[i], lA2[i], lBeta1[i], lBeta2[i], KA, KB, ASAT, S_F,
                      S_M, GAMMA_F, GAMMA_M, UA, UB))
    ##########----------Initial calculation of absolute fitness
    for INDIV in POPUL_F + POPUL_M:
        INDIV.absFitness = calcAbsFitness(INDIV.getSex(), INDIV.getStExpB(), SR, S_F, S_M, GAMMA_F, GAMMA_M)

    ##########----------Run for NUMGEN generations
    ##########  if 'A' gets fixed (freq = 0 or 1), then sim stops after MaxExtraGEN generations
    extraGen, MaxExtraGen = 0, 500
    for gen in range(NUMGEN):
        #####-----determine highest fitness
        maxAbsFit_F = getMaxAbsFit(POPUL_F)
        maxAbsFit_M = getMaxAbsFit(POPUL_M)
        #####-----get gametes using relative fitness
        listFemGam = selectionGamete(POPUL_F, maxAbsFit_F, N)
        listMalGam = selectionGamete(POPUL_M, maxAbsFit_M, N)
        lAlpha1, lA1, lBeta1 = listFemGam[0], listFemGam[1], listFemGam[2]
        lAlpha2, lA2, lBeta2 = listMalGam[0], listMalGam[1], listMalGam[2]
        #####-----create offspring from gametes
        offspringSex = ['f'] * int(N_F) + ['m'] * int(N_M)
        random.shuffle(offspringSex)
        OFFSPRING_F, OFFSPRING_M = [], []
        for i in range(N):
            SEX = offspringSex[i]
            if SEX == 'f':
                tempInd = indiv('f', CONC_D, GA, lAlpha1[i], lAlpha2[i], lA1[i], lA2[i], lBeta1[i], lBeta2[i], KA, KB,
                                ASAT, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB)
                OFFSPRING_F.append(tempInd)
            else:
                tempInd = indiv('m', CONC_D, GA, lAlpha1[i], lAlpha2[i], lA1[i], lA2[i], lBeta1[i], lBeta2[i], KA, KB,
                                ASAT, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB)
                OFFSPRING_M.append(tempInd)
        #####-----new population created from offspring
        POPUL_F = OFFSPRING_F.copy()
        POPUL_M = OFFSPRING_M.copy()

        #####-----Calculate fitness
        for INDIV in POPUL_F + POPUL_M:
            INDIV.absFitness = calcAbsFitness(INDIV.getSex(), INDIV.getStExpB(), SR, S_F, S_M, GAMMA_F, GAMMA_M)

        #####-----Collect data every 10 generations
        if (gen + 1) % 10 == 0:
            vParam = [SIMTYPE, CONC_D, GA, KA, KB, ASAT, SR, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB, N, PROPF, NUMGEN, gen+1]
            vMaxFit = [getMaxAbsFit(POPUL_F), getMaxAbsFit(POPUL_M), getMaxAbsFit(POPUL_F+POPUL_M)]
            vAvgFit = [getAvgAbsFit(POPUL_F), getAvgAbsFit(POPUL_M), getAvgAbsFit(POPUL_F+POPUL_M)]
            vAvgFit_GtA = getAvgAbsFit_GtA(POPUL_F) + getAvgAbsFit_GtA(POPUL_M) + getAvgAbsFit_GtA(POPUL_F+POPUL_M)
            vAvgStExpA = [getAvgExpA(POPUL_F), getAvgExpA(POPUL_M), getAvgExpA(POPUL_F + POPUL_M)]
            vVarStExpA = [getVarExpA(POPUL_F), getVarExpA(POPUL_M), getVarExpA(POPUL_F + POPUL_M)]
            vAvgStExpB = [getAvgStExpB(POPUL_F), getAvgStExpB(POPUL_M), getAvgStExpB(POPUL_F+POPUL_M)]
            vVarStExpB = [getVarStExpB(POPUL_F), getVarStExpB(POPUL_M), getVarStExpB(POPUL_F+POPUL_M)]
            vAvgAlleleExpA_het = getAvgAlleleExpA_Het(POPUL_F) + getAvgAlleleExpA_Het(POPUL_M)
            vFreqA = [getFreqA(POPUL_F), getFreqA(POPUL_M), getFreqA(POPUL_F+POPUL_M)]
            vGtA = getGtFreqA(POPUL_F) + getGtFreqA(POPUL_M) + getGtFreqA(POPUL_F+POPUL_M)
            vAvgAlpha = [getAvgAlpha(POPUL_F), getAvgAlpha(POPUL_M), getAvgAlpha(POPUL_F+POPUL_M)]
            vAvgBeta = [getAvgBeta(POPUL_F), getAvgBeta(POPUL_M), getAvgBeta(POPUL_F+POPUL_M)]
            vCorrel_F, vCorrel_M = getCorrel(POPUL_F), getCorrel(POPUL_M)

            vDistAlpha_F, vDistAlpha_M = getDistAlpha(POPUL_F), getDistAlpha(POPUL_M)
            vDistBeta_F, vDistBeta_M = getDistBeta(POPUL_F), getDistBeta(POPUL_M)
            vDistAlpha, vDistBeta = getDistAlpha(POPUL_F+POPUL_M), getDistBeta(POPUL_M+POPUL_M)
            vDistStExpB_F, vDistStExpB_M = getDistStExpB(POPUL_F), getDistStExpB(POPUL_M)

            vAllData = vParam + vMaxFit + vAvgFit + vAvgFit_GtA + \
                       vAvgStExpA + vVarStExpA + vAvgStExpB + vVarStExpB + vAvgAlleleExpA_het + \
                       vFreqA + vGtA + vAvgAlpha + vAvgBeta + vCorrel_F + vCorrel_M
            dfGeneral.loc[len(dfGeneral)] = vAllData

            vDistData = vParam + vDistAlpha_F + vDistAlpha_M + vDistBeta_F + vDistBeta_M + vDistAlpha + vDistBeta + vDistStExpB_F + vDistStExpB_M
            dfDist.loc[len(dfDist)] = vDistData

        '''
        #####-----Check freq of A, if fixed then stops after MaxExtraGen generations
        freqA = getFreqA(POPUL_F+POPUL_M)
        if freqA == 0 or freqA == 1:
            extraGen += 1
        if extraGen == MaxExtraGen:
            break
        '''
        if gen == NUMGEN-1:
            COLNAME = ['ID', 'SEX', 'GTA', 'A1', 'A2', 'ALPHA1', 'ALPHA2', 'BETA1', 'BETA2', 'concA1', 'concA2', 'stconB', 'absFitness']
            dfIndiv = pd.DataFrame(columns=COLNAME)
            #[self.sex, self.gtA, self.a1, self.a2, self.alpha1, self.alpha2, self.beta1, self.beta2, self.concA1, self.concA2, self.stConB, self.absFitness]
            for i, INDIV in enumerate(POPUL_F+POPUL_M):
                SMY = INDIV.getSmy()
                dfIndiv.loc[len(dfIndiv)] = [i]+SMY

    #####-----return data
    return dfGeneral, dfDist, dfIndiv



####################--------------------Function that runs multiple replicates of simulation
####################    each replicate output file has a different name (increment by an integer)
def runReps(DIR, NAME, REPS, CONC_D, GA, KA, KB, ASAT, SR, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB, N, PROPF, NUMGEN):
    ##########----------Run simulations for REPS number of replicates
    for r in range(REPS):
        print ('Running [{}] rep [{}]'.format(NAME, r + 1))
        ##########----------run simulation
        dfGeneral, dfDist, dfIndiv = runSim01(CONC_D, GA, KA, KB, ASAT, SR, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB, N, PROPF, NUMGEN)
        ##########----------output data
        outGeneral = DIR + NAME + '.' + str(r + 1) + '.general.csv'
        outDist = DIR + NAME + '.' + str(r + 1) + '.dist.csv'
        outIndiv = DIR + NAME + '.' + str(r + 1) + '.indiv.csv'
        dfGeneral.to_csv(path_or_buf=outGeneral, index=False)
        dfDist.to_csv(path_or_buf=outDist, index=False)
        dfIndiv.to_csv(path_or_buf=outIndiv, index=False)


####################--------------------Function that reads parameters from csv file and runs simulations
def runParams(NPROC, PARFILE):
    ##########----------read parameter file
    try:
        PARAM = pd.read_csv(PARFILE)
    except OSError:
        print('ERROR: could not open [{}]'.format(PARFILE))
        sys.exit()
    ##########----------Check that all NAMES for outputing files are unique
    nRow = len(PARAM)
    nUniqName = len(pd.unique(PARAM['NAME']))
    if nUniqName != nRow:
        print('ERROR: one or more NAME are repeated: {}'.format(PARAM['NAME'].tolist()))
        sys.exit()
    ##########----------Run simulation for each set of parameters
    pool = Pool(processes=NPROC)
    for i in range(len(PARAM)):
        #####-----get params for row
        parRow = PARAM.loc[i]
        DIR, NAME, REPS, N, PROPF, NUMGEN = str(parRow['DIR']), str(parRow['NAME']), int(parRow['REPS']), int(
            parRow['N']), float(parRow['PROPF']), int(parRow['NUMGEN'])
        CONC_D, GA, KA, KB, ASAT = float(parRow['CONC_D']), float(parRow['GA']), float(parRow['KA']), float(
            parRow['KB']), float(parRow['ASAT'])
        SR, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB = int(parRow['SR']), float(parRow['S_F']), float(parRow['S_M']), float(
            parRow['GAMMA_F']), float(parRow['GAMMA_M']), float(parRow['UA']), float(parRow['UB'])
        #print(DIR, NAME, REPS, CONC_D, GA, KA, KB, ASAT, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB, N, PROPF, NUMGEN)
        #####-----run simlation if directory exists
        if os.path.isdir(DIR):
            print('Running simulations for [{}]'.format(NAME))
            pool.apply_async(runReps, args=(DIR, NAME, REPS, CONC_D, GA, KA, KB, ASAT, SR, S_F, S_M, GAMMA_F, GAMMA_M, UA, UB, N, PROPF, NUMGEN))
        else:
            print('DIR [{}] does not exist. Simulations for [{}] were not performed.'.format(DIR, NAME))
    #####-----makes sure all processes finish before continung script
    pool.close()
    pool.join()

####################--------------------Read arguments for MAIN
def parseArg():
    ##########----------Initialize parser
    parser = argparse.ArgumentParser()
    ##########----------Read arguments
    parser.add_argument("-n", "--NPROC", help="Number of processes to run simultaneously", type=int, required=True)
    parser.add_argument("-p", "--PARFILE", help="CSV file containing parameters for simulations", type=str, required=True)
    args = parser.parse_args()
    ##########----------Return arguments
    return args.NPROC, args.PARFILE

####################--------------------MAIN
def main():
    ##########----------get arguments from command line
    NPROC, PARFILE = parseArg()
    print("Num of processes to run: [{}]".format(NPROC))
    print("Reading parameters from: [{}]\n".format(PARFILE))
    ##########----------start program
    starttime = time.time()
    #PARFILE = 'D:/Users/eddie/DomSimData/TestMulti/PARAM.csv'
    runParams(NPROC, PARFILE)
    print('That took {} seconds'.format(time.time() - starttime))

if __name__ == "__main__":
    main()


