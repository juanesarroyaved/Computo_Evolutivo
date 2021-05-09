import matplotlib.pyplot as plt
import random
import numpy as np
from numpy.core.fromnumeric import size, cumsum
from numpy.core.numeric import ones, zeros_like

def AG_simple(Nind, Lind, Pc, Pm, Maxgen, Nvar, rango, Mejor, Mejor_cromosoma):
    genotipo = creaprob(Nind, Lind)
    fenotipo = decodifica(genotipo, rango)
    objv = objfun(fenotipo)
    generaciones = 1

    while generaciones < Maxgen:
        aptitud = rankeo(objv, 1)
        nuevo_gen = ruleta(genotipo, fenotipo, aptitud)
        nuevo_gen = xpunto(nuevo_gen, Pc)
        nuevo_gen = muta(nuevo_gen, Pm)
        nuevo_feno = decodifica(nuevo_gen, rango)
        nuevo_objv = objfun(nuevo_feno)
        genotipo = nuevo_gen
        objv = nuevo_objv
        [valor_idx, idx] = max(objv)
        Mejor[generaciones] = valor
        Mejor_cromosoma[generaciones, :] = genotipo[idx, :]
        plt.plot(Mejor, 'ro')
        generaciones += 1

def creaprob(Nind, Lind):
    genotipo = 0.5 > random.random(Nind, Lind)
    
    return genotipo

def decodifica(genotipo, rango):
    Nvar = size(rango, 2)
    [Nind, Lind] = size(genotipo)
    Lvar = Lind / Nvar
    potencias = 2.^(0:Lvar-1)
    
    for i in len(Nind):
        for j in len(Nvar):
            fenotipo[i, j] = sum(potencias*genotipo[i,(j-1)*Lvar+1:j*Lvar])
    
    for i in len(Nvar):
        fenotipo[:,1] = rango[1,i] + ((rango[2,i]-rango[1,i])/(2^Lvar-1))*fenotipo[:,1]

    return fenotipo

def objfun(fenotipo):
    [Nind, Lind] = size(fenotipo)
    for i in len(Nind):
        objv[i, 1] = 21.5 + fenotipo[i,1]*sin(4*pi*fenotipo[i,1]) + fenotipo[i,2]*sin(20*pi*fenotipo[i,2])

    return objv

def rankeo(objv, direccion):
    SP = 2
    [Nind, Nobj] = size(objv)
    
    if direccion == 1:
        [nuevo_obj, posori] = sort(objv)
    else:
        [nuevo_obj, posori] = sort(-1*objv)
    
    apt = 2 - SP + 2*(SP-1)*((1:Nind)-1)/(Nind-1)
    aptitud[posori,1] = apt

    return aptitud

def ruleta(genotipo, fenotipo, aptitud):
    [Nind, aux] = size(aptitud)
    total = sum(aptitud)
    probabilidad = aptitud/total
    acumulada = cumsum(probabilidad)
    for i in len(Nind):
        selecciona = random.random()
        aux = find(acumulada >= selecciona)
        idx[i,1] = aux[1]
    
    nuevo_gen = genotipo[idx,:]
    
    return nuevo_gen

def xpunto(nuevo_gen, Pc):
    [Nind, Lind] = size(nuevo_gen)
    aux_gen = []
    par = Nind % 2
    for i in range(0,len(Nind)-1,2):
        cruza = random.random()
        if cruza <= Pc:
            corte = ceil((Lind-1)*random.random()) #entero mayor
            aux_gen[i,:] = [nuevo_gen[i,1:corte], nuevo_gen[i+1, corte+1:Lind]]
            aux_gen[i+1,:] = [nuevo_gen[i+1,1:corte], nuevo_gen[i, corte+1:Lind]]
        else:
            aux_gen[i,:] = nuevo_gen[i,:]
            aux_gen[i+1,:] = nuevo_gen[i+1,:]
    
    if par == 1:
        aux_gen[Nind,:] = nuevo_gen[Nind,:]
    
    nuevo_gen = aux_gen
    
    return nuevo_gen

def muta(nuevo_gen, Pm):
    [Nind, Lind] = size(nuevo_gen)
    valores = random.randrange(Nind, Lind)
    muta = valores <= Pm
    nuevo_gen = xor(nuevo_gen, muta)

    return nuevo_gen

if __name__ == '__main__':
    Nind = 100
    Lind = 20
    Pc = 0.9
    Pm = 0.01
    Maxgen = 1000
    Nvar = 2
    rango = [-3, 4.1; 12.1,5.8]
    Mejor = np.nan * ones((Maxgen, 1))
    Mejor_cromosoma = zeros_like((Nind, Lind*Nvar))

    AG_simple(Nind, Lind, Pc, Pm, Maxgen, Nvar, rango, Mejor, Mejor_cromosoma)