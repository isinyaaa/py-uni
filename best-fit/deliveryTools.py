import pandas as pd

def importData(name):
    data = pd.read_csv(name)
    df = pd.DataFrame(data, columns=['length','height','width','weight'])

    return [item(i) for i in df.values.tolist()]

# organizamos os itens pensando em na diferença entre suas dimensões,
# já que itens mais "quadrados" provavelmente cabem melhor
# para medir quão quadrados eles são eu comparo a derivada da área como função do perímetro
# e quanto menor ela for, mais próximo o item é de ser quadrado
def sortItems(data):
    return sorted(data, key=lambda i: abs(i.Len+i.Wid-2*i.Wid))

from deliveryObjects import *
import itertools as it

# gerar possibilidades
def genPossibilities(back,veh,data,trips,lastitemset):
    
    # contador de itens da rodada atual
    run = len(data)
    # enquanto tivermos elementos ou rodadas
    while run > 0:

        # aqui guardamos as possibilidades da rodada atual
        possibilities = sorted(it.combinations(data, run), key=lambda c: sum([i.W for i in c]))

        # print para mostrar o que está acontecendo
        print(f"conjuntos de {run} itens")
        for j in possibilities:
            print([k.id for k in j])

        # loopamos as possibilidades
        for p in possibilities:
            # e então loopamos cada veículo, começando pelo menor
            i = 0
            while i < len(veh):
                # se a coleção de itens couber aqui
                if veh[i].fit(p):
                    # pegamos os itens possíveis
                    itemset = set([k.id for k in p])
                    # salvamos na lista de possíveis viagens
                    while i < len(veh):
                        trips[i].append(itemset)
                        i += 1

                    # removemos dos dados
                    for j in p:
                        data.remove(j)
                    
                    if len(data) > 0:
                        # vamos para a próxima possibilidade
                        genPossibilities(back,veh,data,trips,itemset)
                    else:
                        # e caso tenha zerado
                        run = 0 # para finalizar o último loop
                        break
                i += 1

        run -= 1

    # backtracking
    for b in back:
        if b.id in lastitemset:
            data.append(b)
    
# peso total dos itens de uma viagem
def totalWeight(itemlist,trip):
    return sum([j.W for j in filter(lambda i: i.id in itemlist, trip)])

# listagem das possibilidades itens
def listPossibilities(comp,dataset):
    # cópia dos veículos
    veh = [i.clone() for i in comp.vehicles]
    # possibilidades de fitting para cada veículos
    trips = [list() for i in range(len(veh))]
    # itens (para uso)
    data = [i.clone() for i in dataset]
    # backup para o backtracking
    back = [i.clone() for i in dataset]

    # primeira chamada da recursão
    genPossibilities(back,veh,data,trips,set())
    print("antes de arrumar")
    print(trips)

    # pegamos a "melhor" dentre todas as possíveis (aquela com o maior número de objetos sendo carregados)
    trips = [sorted(v, key = lambda k: len(k)) for v in trips]
    # vamos devolver a lista de viagens com [ # veículo -> [{itens}, peso total] ]
    for i in range(len(trips)):
        if len(trips[i]) != 0:
            best = trips[i][len(trips[i]) - 1]
            trips[i] = [list(best), totalWeight(best,back)]
        else:
            trips[i] = [None, 0]
    
    print("pronto para ser impresso")
    print(trips)

    return trips

