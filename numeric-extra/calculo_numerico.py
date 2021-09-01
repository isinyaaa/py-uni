from math import * # usaremos a biblioteca padrão math

norm1 = lambda v : sum(map(fabs, v)) # norma_1

# multiplicação de matriz por vetor
matVecMul = lambda v, w : [sum([v[i][k] * w[k] for k in range(0, len(w))]) for i in range(0,len(v))]

# multiplicação de matriz por escalar
scalarMult = lambda v, c : [scalarMult(i, c) if type(i) is list else c*i for i in v]

# soma de matrizes
matSum = lambda v, w : [matSum(v[i], w[i]) if type(v[i]) is list else v[i] + w[i] for i in range(0, len(v))]

# calcula a matriz A para a primeira tarefa
def AMat(v):
  result = []
  for i in range(0,len(v)):
    element = []
    for j in range(0,len(v)):
      if j + 1 in v[i]:
        element.append(1/len(v[j]))
      else:
        element.append(0);
    result.append(element)
  return result

# calcula a matriz S para a primeira tarefa
SMat = lambda v : [[1/len(v) for j in range(0,len(v))] for i in range(0,len(v))]

# função a ser iterada para a primeira tarefa
# x_i = (1 - m)A x_{i-1} + m x_0
iterOp = lambda v, w, z, c : matSum(matVecMul(scalarMult(v, (1 - c)),w), scalarMult(z, c))

def main_tarefa1():
  x0 = [1/8 for i in range(0,8)]

  m = 0.15

  matPages = [{2,3},  # pagina 1
              {3,4},  # pagina 2
              {6,7},  # pagina 3
              {5,6},  # pagina 4
              {6},    # pagina 5
              {7},    # pagina 6
              {8},    # pagina 7
              {1,3}   # pagina 8
              ]

  A = AMat(matPages)
  S = SMat(matPages)

  M = matSum(scalarMult(A, (1 - m)), scalarMult(S, m))

  # x_{l - 1}
  previous = iterOp(M,x0,x0,m)
  # x_l
  current = iterOp(M,previous,x0,m)

  # norma da comparação entre os valores de x's sucessivos
  nval = norm1(matSum(current,scalarMult(previous, -1)))

  lowerBound = 10**(-5)
  while nval >= lowerBound:
    previous = current.copy()
    current = iterOp(M,previous,x0,m)

    nval = norm1(matSum(previous,scalarMult(current, -1)))

  # vamos parear cada peso com a página correspondente
  zippedList = list(zip(range(1, len(current) + 1), current))
  # depois ordenamos e invertemos
  sortedList = sorted(zippedList, key=lambda weight : weight[1])[::-1]
  print("Tarefa 1")
  # agora vamos printar as páginas em ordem
  count = 1
  for i in sortedList:
    print(f"#{count} página {i[0]} com peso %.3f" %i[1])
    count += 1

main_tarefa1()

# função a ser iterada para a segunda tarefa
# x_i = (1 - m)y_{i - 1} + m x_0
iterOp2 = lambda v, w, c : matSum(scalarMult(v, (1 - c)), scalarMult(w, c))

def main_tarefa2():
  V = [] # matriz de pesos A alterada
  L = [] # linha onde o vertice i se encontra
  C = [] # coluna onde o vértice i se encontra
  total = 230
  m = 0.15
  tribes = 20

  for i in range(1, tribes + 1): # vamos ter 20 tribos
    for j in range(0, i + 1): # cada tribo tem i + 1 participantes
      cacique = i*(i + 1)//2
      if j == 0: # o primeiro será o cacique
        for k in range(1, tribes + 1): # que se liga a todos os outros caciques
          # nenhum índio/cacique se liga a si mesmo, portanto
          if i == k:
            continue
          # todo cacique terá número de índios na tribo + número de caciques - 2 em seus pesos, portanto
          V.append(1/(i + 1 + tribes - 2))
          # a coluna corresponde àquela do cacique atual, portanto
          C.append(cacique)
          # a linha corresponde ao cacique do loop em k
          cacique2 = k*(k + 1)//2
          L.append(cacique2)

        # mas cada cacique também se liga a todos os outros índios da tribo, como um índio qualquer

      for k in range(0, i + 1):
        # nenhum índio se liga a si mesmo, portanto
        if j == k:
          continue
        # todo índio terá número de índios na tribo - 1 em seus pesos, portanto
        if j == 0:
          V.append(1/(i + 1 + tribes - 2))
        else:
          V.append(1/(i + 1 - 1))
        # a coluna corresponde àquela do cacique + índio, portanto
        C.append(cacique + j)
        # a linha corresponde àquela do cacique + índio2, portanto
        L.append(cacique + k)

  x0 = [1/total for i in range(0,total)] # x0 vai ser um vetor com elementos 1/230
  # x_{l - 1}
  previous = x0.copy()

  y = [0 for i in range(0,total)] # y0 será um vetor nulo
  # x_l
  current = iterOp2(y, x0, m) # e o elemento x1 será (1 - m) y0 + m x0

  # norma da comparação entre os valores de x's sucessivos
  nval = norm1(matSum(previous, scalarMult(current, -1)))

  lowerBound = 10**(-5)
  while nval >= lowerBound:
    for i in range(0, len(V)):
      y[L[i] - 1] += V[i] * previous[C[i] - 1] # y_i = A x_i = y_{i-1} + V x_{i - 1}

    current = iterOp2(y, x0, m)
    # previous será copiado após o cálculo para que, na primeira iteração, utilizemos x0
    previous = current.copy()

    nval = norm1(matSum(previous, scalarMult(current, -1)))

  # lista sem duplicados
  dedList = []
  for i in range(1, tribes + 1):
    # cada cacique será guardado pelo seu número de página normalmente
    dedList.append([i*(i + 1)//2, current[i*(i + 1)//2 - 1]])
    # e em seguida um índio representativo (o primeiro da tribo), com [seu número de página, total de índios na tribo]
    dedList.append([[i*(i + 1)//2 + 1, i], current[i*(i + 1)//2]])

  # agora ordenamos a lista como na tarefa 1
  sortedList = sorted(dedList, key=lambda weight : weight[1])[::-1]

  print("Tarefa 2")
  count = 0
  prevcount = count
  for i in sortedList:
    if type(i[0]) is not list:
      count += 1
      print("#%3d" %count, 7*" ", "página  %3d" %i[0], 6*" ", " com peso %.4f" %i[1],sep="")
    elif i[0][1] == 1:
      count += 1
      print("#%3d" %count, 7*" ", "página  %3d" %i[0][0], 6*" "," com peso %.4f" %i[1],sep="")
    else:
      count += 1
      prevcount = count
      count += i[0][1] - 1
      print("#%3d à %3d" %(prevcount,count), " páginas %3d à %3d" %(i[0][0], i[0][0]+i[0][1]-1)," com peso %.4f" %i[1],sep="")

print()
main_tarefa2()
