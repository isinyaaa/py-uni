import numpy as np
import math


def GivensRotation(alphas, betas, gammas, deltas, i):
    # aplica a rotação de givens com as dicas do enunciado
    j = i + 1
    (sk, ck) = SkCk(alphas[i], betas[i])

    # muda sempre no maximo 3 colunas
    for k in range(3):
        if k == 2 and i != len(betas) - 1:
            prev_delta = deltas[i]
            prev_gamma = gammas[j]
            gammas[j] = sk*prev_delta + ck*prev_gamma
            deltas[i] = ck*prev_delta - sk*prev_gamma

        elif k == 1:
            prev_gamma = gammas[i]
            prev_alpha = alphas[j]
            gammas[i] = ck*prev_gamma - sk*prev_alpha
            alphas[j] = sk*prev_gamma + ck*prev_alpha


        if k == 0:
            prev_alpha = alphas[i]
            prev_beta = betas[i]
            betas[i] = sk*prev_alpha + ck*prev_beta
            alphas[i] = ck*prev_alpha - sk*prev_beta

    return sk, ck

def SkCk(a, b):
    # calcula sk e ck da mesma forma que os exemplos no enunciado
    dem = math.sqrt(a**2 + b**2) # diagonal
    # projeções
    ck = a/dem
    sk = -b/dem
    return (sk, ck) # vetor

def sign(n):
    # função sinal
    return -1 if n < 0 else 1

def RQ(alphas, betas, gammas, vector_ck, vector_sk):
    # multiplica RQ

    for i in range(len(alphas) - 1):
        j = i + 1
        sk, ck = vector_sk[i], vector_ck[i]

        alpha_antigo = alphas[i]
        alphas[i] = - sk*gammas[i] + ck*alpha_antigo
        #

        alpha_antigo = alphas[j]
        beta_antigo = betas[i]
        betas[i] = ck*beta_antigo - sk*alpha_antigo
        alphas[j] = ck*alpha_antigo + sk*beta_antigo

    for i in range(len(gammas)):
        gammas[i] = betas[i]

def VQ(V, vector_ck, vector_sk):
    # multiplica VQ
    n = len(V)
    for i in range(n-1):
        j =i+1
        sk, ck = vector_sk[i], vector_ck[i]

        for k in range(n):
            Vki_antigo = V[k, i]
            Vkj_antigo = V[k, j]
            V[k, j] = ck*Vkj_antigo + sk*Vki_antigo
            V[k, i] = ck*Vki_antigo - sk*Vkj_antigo

def QR(eps, wilkinson, alphas, betas, gammas, deltas, V):
    # algoritmo iterativo qr para achar autovetores e autovalores

    n = len(alphas)
    #V = V0 #np.identity(n)
    k = 0

    for m in range(len(betas), 0, -1):
        while abs(betas[m-1]) > eps:
            # reseta as variáveis
            alphas_uk = np.zeros(n)
            cks = np.zeros(n-1)
            sks = np.zeros(n-1)

            if not wilkinson or k == 0: # não usa wilkinson
                uk = 0
            else:
                dk = (alphas[m-1] - alphas[m])/2
                uk = alphas[m] + dk - sign(dk)*(dk**2 + betas[m-1]**2)**(1/2)

            # deslocamento
            for i in range(n):
                alphas_uk[i] = alphas[i] - uk

            for i in range(n-1):
                (sks[i] , cks[i]) = GivensRotation(alphas_uk, betas, gammas, deltas, i)

            RQ(alphas_uk, betas, gammas, cks, sks)

            VQ(V, cks, sks)

            # anti-deslocamento
            for i in range(n):
                alphas[i] = alphas_uk[i] + uk

            deltas = np.zeros(n-2)

            if (betas[m-1]**2)**0.5 < eps:
                betas[m-1] = 0  # zera beta se for muito pequeno
            k += 1

    return k, V


def gen_a_i(A, i):
    a_i = np.array(A[i + 1:, i], copy=True)
    return a_i


def w_i(A, i):
    # calculamos a_i
    a_i = gen_a_i(A, i)
    # pegamos o sinal a_i
    delta = -1 if a_i[0] < 0 else 1
    # geramos o vetor unitário e_i
    e_i = np.zeros(np.size(a_i))
    e_i[0] = 1
    # pela fórmula dada no enunciado, temos
    return a_i + delta * np.linalg.norm(a_i) * e_i


def sandwich(A, i, w):
    # produto do vetor w com as colunas de A
    H_wbar_left(A, i, w)
    # produto do vetor w com as linhas de A
    H_wbar_right(A, i, w)


def H_wbar_right(A, i, w):
    # vamos da (i+1)-ésima linha até o final da matriz
    for j in range(i + 1, np.size(A[:, 0])):
        # x é um corte da linha atual partindo da posição i
        x = A[j, i + 1:]
        # H · x
        new_x = x - 2*(np.inner(w, x)/np.inner(w, w)) * w
        # e atualizamos os valores em A
        for k in range(i + 1, np.size(A[j, :])):
            if k - (i + 1) < 0:
                A[j, k] = A[k, j]
            else:
                n = new_x[k - (i + 1)]
                A[j, k] = n

def H_wbar_left(A, i, w):
    # idem à rotina acima, exceto que agora olhando as colunas {i, n}
    for j in range(i, np.size(A[0, :])):
        x = A[i + 1:, j]
        new_x = x - 2*(np.inner(w, x)/np.inner(w, w)) * w
        for k in range(np.size(A[:, j])):
            if k - (i + 1) < 0:
                A[k, j] = A[j, k]
            else:
                n = new_x[k - (i + 1)]
                A[k, j] = n

class Bar:
    def __init__(self, i, j, theta, length, fixed_edge = False):
        self.i = i
        self.j = j
        self.theta = theta
        self.length = length
        self.fixed_edge = fixed_edge


def get_bar_K_ij(bar, Area, E):
   C = math.cos(bar.theta)
   S = math.sin(bar.theta)
   stiffness_vec = np.array([-C, -S, C, S], dtype=np.float64)
   return Area * E/bar.length * np.outer(stiffness_vec, stiffness_vec)


def main():
    while True:
        question = input("Digite alguma tarefa que gostaria de fazer (a, b ou c) ou qualquer outra coisa para sair: ").lower()
        if question not in 'abc':
            break

        if question in 'ab':

            A = []
            if question == 'a':
                with open('./input-a') as f:
                    lines = f.readlines()
                    shape = int(lines[0])
                    for i in range(shape):
                        A.append(list(map(float, lines[i + 1].split())))

                preview_ev = [-1, -2, 2, 7]
            else:
                shape = 0
                with open('./input-b') as f:
                    lines = f.readlines()
                    shape = int(lines[0])
                    for i in range(shape):
                        A.append(list(map(float, lines[i + 1].split())))

                formula = lambda i, n : (1-math.cos((2*i-1)*math.pi/(2*n+1)))**(-1)/2
                preview_ev = np.array([formula(i, shape) for i in range(1, 21)])

            A = np.array(A, dtype=np.float64)

            original_mat = A.copy()

            # algoritmo básico de tridiagonalização +
            # processamento de H^T
            size = np.size(A[0, :])
            H = np.eye(size)
            for i in range(size-2):
                w = w_i(A, i)
                sandwich(A, i, w)
                # Multiplicação por H^T não é sujeita a redução dimensional
                Hw = np.zeros(size - 1)
                for j in range(len(w)):
                    Hw[size - 2 - j] = w[len(w) - 1 - j]
                H_wbar_right(H, 0, Hw)

            # aqui isolamos a diagonal superior e sub-diagonal de A tridiagonalizado
            alphas = np.zeros(size)
            betas = np.zeros(size - 1)

            for i in range(size):
                alphas[i] = A[i, i]
                if i < size - 1:
                    betas[i] = A[i + 1, i]

            # agora, basta replicarmos \beta em \gamma já que a matriz é simétrica
            gammas = betas.copy()
            deltas = np.zeros(size - 2)
            epsilon = 1e-6
            k, V = QR(epsilon, True, alphas, betas, gammas, deltas, H.copy())

            V = np.transpose(V)

            print(f"Matriz de escolha da tarefa {question.upper()}:")
            print(original_mat)
            print()
            print("Autovalores calculados:")
            print(alphas)

            print("Autovalores previstos:")
            print(preview_ev)

            print()

            print("Autovetores:")
            print(V)

            print()
            print("Av == lambda v")
            for value, vector in zip(alphas, V):
                print(np.matmul(original_mat, vector), " == ", value * vector)

            print()
            print("AA^T=")
            I = np.matmul(V, np.transpose(V))
            # limpar o array de elementos menores que nosso epsilon
            I[I < epsilon] = 0
            print(I)

        elif question == 'c':
            # inicializamos as variáveis
            edge_data = []
            non_fixed_bars = 0

            rho = 0
            Area = 0
            E = 0

            # e aqui vamos achar seus valores reais
            with open('./input-c') as f:
                lines = f.readlines()
                shape_info = lines[0].split()
                const_info = lines[1].split()

                vertex_count = int(shape_info[0])
                non_fixed_vertex_count = int(shape_info[1])
                edge_count = int(shape_info[2])

                rho = float(const_info[0])
                Area = float(const_info[1])
                E = float(const_info[2])

                def edge_input(a, b, theta, length):
                    return int(a), int(b), float(theta), float(length)

                for i in range(edge_count):
                    data = lines[2 + i].split()
                    data = edge_input(*data)
                    if data[1] > non_fixed_vertex_count:
                        edge_data.append(Bar(*data, True))
                    else:
                        edge_data.append(Bar(*data))
                        non_fixed_bars += 1

            K = np.zeros([non_fixed_bars, non_fixed_bars])
            #M = np.zeros([non_fixed_bars, non_fixed_bars])

            # pega o índice correto para K da estrutura
            def get_index(n, bar):
                if n % 4 == 0:
                    index = 2 * bar.i - 1
                elif n % 4 == 1:
                    index = 2 * bar.i
                elif n % 4 == 2:
                    index = 2 * bar.j - 1
                else:
                    index = 2 * bar.j
                # começam em 1
                return index - 1


            # gera K da estrutura
            for i, bar in enumerate(edge_data):
                K_ij = get_bar_K_ij(bar, Area, E)

                for rows in range(np.size(K_ij[:, 0])):
                    for columns in range(np.size(K_ij[0, :])):
                        k_i, k_j = get_index(rows, bar), get_index(columns, bar)
                        # para os componentes relativos aos vértices 13 e 14
                        if bar.fixed_edge and max(k_i, k_j) >= non_fixed_bars:
                            continue
                        num = K_ij[rows, columns]
                        K[k_i, k_j] += num

                # atualiza M
                #mass = 0.5 * rho * Area * bar.length

                #M[2 * bar.i - 2, 2 * bar.i - 2] += mass
                #M[2 * bar.i - 1, 2 * bar.i - 1] += mass

                #if bar.j > non_fixed_vertex_count:
                    #continue

                #M[2 * bar.j - 2, 2 * bar.j - 2] += mass
                #M[2 * bar.j - 1, 2 * bar.j - 1] += mass

            # abaixo o mesmo algoritmo das questões anteriores
            epsilon = 1e-6

            original_mat = K.copy()

            size = np.size(K[0, :])
            H = np.eye(size)
            for i in range(size-2):
                w = w_i(K, i)
                sandwich(K, i, w)
                Hw = np.zeros(size - 1)
                for j in range(len(w)):
                    Hw[size - 2 - j] = w[len(w) - 1 - j]
                H_wbar_right(H, 0, Hw)

            alphas = np.zeros(size)
            betas = np.zeros(size - 1)

            for i in range(size):
                alphas[i] = K[i, i]
                if i < size - 1:
                    betas[i] = K[i + 1, i]

            gammas = betas.copy()
            deltas = np.zeros(size - 2)
            k, V = QR(epsilon, True, alphas, betas, gammas, deltas, H.copy())

            V = np.transpose(V)

            #print(f"Matriz de escolha da tarefa {question.upper()}:")
            #print(original_mat)
            #print()
            #print("Autovalores calculados:")
            #print(alphas)

            #print("Autovetores:")
            #print(V)

            #print()
            print("Av == lambda v")
            for value, vector in zip(alphas, V):
                print(np.matmul(original_mat, vector), " == ", value * vector)

            #print()
            print("AA^T=")
            I = np.matmul(V, np.transpose(V))
            # limpar o array de elementos menores que nosso epsilon
            I[I < epsilon] = 0
            print(I)

            # raíz de w^2 gerada pelo algoritmo
            for i in range(np.size(alphas)):
                alphas[i] = math.sqrt(alphas[i])

            # ordenação dos vetores para encontrarmos as menores frequências
            idx_ordered = np.argsort(alphas)
            alphas = alphas[idx_ordered]
            V = V[idx_ordered, :]

            print("5 menores frequências de vibração")
            print(alphas[:5])

            print()

            vibe_modes = []
            for i in range(5):
                vibe_modes.append(np.matmul(H, V[i]))
            vibe_modes = np.array(vibe_modes)

            print("Modos de vibração correspondentes")
            print(vibe_modes)


if __name__ == "__main__":
    main()