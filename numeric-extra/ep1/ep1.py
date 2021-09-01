import numpy as np
import matplotlib.pyplot

def GivensRotation(alphas, betas, gammas, deltas, i):
    # aplica a rotação de givens com as dicas do enunciado
    j = i + 1
    (sk, ck) = SkCk(alphas[i], betas[i])

    # muda sempre no maximo 3 colunas
    for k in range(3):
        if k == 2 and i != len(betas) - 1:
            delta_antigo = deltas[i]
            gamma_antigo = gammas[j]
            gammas[j] = sk*delta_antigo + ck*gamma_antigo
            deltas[i] = ck*delta_antigo - sk*gamma_antigo

        elif k == 1:
            gamma_antigo = gammas[i]
            alpha_antigo = alphas[j]
            gammas[i] = ck*gamma_antigo - sk*alpha_antigo
            alphas[j] = sk*gamma_antigo + ck*alpha_antigo


        if k == 0:
            alpha_antigo = alphas[i]
            beta_antigo = betas[i]
            betas[i] = sk*alpha_antigo + ck*beta_antigo
            alphas[i] = ck*alpha_antigo - sk*beta_antigo

    return sk, ck

def SkCk(a, b):
    # calcula sk e ck da mesma forma que os exemplos no enunciado
    dem = (a*a + b*b)**(1/2) # diagonal
    # projeções
    ck = a/dem
    sk = -b/dem
    return (sk, ck) # vetor

def sign(n):
    # função sinal
    if n >= 0:
        return 1
    else:
        return -1

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



def QR(eps, wilkinson, alphas, betas, gammas, deltas):
    # algoritmo iterativo qr para achar autovetores e autovalores

    n = len(alphas)
    V = np.identity(n)
    k = 0

    for m in range(len(betas), 0, -1):
        while abs(betas[m-1]) >= eps:
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
                (sks[i] , cks[i]) = GivensRotation(alphas_uk , betas, gammas, deltas, i)

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


first_time = True # entra no while
user = ''
while user in ['A', 'B', 'C'] or first_time == True:
    user = input('Escolha a tarefa A, B ou C. Outro botão sai do loop \n')
    first_time = False
    user = user.upper()

    if user == 'A':
        possible_n = [4, 8, 16, 32]
        epsilon = 1e-6
        for n in possible_n:
            second = False
            wilkinson = False
            while not second:
                alphas = n*[2]
                betas = (n-1)*[-1]
                gammas = list(betas)

                deltas = (len(alphas) - 2)*[0]

                k, V = QR(epsilon, wilkinson, alphas, betas, gammas, deltas)

                print('n =', n, 'k =', k, 'Wilkinson =', wilkinson)

                if wilkinson:
                    second = True
                if not wilkinson:
                    wilkinson = True

                idx_ordered = np.argsort(alphas)

                alphas = np.array(alphas)[idx_ordered]
                V = V[:, idx_ordered]

                eigenvectors = []
                for j in range(1, n+1):
                    i = np.arange(1, n+1)
                    correct_values = np.sin((i*j*np.pi)/(n+1))
                    eigenvectors.append(correct_values)

                eigenvectors = np.array(eigenvectors)
                for i in range(len(eigenvectors)):
                    eigenvectors[:, i] = V[0, i]*eigenvectors[:, i]/eigenvectors[0, i]
                print('Maior diferença para os autovetores corretos =')
                print(np.max(abs(V-eigenvectors)))

                print('Maior diferença para os autovalores corretos = ')
                diffs = []
                for j in range(1, n+1):
                    diffs.append(2*(1-np.cos((np.pi*j)/(1+n)))- alphas[j-1])
                print(max(diffs))
                print()




    if user == 'B' or user == 'C':
        initial_states = [[-2, -3, -1, -3, -1],
                [1, 10, -4, 3, -2],
                [0.18933466, -0.39110521, 0.55766111,-0.588202, 0.39271058]]
        wilkinson = True
        epsilon = 1e-15
        mass = 2

        if user == 'B':
            ki = [42, 44, 46, 48, 50, 52]
            n = 5
        else:
            ki = [38, 40, 38, 40, 38, 40, 38, 40, 38, 40, 38]
            n = 10


        for jj in range(3):
            x_init = list(initial_states[jj])

            if user == 'C':
                if x_init == [0.18933466, -0.39110521, 0.55766111,-0.588202, 0.39271058]:
                    x_init = [0.12524522, -0.22901235, 0.32440201, -0.38597125, 0.42149319,
                            -0.42149319, 0.38597125, -0.32440201, 0.22901235, -0.12524522]
                else:
                    for i in range(5):
                        x_init.append(x_init[i])

            alphas = []
            betas = []
            deltas = (n-2)*[0]

            for i in range(n-1):
                betas.append(-ki[i+1]/mass)

            for i in range(n):
                alphas.append((ki[i] + ki[i+1])/mass)

            gammas = list(betas)

            k, V = QR(epsilon, wilkinson, alphas, betas, gammas, deltas)
            y_t = V.T.dot(x_init)

            freqs = [alphas[i]**0.5 for i in range(n)]

            print('Item', user, ' X(0) =', x_init)
            print('Alphas obtidos, que são o quadrado das frequências angulares =', alphas)

            x_t = np.array([v*y_t for v in V])

            t = np.array(range(10*400))/400

            fig, axs = matplotlib.pyplot.subplots(len(x_t))
            soma_total = []
            for i in range(len(x_t)):
                soma_total.append(0)
                for j in range(len(x_t[0])):
                    soma_total[i] += x_t[i, j]*np.cos(freqs[j]*t)

            for i in range(n):
                axs[i].plot(t, soma_total[i])

            strings = []
            for i in range(n):
                strings.append('')
                for j in range(len(x_t[0])):
                    strings[i] += str(round(x_t[i, j], 4))
                    strings[i] += 'cos(' + str(round(alphas[j]**(1/2), 3)) + 't) + '
                print(f'x_{i+1} =', strings[i-1][:-3])

            title = 'Tarefa ' + user + '   X(0) =' + str(x_init)
            fig.suptitle(title)
            fig.savefig(f"{user}_{jj}")