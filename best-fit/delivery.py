from deliveryTools import *

# caso executemos esse arquivo, realizamos um pequeno teste
if __name__ == "__main__":
    # inicializo os veículos de cada empresa manualmente
    LalaCars = list()

    LalaCars.append(vehicle(35,40,30,20))
    LalaCars.append(vehicle(188,133,108,500))
    LalaCars.append(vehicle(300,180,200,1500))
    LalaCars[0].id = "Moto"
    LalaCars[1].id = "Fiorino"
    LalaCars[2].id = "Carreto"

    Lala = company(0.001, LalaCars)

    OgiCars = list()

    OgiCars.append(vehicle(52,36,52,20))
    OgiCars.append(vehicle(125,80,60,200))
    OgiCars[0].id = "Moto"
    OgiCars[1].id = "SUV"

    Ogi = company(0.001, OgiCars)

    currDS = importData(r"./meow.csv")

    print(currDS[0].W)
    for i in range(len(currDS)):
        currDS[i].id = i
    
    currDS = sortItems(currDS)
    for i in currDS:
        print(i)
    
    print()

    print("Empresa: Ogi")
    lpOgi = listPossibilities(Ogi,currDS)
    for i in range(len(lpOgi)):
        w = lpOgi[i][1]
        if w == 0:
            continue
        print(f"Para os itens dados, temos R${10000*Ogi.estimate(i,w)//1/100}/km no carro {Ogi.vehicles[i].id}, com {10000*(w/Ogi.vehicles[i].W)//1/100}% de ocupação.\n")

    print("Empresa: Lala")
    lpLala = listPossibilities(Lala,currDS)
    for i in range(len(lpLala)):
        w = lpLala[i][1]
        if w == 0:
            continue
        print(f"Para os itens dados, temos R${10000*Lala.estimate(i,w)//1/100}/km no carro {Lala.vehicles[i].id}, com {10000*(w/Lala.vehicles[i].W)//1/100}% de ocupação.\n")

