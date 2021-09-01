# classe para um objeto qualquer que vamos estar lidando
class thing:
    def __init__(self,Len,H,Wid,W):
        self.Len = Len
        self.H = H
        self.Wid = Wid
        self.W = W

# classe para um veículo
class vehicle(thing):
    def __init__(self,maxLen,maxH,maxWid,maxW):
        thing.__init__(self,maxLen,maxH,maxWid,maxW)
        self.inside = list()
        self.id = ""
    
    def clone(self):
        veh = vehicle(self.Len,self.H,self.Wid,self.W)
        veh.id = self.id
        return veh
    
    def fit(self,items):
        # se passarem do peso não há nada ser feito
        if self.W < sum([i.W for i in items]):
            return False

        # então nós organizamos por área e vamos fazer height fitting
        items = sorted(items, key=lambda i: i.Len*i.Wid)

        # caso o maior dos objetos não caiba em algum dos quesitos não há nada a ser feito
        if items[0].Len > self.Len*self.Wid/items[0].Wid:
            return False

        # lista com os objetos posicionados
        pos = list()

        while len(items) != 0:
            pos.append([])
            # altura disponível
            rmH = self.H

            row = list()
            # primeiro vamos fazer height filling usando a área
            for i in range(len(items)):
                j = i - len(row)
                if items[j].H < rmH:
                    rmH -= items[j].H
                    row.append(items.pop(j))
                else:
                    break

            if len(pos) == 0:
                return False

            # agora fazemos row filling usando largura
            row = sorted(row, key=lambda i: i.Wid)
            rmLen = self.Len
            for i in range(len(row)):
                if row[i].Len < rmLen:
                    rmLen -= row[i].Len
                else:
                    break
        
        return True


# classe para objetos que serão transportados
class item(thing):
    def __init__(self,lt):
        thing.__init__(self,lt[0],lt[1],lt[2],lt[3])
        self.id = -1

    def __str__(self):
        return f"len: {self.Len}, h: {self.H}, wid: {self.Wid}, w {self.W}"

    def clone(self):
        i = item([self.Len,self.H,self.Wid,self.W])
        i.id = self.id
        return i

# classe para cada empresa de transporte
class company:
    def __init__(self,baseValue,vehicles):
        self.vehicles = vehicles

        # suponha que cada carro é um vetor em R4
        # então podemos tomar seu módulo (estimado a partir da sua capacidade)
        # e estabelecer valores (possivelmente) relativos entre os veículos
        self.fares = list()
        for v in vehicles:
            self.fares.append((baseValue)*(v.Len**2 + v.Wid**2 + v.H**2 + v.W**2)**(1/2))

    def estimate(self,which,weight):
        return self.fares[which]*(1 + weight/self.vehicles[which].W)


