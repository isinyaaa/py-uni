class response:
    def __init__(self):
        self.result = -1;
        self.status = None

class fibo:    

    # get the n-th fibonacci number
    @classmethod
    def get(cls,n):
        res = response()
        if not (type(n) is int):
            try:
                n = int(n)
            except:
                res.status = "Invalid input"
                return response
        if n >= 0:
            res.result = cls.__fibo(n)
            res.status = "Success!"
        else:
            res.status = "Invalid input"
        return res
    
    @staticmethod
    def __fibo(n):
        a = b = 0
        c = 1
        for i in range(n-1):
            a = b
            b = c
            c = a + b
        return b

if __name__ == "__main__":
    ans = input("qual número você deseja? ")

    res = fibo.get(ans)
    print(f"n-th fibonacci number: {res.result}\nprogram status: {res.status}")
