import math
def bob_func(k,top,bot=2):
    l=k-1
    x=0
    for i in range(bot,top+1):
        x=x+1./math.pow(i,l)
    if l<=0:
        print('Infinity')
    else:
#        print x
        x=x+(-1./(1.-l))*math.pow(top,1.-l)
        print x
        return x
    
