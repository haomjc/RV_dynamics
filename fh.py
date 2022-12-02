def fh(x):
    #b=1*10**(-6)
    b=1/2*10**(-6)
    if x > b:
        fhx = x - b
    elif abs(x) < b:
        fhx = 0
    else:
        fhx = x + b
    
    return fhx


