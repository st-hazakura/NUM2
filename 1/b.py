import numpy as np

def nejv_nd_prvek(matice):
    pocet_radku = len(matice)
    max_ned_prvek, sum_nediag_prvku, radek, sloupec = 0, 0.0, 0, 0

    for i in range(pocet_radku - 1):
        for j in range(i+1, pocet_radku):
            sum_nediag_prvku = sum_nediag_prvku + matice[i,j]**2 #abs(matice[i,j])  # matice[i,j]**2 ???
            if abs(matice[i,j]) > max_ned_prvek:
                max_ned_prvek = abs(matice[i,j])
                radek = i
                sloupec = j
    return radek, sloupec, sum_nediag_prvku

def vypoc_uhlu(matice, radek, sloupec):
    w = - matice[radek, sloupec]
    m = ( matice[radek,radek] - matice[sloupec, sloupec] ) / 2.0
    v = np.sqrt(w**2.0 + m**2.0)

    c = np.sqrt(( v + abs(m) ) / (2.0 * v))
    s = (np.sign(m) * w) / (2.0 * v * c)
    
    return c, s

def rotace(matice, radek, sloupec, s, c):
    a_ss= matice[sloupec, sloupec]
    a_rr = matice[radek,radek]
    a_rs = matice[radek, sloupec]
    
    for i in range (len(matice)):
        if i!= radek and i!= sloupec:
            matice[i,radek] = matice[radek, i] = c * matice[i, radek] - s * matice[i, sloupec]
            matice[i,sloupec] = matice[sloupec, i] = c * matice[i, sloupec] + s * matice[i, radek]
    
    matice[radek, radek] = c**2 * a_rr + s**2 * a_ss - 2.0 * c * s * a_rs
    matice[sloupec, sloupec] = c**2 * a_ss + s**2 * a_rr + 2.0 *c * s * a_rs
    matice[radek, sloupec] = 0.0
    matice[sloupec,radek] = 0.0
    return matice

def vytv_aktualiz_t(pocet_prvku, radek, sloupec, s, c):
    t = np.eye(pocet_prvku, dtype= 'f' )
    t[radek,radek] = t[sloupec,sloupec] = c
    t[radek,sloupec]  = s
    t[sloupec, radek] = -s    
    return t

def aktualiz_q(q,t):
    if q is None:
        return t
    else:
        return np.dot(q,t)
    
    
def Jakobian(matice, presnost, max_iter):
    
    q = None
    pocet_prvku = len(matice)
    
    for i in range(max_iter):
        
        radek, sloupec, sum_nediag_prvku = nejv_nd_prvek(matice)
        if sum_nediag_prvku*2.0 <= presnost:
            vlastni_cisla = np.diag(matice)
            vlastni_vektory = q
            return vlastni_cisla, vlastni_vektory, i
        
        c,s = vypoc_uhlu(matice, radek, sloupec)
        matice = rotace(matice, radek, sloupec, s, c)
        t = vytv_aktualiz_t(pocet_prvku, radek, sloupec, s, c)
        q = aktualiz_q(q, t)
        
    vlastni_cisla = np.diag(matice)
    vlastni_vektory = q
        
    return vlastni_cisla, vlastni_vektory, i

def main():
    matice = np.array([
        [ 3, -1,  1], 
        [-1,  5, -1], 
        [ 1, -1,  3]], dtype= 'f' )

    presnost = 1e-6
    max_iter = 1000

    vlastni_cisla, vlastni_vektory, pocet_iteraci = Jakobian(matice, presnost, max_iter)
    # VC = [round(i) for i in vlastni_cisla ]
    print( f"VC{ vlastni_cisla}" )
    print(vlastni_vektory)
    
if __name__ == "__main__":
    main()