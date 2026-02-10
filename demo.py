from smith import *
def parse_p(s):
    parts = s.strip().split()
    clist = []
    for p in parts: clist.append(complex(p))
    return poly(clist)
def parse_r(s):
    if '/' in s:
        bits = s.split('/')
        return rtnlf(parse_p(bits[0]), parse_p(bits[1]))
    return rtnlf(parse_p(s))
def main():
    print("1. Integers  2. Gaussian Integers  3. Polynomials  4. Proper Rational Functions")
    c = input("Choice: ")
    raw = input("Enter matrix (Rows ';', Elements ','):\n")
    row_strings = raw.split(';')
    A = []
    for rs in row_strings:
        elements = rs.split(',')
        row_data = []
        for e in elements:
            e = e.strip()
            if c == "1": row_data.append(int(e))
            elif c == "2":
                tmp = complex(e)
                row_data.append(gint(tmp.real, tmp.imag))
            elif c == "3": row_data.append(parse_p(e))
            elif c == "4": row_data.append(parse_r(e))
        A.append(row_data)
    if c == "1": d = intd()
    elif c == "2": d = gssd()
    elif c == "3": d = pld()
    else: d = rtnd()
    res = snf(A, d)
    print("\nResult:")
    for row in res:
        line = ""
        for item in row:
            if line != "": line += " , "
            line += str(item)
        print("[" + line + "]")
if __name__ == "__main__":
    main()