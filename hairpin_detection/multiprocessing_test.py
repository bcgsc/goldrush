from multiprocessing import Pool

def f(x):
    print(x)
    return "DONE"

if __name__ == '__main__':
    with open("/dev/stdin", 'r') as fin:
        with Pool(5) as p:
            result = p.imap(f, (line.strip() for line in fin))
            for r in result:
                print(r)

