import math

#Lenos Trigeorgis
def trigeorgisTree(K, T, S, sig, r, N, PorC):
    dnt = T / N;
    u = r - 0.5 * sig * sig;
    x = math.sqrt(sig * sig * dnt + (u * dnt) ** 2);
    dxu = math.exp(x);



    dxd = 1 / dxu;
    pu = 0.5 + 0.5 * (u * dnt / x);
    pd = 1 - pu;
    disc = math.exp(-r * dnt);

    St = [0] * (N + 1)
    C = [0] * (N + 1)

    St[0] = S * dxd ** N;

    for j in xrange(1, N + 1):
        St[j] = St[j - 1] * dxu / dxd;

    for j in xrange(1, N + 1):
        if PorC == 'p':
            C[j] = max(K - St[j], 0);
        elif PorC == 'c':
            C[j] = max(St[j] - K, 0);

    for i in xrange(N, 0, -1):
        for j in xrange(0, i):
            C[j] = disc * (pu * C[j + 1] + pd * C[j]);

    return C[0]


def CRRTree(K, T, S, sig, r, N, PorC):
    dnt = T / N;
    dxu = math.exp(sig * math.sqrt(dnt));
    dxd = math.exp(-sig * math.sqrt(dnt));
    pu = ((math.exp(r * dnt)) - dxd) / (dxu - dxd);
    pd = 1 - pu;
    disc = math.exp(-r * dnt);

    St = [0] * (N + 1)
    C = [0] * (N + 1)

    St[0] = S * dxd ** N;

    for j in xrange(1, N + 1):
        St[j] = St[j - 1] * dxu / dxd;

    for j in xrange(1, N + 1):
        if PorC == 'p':
            C[j] = max(K - St[j], 0);
        elif PorC == 'c':
            C[j] = max(St[j] - K, 0);

    for i in xrange(N, 0, -1):
        for j in xrange(0, i):
            C[j] = disc * (pu * C[j + 1] + pd * C[j]);

    return C[0]


def JRTree(K, T, S, sig, r, N, PorC):
    dnt = T / N;
    dxu = math.exp((r - (sig ** 2 / 2)) * dnt + sig * math.sqrt(dnt));
    dxd = math.exp((r - (sig ** 2 / 2)) * dnt - sig * math.sqrt(dnt));
    pu = 0.5;
    pd = 1 - pu;
    disc = math.exp(-r * dnt);

    St = [0] * (N + 1)
    C = [0] * (N + 1)

    St[0] = S * dxd ** N;

    for j in xrange(1, N + 1):
        St[j] = St[j - 1] * dxu / dxd;

    for j in xrange(1, N + 1):
        if PorC == 'p':
            C[j] = max(K - St[j], 0);
        elif PorC == 'c':
            C[j] = max(St[j] - K, 0);

    for i in xrange(N, 0, -1):
        for j in xrange(0, i):
            C[j] = disc * (pu * C[j + 1] + pd * C[j]);

    return C[0]
