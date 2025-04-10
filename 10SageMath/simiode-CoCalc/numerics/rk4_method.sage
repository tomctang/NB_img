def rk4_method(f,ic,stepsize,finalT):
    #Basic RK4 Method implementation. If number of steps from initial to final time is not an integer, last step is incremental of proper length.
    #Also, can solve backwards in time. Uses transformation t->-t, so u'(t) = f(t,u(t)) with u(a) = b becomes v'(t) = -f(t,v(t)) with v(-a) = b under time reversal.

    if finalT > ic[0]: #Solve forward in time
        F = 1
        t0 = ic[0]
        u0 = ic[1]
        tf = finalT
        f2 = lambda t,u: f(t,u)
    else:              #Solving backward in time
        F = -1
        t0 = -ic[0]
        u0 = ic[1]
        tf = -finalT
        f2 = lambda t,u: -f(-t,u)

    #Make sure everyone is a float
    t0 = t0.n(); u0 = u0.n(); tf = tf.n();

    N = floor((tf-t0)/stepsize) #Round down
    T1 = t0+N*stepsize #Corresponding time
    Nf = N+1 #Number of steps we'll take, unless T1 = finalT
    if abs(T1-tf) < 1.0e-8: #To some tolerance
        Nf = N
    sol = matrix(SR,Nf+1,2) #Solution time and values
    sol[0,0] = t0 #Initial time
    sol[0,1] = u0 #Initiai value

    #March out in time
    for i in range(1,Nf):
        sol[i,0] = sol[i-1,0]+stepsize
        m1 = f2(sol[i-1,0],sol[i-1,1])
        m2 = f2(sol[i-1,0]+0.5*stepsize,sol[i-1,1]+0.5*stepsize*m1)
        m3 = f2(sol[i-1,0]+0.5*stepsize,sol[i-1,1]+0.5*stepsize*m2)
        m4 = f2(sol[i-1,0]+stepsize,sol[i-1,1]+stepsize*m3)
        m = (m1+2*m2+2*m3+m4)/6.0
        sol[i,1] = sol[i-1,1] + stepsize*m

    #Last step, might be special.
    sol[Nf,0] = tf
    step2 = tf-sol[Nf-1,0]
    m1 = f2(sol[Nf-1,0],sol[Nf-1,1])
    m2 = f2(sol[Nf-1,0]+0.5*step2,sol[Nf-1,1]+0.5*step2*m1)
    m3 = f2(sol[Nf-1,0]+0.5*step2,sol[Nf-1,1]+0.5*step2*m2)
    m4 = f2(sol[Nf-1,0]+step2,sol[Nf-1,1]+step2*m3)
    m = (m1+2*m2+2*m3+m4)/6.0
    sol[Nf,1] = sol[Nf-1,1] + step2*m

    #If necessary, negate t values and return
    if F == -1:
        for i in range(0,Nf+1):
            sol[i,0] = -sol[i,0]

    return sol
