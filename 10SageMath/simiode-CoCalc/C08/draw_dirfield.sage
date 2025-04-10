def draw_dirfield(f, rang, *args):
    #Usage: draw_dirfield(f, rang, *args)
    #Draws direction field for u'=f(t,u). Here rang is a list of the form
    #[tlow, thigh, ulow, uhigh]. If args is defined then ic = args[0]
    #and ic is an optional list of initial conditions for solution curves,
    #of the form ic = [[t0,u0],[t1,u1],...,[tn,un]]

    var('t u')

    #Pick of low, high in each variable, make sure all are floating point.
    tlow = rang[0].n(); thigh = rang[1].n();
    ulow = rang[2].n(); uhigh = rang[3].n();

    #Plot slope field itself with built-in Sage command
    p = plot_slope_field(f, (t,tlow,thigh), (u,ulow,uhigh));

    #Load in homegrown rk4 solve to sketch solution curves
    #Housman: commented out original line and replaced 
    load('../numerics/rk4_method.sage')
    #load("./rk4_method.sage")

    if len(args)>0:
        ic = args[0]
        N = len(ic)

        for i in range(0,N):
            t0 = ic[i][0]; u0 = ic[i][1];
            t0 = t0.n(); u0 = u0.n();
            h = (thigh-tlow)/40.0;
            rk4_results = rk4_method(f,[t0,u0],h,thigh)
            p2 = line(rk4_results,rgbcolor=[1,0,0]);
            p = p + p2

            rk4_results = rk4_method(f,[t0,u0],h,tlow)
            p3 = line(rk4_results,rgbcolor=[1,0,0]);
            p = p + p3

    show(p,xmin=tlow,xmax=thigh,ymin=ulow,ymax=uhigh)

    return
