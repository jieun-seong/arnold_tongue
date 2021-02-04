main()

function main()

    PREC = 128
    setprecision(PREC)
    zero = BigFloat(0)
    one = BigFloat(1)
    two = BigFloat(2)
    three = BigFloat(3)
    five = BigFloat(5)
    ten = BigFloat(10)

    function get_wt_array(L)
        wt_array = collect(1:1:L)
        wt_array = wt_array/BigFloat(L)
        wt_array = wt_array.*(wt_array.-1).*(-one)
        wt_array = wt_array.^(-1).*(-one)
        wt_array = exp.(wt_array)
        wt_array[L] = zero
        return wt_array
    end
    
    function map(x)
        temp = x+omega-lambda*sin(two*pi*x)/(two*pi)
        temp = temp-floor(temp)
    end
    
    function rotnum(omegahere)
        x = initx
        omega = omegahere
        for i in 1:N_t
            orbit_t[i] = step(x)
            x = map(x)
        end
        return dot(orbit_t, wt_array_t[m_t])/tot_wt_t[m_t] - target
    end
    
    function dfdx(x) return (one-lambda*cos(two*pi*x)) end

    function ddfdxx(x) return (two*pi*lambda*sin(two*pi*x)) end

    function dfdl(x) return (-sin(two*pi*x)/(two*pi)) end

    function ddfdll(x) return zero end

    function ddfdldx(x) return (-cos(two*pi*x)) end

    function dfdw(x) return one end

    function ddfdww(x) return zero end

    function ddfdwdx(x) return zero end

    function ddfdwdl(x) return zero end
    
    ITMAX = 100
    if (PREC==64)
        EPS = one/BigFloat(10^16)
    elseif (PREC==128)
        EPS = one/BigFloat(10^32)
    elseif (PREC==256)
        EPS = one/BigFloat(10^64)
    end #if
    
    function zbrent(f, x1, x2, tol)
        a=x1
        b=x2
        #c,d,e,min1,min2,fc,p,q,r,s,tol1,xm;
        fa=f(a)
        fb=f(b)
        if (fb*fa>zero) throw(ErrorException("Root must be bracketed in ZBRENT")) end
        fc=fb
        for iter in 1:ITMAX
            if (fb*fc>zero)
                c=a
                fc=fa
                e=d=b-a
            end
            if (abs(fc)<abs(fb))
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            end
            tol1=two*EPS*abs(b)+(tol/two)
            xm=(c-b)/two
            if (abs(xm)<=tol1||fb==zero)
              return b
            end
            if (abs(e)>=tol1 && abs(fa)>abs(fb))
                s=fb/fa;
                if (a==c)
                    p=two*xm*s
                    q=one-s
                else
                    q=fa/fc
                    r=fb/fc
                    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
                    q=(q-one)*(r-one)*(s-one)
                end
                if(p>zero)
                    q=-q
                end
                p=ABS(p)
                min1=three*xm*q-abs(tol1*q)
                min2=abs(e*q)
                if ((two*p)<(min1<min2?min1:min2))
                    e=d
                    d=p/q
                else
                    d=xm
                    e=d
                end
            else
                d=xm
                e=d
            end
            a=b
            fa=fb
            if (abs(d)>tol1)
                b=b+d
            else
                b=b+(xm>zero?abs(tol1):-abs(tol1))
            end
            fb=f(b)
        end
        throw(ErrorException("Maximum number of iterations exceeded in ZBRENT"))
        return zero
    end #zbrent

    println("#1:   #2:   #3:   #4:")










    for N in 1:1000

        wt_array = get_wt_array(2*N+1)
        tot_wt = sum(wt_array)

        omega = (sqrt(five)-one)/two
        x(p) = sin(omega*p)
        #x(p) = p+omega+BigFloat(15)/BigFloat(10)*sin(two*pi*p) # should be x(p+1) = x(p) + 1, not equal to x(p)
        n = collect(-BigFloat(N):one:BigFloat(N)) # need to be careful when taking mod -> discontinuity
        x_n_array = x.(n)

        for k in -2:2

            print(N, "  ", k, "  ")

            #Classical
            x_k = 1/(two*BigFloat(N)+one)*sum(x_n_array.*exp.(-two*pi*im*BigFloat(k)*omega*n))
            print(abs(x_k), "  ")

            #DSY with BigFloat
            x_k = sum(wt_array.*x_n_array.*exp.(-two*pi*im*BigFloat(k)*omega*n))/tot_wt
            println(abs(x_k))

        end #for k

    end #for N

end #function main