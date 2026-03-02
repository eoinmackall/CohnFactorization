using Oscar
using Base.Threads

function polynomial_bases_count_threaded(n::Integer; verbose::Bool=false)
    R, x = polynomial_ring(ZZ, "x")
    
    success_count = Atomic{Int}(0)
    test_count = Atomic{Int}(0)
    
    print_lock = ReentrantLock()
    
    poly_buffers = [R() for _ in 1:nthreads()]

    limit = isqrt(n)

    try
        chunk_size = max(1, cld(limit - 1, nthreads()))
        chunks = collect(Iterators.partition(2:limit, chunk_size))
    
        @threads for chunk in chunks            
            
            f = R() 
            for z in chunk
                if verbose
                    atomic_add!(test_count, 1)
                end

                zero!(f) 
                
                temp_n = n
                i = 0
                g = 0
                
                while temp_n > 0
                    temp_n, c = divrem(temp_n, z)
                    setcoeff!(f, i, c)
                    if g != 1
                        g = gcd(g, c)
                    end
                    i += 1
                end

                if g > 1
                    atomic_add!(success_count, 1)
                    if verbose
                        lock(print_lock) do
                            println(divisor)
                            println("Found at z = ", z)
                        end
                    end
                    continue
                end           

                factors = Nemo.factor(f)
                
                for (divisor, mult) in factors
                    if degree(divisor) < degree(f)
                        if verbose
                            lock(print_lock) do
                                println(divisor)
                                println("Found at z = ", z)
                            end
                        end
                        atomic_add!(success_count, 1)
                        break
                    end
                end
            end
        end
    finally
        if verbose
            println("Test ended on iteration: ", test_count[])
        end
    end
    return success_count[]
end

function factor_polynomial_bases(n::Integer; verbose::Bool=false)
    R, x = polynomial_ring(ZZ, "x")

    limit = isqrt(n)

    chunk_size = max(1, cld(limit - 1, nthreads()))
    chunks = collect(Iterators.partition(2:limit, chunk_size))
    
    found = Atomic{Bool}(false)
    fac_result = Ref{typeof(ZZ(n))}() 
    fac_lock = ReentrantLock()

    @threads for chunk in chunks            
        
        f = R() 
        for z in chunk
            if found[]
                break
            end

            if verbose
                atomic_add!(test_count, 1)
            end

            zero!(f) 
            
            temp_n = n
            i = 0
            g = 0
            
            while temp_n > 0
                temp_n, c = divrem(temp_n, z)
                setcoeff!(f, i, c)
                if g != 1
                    g = gcd(g, c)
                end
                i += 1
            end

            if g > 1
                lock(fac_lock) do
                    if !found[]
                        fac_result[] = ZZ(g)
                        found[] = true
                    end
                end
                break
            end           

            factors = Nemo.factor(f)
            
            for (divisor, mult) in factors
                if degree(divisor) < degree(f)
                    lock(fac_lock) do
                        if !found[]
                            fac_result[] = ZZ(evaluate(divisor, z))
                            found[] = true
                        end
                    end
                    break
                end
            end
        end
    end
    if found[]
        return (true, fac_result[])
    else
        return (false, ZZ(n))
    end
end
