mutable struct ParallelPlan

    pids           :: AbstractArray{Integer, 1}
    tmd_slave_pids :: AbstractArray{Integer, 1}

    #
    # `visible_rng` is the range slave can view the master variable 
    # The range index is the index of master variable. It is used to
    # create a view of master variable.
    #
    # `update_rng` is the range slave should update the master variable.
    # The range index is the index of slave variable. It is used to
    # create a view of slave variable.
    #

    visible_rngs_T :: Array{UnitRange}
    visible_rngs_V :: Array{UnitRange}
    update_rngs_T  :: Array{UnitRange}


    function ParallelPlan(
        Ny      :: Int64;
        pids    :: AbstractArray{Int64, 1} = workers(),
    )
            
        if length(pids) >= 2
            tmd_slave_pids = pids[2:end]
        elseif length(pids) == 1
            tmd_slave_pids = pids
        else
            throw(ErrorException("No available workers!"))
        end

        visible_rngs_T, visible_rngs_V, update_rngs_T = calParallizationRange(N=Ny, P=length(tmd_slave_pids))

        return new(
            pids,
            tmd_slave_pids,
            visible_rngs_T,
            visible_rngs_V,
            update_rngs_T,
        ) 

    end

end

function calParallizationRange(;
    N :: Integer,     # Total T grids
    P :: Integer,     # Number of procs
)
    if N < P
        throw(ErrorException("Condition must be satisfied: N >= max(1, L) * P"))
    end

    n̄ = floor(Integer, N / P)
    R = N - n̄ * P


    visible_rngs_T = Array{Union{UnitRange, Nothing}}(undef, P)
    visible_rngs_V = Array{Union{UnitRange, Nothing}}(undef, P)
    update_rngs_T  = Array{Union{UnitRange, Nothing}}(undef, P)
    


    if P == 1   

        visible_rngs_T[1] = 1:N
        visible_rngs_V[1] = 1:N+1
        update_rngs_T[1]  = 1:N

    else 
   
        # `cnt` represent the starting update grid 
        cnt = 1
        for p = 1:P

            m = (p <= R) ? n̄ + 1 : n̄  # assigned grids

            if p == 1

                visible_rngs_T[p] = cnt:cnt+(m-1)+L
                visible_rngs_V[p] = cnt:cnt+(m-1)+L+1
                update_rngs_T[p]  = 1:m

            elseif p == P
 
                visible_rngs_T[p] = cnt-L:cnt+(m-1)
                visible_rngs_V[p] = cnt-L:cnt+(m-1)+1
                update_rngs_T[p]  = L+1:L+m
           
            else

                visible_rngs_T[p] = cnt-L:cnt+(m-1)+L
                visible_rngs_V[p] = cnt-L:cnt+(m-1)+L+1
                update_rngs_T[p]  = L+1:L+m

            end
            cnt += m

        end

        if (cnt-1) != N

            throw(ErrorException("Error: Range calculation wrong."))

        end
    end

    return visible_rngs_T, visible_rngs_V, update_rngs_T

end

