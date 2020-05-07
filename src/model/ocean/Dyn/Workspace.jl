mutable struct Workspace

    T :: Array
    U :: Array
    V :: Array

    sT :: Array
    sU :: Array
    sV :: Array


    ptr :: Dict


    function Workspace(;
        Nx :: Int64,
        Ny :: Int64,
        Nz :: Int64,
        T :: Int64=0,
        U :: Int64=0,
        V :: Int64=0,
        sT :: Int64=0,
        sU :: Int64=0,
        sV :: Int64=0,
    )


        _T = []
        _U = []
        _V = []

        _sT = []
        _sU = []
        _sV = []

        for i=1:T
            push!(_T, zeros(Float64, Nx, Ny, Nz))
        end 

        for i=1:U
            push!(_U, zeros(Float64, Nx, Ny, Nz))
        end 

        for i=1:V
            push!(_V, zeros(Float64, Nx, Ny+1, Nz))
        end 

        for i=1:sT
            push!(_sT, zeros(Float64, Nx, Ny))
        end 

        for i=1:sU
            push!(_sU, zeros(Float64, Nx, Ny))
        end 

        for i=1:sV
            push!(_sV, zeros(Float64, Nx, Ny+1))
        end 


        ptr = Dict(
            :T => 1,
            :U => 1,
            :V => 1,
            :sT => 1,
            :sU => 1,
            :sV => 1,
        )


        return new(
            _T, _U, _V,
            _sT, _sU, _sV,
            ptr,
        )

    end
    

end

function getSpace!(
    wksp :: Workspace,
    grid :: Symbol;
    flat :: Bool = false
)
    i = wksp.ptr[grid]
    list = getfield(wksp, grid)
    wksp.ptr[grid] += 1

    if i > length(list)
        throw(ErrorException("Running out of workspace of " * string(grid)))
    end

    return ( flat ) ? view(list[i], :) : list[i]
end

function reset!(
    wksp :: Workspace,
    grid :: Symbol=:ALL,
)
    if grid == :ALL
        for k in keys(wksp.ptr)
            wksp.ptr[k] = 1
        end
    else
        wksp.ptr[grid] = 1
    end

end

