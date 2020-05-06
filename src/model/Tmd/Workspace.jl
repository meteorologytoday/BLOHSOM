mutable struct Workspace

    Nx :: Int64
    Ny :: Int64
    Nz :: Int64

    T :: Array
    U :: Array
    V :: Array
    W :: Array

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
        W :: Int64=0,
        sT :: Int64=0,
        sU :: Int64=0,
        sV :: Int64=0,
    )


        _T = []
        _U = []
        _V = []
        _W = []

        _sT = []
        _sU = []
        _sV = []

        ptr = Dict(
            :T => 1,
            :U => 1,
            :V => 1,
            :W => 1,
            :sT => 1,
            :sU => 1,
            :sV => 1,
        )

        return new(
            Nx, Ny, Nz,
            _T, _U, _V, _W,
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

#    println(string(grid), " => ", i, "; length(list) = ", length(list))

    if i > length(list)
        println("Running out of workspace of " * string(grid) * ", create new...")
        push!(list, genEmptyGrid(wksp, Float64, grid))
    end
    
    wksp.ptr[grid] += 1

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

function genEmptyGrid(
    wksp  :: Workspace,
    dtype :: DataType,
    grid  :: Symbol,
)
    Nx, Ny, Nz = wksp.Nx, wksp.Ny, wksp.Nz
    dim = Dict(
        :T => (Nz,   Nx, Ny  ),
        :U => (Nz,   Nx, Ny  ),
        :V => (Nz,   Nx, Ny+1),
        :W => (Nz+1, Nx, Ny  ),
        :sT => ( Nx, Ny  ),
        :sU => ( Nx, Ny  ),
        :sV => ( Nx, Ny+1),
    )[grid]

    return zeros(dtype, dim...)

end
