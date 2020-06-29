function allocate(
    dtype :: DataType,
    dims...
)

    return SharedArray{dtype}(dims...)

end
