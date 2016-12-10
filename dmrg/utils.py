def index_map(array):
    d = {}
    for index, value in enumerate(array):
        d.setdefault(value, []).append(index)
    return d

