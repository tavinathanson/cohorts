import copy

# courtesy of https://stackoverflow.com/a/8714242/3457743

DictProxyType = type(object.__dict__)

def make_hash(o):

  """
  Makes a hash from a dictionary, list, tuple or set to any level, that 
  contains only other hashable types (including any lists, tuples, sets, and
  dictionaries). In the case where other kinds of objects (like classes) need 
  to be hashed, pass in a collection of object attributes that are pertinent. 
  For example, a class can be hashed in this fashion:

    make_hash([cls.__dict__, cls.__name__])

  A function can be hashed like so:

    make_hash([fn.__dict__, fn.__code__])
  """

  if type(o) == DictProxyType:
    o2 = {}
    for k, v in o.items():
      if not k.startswith("__"):
        o2[k] = v
    o = o2  

  if isinstance(o, (set, tuple, list)):

    return tuple([make_hash(e) for e in o])    

  elif not isinstance(o, dict):

    return hash(o)

  new_o = copy.deepcopy(o)
  for k, v in new_o.items():
    new_o[k] = make_hash(v)

  return hash(tuple(frozenset(sorted(new_o.items()))))
