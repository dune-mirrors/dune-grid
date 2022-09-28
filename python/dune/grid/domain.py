from ._grids import CartesianDomain

def cartesian(lower, upper, division, **parameters):
    return CartesianDomain(lower,upper,division,True,**parameters)
