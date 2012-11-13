A simple extension to the Event and MultiEvent classes defined in Survival.py.
Defines __Event and __MulitEvent classes through the Python C-Api in order to
allow for the registration of the defined classes with the Numpy system enambling
the creation of ufuncs to replace the logp and random methods also defined in
Survival.py such that broadcasting is still supported.
