from distutils.core import setup, Extension
setup(name="__Survival", version="1.0",
      ext_modules=[Extension("__Survival", ["__Survival.c"])])