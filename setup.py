from distutils.core import setup, Extension

extension_mod = Extension("minhash", ["minhashmodule.c", "minhash.c"])

setup(name = "minhash", ext_modules=[extension_mod])
