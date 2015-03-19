MinHash
=======

A simple and fast MinHash implementation in C, with a Python wrapper.

Performs only a single pass (one hash function) over all tiles of given size in the string.

Returns: minimum hash value (uint32), position of tile with min hash value


Install
-------

C: gcc minhash.c -o minhash

Python: python setup.py install


Usage
-----

Command line: minhash <seq> <tile size> <seed>

Python:
  import minhash
  hash_val, kmer_pos = minhash.minhash("stringy", 3, random.randint(0,2**32-1))
