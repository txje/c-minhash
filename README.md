MinHash
=======

A simple and fast MinHash implementation in C, with a Python wrapper.


Install
-------

C:

    gcc minhash.c -o minhash
    gcc mh_fasta.c -o mh_fasta

Python:

    python setup.py install


Usage
-----

Command line:

    minhash <seq> <tile size> <seed>
    mh_fasta <query_fasta> <target_fasta> <k> <h> <seed> <threshold> > <output>

Python:

    import minhash
    hash_val, kmer_pos = minhash.minhash("stringy", 3, random.randint(0,2**32-1))


minhash.c
---------

Performs only a single pass (one hash function) over all tiles of given size in the string.

Returns: minimum hash value (uint32), position of tile with min hash value


mh_fasta.c
----------

Parameters: &lt;query_fasta> &lt;target_fasta> &lt;k> &lt;h> &lt;seed> &lt;threshold>

Query and target files may be the same, in which case it will do pairwise comparison

Results should be piped to an output file, and are in the form:
&lt;query_idx>,&lt;query_reverse?>,&lt;target_idx>,&lt;# matches>,&lt;offset>

Where query and target are indices into their respective FASTA, query_reverse is 1(reverse complemented) or 0, # matches is the integer total of matched k-mers (out of a maximum of h), and the offset is the computed average offset between matched k-mers.

Values >16 for k are not allowed because I'm packing the k-mer into a uint32. This may change in the future.
