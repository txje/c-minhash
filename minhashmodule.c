#include "minhash.h"
#include <Python.h>

static PyObject* minhash_wrapper(PyObject * self, PyObject * args) {
  char* s;
  unsigned int k;
  unsigned int seed;
  Min result;
  PyObject* ret = PyList_New(2);

  if (!PyArg_ParseTuple(args, "sII", &s, &k, &seed)) {
    return NULL;
  }

  result = minhash(s, k, seed);

  PyObject* hash_val = PyLong_FromUnsignedLong(result.hash);
  PyObject* kmer_pos = PyLong_FromUnsignedLong(result.pos);
  PyList_SetItem(ret, 0, hash_val);
  PyList_SetItem(ret, 1, kmer_pos);

  return ret;
}

static PyMethodDef MinHashMethods[] = {
  { "minhash", minhash_wrapper, METH_VARARGS, "MinHash" },
  { NULL, NULL, 0, NULL }
};

DL_EXPORT(void) initminhash(void) {
  Py_InitModule("minhash", MinHashMethods);
}
