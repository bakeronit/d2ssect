// The code for this module is based heavily on the code example count_in_file.cc given as part of the jellyfish source distribution
//
// modified from https://github.com/gmarcais/Jellyfish/blob/master/examples/count_in_file/count_in_file.cc
//
// Given a number of jellyfish databases, it retrieves all the k-mers
// (in the union of the database) the number of occurrences in each of
// the database and uses this for the d2s calculation

// LIMITATION: all the databases must have been created with the same size
// parameter to Jellyfish. It is advised to use the --disk switch of
// 'jellyfish count' to ensure that all the database are created with the
// same size. For example: jellyfish count -s 100M --disk ...

// Some details
// ============

// The Jellyfish database are sorted list of k-mers. The sorting order is
// based on the hash function (a random binary matrix). The program works
// similarly to how merging sorted lists is done, but instead of summing
// the count of each k-mers, it displays the count in each of the sorted
// list.

// Because all the files must be sorted according to the same criteria,
// the size of the hash (and therefore the hash function) must be the
// same for all the input files. See the LIMITATION paragraph above to
// run Jellyfish properly. Note that the size passed to 'jellyfish count'
// is NOT a limit on the number k-mers in the file (only a limit on the
// number of mers in memory).


#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <map>

#include "jellyfish/err.hpp"
#include "jellyfish/misc.hpp"
#include "jellyfish/mer_heap.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/rectangular_binary_matrix.hpp"
#include "jellyfish/cpp_array.hpp"

namespace err = jellyfish::err;

using jellyfish::file_header;
using jellyfish::RectangularBinaryMatrix;
using jellyfish::mer_dna;
using jellyfish::cpp_array;
typedef std::unique_ptr<binary_reader>           binary_reader_ptr;
typedef std::unique_ptr<text_reader>             text_reader_ptr;

struct file_info {
  std::ifstream is;
  file_header   header;

  file_info(const char* path) :
    is(path),
    header(is)
  { }
};


struct common_info {
  unsigned int            key_len;
  size_t                  max_reprobe_offset;
  size_t                  size;
  unsigned int            out_counter_len;
  std::string             format;
  RectangularBinaryMatrix matrix;

  common_info(RectangularBinaryMatrix&& m) : matrix(std::move(m))
  { }
};


common_info read_headers(int argc, char* input_files[], cpp_array<file_info>& files) {
  // Read first file
  files.init(0, input_files[0]);
  if(!files[0].is.good())
    err::die(err::msg() << "Failed to open input file '" << input_files[0] << "'");

  file_header& h = files[0].header;
  common_info res(h.matrix());
  res.key_len            = h.key_len();
  res.max_reprobe_offset = h.max_reprobe_offset();
  res.size               = h.size();
  res.format = h.format();
  size_t reprobes[h.max_reprobe() + 1];
  h.get_reprobes(reprobes);
  res.out_counter_len = h.counter_len();

  // Other files must match
  for(int i = 1; i < argc; i++) {
    files.init(i, input_files[i]);
    file_header& nh = files[i].header;
    if(!files[i].is.good())
      err::die(err::msg() << "Failed to open input file '" << input_files[i] << "'");
    if(res.format != nh.format())
      err::die(err::msg() << "Can't compare files with different formats (" << res.format << ", " << nh.format() << ")");
    if(res.key_len != nh.key_len())
      err::die(err::msg() << "Can't compare hashes of different key lengths (" << res.key_len << ", " << nh.key_len() << ")");
    if(res.max_reprobe_offset != nh.max_reprobe_offset())
      err::die("Can't compare hashes with different reprobing strategies");
    if(res.size != nh.size())
      err::die(err::msg() << "Can't compare hash with different size (" << res.size << ", " << nh.size() << ")");
    if(res.matrix != nh.matrix())
      err::die("Can't compare hash with different hash function");
  }

  return res;
}




static PyObject *JellyfishError;
// static PyObject * jellyfish_system(PyObject *self, PyObject *args);

static PyObject * jellyfish_d2s(PyObject *self, PyObject *args);

static PyMethodDef JellyfishMethods[] = {
    {"d2s",  jellyfish_d2s, METH_VARARGS,
     "Calculate self d2s score"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


static struct PyModuleDef jellyfishmodule = {
    PyModuleDef_HEAD_INIT,
    "jellyfish",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    JellyfishMethods
};




PyMODINIT_FUNC
PyInit_jellyfish(void)
{
    PyObject *m;

    m = PyModule_Create(&jellyfishmodule);
    if (m == NULL)
        return NULL;

    JellyfishError = PyErr_NewException("jellyfish.error", NULL, NULL);
    Py_XINCREF(JellyfishError);
    if (PyModule_AddObject(m, "error", JellyfishError) < 0) {
        Py_XDECREF(JellyfishError);
        Py_CLEAR(JellyfishError);
        Py_DECREF(m);
        return NULL;
    }

    return PyModule_Create(&jellyfishmodule);
}

// ----------
typedef std::map<char, double> BaseFrequencyMap;


template<typename reader_type>
double d2s_raw(cpp_array<file_info>& files, long n_seq1, long n_seq2, long total_len1, long total_len2, BaseFrequencyMap cf1, BaseFrequencyMap cf2) {
  cpp_array<reader_type> readers(files.size());
  typedef jellyfish::mer_heap::heap<mer_dna, reader_type> heap_type;
  typedef typename heap_type::const_item_t heap_item;
  heap_type heap(files.size());

  // Prime heap
  for(size_t i = 0; i < files.size(); ++i) {
    readers.init(i, files[i].is, &files[i].header);
    if(readers[i].next())
      heap.push(readers[i]);
  }

  heap_item          head      = heap.head();
  mer_dna            key;
  const int          num_files = files.size();
  const reader_type* base      = &readers[0];
  uint64_t           counts[num_files];

  std::string keystring;

  double p_kmer1,norm_kmc1;
  double p_kmer2,norm_kmc2;

  double d2s = 0;
  double d2s_self1 = 0;
  double d2s_self2 = 0;

  while(heap.is_not_empty()) {
    key = head->key_;
    memset(counts, '\0', sizeof(uint64_t) * num_files);
    do {
      counts[head->it_ - base] = head->val_;
      heap.pop();
      if(head->it_->next())
        heap.push(*head->it_);
      head = heap.head();
    } while(head->key_ == key && heap.is_not_empty());

    p_kmer1=1;
    p_kmer2=1;
    keystring = key.to_str();
    for(int ki = 0;ki < key.k();ki++){
    	p_kmer1*=cf1[keystring[ki]];
    	p_kmer2*=cf1[keystring[ki]];
    }
    p_kmer1*=(total_len1 - (n_seq1*(key.k()-1)));
    p_kmer2*=(total_len2 - (n_seq2*(key.k()-1)));

    if ( counts[0]> 0){
	    norm_kmc1 = counts[0] - p_kmer1;
	    d2s_self1 += norm_kmc1 * norm_kmc1 / sqrt((2*norm_kmc1*norm_kmc1));

	}
	if ( counts[1]>0){
	    norm_kmc2 = counts[1] - p_kmer2;
	    d2s_self2 += norm_kmc2 * norm_kmc2 / sqrt((2*norm_kmc2*norm_kmc2));	    
	}
	if ( counts[0]>0 && counts[1] > 0){
	    d2s += norm_kmc1 * norm_kmc2 / sqrt((norm_kmc1*norm_kmc1  + norm_kmc2*norm_kmc2));
	 }

  }

  d2s = abs(log(d2s/sqrt(d2s_self1*d2s_self2)));

  return d2s;
}


static PyObject *
jellyfish_d2s(PyObject *self, PyObject *args)
{
    char *jf1;
    char *jf2;
    long n_seq1;long n_seq2;
    long total_len1;long total_len2;
    PyObject* cf1_po;
    PyObject* cf2_po;
    double d2s;
    if (!PyArg_ParseTuple(args, "ssllllO!O!", &jf1,&jf2,&n_seq1,&n_seq2,&total_len1,&total_len2,  &PyDict_Type, &cf1_po,&PyDict_Type, &cf2_po))
        return NULL;

    char *file_paths[] = {jf1,jf2};

	BaseFrequencyMap cf1;
	BaseFrequencyMap cf2;

    cf1['C'] = PyFloat_AsDouble(PyDict_GetItemString(cf1_po, "C"));
    cf1['G'] = PyFloat_AsDouble(PyDict_GetItemString(cf1_po, "G"));
    cf1['T'] = PyFloat_AsDouble(PyDict_GetItemString(cf1_po, "T"));
    cf1['A'] = PyFloat_AsDouble(PyDict_GetItemString(cf1_po, "A"));

    cf2['C'] = PyFloat_AsDouble(PyDict_GetItemString(cf2_po, "C"));
    cf2['G'] = PyFloat_AsDouble(PyDict_GetItemString(cf2_po, "G"));
    cf2['T'] = PyFloat_AsDouble(PyDict_GetItemString(cf2_po, "T"));
    cf2['A'] = PyFloat_AsDouble(PyDict_GetItemString(cf2_po, "A"));

	// Read the header of each input files and do sanity checks.
  	cpp_array<file_info> files(2);
	common_info cinfo = read_headers(2, file_paths, files);
	mer_dna::k(cinfo.key_len / 2);

  if(cinfo.format == binary_dumper::format)
    d2s = d2s_raw<binary_reader>(files,n_seq1,n_seq2,total_len1,total_len2,cf1,cf2);
  else if(cinfo.format == text_dumper::format)
    d2s = d2s_raw<text_reader>(files,n_seq1,n_seq2,total_len1,total_len2,cf1,cf2);
  else
    err::die(err::msg() << "Format '" << cinfo.format << "' not supported\n");


  return PyFloat_FromDouble(d2s);
}


