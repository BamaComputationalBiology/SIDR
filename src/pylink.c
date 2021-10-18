

// https://docs.python.org/3/extending/embedding.html

#include "pylink.h"

int run_analysis(PARAM *param) {
//	printf("inside run_analysis pylink.c\n");
        PyObject *function = PyObject_GetAttrString(param->module, "run_analysis");

        if (function && PyCallable_Check(function)) {

            PyObject *args = PyTuple_New(1);
            PyObject *kmer_list = PyList_New(param->n_kmers);

            if (!args || !kmer_list) {
                printf("Cannot create analysis args\n");
                Py_Finalize();
                return -1;
            }

            for(int n = 0; n < param->n_kmers; n++) {
                PyList_SetItem(kmer_list, n, Py_BuildValue("i", param->kmer_list[n]));
            }

            PyTuple_SetItem(args, 0, kmer_list);

            PyObject *ret = PyObject_CallObject(function, args);
            Py_DECREF(args);

            if (ret != NULL){
                Py_Finalize();
                return 0;
            }

            else {
                PyErr_Print();
                printf("Analysis call failed\n");
                Py_Finalize();
                return -1;
            }
        }

        else {
            if (PyErr_Occurred()) PyErr_Print();
            printf("Cannot find analysis function\n");
            Py_Finalize();
            return -1;
        }
}



int export_seqdata(PARAM *param, SEQDATA *seqdata) {
//printf("inside export seqdata in pylink.c\n");
    PyObject *function = PyObject_GetAttrString(param->module, "import_seqdata");

    if (function && PyCallable_Check(function)) {
        PyObject *args = PyTuple_New(7);
        PyObject *dist_list = PyList_New(param->n_kmers);

        if (!args || !dist_list) {
            printf("Cannot create export args\n");
            return -1;
        }

        for(uint32_t n = 0; n < param->n_kmers; n++) {
            PyObject *dist = PyList_New(0);
            for(int i = 0; i < param->max_kmer_freq; i++){
                PyList_Append(dist, Py_BuildValue("i", get_kpoint(seqdata, n, i)));
            }
            PyList_SetItem(dist_list, n, dist);
        }
	printf("after for loop");
        PyTuple_SetItem(args, 0, Py_BuildValue("s", get_name(seqdata)));
        PyTuple_SetItem(args, 1, Py_BuildValue("I", get_length(seqdata)));
        PyTuple_SetItem(args, 2, Py_BuildValue("f", get_gc(seqdata)));
        PyTuple_SetItem(args, 3, Py_BuildValue("f", get_cov(seqdata)));
        PyTuple_SetItem(args, 4, Py_BuildValue("i", get_blast(seqdata)));
        PyTuple_SetItem(args, 5, Py_BuildValue("i", get_tax(seqdata)));
        PyTuple_SetItem(args, 6, dist_list);

        PyObject *ret = PyObject_CallObject(function, args);
        Py_DECREF(args);
        freeSEQDATA(seqdata, param->n_kmers);

        if (ret != NULL){
            Py_DECREF(ret);
            printf("exporting data \n");
            return 0;
        }

        else {
            PyErr_Print();
            printf("Export call failed\n");
            return -1;
        }
    }

    else {
        if (PyErr_Occurred()) PyErr_Print();
        printf("Cannot find export function\n");
        return -1;
    }
}
