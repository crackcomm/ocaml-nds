// OCaml bindings to multicore ndsort algorithm.
// There were no modifications in the algorithm.
//
// Source: https://github.com/juanjonrg/nds
// SPDX-License-Identifier: Apache-2.0

#define _GNU_SOURCE
#define CAML_INTERNALS

#include <caml/memory.h>
#include <caml/bigarray.h>

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>

int n, m;
const float *population;
int *sorted_pop;
int ***all_fronts;
int **all_front_counts;
int *ranks;

int comp(const int a_idx, const int b_idx, const int j)
{
    float a = population[a_idx * m + j];
    float b = population[b_idx * m + j];
    if (a < b)
        return -1;
    else if (a > b)
        return 1;
    return 0;
}

int compare(const void *a_in, const void *b_in, void *thunk_in)
{
    int j = *((int *)thunk_in);
    int a_idx = *(int *)a_in;
    int b_idx = *(int *)b_in;
    int result = comp(a_idx, b_idx, j);

    if (result == 0)
    {
        for (int i = 1; i < m; i++)
        {
            result = comp(a_idx, b_idx, (j + i) % m);
            if (result != 0)
            {
                return result;
            }
        }
    }
    return result;
}

int is_dominated(int a, int b)
{
    bool equal = true;
    for (int j = 0; j < m; j++)
    {
        if (population[a * m + j] < population[b * m + j])
        {
            return 0;
        }
        else if (equal && population[a * m + j] > population[b * m + j])
        {
            equal = false;
        }
    }
    return !equal;
}

int add_to_front(int idx, int f, int last_front, int **fronts, int *front_counts)
{
    int count = front_counts[f];
    if (front_counts[f] == 0)
    {
        // Initialize front
        fronts[f] = malloc(n * sizeof(int));
    }
    fronts[f][count] = idx;
    ranks[idx] = f;
    front_counts[f]++;

    if (f > last_front)
    {
        return f;
    }
    return last_front;
}

static void *threaded_nds(void *thread_id)
{
    int j = (long)thread_id;

    // Initalize population indexes.
    for (int i = 0; i < n; i++)
    {
        sorted_pop[j * n + i] = i;
    }

    // Sort population by objective j.
    qsort_r(&sorted_pop[j * n], n, sizeof(int), compare, &j);

    // Initialize front structures.
    all_fronts[j] = malloc(n * sizeof(int *));
    all_front_counts[j] = calloc(n, sizeof(int));
    int **fronts = all_fronts[j];
    int *front_counts = all_front_counts[j];
    int last_front = 0;

    // Lets rank!
    for (int i = 0; i < n; i++)
    {
        int idx = sorted_pop[j * n + i];
        int rank = ranks[idx];

        if (rank < 0)
        {
            // Individual not ranked!
            bool check = true;
            for (int x = 0; x <= last_front; x++)
            {
                check = false;
                for (int y = 0; y < front_counts[x]; y++)
                {
                    int member = fronts[x][y];
                    check = is_dominated(idx, member);
                    if (check)
                    {
                        // Not on this front!
                        break;
                    }
                }
                if (!check)
                {
                    // On this front!
                    last_front = add_to_front(idx, x, last_front, fronts, front_counts);
                    break;
                }
            }
            if (check)
            {
                // Dominated by the last front, on a new front!
                last_front = add_to_front(idx, last_front + 1, last_front, fronts, front_counts);
            }
        }
        else
        {
            // Individual ranked by other thread, add it to our fronts!
            last_front = add_to_front(idx, rank, last_front, fronts, front_counts);
        }
    }

    for (int i = 0; i <= last_front; i++)
    {
        free(fronts[i]);
    }
    free(fronts);

    return NULL;
}

#define checkPth(rc) __checkPth(rc, __FILE__, __LINE__)
static inline void __checkPth(int rc, const char *file, const int line)
{
    if (rc != 0)
    {
        fprintf(stderr, "PTH error at %s:%i: %d\n", file, line, rc);
        exit(1);
    }
}

void ndsort()
{
    sorted_pop = malloc(n * m * sizeof(int));
    all_fronts = malloc(m * sizeof(int **));
    all_front_counts = malloc(m * sizeof(int *));

    for (int i = 0; i < n; i++)
    {
        ranks[i] = -1;
    }

    int num_threads = m;
    pthread_t threads[num_threads];

    for (int i = 0; i < num_threads; i++)
    {
        checkPth(pthread_create(&threads[i], NULL, threaded_nds, (void *)(intptr_t)i));
    }
    for (int i = 0; i < num_threads; i++)
    {
        checkPth(pthread_join(threads[i], NULL));
    }

    for (int i = 0; i < num_threads; i++)
    {
        free(all_front_counts[i]);
    }
    free(sorted_pop);
    free(all_fronts);
    free(all_front_counts);
}

value ml_ndsort(value vb, value vr)
{
    CAMLparam2(vb, vr);
    struct caml_ba_array *b = Caml_ba_array_val(vb);
    struct caml_ba_array *r = Caml_ba_array_val(vr);
    assert(b->num_dims == 2);
    assert(r->num_dims == 1);
    assert(sizeof(float) == caml_ba_element_size[b->flags & CAML_BA_KIND_MASK]);
    assert(sizeof(int) == caml_ba_element_size[r->flags & CAML_BA_KIND_MASK]);
    assert(r->dim[0] == b->dim[0]);
    n = b->dim[0];
    m = b->dim[1];

    population = (float *)b->data;
    ranks = (int *)r->data;

    ndsort();

    population = NULL;
    ranks = NULL;

    CAMLreturn(Val_unit);
}
