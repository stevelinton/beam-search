#include "beam.h"
#include <omp.h>
#include <stdio.h>

#define cpu_relax() asm volatile("pause\n" : : : "memory")

typedef struct s_hashtab {
  uint32_t *fitness; 
  char *data;
  size_t data_size;
  size_t tabsize;
  uint32_t (*fitness_func)(const char *);
  bool (*equal)(const char *, const char *);
  uint64_t (*hash)(const char *);
  uint64_t nprobes;
  void (*print_item)(const char *);
} * hashtab;

static hashtab new_ht(size_t data_size, size_t tabsize,
                      uint32_t (*fitness_func)(const char *),
                      bool (*equal)(const char *, const char *),
                      uint64_t (*hash)(const char *), uint64_t nprobes,
                      void (*print_item)(const char *)) {
  hashtab h = (hashtab)malloc(sizeof(struct s_hashtab));
  if (tabsize < 17)
    tabsize = 17;
  h->fitness = (uint32_t *)calloc(4, tabsize);
  h->data = calloc(data_size, tabsize);
  h->data_size = data_size;
  h->tabsize = tabsize;
  h->fitness_func = fitness_func;
  h->equal = equal;
  h->hash = hash;
  h->nprobes = nprobes;
  h->print_item = print_item;
  return h;
}

static void free_ht(hashtab h) {
  free(h->fitness);
  free(h->data);
  free(h);
}

#define IN_USE 0xFFFFFFFF

uint32_t get_control(hashtab h, uint64_t k, uint32_t fit) {
  while (1) {
    uint32_t nfit = __sync_val_compare_and_swap(&(h->fitness[k]), fit, IN_USE);
    if (nfit == fit) {
      //            printf("Locked %li %i %i\n",k, omp_get_thread_num(), fit);
      return fit;
    }
    if (nfit != IN_USE)
      fit = nfit;
    cpu_relax();
  }
}

static bool earlystop;

static void ht_probe(hashtab h, const char *item) {
  uint64_t key = h->hash(item);
  uint64_t key1 = 13 - key % 13;
  uint32_t myfit = h->fitness_func(item);
  if (myfit == stop_fitness)
      earlystop = true;
  void *tmp_item[2] = {alloca(h->data_size), alloca(h->data_size)};
  int nexttmp = 0;
  bool havelock = false;
  //    printf("probing ");
  // h->print_item(item);
  for (int i = 0; i < h->nprobes; i++) {
    key %= h->tabsize;
    uint32_t fit;
    fit = h->fitness[key];
    while (fit == IN_USE) {
      cpu_relax();
      fit = h->fitness[key];
    }
    if (!fit) {
        fit = get_control(h, key, 0);
      havelock = true;
      if (!fit) {
        memcpy(h->data + h->data_size * key, item, h->data_size);
        __sync_synchronize();
        h->fitness[key] = myfit;
        //                printf("Unlocked %li %i %i\n",key,
        //                omp_get_thread_num(), myfit);
        // printf(" into empty slot %i\n",i);
        return;
      }
    }
    if (fit < myfit) {
      if (!havelock) {
          fit = get_control(h, key, fit);
        havelock = true;
      }
      if (fit < myfit) {
        __sync_synchronize();
        memcpy(tmp_item[nexttmp], h->data + h->data_size * key, h->data_size);
        memcpy(h->data + h->data_size * key, item, h->data_size);
        __sync_synchronize();
        h->fitness[key] = myfit;
        //                printf("Unlocked %li %i %i\n",key,
        //                omp_get_thread_num(), myfit);
        havelock = false;
        myfit = fit;
        item = tmp_item[nexttmp];
        nexttmp ^= 1;
        // printf(" swapped %i ",i);
        // h->print_item(item);
        key += key1;
        continue;
      }
    }
    if (fit == myfit) {
      if (!havelock) {
          fit = get_control(h, key, fit);
        havelock = true;
      }
      if (fit == myfit) {
        __sync_synchronize();
        if (h->equal(item, h->data + h->data_size * key)) {
          // printf(" dup %i\n",i);
          h->fitness[key] = fit;
          // printf("Unlocked %li %i %i\n",key, omp_get_thread_num(), fit);
          return;
        } else {
          h->fitness[key] = fit;
          // printf("Unlocked %li %i %i\n",key, omp_get_thread_num(), fit);
          havelock = false;
        }
      }
    }
    if (havelock) {
      h->fitness[key] = fit;
      havelock = false;
      // printf("Unlocked %li %i %i\n",key, omp_get_thread_num(), fit);
    }
    key += key1;
  }
  // printf(" ran out\n");
}

static void probe_multi(hashtab h, const char *items, int nitems) {
  for (int j = 0; j < nitems; j++) {
    ht_probe(h, items + j * h->data_size);
  }
}

static void visit(const char *item, void *context) {
  ht_probe((hashtab)context, item);
}

static hashtab nextgen(const hashtab h,
                       void visit_children(const char *,
                                           void (*visit)(const char *, void *),
                                           void *),
                       int beamsize) {
  hashtab newtab = new_ht(h->data_size, beamsize, h->fitness_func, h->equal,
                          h->hash, h->nprobes, h->print_item);
#pragma omp parallel for
  for (int i = 0; i < h->tabsize; i++) {
    if (h->fitness[i] != 0) {
      visit_children((char *)(h->data + h->data_size * i), visit, newtab);
    }
  }
  return newtab;
}

char *
beam_search(const char *seeds, int nseeds,
            void visit_children(const char *,
                                void (*visit)(const char *, void *), void *),
            int beamsize, int ngens, size_t data_size,
            uint32_t fitness_func(const char *),
            bool equal(const char *, const char *), uint64_t hash(const char *),
            int nprobes, void print_item(const char *), int *nresults) {
  hashtab current = new_ht(data_size, beamsize, fitness_func, equal, hash,
                           nprobes, print_item);
  probe_multi(current, seeds, nseeds);
  earlystop = false;
  for (int i = 0; i < ngens; i++) {
    printf("GENERATION %i\n", i);
    hashtab next = nextgen(current, visit_children, beamsize);
    free_ht(current);
    current = next;
    if (earlystop)
        break;
  }
  char *results = malloc(data_size * current->tabsize);
  int nres = 0;
  for (int i = 0; i < current->tabsize; i++) {
    if (current->fitness[i]) {
      memcpy(results + nres * data_size, current->data + i * data_size,
             data_size);
      nres++;
    }
  }
  *nresults = nres;
  free_ht(current);
  return results;
}
