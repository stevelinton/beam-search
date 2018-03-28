#include "beam.h"
#include <stdio.h>

#define NB 10 
#define maxLen 64
typedef uint16_t elt;

typedef struct s_code {
    int len;
    uint32_t fitness;
    elt code[maxLen];
    char mask[(1<<NB)]; 
} *code;

#define data_size sizeof(struct s_code)

static uint32_t fitness(const char *cv) {
    code c = (code)cv;
    return c->fitness;
}

static void visit_children(const char *parent, void visit(const char *, void *), void *context) {
    int ct = 0;
    code c = (code)parent;
    int l = c->len;
    char ch[data_size];
    code child = (code)ch;
    elt x;
    for (x = 0; x < (1 << NB); x++) {
        memcpy(ch, parent, data_size);
        child->code[l] = x;
        child->len++;
        for (int i = 0; i < l; i++) {
            elt y = x^child->code[i];
            if (!child->mask[y]) {
                child->fitness++;
                child->mask[y] = 1;
            }            
        }
        visit(ch, context);
    }
}

static bool equal(const char *a1, const char *a2) {
    return (0 == strncmp(a1,a2,data_size));
}

#define fnvp 1099511628211ULL
#define fnvob 14695981039346656037ULL

static uint64_t hash( const char *c) {
    uint64_t h = fnvob;
    const code cc = (code)c;
    for (int i = 0; i < sizeof(elt)*(cc->len); i++) {
        h = (h*fnvp) ^ ((char *)(cc->code))[i];
    }
    return h;    
}

static void print_code(const char *i) {
    const code c = (code) i;
    printf("<code");
    for (int j = 0; j < c->len; j++)
        printf(" %i",c->code[j]);
    printf(">");
}


int main(int argc, char **argv) {
    int len = atoi(argv[1]);
    if (len > maxLen)
        exit(EXIT_FAILURE);
    int beamsize = 10000;
    if (argc >= 3)
        beamsize = atoi(argv[2]);
    int nprobes = 3;
    char * seed = malloc(data_size);
    memset(seed, 0, data_size);
    ((code)seed)->len = NB+1;
    ((code)seed)->fitness = 1+(NB*(NB+1))/2;    
    ((code)seed)->code[0] = 0;
    ((code)seed)->mask[0] = 1;
    for(int i = 0; i < NB; i++) {
        ((code)seed)->code[i+1] = 1 << i;
        ((code)seed)->mask[1 << i] = 1;
        for (int j = 0; j < i; j++)
            ((code)seed)->mask[(1 << i) | (1 << j)] = 1;
    }
    int nresults;
    char * results = beam_search(seed,1,visit_children, beamsize, len-NB-1,
                                 data_size,  fitness, equal, hash, nprobes, print_code, &nresults);
    int maxfitness = 0;
    const char * bestcode = NULL;
    int *fitcounts = calloc(sizeof(int),(1<<NB)+1);
    int count = 0;
    for (int i = 0; i < nresults; i++) {
        const char *c = results + i*data_size;
        int f = fitness(c);
        if (f == stop_fitness)
            f = 1024;
        if (f > maxfitness) {
            maxfitness = f;
            bestcode = c;
        }
        fitcounts[f]++;
        count++;
    }
    printf("%i solutions found", count);
    if (count) {
        printf(", best has fitness %i ", maxfitness);
        print_code(bestcode);
        printf("\nfitness counts:\n");
        for (int i = 0; i <= 1<<NB; i++) {
            if (fitcounts[i]) {
                printf("%i %i\n",i,fitcounts[i]);
            }
        }
    }
    exit(EXIT_SUCCESS);
}
