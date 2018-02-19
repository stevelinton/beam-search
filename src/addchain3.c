#include "beam.h"
#include <stdio.h>

#define maxP 1024
#define maxLen 128
typedef uint16_t elt;

elt code[maxLen];
bool targets[maxP];
int codelen;

typedef struct s_chain {
    int len;
    uint32_t fitness;
    elt chain[maxLen];
} *chain;

static int P;
#define data_size sizeof(struct s_chain)

static uint32_t fitness(const char *cv) {
    chain c = (chain)cv;
    return c->fitness;
}

static void visit_children(const char *parent, void (*visit)(const char *, void *), void *context) {
    chain c = (chain)parent;
    int l = c->len;
    char ch[data_size];
    chain child = (chain)ch;
    for (int i =1; i < l; i++)
        for (int j = 1; j <= i; j++) {
            int k = (c->chain[i] + c->chain[j]) % P;
            bool isnew = true;
            for (int x = 0; x < l; x++)
                if (c->chain[x] == k) {
                    isnew = false;
                    break;
                }
            if (isnew) {
                memcpy(child, parent, data_size);
                child->chain[child->len++] = k;
                if (targets[k])
                    child->fitness++;
                if (child->fitness == codelen)
                    child->fitness = stop_fitness;
                visit(ch, context);
            }
        }
}

static bool equal(const char *a1, const char *a2) {
    return (0 == strncmp(a1,a2,data_size));
}

#define fnvp 1099511628211ULL
#define fnvob 14695981039346656037ULL

static uint64_t hash( const char *c) {
    uint64_t h = fnvob;
    const chain cc = (chain)c;
    for (int i = 0; i < sizeof(elt)*(cc->len); i++) {
        h = (h*fnvp) ^ ((char *)(cc->chain))[i];
    }
    return h;    
}

static void print_chain(const char *i) {
    const chain c = (chain) i;
    printf("<chain");
    for (int j = 0; j < c->len; j++)
        printf(" %i",c->chain[j]);
    printf(">");
}


int main(int argc, char **argv) {
    P = atoi(argv[1]);
    int len = atoi(argv[2]);
    if (P > maxP || len > maxLen)
        exit(EXIT_FAILURE);
    int beamsize = 10000;
    if (argc >= 4)
        beamsize = atoi(argv[3]);
    int nprobes = 3;
    if (argc >= 5)
        nprobes = atoi(argv[4]);
    scanf("%i",&codelen);
    for (int i = 0; i < P; i++)
        targets[i] = false;
    for (int i = 0; i < codelen; i++) {
        int x;
        scanf("%i",&x);
        code[i] = x;
    }
    int a,b;
    char * seed = calloc(1,data_size);
    ((chain)seed)->len = 2;
    ((chain)seed)->chain[0] = 0;
    ((chain)seed)->chain[1] = 1;
    ((chain)seed)->fitness = 2;
    for (int z = 0; z < codelen; z++)
        {
            for( int o = 0; o < codelen; o++) {
                if (o == z) continue;
                a = (P - code[z]) % P ;
                // printf("\nz: %i o:%i code[z]:%i code[o]:%i a: %i",z,o,code[z],code[o],a);
                int c = (P + code[o] - code[z]) % P;
                for (b = 1; b < P; b++)
                    if ((b*c) % P == 1)
                        break;
                memset(targets, 0, sizeof(targets));
                for (int i = 0; i < codelen; i++) {
                    int x = (b*(a + code[i])) % P;
                    targets[x] = true;
                }
                int nresults;
                char * results = beam_search(seed,1,visit_children, beamsize, len-2,
                                             data_size,  fitness, equal, hash, nprobes,
                                             print_chain, &nresults);
                int maxfitness = 0;
                const char * bestchain = NULL;
                int *fitcounts = calloc(sizeof(int),P+1);
                int count = 0;
                for (int i = 0; i < nresults; i++) {
                    const char *c = results + i*data_size;
                    int f = fitness(c);
                    if (f == stop_fitness)
                        f = P;
                    if (f > maxfitness) {
                        maxfitness = f;
                        bestchain = c;
                    }
                    fitcounts[f]++;
                    count++;
                }
                if (maxfitness >= codelen) {
                    printf("Code: ");
                    for (int i = 0; i < P; i++)
                        if (targets[i])
                            printf("%i ",i);
                    printf("\nSolution: ");
                    print_chain(bestchain);
                    printf("\n");
                }
                free(results);
        }
    }
    exit(EXIT_SUCCESS);
    
}
