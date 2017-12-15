#include "beam.h"
#include <stdio.h>

#define maxP 512
#define maxLen 64
typedef uint16_t elt;

typedef struct s_code {
    int len;
    uint32_t fitness;
    elt code[maxLen];
    char mask[maxP]; // 1 in code, 2 reachable in AS
} *code;

static int P;
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
    for (int k = 2; k < P; k++) {
        if (c->mask[k] != 1) {
            memcpy(ch, parent, data_size);
            child->code[child->len++] = k;
            if (child->mask[k] == (char)0)
                child->fitness++;
            child->mask[k] = 1;
            for (int a = 0; a < child->len-1; a++) {
                int x= child->code[a];
                int c = (P + k +k -x) %P;
                if (child->mask[c] == (char)0) {
                    child->fitness ++;
                    child->mask[c] = 2;
                }                
                for (int b = 0; b <= a; b++) {
                    int y = child->code[b];
                     c = (P + x + y -k) % P;
                    if (child->mask[c] == (char)0) {
                        child->fitness ++;
                        child->mask[c] = 2;
                    }
                    c = (P+x+k-y) %P;
                    if (child->mask[c] == (char)0) {
                        child->fitness ++;
                        child->mask[c] = 2;
                    }
                    c = (P+y+k-x) % P;
                    if (child->mask[c] == (char)0) {
                        child->fitness ++;
                        child->mask[c] = 2;
                    }
                    
                }
            }
            if (child->fitness == P)
                child->fitness=stop_fitness;
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
    char * seed = calloc(1,data_size);
    ((code)seed)->len = 2;
    ((code)seed)->fitness = 4;
    ((code)seed)->code[0] = 0;
    ((code)seed)->code[1] = 1;
    ((code)seed)->mask[0] = 1;
    ((code)seed)->mask[1] = 1;
    ((code)seed)->mask[P-1] = 2;
    ((code)seed)->mask[2] = 2;
    int nresults;
    char * results = beam_search(seed,1,visit_children, beamsize, len-2,
                                 data_size,  fitness, equal, hash, nprobes, print_code, &nresults);
    int maxfitness = 0;
    const char * bestcode = NULL;
    int *fitcounts = calloc(sizeof(int),P+1);
    int count = 0;
    for (int i = 0; i < nresults; i++) {
        const char *c = results + i*data_size;
        int f = fitness(c);
        if (f == stop_fitness)
            f = P;
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
        for (int i = 0; i <= P; i++) {
            if (fitcounts[i]) {
                printf("%i %i\n",i,fitcounts[i]);
            }
        }
    }
    exit(EXIT_SUCCESS);
}
