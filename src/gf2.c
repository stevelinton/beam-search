#include "beam.h"
#include <stdio.h>
#include <stdlib.h>

#define maxLines 8
typedef uint8_t halfline;
typedef uint16_t vector;
typedef struct { halfline a;
    halfline b;} line;

#define maxDim 12

typedef struct {
    uint8_t dim;
    uint8_t pivs[maxDim];
    vector ech[maxDim];
} space;



/*

THis is the set of 4 vectors we have to hit [i,j] is a_i tensor b_j

[[0,0],[1,2]]
[[0,1],[1,3]]
[[2,0],[3,2]]
[[2,1],[3,3]]

0041
0082
4100
8200

As vectors in hex, generators of the 8 space target
*/
typedef struct s_lines {
    int len;
    uint32_t fitness;
    line lines[maxLines];
    space lspace;
    space sumspace;
} *soln;

#define data_size sizeof(struct s_lines)

static void print_soln(const char *i) {
    const soln c = (soln) i;
    printf("<lines");
    for (int j = 0; j < c->len; j++)
        printf(" %x.%x",(int)c->lines[j].a, (int)c->lines[j].b);
    printf(">");
}


static uint32_t fitness(const char *cv) {
     soln c = (soln)cv;
    return c->fitness;
}

vector tensor(line l) {
    vector t = 0;
    for (int i = 0; i < 4; i++)
        if (l.a & (1L << i))
            t |= (((vector)l.b) << (4*i));
    return t;
}

vector clean( space *s, vector v) {
    for (int i = 0; i < s->dim; i++) {
        if (v & (1 << s->pivs[i]))
            v ^= s->ech[i];
    }
    return v;
}

void extend( space *s, vector v) {
    for (int i = 0; i <16; i++)
        if (v & (1L <<i)) {
            s->ech[s->dim] = v;
            s->pivs[s->dim++] = i;
            /* for (int j = 0; j +1 < s->dim; j++) */
            /*     if (s->pivs[j] == i) { */
            /*         printf("%i %i %i\n",i, s->dim, j);                     */
            /*     } */
            return;
        }
}

// return 0 -- line already in lspace, 1 -- intersection dimn increased
// 2 -- i12 dimn inscreaed but not i8, 3 -- none of the above

int newline(soln sol, line l) {
    vector v = tensor(l);
    v = clean(&(sol->lspace), v);
    if (!v)
        return 0;
    extend(&(sol->lspace), v);
    sol->lines[sol->len++] = l;
    v = clean(&(sol->sumspace), v);
    if (!v)
        return 1;
    extend(&(sol->sumspace),v);
    return 2;
}

line fix[] = {{1,1},{2,4},{1,2},{2,8},{4,1},{8,4},{4,2},{8,8}};

static void visit_children(const char *parent, void visit(const char *, void *), void *context) {
    int ct = 0;
    soln c = (soln)parent;
    //print_soln(parent);
    char ch[data_size];
    soln child = (soln)ch;
    int start = 1;
    if (c->len) {
        line ll = c->lines[c->len-1];
        start += ll.a + 16 *ll.b+1;
    }
    for (int i = start; i < (1<<8); i++)
        {
            line l;
            l.a = i &0x0F ;
            //            l.b = transtab[i];
            l.b = i >>4;

            // l = fix[c->len];

            if (l.a && l.b) {
                memcpy(ch, parent, data_size);
                int status = newline(child, l);
                //                printf("%i",status);
                if (status) {
                    if (status == 1) {
                        child->fitness ++;
                        //                            print_soln(ch);
                    }
                    visit(ch, context);
                }
            }
        }
}

static bool equal(const char *a1, const char *a2) {
    soln s1 = (soln)a1;
    soln s2 = (soln)a2;
    return s1->len == s2->len &&
        0 == strncmp((char *)&(s1->lines),(char *)&(s2->lines),sizeof(line)*s1->len);
}

#define fnvp 1099511628211ULL
#define fnvob 14695981039346656037ULL

static uint64_t hash( const char *c) {
    uint64_t h = fnvob;
    const soln cc = (soln)c;
    for (int i = 0; i < sizeof(line)*(cc->len); i++) {
        h = (h*fnvp) ^ ((char *)(cc->lines))[i];
    }
    return h;    
}



int main(int argc, char **argv) {
    int beamsize;
    char * seed = malloc(data_size);
    beamsize = atoi(argv[1]);
    memset(seed, 0, data_size);
    ((soln)seed)->len = 0;
    ((soln)seed)->fitness = 1;
    ((soln)seed)->sumspace.ech[0] = 0x0041;
    ((soln)seed)->sumspace.pivs[0] = 0;
    ((soln)seed)->sumspace.ech[1] = 0x0082;
    ((soln)seed)->sumspace.pivs[1] = 1;
    ((soln)seed)->sumspace.ech[2] = 0x4100;
    ((soln)seed)->sumspace.pivs[2] = 8;
    ((soln)seed)->sumspace.ech[3] = 0x8200;
    ((soln)seed)->sumspace.pivs[3] = 9;
    ((soln)seed)->sumspace.dim = 4;
    
    int nresults;
    char * results = beam_search(seed,1,visit_children, beamsize, 6,
                                 data_size,  fitness, equal, hash, 3, print_soln, &nresults);
    int maxfitness = 0;
    const char * bestsoln = NULL;
    int count = 0;
    for (int i = 0; i < nresults; i++) {
        const char *c = results + i*data_size;
        uint32_t f = ((soln)c)->fitness;
        if (f == stop_fitness)
            f = 4;
        if (f > maxfitness) {
            maxfitness = f;
            bestsoln = c;
        }
        count++;
    }
    printf("%i solutions found", count);
    if (count) {
        printf(", best has fitness %i ", maxfitness);
        print_soln(bestsoln);
    }
    exit(EXIT_SUCCESS);
}
