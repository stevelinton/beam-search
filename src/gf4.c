#include "beam.h"
#include <stdio.h>
#include <stdlib.h>

#define maxLines 21
typedef uint8_t halfline;
typedef uint64_t vector;
typedef struct { halfline a;
    halfline b;} line;

#define maxDim 34

typedef struct {
    uint8_t dim;
    uint8_t pivs[maxDim];
    vector ech[maxDim];
} space;

uint8_t transtab[256];

void fill_transtab(void) {
    for (int i = 0; i < 256; i++)
        transtab[i] = (i & 0xC0) | ((i & 0x0c) << 2) | ((i & 0x30)>>2) | (i & 0x03);
}

/*

THis is the set of 8 vectors we have to hit [i,j] is a_i tensor b_j

[ [ [ 0, 0 ], [ 1, 1 ], [ 2, 4 ], [ 3, 5 ] ], 
  [ [ 0, 1 ], [ 1, 1 ], [ 1, 0 ], [ 2, 5 ], [ 3, 5 ], [ 3, 4 ] ], 
  [ [ 0, 2 ], [ 1, 3 ], [ 2, 6 ], [ 3, 7 ] ], 
  [ [ 0, 3 ], [ 1, 3 ], [ 1, 2 ], [ 2, 7 ], [ 3, 7 ], [ 3, 6 ] ], 
  [ [ 4, 0 ], [ 5, 1 ], [ 6, 4 ], [ 7, 5 ] ], 
  [ [ 4, 1 ], [ 5, 1 ], [ 5, 0 ], [ 6, 5 ], [ 7, 5 ], [ 7, 4 ] ], 
  [ [ 4, 2 ], [ 5, 3 ], [ 6, 6 ], [ 7, 7 ] ], 
  [ [ 4, 3 ], [ 5, 3 ], [ 5, 2 ], [ 6, 7 ], [ 7, 7 ], [ 7, 6 ] ] ]

[ "20100201", "30200302", "80400804", "C0800C08", "2010020100000000", "3020030200000000",
  "8040080400000000", "C0800C0800000000" ]

Per Richards suggestion we look also at the even projectsions of this
[[0,0],[2,4]]
[[0,2],[2,6]]
[[4,0],[6,4]]
[[4,2],[6,6]]

100001
400004
10000100000000
40000400000000


As vectors in hex, generators of the 8 space target
*/
typedef struct s_lines {
    int len;
    uint32_t fitness;
    line lines[maxLines];
    space lspace;
    space sumspace;
    space sum12space;
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
    for (int i = 0; i < 8; i++)
        if (l.a & (1L << i))
            t |= (((uint64_t)l.b) << (8*i));
    return t;
}

vector clean( space *s, vector v) {
    for (int i = 0; i < s->dim; i++) {
        if (v & (1L << s->pivs[i]))
            v ^= s->ech[i];
    }
    return v;
}

void extend( space *s, vector v) {
    for (int i = 0; i < 64; i++)
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
    v = clean(&(sol->sum12space),v);
    if (!v)
        return 2;
    extend(&(sol->sum12space),v);
    return 3;
}

static void visit_children(const char *parent, void visit(const char *, void *), void *context) {
    int ct = 0;
    soln c = (soln)parent;
    //print_soln(parent);
    char ch[data_size];
    soln child = (soln)ch;
    int start = 1;
    if (c->len) {
        line ll = c->lines[c->len-1];
        start += ll.a * 256 + ll.b+1;
    }
    for (int i = start; i < (1<<16); i++)
        {
            line l;
            l.a = i &0xFF ;
            //            l.b = transtab[i];
            l.b = i >>8;
            memcpy(ch, parent, data_size);
            int status = newline(child, l);
            if (status) {
                if (status == 1) {
                    child->fitness += 256;
                }
                if (status == 2) {
                    child->fitness++;
                }
                visit(ch, context);
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
    fill_transtab();
    memset(seed, 0, data_size);
    ((soln)seed)->len = 0;
    ((soln)seed)->fitness = 1;
    ((soln)seed)->sumspace.ech[0] = 0x20100201UL;
    ((soln)seed)->sumspace.pivs[0] = 0;
    ((soln)seed)->sumspace.ech[1] = 0x30200302UL;
    ((soln)seed)->sumspace.pivs[1] = 1;
    ((soln)seed)->sumspace.ech[2] = 0x80400804UL;
    ((soln)seed)->sumspace.pivs[2] = 2;
    ((soln)seed)->sumspace.ech[3] = 0xC0800C08UL;
    ((soln)seed)->sumspace.pivs[3] = 3;
    ((soln)seed)->sumspace.ech[4] = 0x2010020100000000UL;
    ((soln)seed)->sumspace.pivs[4] = 32;
    ((soln)seed)->sumspace.ech[5] = 0x3020030200000000UL;
    ((soln)seed)->sumspace.pivs[5] = 33;
    ((soln)seed)->sumspace.ech[6] = 0x8040080400000000UL;
    ((soln)seed)->sumspace.pivs[6] = 34;
    ((soln)seed)->sumspace.ech[7] = 0xC0800C0800000000UL;
    ((soln)seed)->sumspace.pivs[7] = 35;
    ((soln)seed)->sumspace.dim = 8;
    memcpy(&((soln)seed)->sum12space, &((soln)seed)->sumspace, sizeof(((soln)seed)->sumspace));
    ((soln)seed)->sum12space.ech[8] = 0x20000200UL;
    ((soln)seed)->sum12space.pivs[8] = 10;
    ((soln)seed)->sum12space.ech[8] = 0x80000800UL;
    ((soln)seed)->sum12space.pivs[8] = 12;
    ((soln)seed)->sum12space.ech[8] = 0x2000020000000000UL;
    ((soln)seed)->sum12space.pivs[8] = 42;
    ((soln)seed)->sum12space.ech[8] = 0x8000080000000000UL;
    ((soln)seed)->sum12space.pivs[8] = 44;                                            
    ((soln)seed)->sum12space.dim = 12;
    
    int nresults;
    char * results = beam_search(seed,1,visit_children, beamsize, 21,
                                 data_size,  fitness, equal, hash, 3, print_soln, &nresults);
    int maxfitness = 0;
    const char * bestsoln = NULL;
    int count = 0;
    for (int i = 0; i < nresults; i++) {
        const char *c = results + i*data_size;
        uint32_t f = ((soln)c)->fitness;
        if (f == stop_fitness)
            f = 8;
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
