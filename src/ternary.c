#include "beam.h"
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#define MAXMOVE 16
#define ANDN 1

typedef uint16_t regs_t;
typedef uint16_t res_t;
typedef uint8_t nstates_t;
#define NREGS 8*sizeof(regs_t)
#define MAXSTATES (1<<(8*sizeof(nstates_t)))
#define FAIL (MAXSTATES -1)


typedef struct {
    uint8_t arity;
    uint8_t r1;
    uint8_t r2;
    uint8_t r3;
    uint8_t op;
    uint8_t drop;
} move;

typedef struct {
    regs_t regs;
    res_t res;
} state;

typedef struct s_node {
    nstates_t s; // number of states
    uint8_t r; // number of registers
#ifdef TRACKMOVES
    uint8_t nmoves; // number of moves to get here
    move moves[MAXMOVE]; // record of those moves
#endif
    state states[];
} node;

typedef struct {
    uint8_t len;
    nstates_t size;
    state codewords[];
} coding;

static size_t data_size;

static fitness_t valreg, valstate, valh;



uint8_t ternary(uint8_t op, uint8_t in1, uint8_t in2, uint8_t in3) {
    return (op >>  (in1 | in2 | in3)) & 1;
}

uint8_t binary(uint8_t op, uint8_t in1, uint8_t in2) {
    return (op >> (in1 | in2)) &1;
}

const static regs_t mask1[16] = {
                                   32767, 16383, 8191, 4095, 
                                   2047, 1023, 511, 255, 
                                   127, 63, 31, 15, 
                                   7, 3, 1, 0, 
};
const static regs_t mask2[17] = {
                                   0, 32768, 49152, 57344, 
                                   61440, 63488, 64512, 65024, 
                                   65280, 65408, 65472, 65504, 
                                   65520, 65528, 65532, 65534, 65535
};

void drop(state *st, uint8_t pos) {
    regs_t l = st->regs & mask1[pos];
    regs_t r = st->regs & mask2[pos];
    st->regs = r | (l << 1);
}

uint8_t reg_extract_topos(regs_t r, uint8_t pos, uint8_t to) {
    return ((r << pos)>> (NREGS-1-to)) & (1 << to);
}

uint8_t reg_extract(regs_t r, uint8_t pos) {
    return reg_extract_topos(r,pos,0);
}

void reg_set(regs_t *r, uint8_t pos, uint8_t val) {
    *r |= ((val << (NREGS-1)) >> pos);
}


uint64_t hamming_states(node *n) {
    uint64_t score = 0;
    for (int i = 0; i < n->s-1; i++) {
        res_t r = n->states[i].res;
        regs_t rg = n->states[i].regs;
        for (int j = i+1; j < n->s; j++) {
            if (n->states[j].res == r) {
                score += __builtin_popcount(rg ^ n->states[j].regs);
            }
        }
    }
    return score;
}


static uint32_t fitness(const char *cv) {
    node *n = (node *)cv;
    fitness_t f =  1000000 - valstate*n->s - valreg*n->r;
    if (valh)
        f -= valh*hamming_states(n);
    return f;
}

void apply1(state *st, const move *m, uint8_t nextreg) {
    uint8_t in1, in2, in3;
    switch(m->arity) {
    case 1:
        in1 = reg_extract(st->regs, m->r1);
        reg_set(&(st->regs), nextreg, !in1);
        if (m->drop & 1)
            drop(st, m->r1);
        break;

    case 2:
        in1 = reg_extract_topos(st->regs, m->r1,1);
        in2 = reg_extract(st->regs, m->r2);
        reg_set(&(st->regs), nextreg, binary(m->op, in1, in2));
        if (m->drop & 2)
            drop(st, m->r2);
        if (m->drop & 1)
            drop(st, m->r1);
        break;

    case 3:
        in1 = reg_extract_topos(st->regs, m->r1,2);
        in2 = reg_extract_topos(st->regs, m->r2,1);
        in3 = reg_extract(st->regs, m->r3);
        reg_set(&(st->regs), nextreg, ternary(m->op, in1, in2, in3));
        if (m->drop &4)
            drop(st, m->r3);        
        if (m->drop & 2)
            drop(st, m->r2);
        if (m->drop & 1)
            drop(st, m->r1);
        break;
    }
    return;
}

const static uint8_t del2[4] = {1,0,0,-1};
const static uint8_t del3[8] = {1,0,0,-1,0,-1,-1,-2};



// returns number of states remaining, or FAIL
nstates_t sort_and_merge_states (state *states, nstates_t nstates) {
    if (nstates == 0)
        return 0;    
    nstates_t newstates = 1;
    for (int i = 1; i < nstates; i++) {
        state s = states[i];
        regs_t x = s.regs;
        int lo = 0;
        int hi = newstates-1;
        int merged = 0;
        while (hi >= lo) {
            int mid = (lo + hi)/2;
            regs_t y = states[mid].regs;
            if (x == y) {
                if (s.res != states[mid].res)
                    return FAIL;
                merged = 1;
                break;
            }
            if (x > y)
                lo = mid+1;
            else
                hi = mid-1;
        }
        if (!merged) {
            for (int j = newstates-1; j > hi; j--) {
                states[j+1] = states[j];
            }
            newstates++;
            states[hi+1] = s;
        }
    }
    return newstates;
}


static bool apply(node *c, const move *m) {
    if (c->r >= NREGS) {
        printf("Register overflow\n");
        return false;
    }
    for (int i = 0; i < c->s; i++)
        apply1(&(c->states[i]), m, c->r);
    // adjust number of registers
    switch(m->arity) {
    case 1:
        c->r += 1-m->drop;
        break;
    case 2:
        c ->r += del2[m->drop];
        break;
    case 3:
        c->r += del3[m->drop];
        break;
    }
    if (m->arity > 1 && m->drop) {
        c->s = sort_and_merge_states(c->states, c->s);
        if (c->s == FAIL)
            return false;
    }
#ifdef TRACKMOVES
    // record the move
    c->moves[c->nmoves++] = *m;
#endif
    return true;
}


static inline node *ind(node *arr, int i) {
    return (node *)(((char *)arr)+i*data_size);
}

static inline uint8_t makefrom(node *o, int i, int j, int d) {
    node *n = ind(o,i);
    memcpy(n,ind(o,j),data_size);
    state *ss = n->states;
    nstates_t ns = n->s;
    for (int i = 0; i < ns; i++) {
        drop(ss+i, d);
    }
    n->r--;
    nstates_t x = sort_and_merge_states(ss, ns);
    if (x != FAIL) {
        n->s = x;
#ifdef TRACKMOVES
        n->moves[n->nmoves-1].drop = i;
#endif
        return 1 << i;
    }
    return 0;
}

#ifndef BINARY

static uint8_t apply8(const node *c, const move *m, node *o) {
    uint8_t ok = 0;
    assert(m->arity == 3 && m->drop == 0);
    memcpy(o,c,data_size);
    if (!apply(o,m)) {
        // can only be register overflow
        return 0;
    }
    ok = 1;    
    ok |= makefrom(o,1,0,m->r1);
    ok |= makefrom(o,2,0,m->r2);
    ok |= makefrom(o,4,0,m->r3);
    if ((ok & 3) == 3)
        ok |= makefrom(o,3,2,m->r1);
    if ((ok & 5) == 5)
        ok |= makefrom(o,5,4,m->r1);
    if ((ok & 6) == 6)
        ok |= makefrom(o,6,4,m->r2);
    if ((ok & 104) == 104)
        ok |= makefrom(o,7,6,m->r1);
    return ok;
}

#endif

static uint8_t apply4(const node *c, const move *m, node *o) {
    uint8_t ok = 0;
    assert(m->arity == 2 && m->drop == 0);
    memcpy(o,c,data_size);
    if (!apply(o,m)) {
        // can only be register overflow
        return 0;
    }
    ok = 1;    
    ok |= makefrom(o,1,0,m->r1);
    ok |= makefrom(o,2,0,m->r2);
    if ((ok & 3) == 3)
        ok |= makefrom(o,3,2,m->r1);
    return ok;
}    

// must take 000->0 so LSB is zero
//
// Reps under permutation of inputs
// Unary functions omitted 

uint8_t TernaryOps [] = {
#ifdef ARM
                         226, 216, 202, 172, 228, 184 //BIF, BIT, BSL
#else
    2, 4, 6, 8, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 38, 40, 42, 44, 46, 50, 52,
    54, 56, 58, 62, 64, 66, 70, 72, 74, 76, 78, 82, 84, 86, 88, 92, 94, 96, 98, 100,
    104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134,
    138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 162, 164, 166, 168, 172,
    174, 176, 178, 180, 182, 184, 186, 188, 190, 194, 196, 198, 200, 202, 206, 208,
    210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 242, 
    244, 246, 248, 254
#ifdef ALLTERN
    ,1, 7, 9, 11, 13, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 53, 
  55, 57, 59, 61, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 87, 89, 91, 93, 97, 99, 101, 
  103, 105, 107, 109, 111, 113, 115, 117, 121, 123, 125, 127, 129, 131, 133, 135, 137, 
  139, 141, 143, 145, 147, 149, 151, 155, 157, 159, 161, 163, 167, 169, 171, 173, 177, 
  179, 181, 183, 185, 189, 191, 193, 197, 199, 201, 203, 205, 209, 211, 213, 215, 217, 
  219, 223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 247, 249, 251, 253 
#endif
#endif
};

uint8_t BinaryOps []  = {
#if ANDN || defined(ARM)                         
    2,4,
#endif
    6,8,14
#ifdef ARM
    ,11,13 // ORNOT
#endif
        };

const int nterns =sizeof(TernaryOps);
const int nbins =sizeof(BinaryOps);

void print_coding_inner(int nstates, int len, state *states) {
    int i = 0;
    int lastres;
    int started = 0;
    while (i < nstates) {
        if (started) {
            if (states[i].res != lastres) {
                lastres = states[i].res;
                printf(") %i(", lastres);
            } else
                printf(", ");
        } else {
            lastres = states[i].res;
            started = 1;
            printf("%i(",lastres);
        }
        for (int j = len-1; j >= 0; j--)
            printf("%c", reg_extract(states[i].regs, j) ? '1':'0');
        i++;            
    }
    printf(")");
}

void print_coding(coding *c) {
    print_coding_inner(c->size, c->len, c->codewords);    
}


void print_move(const move *m) {
    switch(m->arity) {
    case 1:
        printf("not(%i);",(int)m->r1);
        break;
    case 2:
        printf("b0x%x(%i,%i);",
               (int)m->op, (int)m->r1, (int)m->r2);
        break;
    case 3:
        printf("t0x%x(%i,%i,%i);",
               (int)m->op, (int)m->r1, (int)m->r2, (int)m->r3);
        break;
    }
    if (m->drop != 0) {
        int started = 0;
        printf(" drop(");
        if (m->drop & 1) {
            printf("%i",(int)m->r1);
            started = 1;
        }
        if (m->drop &2) {
            if (started)
                printf(",");
            printf("%i",(int)m->r2);
            started = 1;
        }
        if (m->drop &4) {
            if (started)
                printf(",");
            printf("%i",(int)m->r3);
        }
        printf(");");
    }
}

static void print_node(const char *np) {
    const node *n = (node *) np;
    state states[n->s];
    memcpy(states, n->states, sizeof(state)*n->s);
    for (int i = 1; i < n->s; i++) {
        state s = states[i];
       int j = i-1;
       while (j >= 0 && (states[j].res > s.res ||
                         (states[j].res == s.res && states[j].regs > s.regs))) {
            states[j+1] = states[j];
            j--;
        }
        states[j+1] = s;
    }
    print_coding_inner(n->s, n->r, states);
#ifdef TRACKMOVES
    if (n->nmoves) {
        printf(" via");
        for (int i = 0; i < n->nmoves; i++) {
            printf(" ");
            print_move(n->moves + i);
        }
    }
#endif
}




static void visit_children(const char *parent, void visit(const char *, void *), void *context) {
    node *n = (node *)parent;
#ifdef DEBUG
    printf("VC ");
    print_node((const char *)n);
    printf("\n"); 
#endif
    node *ch = malloc(data_size);
    node *children = malloc(data_size*8);
    move m;
    for (int i = 0; i < n->r; i++) {
        m.r1 = i;
#if defined(BINARY) || defined(ARM)        
        m.arity = 1;
        m.op = 1;
        for (int drop = 0; drop < 2; drop++) {
            m.drop = drop;
            memcpy(ch, n, data_size);
            if (apply(ch, &m)) {
#ifdef DEBUG
                printf("MOVE ");
                print_move (&m);                            
                printf(" CHILD ");
                print_node((const char *)ch);
                printf("\n");
#endif
                (*visit)((char *)ch, context);
            }
        }
#endif
        for (int j = i+1; j < n->r; j++) {
            m.arity = 2;
            m.drop = 0;
            m.r2 = j;
            for (int op = 0; op < nbins; op++) {
                m.op = BinaryOps[op];                
                uint8_t cases = apply4(n,&m,children);
                for (int drop =0; drop < 4; drop ++) {
                    if (cases & (1 << drop)) {
                        const char *child = ((const char *)children) + data_size*drop;
#ifdef DEBUG
                        m.drop = drop;
                        printf("MOVE ");
                        print_move (&m);                            
                        printf(" CHILD ");
                        print_node(child);
                        printf("\n");
                        m.drop = 0;
#endif
                        (*visit)(child, context);
                    }
                }
                
            }
#ifndef BINARY
            m.arity = 3;
            m.drop = 0;
            for (int k = j+1; k < n->r; k++) {
                m.r3 = k;
                for (int op = 0; op < nterns; op++) {
                    m.op = TernaryOps[op];
                    uint8_t cases = apply8(n,&m,children);
                    for (int drop =0; drop < 8; drop ++) {
                        if (cases & (1 << drop)) {
                            const char *child = ((const char *)children) + data_size*drop;
#ifdef DEBUG
                            m.drop = drop;
                            printf("MOVE ");
                            print_move (&m);                            
                            printf(" CHILD ");
                            print_node(child);
                            printf("\n");
                            m.drop = 0;
#endif
                            (*visit)(child, context);
                        }
                    }
                }
            }
#endif
        }
    }
    free(children);
    free(ch);
}

static bool equal(const char *a1, const char *a2) {
    node *n1 = (node *)a1;
    node *n2 = (node *)a2;
    return n1->r == n2->r &&
        n1->s == n2->s &&
        !memcmp((void *)n1->states, (void *)n2->states, sizeof(state)*n1->s);
}

#define fnvp 1099511628211ULL
#define fnvob 14695981039346656037ULL

// could shift to a faster hash
static uint64_t hash( const char *c) {
    uint64_t h = fnvob;
    node *n = (node *)c;
    h = (h*fnvp) ^ n->r;
    h = (h*fnvp) ^ n->s;    
    for (int i = 0; i < n->s; i++) {
        h = (h*fnvp) ^ n->states[i].regs;
        h = (h*fnvp) ^ n->states[i].res;
    }

    return h;    
}

int fgetc1(FILE *f) {
    int c;
    while (isspace(c = fgetc(f)))
        ;
    return c;
}

coding *read_coding(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f)
        return NULL;
    coding *c = calloc(sizeof(coding) + sizeof(state)*256,1);
    uint8_t len = 0;
    nstates_t size = 0;
    res_t res;
    while (1) {
        int x;
        if (1 != fscanf(f,"%i", &x))
            break;
        res = x;
        int ch = fgetc1(f);
        if (ch != '(') {
            return NULL;
        }
        regs_t wd;
        uint8_t thislen;
        while (1) {
            wd = 0;
            thislen = 0;
            while (1) {
                ch = fgetc1(f);
                if (ch == '0' || ch == '1') {
                    wd >>= 1;
                    if (ch == '1')
                        wd |= 0x8000;
                    thislen++;
                } else if (ch == ',' || ch == ')') {
                    c->codewords[size].regs = wd;
                    c->codewords[size].res = res;
                    size++;
                    if (len == 0) {
                        len = thislen;                        
                    } else if (len != thislen) {
                        return NULL;
                    }
                    break;
                } else
                    return NULL;
            }
            if (ch == ')') {
                break;
            }
        }
    }
    c->len = len;
    c->size = size;
    fclose(f);
    return c;
}

int read_params(const char *fn, int *P, int *steps, fitness_t *valreg, fitness_t *valstate, fitness_t *valh,
                size_t *beamsize, fitness_t *maxval) {
    FILE *f = fopen(fn, "r");
    int a,b,c,d,e,g,v;
    if (!f || 7 != fscanf(f, "%i%i%i%i%i%i%i", &a, &b, &c, &d, &v, &e, &g))
        return 0;
    *P = a;
    *steps = b;
    *valreg = c;
    *valstate = d;
    *valh = v;
    *beamsize = e;
    *maxval = g;
    fclose(f);
    return 1;
}

static node *make_seed(coding *b, coding *c, int P) {
    if (b->size*c->size > MAXSTATES) {
        printf("Too many states\n");
        exit(EXIT_FAILURE);
    }
    node* seed = calloc(data_size,1);
    seed->s = b->size*c->size;
    seed->r = b->len+c->len;
    nstates_t k = 0;
    for (int i = 0; i < b->size; i++) {
        regs_t r = b->codewords[i].regs >> c->len;
        res_t s = b->codewords[i].res;
        for (int j = 0; j < c->size; j++) {
            seed->states[k].regs  = r | c->codewords[j].regs;
            seed->states[k].res  = (s+c->codewords[j].res) % P;
            k++;
        }
    }
    seed->s = sort_and_merge_states(seed->states, seed->s);
    if (seed->s == FAIL) {
        printf("Contradiction found in seed state\n");
        exit(EXIT_FAILURE);
    }
    return seed;
}
    

int main(int argc, char **argv) {
    size_t beamsize;
    int nprobes = 4;
    size_t nresults;
    int steps;
    fitness_t maxval;
    int P;
    if (argc < 4) {
        printf("Usage: ternary <b-code> <c-code> <params>\n");
        exit(EXIT_FAILURE);
    }
    coding *b = read_coding(argv[1]);
    if(b) {
        printf("B coding:");
        print_coding(b);
        printf("\n");
    }
    coding *c = read_coding(argv[2]);
    if (c) {
        printf("C coding:");
        print_coding(c);
        printf("\n");
    }
    if (!b || !c || !read_params(argv[3], &P, &steps, &valreg, &valstate, &valh, &beamsize, &maxval)) {
        printf("Error reading files\n");
        exit(EXIT_FAILURE);
    }
#ifdef TRACKMOVES
    if (steps > MAXMOVE) {
        printf("No room to record that many steps -- recompile with bigger MAXMOVE\n");
        exit(EXIT_FAILURE);
    }
#endif
    data_size = sizeof(node) + sizeof(state)*b->size*c->size;
    printf("Parameters: %i %u %u %i %i %lu %u\n", P, steps, valreg, valstate, valh, beamsize,maxval);
    node *seed = make_seed(b,c,P);
    printf("Starting search at ");
    print_node((char *)seed);
    printf("\n");
    char * results = beam_search((char *)seed,1,visit_children, beamsize, steps,
                                 data_size,  fitness, equal, hash, nprobes, print_node, &nresults);
    for (int i = 0; i < nresults; i++) {
        const char *n = results + i*data_size;
        int f = fitness(n);
        if (f >= 1000000 - maxval) {
            print_node(n);
            printf("\n");            
        }
    }
    exit(EXIT_SUCCESS);
}
