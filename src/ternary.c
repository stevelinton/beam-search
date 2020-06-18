#include "beam.h"
#include <stdio.h>
#include <ctype.h>

#define MAXMOVE 16
// #define USE_SAM

typedef struct {
    bool isBin;
    uint8_t r1;
    uint8_t r2;
    uint8_t r3;
    uint8_t op;
    uint8_t drop;
} move;

typedef struct {
    uint16_t regs;
    uint8_t res;
} state;

typedef struct s_node {
    uint8_t s; // number of states
    uint8_t r; // number of registers
#ifdef TRACKMOVES
    uint8_t nmoves; // number of moves to get here
    move moves[MAXMOVE]; // record of those moves
#endif
    state states[];
} node;

typedef struct {
    int len;
    int size;
    state codewords[];
} coding;

static int data_size;


static node *make_seed(coding *b, coding *c, int P) {
    node* seed = calloc(data_size,1);
    ((node *)seed)->s = b->size*c->size;
    ((node *)seed)->r = b->len+c->len;
    int k = 0;
    for (int i = 0; i < b->size; i++) {
        uint16_t r = b->codewords[i].regs << c->len;
        uint8_t s = b->codewords[i].res;
        for (int j = 0; j < c->size; j++) {
            ((node *)seed)->states[k].regs  = r | c->codewords[j].regs;
            ((node *)seed)->states[k].res  = (s+c->codewords[j].res) % P;
            k++;
        }
    }
    return seed;
}

static int valreg, valstate;

static uint32_t fitness(const char *cv) {
    node *n = (node *)cv;
    return 1000000 - valstate*n->s - valreg*n->r;
}


int ternary(uint8_t op, int in1, int in2, int in3) {
    return (op >> (in1 *4 + in2 * 2 + in3)) & 1;
}

int binary(uint8_t op, int in1, int in2) {
    return (op >> (in1 *2 + in2)) &1;
}

const static uint16_t mask1[16] = {0, 1, 3, 7, 15, 31, 63, 127,
                      255, 511, 1023, 2047, 4095, 8191, 16383, 32767};
const static uint16_t mask2[16] = {65534, 65532, 65528, 65520, 65504,
                                   65472, 65408, 65280, 65024, 64512, 63488, 61440, 57344, 49152, 32768, 0 };                                   
void drop(state *st, uint8_t bit) {
    uint16_t l = st->regs & mask1[bit];
    uint16_t r = st->regs & mask2[bit];
    st->regs = l | (r >> 1);
}



void apply1(state *st, const move *m, uint8_t nextreg) {
    int in1 = (st->regs >> m->r1) & 1;
    int in2 = (st->regs >> m->r2) & 1;
    if (m->isBin) {
        st->regs |= (binary(m->op, in1, in2) << nextreg);
    } else {
        int in3 = (st->regs >> m->r3) & 1;
        st->regs |= (ternary(m->op, in1, in2, in3) << nextreg);
        if (m->drop &4)
            drop(st, m->r3);        
    }
    if (m->drop & 2)
        drop(st, m->r2);
    if (m->drop & 1)
        drop(st, m->r1);
    return;
}

const static int del2[4] = {1,0,0,-1};
const static int del3[8] = {1,0,0,-1,0,-1,-1,-2};

// returns number of states remaining, or -1 if contradiction found 
int sort_and_merge_states (state *states, int nstates) {
    if (nstates == 0)
        return 0;
    int newstates = 1;
    for (int i = 1; i < nstates; i++) {
        state s = states[i];
        uint16_t x = s.regs;
        int lo = 0;
        int hi = newstates-1;
        int merged = 0;
        while (hi >= lo) {
            int mid = (lo + hi)/2;
            uint16_t y = states[mid].regs;
            if (x == y) {
                if (s.res != states[mid].res)
                    return -1;
                merged = 1;
                break;
            }
            if (x > y)
                lo = mid+1;
            else if (x < y)
                hi = mid -1;
        }
        if (!merged) {
            memmove(states+hi+2,states+hi+1, sizeof(state)*(newstates -hi-1));
            newstates++;
            states[hi+1] = s;
        }
    }
    return newstates;
}

static void sortstates(state *states, int nstates) {
    // insertion sort the states
    for (int i = 1; i < nstates; i++) {
        state s = states[i];
        uint16_t x = s.regs;
        int j;
        for (j = i-1; j >= 0; j--) {
            uint16_t y = states[j].regs;
            if (y > x) {
                states[j+1] = states[j];
            } else break;
        }
        states[j+1] = s;
    }
}

static bool apply(node *c, const move *m) {
        for (int i = 0; i < c->s; i++)
            apply1(&(c->states[i]), m, c->r);
        // adjust number of registers
        if (m->isBin)
            c ->r += del2[m->drop];
        else
            c->r += del3[m->drop];
#ifdef USE_SAM
        c->s = sort_and_merge_states(c->states, c->s);
        if (c->s == (uint8_t)(-1))
            return false;
#else
        sortstates(c->states, c->s);
        int j = 0;
        for (int i = 1; i < c->s; i++) {
            if (c->states[i].regs == c->states[j].regs) {
                if (c->states[i].res != c->states[j].res) {
                    return false;
                }
            } else {
                c->states[++j] = c->states[i];
            }
        }
        c->s = j+1;
#endif
#ifdef TRACKMOVES
        // record the move
        c->moves[c->nmoves++] = *m;
#endif
        return true;
}

// must take 000->0 so LSB is zero
//
// Reps under permutation of inputs
// Unary functions omitted 

uint8_t TernaryOps [] = {
    2, 4, 6, 8, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 38, 40, 42, 44, 46, 50, 52, 54, 56, 58, 62, 64, 66, 70, 72, 74, 76, 78, 82, 
  84, 86, 88, 92, 94, 96, 98, 100, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 162, 164, 
  166, 168, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 194, 196, 198, 200, 202, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 242, 
  244, 246, 248, 254  };

uint8_t BinaryOps []  = {
    2,4,6,8,14
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
            printf("%c", (states[i].regs & (1 <<j)) ? '1':'0');
        i++;            
    }
    printf(")");
}

void print_coding(coding *c) {
    print_coding_inner(c->size, c->len, c->codewords);    
}

static void print_move(const move *m) {
    if (m->isBin) {
        printf("b0x%x(%i,%i);",
               (int)m->op, (int)m->r1, (int)m->r2);
    } else {
        printf("t0x%x(%i,%i,%i);",
               (int)m->op, (int)m->r1, (int)m->r2, (int)m->r3);        
    }
    if (m->drop != 0) {
        int started = 0;
        printf(" drop(");
        if (m->drop & 1) {
            printf("%i",m->r1);
            started = 1;
        }
        if (m->drop &2) {
            if (started)
                printf(",");
            printf("%i",m->r2);
            started = 1;
        }
        if (m->drop &4) {
            if (started)
                printf(",");
            printf("%i",m->r3);
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
        while (states[j].res > s.res ||
               (states[j].res == s.res && states[j].regs > s.regs)) {
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
    /*    printf("VC ");
    print_node((const char *)n);
    printf("\n"); */
    node *ch = malloc(data_size);
    move m;
    for (int i = 0; i < n->r-2; i++) {
        m.r1 = i;
        for (int j = i+1; j < n->r; j++) {
            m.r2 = j;
            m.isBin = true;
            for (int op = 0; op < nbins; op++) {
                m.op = BinaryOps[op];
                for (int drop = 0; drop < 4; drop ++) {
                    m.drop = drop;
                    memcpy(ch, n, data_size);
                    if (apply(ch, &m)) {
                        /*            printf("MOVE ");
                              print_move (&m);                            
                              printf(" CHILD ");
                              print_node((const char *)ch);
                              printf("\n");  */
                        (*visit)((char *)ch, context);
                        }
                }
            }
            m.isBin = false;
            for (int k = j+1; k < n->r; k++) {
                m.r3 = k;
                for (int op = 0; op < nterns; op++) {
                    m.op = TernaryOps[op];
                    for (int drop =0; drop < 8; drop ++) {
                        m.drop = drop;
                        memcpy(ch, n, data_size);
                        if (apply(ch, &m)) {
                            /*                                printf("MOVE ");
                                  print_move (&m);                            
                                  printf(" CHILD ");
                                  print_node((const char *)ch);
                                  printf("\n");   */
                            (*visit)((char *)ch, context);
                        }
                    }
                }
            }
        }
    }
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

static uint64_t hash( const char *c) {
    uint64_t h = fnvob;
    node *n = (node *)c;
    h = (h*fnvp) ^ n->r;
    h = (h*fnvp) ^ n->s;    
    for (int i = 0; i < n->s; i++) {
        h = (h*fnvp) ^ n->states[i].regs ^n->states[i].res;
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
    int len = 0;
    int size = 0;
    int res;
    while (1) {
        if (1 != fscanf(f,"%i", &res))
            break;
        int ch = fgetc1(f);
        if (ch != '(') {
            return NULL;
        }
        uint16_t wd;
        int thislen;
        while (1) {
            wd = 0;
            thislen = 0;
            while (1) {
                ch = fgetc1(f);
                if (ch == '0' || ch == '1') {
                    wd <<= 1;
                    if (ch == '1')
                        wd++;
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

int read_params(const char *fn, int *P, int *steps,  int *valreg, int *valstate, int *beamsize, int *maxval) {
    FILE *f = fopen(fn, "r");
    if (!f || 6 != fscanf(f, "%i%i%i%i%i%i", P, steps, valreg, valstate, beamsize, maxval))
        return 0;
    fclose(f);
    return 1;
}

int main(int argc, char **argv) {
    int beamsize;
    int nprobes = 4;
    int nresults;
    int steps;
    int maxval;
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
    if (!b || !c || !read_params(argv[3], &P, &steps, &valreg, &valstate, &beamsize, &maxval)) {
        printf("Error reading files\n");
        exit(EXIT_FAILURE);
    }
    data_size = sizeof(node) + sizeof(state)*b->size*c->size;
    printf("Parameters: %i %i %i %i %i %i\n", P, steps, valreg, valstate, beamsize,maxval);
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
