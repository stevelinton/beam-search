#include "beam.h"
#include <stdio.h>

#define P 5

#define MAXMOVE 16

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
    uint8_t s;
    uint8_t r;
    uint8_t nmoves;
    move moves[MAXMOVE];
    state states[P*P];
} node;

#define data_size sizeof(struct s_node)

static uint32_t fitness(const char *cv) {
    node *n = (node *)cv;
    if (n->s == P)
        return stop_fitness;
    return 10000 - 10*n->s - n->r;
}


static int ternary(uint8_t op, int in1, int in2, int in3) {
    return (op & (1 << (in1 *4 + in2 * 2 + in3))) ? 1 : 0;
}

static int binary(uint8_t op, int in1, int in2) {
    return (op & (1 << (in1 *2 + in2))) ? 1 : 0;
}

void drop(state *st, uint8_t bit) {
    uint16_t l = st->regs & (1 << bit -1);
    uint16_t r = st->regs & ~(2 << bit -1);
    st->regs = l | (r >> 1);
}



static void apply1(state *st, move *m, uint8_t nextreg) {
    int in1 = (st->regs & (1 << m->r1)) >> m->r1;
    int in2 = (st->regs & (1 << m->r2)) >> m->r2;
    if (m->isBin) {
        st->regs |= (binary(m->op, in1, in2) << nextreg);
    } else {
        int in3 = (st->regs & (1 << m->r3)) >> m->r3;
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

static bool apply(node *c, move *m) {
        for (int i = 0; i < c->s; i++)
            apply1(&(c->states[i]), m, c->r);
        // adjust number of registers
        if (m->isBin)
            c ->r += 1 - (m->drop &1) -((m->drop & 2) >> 1);
        else
            c->r += 1 - (m->drop & 1) - ((m->drop & 2) >>1) - ((m->drop & 4)>>2);
        // insertion sort the states
        for (int i = 1; i < c->s; i++) {
            state s = c->states[i];
            uint16_t x = s.regs;
            int j;
            for (j = i-1; j >= 0; j--) {
                uint16_t y = c->states[j].regs;
                if (y > x) {
                    c->states[j+1] = c->states[j];
                } else break;
            }
            c->states[j+1] = s;
        }
        // check for duplicate states, or incompatible states
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

        // record the move
        c->moves[c->nmoves++] = *m;
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

static void print_node(const char *i) {
    const node *n = (node *) i;
    printf("<node %i %i:",(int)n->r, (int) n->s);
    for (int j = 0; j < n->s; j++)
        printf(" 0x%hx->%i",n->states[j].regs, (int)n->states[j].res);
    printf(">");
}

static void print_move(move *m) {
    if (m->isBin) {
        printf("<move binop 0x%x regs %i %i drop 0x%x>",
               (int)m->op, (int)m->r1, (int)m->r2, (int)m->drop);
    } else {
        printf("<move ternop 0x%x regs %i %i %i drop 0x%x>",
               (int)m->op, (int)m->r1, (int)m->r2, (int)m->r3, (int)m->drop);
    }
}


static void visit_children(const char *parent, void visit(const char *, void *), void *context) {
    node *n = (node *)parent;
    /*    printf("VC ");
    print_node((const char *)n);
    printf("\n"); */
    node *ch = malloc(sizeof(node));
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
                    memcpy(ch, n, sizeof(node));
                        if (apply(ch, &m)) {
                            /*    printf("MOVE ");
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
                        memcpy(ch, n, sizeof(node));
                        if (apply(ch, &m)) {
                            /*    printf("MOVE ");
                                  print_move (&m);                            
                                  printf(" CHILD ");
                                  print_node((const char *)ch);
                                  printf("\n");  */
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
    return (0 == strncmp(a1,a2,data_size));
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

#if (P==3)
#define codelen 2
int codes[P] = {3,1,2};
#else
#define codelen 4
//int codes[P] = {0,1,2,3,4};
int codes[P] = {0,3,6,9,12};
#endif
           

int main(int argc, char **argv) {
    int beamsize = 50000;
    int nprobes = 4;
    char * seed = malloc(data_size);
    memset(seed, 0, data_size);
    ((node *)seed)->s = P*P;
    ((node *)seed)->r = 2*codelen;
    ((node *)seed)->nmoves = 0;
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++) {
            ((node *)seed)->states[P*i+j].regs  = (codes[i]<<codelen) | codes[j];
            ((node *)seed)->states[P*i+j].res  = (i+j)%P;
        }
    int nresults;
    char * results = beam_search(seed,1,visit_children, beamsize, MAXMOVE,
                                 data_size,  fitness, equal, hash, nprobes, print_node, &nresults);
    const char *winner  = NULL;
    const char *bestnode = NULL;
    int bestfit = 0;
    for (int i = 0; i < nresults; i++) {
        const char *c = results + i*data_size;
        int f = fitness(c);
        if (f == stop_fitness) {
            winner = c;
            break;
        } else if (f > bestfit) {
            bestnode = c;
            bestfit = f;
        }
    }
    if (!winner) {
        printf("No solution found\n");
        printf("Final beam size %i\nBest fitness %i\n",nresults, bestfit);
        print_node(bestnode);
        printf("\n");
        exit(EXIT_FAILURE);
    }
    printf("solution found\n");
    node *n = malloc(sizeof(node));
    memcpy(n, seed, sizeof(node));
    print_node(seed);
    printf("\n");
    for (int i = 0; i < ((node *)winner)->nmoves; i++) {
        move *m = ((node *)winner)->moves + i;
        print_move(m);
        printf(" -> ");
        apply(n,m);
        print_node((char *)n);
        printf("\n");
    }
        
    
    exit(EXIT_SUCCESS);
}
