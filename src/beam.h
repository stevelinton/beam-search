#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Arrays of objects being searched on are represented as char * pointers (mainly const, since they are
generally only written to while being created.  Such an array needs a separate integer variable to record its size.

Parameters: seed and nseeds are the starting objects
            visit_children is called by the search with a parent object a visit function and a context (void *)
                        it should call the visit function on each child of the parent, passing the context as
                        the second argument.
            beamsize is the number of objects held in each generation
            ngens is the number of deearch generations to run.
            data_size is the size in bytes of an object
            fitness should return the fitness of an object. High fitness objects are preferred. Fitness 0 and 0xFFFFFFFF
                         should not be used
            equal   tests two objects for being equal. duplicates are discarded
            hash is a has function on objects. It should respect equality.
            nprobes is the number of times to try inserting each object into the hash table before 
                      giving up in general a small value effectively makes the search "more random" as low fitness
                      objects have more chance to survive.
            print_item is used for debugging
            nresults is used to indicate the number of results being returned (usually the beamsize)

the return value is the array of returned objects.
   
*/

#define stop_fitness 0xFFFFFFFE // any object with this fitness is considered "perfect"
// when it is discovered the search is terminated at the end of the current generation
// if you don't want this, just never use this fitness value.a

extern char *beam_search(
    const char *seeds, int nseeds,
    void visit_children(const char *, void (*)(const char *, void *), void *),
    int beamsize, int ngens, size_t data_size, uint32_t fitness(const char *),
    bool equal(const char *, const char *), uint64_t hash(const char *),
    int nprobes, void print_item(const char *), int *nresults);
