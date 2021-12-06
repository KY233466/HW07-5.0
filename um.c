
/**************************************************************
 *                     um.c
 * 
 *     Assignment: Homework 6 - Universal Machine Program
 *     Authors: Katie Yang (zyang11) and Pamela Melgar (pmelga01)
 *     Date: November 24, 2021
 *
 *     Purpose: This C file will hold the main driver for our Universal
 *              MachineProgram (HW6). 
 *     
 *     Success Output:
 *              
 *     Failure output:
 *              1. 
 *              2. 
 *                  
 **************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "assert.h"
#include "uarray.h"
#include "seq.h"
#include "bitpack.h"
#include <sys/stat.h>

#define BYTESIZE 8


/* Instruction retrieval */
typedef enum Um_opcode {
        CMOV = 0, SLOAD, SSTORE, ADD, MUL, DIV,
        NAND, HALT, ACTIVATE, INACTIVATE, OUT, IN, LOADP, LV
} Um_opcode;


struct Info
{
    uint32_t op;
    uint32_t rA;
    uint32_t rB;
    uint32_t rC;
    uint32_t value;
};


typedef struct Info *Info;

#define num_registers 8
#define two_pow_32 4294967296
uint32_t MIN = 0;   
uint32_t MAX = 255;
uint32_t registers_mask = 511;
uint32_t op13_mask = 268435455;

uint32_t seg_size;
uint32_t seg_capacity;

///////////////////////////////////////////////////////////////////////////////////////////

struct Array_T
{
    uint32_t *array;
    uint32_t size;
};

typedef struct Array_T Array;


struct Memory
{
    Array *segments;
    Seq_T map_queue;
};

typedef struct Memory Memory;

///////////////////////////////////////////////////////////////////////////////////////////
/*Register Manager*/
#define num_registers 8

static inline uint32_t get_register_value(uint32_t *all_registers, uint32_t num_register);

/* Memory Manager */
Memory initialize_memory();
//uint32_t memorylength(Memory memory);         not used
uint32_t segmentlength(Memory *memory, uint32_t segment_index);
void add_to_seg0(Memory *memory, uint32_t word);
uint32_t get_word(Memory *memory, uint32_t segment_index, 
                  uint32_t word_in_segment);
void set_word(Memory *memory, uint32_t segment_index, uint32_t word_index, 
              uint32_t word);
uint32_t map_segment(Memory *memory, uint32_t num_words);
void unmap_segment(Memory *memory, uint32_t segment_index);
void duplicate_segment(Memory *memory, uint32_t segment_to_copy);
void free_segments(Memory *memory);
//////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////


/* Instruction retrieval */
Info get_Info(uint32_t instruction);

static inline void instruction_executer(Info info, Memory *all_segments,
                          uint32_t *all_registers, uint32_t *counter);
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*Making Helper functions for our ARRAYLIST*/

static inline uint32_t array_size(Array array);
static inline uint32_t element_at(Array array, uint32_t index);
static inline void replace_at(Array *array, uint32_t insert_value, uint32_t old_value_index);
static inline void push_at_back(Array *array, uint32_t insert_value);

static inline Array* array_at(Array *segments, uint32_t segment_index);
static inline void seg_push_at_back(Memory *all_segments, Array *insert_seg);
static inline void replace_array(Memory *all_segments, Array *insert_array, uint32_t old_array_index);

//////////////////////////////////////////////////////////////////////////////

    int main(int argc, char *argv[])
{
    /*Progam can only run with two arguments [execute] and [machinecode_file]*/
    
    if (argc != 2) {
        fprintf(stderr, "Usage: ./um [filename] < [input] > [output]\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(argv[1], "r");
    assert(fp != NULL);
    struct stat sb;

    const char *filename = argv[1];

    if (stat(filename, &sb) == -1) {
        perror("stat");
        exit(EXIT_FAILURE); // TODO what should the behavior be
    }

    int seg_0_size = sb.st_size / 4;

    // printf("seg_0_size is %d\n", seg_0_size);

////////////////////////////////////////////////////////////////////////////
    /* EXECUTION.C MODULE  */
    //init memory

    Memory all_segments;
    seg_size = 0;
    seg_capacity = 500;

    all_segments.segments = malloc((sizeof(Array)) * seg_capacity);
    all_segments.map_queue = Seq_new(30);

    Array segment0;
    // Array *segment0;
    segment0.array = malloc(sizeof(uint32_t) * seg_0_size);
    segment0.size = 0;

    // all_segments.segments[0] = segment0;

    seg_push_at_back(&all_segments, &segment0);
    // printf("seg 0 added\n");

    //init registers
    uint32_t all_registers[8] = { 0 };

/////////////////////////////////////// READ FILE ///////////////////////////
    uint32_t byte = 0; 
    uint32_t word = 0;
    int c;
    int counter = 0;

    /*Keep reading the file and parsing through "words" until EOF*/
    /*Populate the sequence every time we create a new word*/
    while ((c = fgetc(fp)) != EOF) {
        word = Bitpack_newu(word, 8, 24, c);
        counter++;
        for (int lsb = 2 * BYTESIZE; lsb >= 0; lsb = lsb - BYTESIZE) {
            byte = fgetc(fp);
            word = Bitpack_newu(word, BYTESIZE, lsb, byte);
        }
        
        /* Populate sequence */
        add_to_seg0(&all_segments, word);
    }

   fclose(fp);

//    printf("seg 0 populated\n");

   /////////////////////////////////////// execution /////////////////////

    uint32_t program_counter = 0; /*start program counter at beginning*/
    
    /*Run through all instructions in segment0, note that counter may loop*/
    while (program_counter < (segmentlength(&all_segments, 0)) ) {
        uint32_t instruction = get_word(&all_segments, 0, program_counter);
        Info info = get_Info(instruction);
        program_counter++;

        // printf("execute program_counter at %u\n", program_counter);

        /*program executer is passed in to update (in loop) if needed*/
        instruction_executer(info, &all_segments, all_registers, 
                             &program_counter);
    }

    free_segments(&all_segments);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return EXIT_FAILURE; /*failure mode, out of bounds of $m[0]*/
}


static inline void seg_push_at_back(Memory *all_segments, Array *insert_seg)
{
    if (seg_size != seg_capacity) {
        // printf("size %d capacity %d\n", seg_size, seg_capacity);
        all_segments->segments[seg_size] = *insert_seg;
        seg_size++;
    }
    else {
        // printf("----expand triggered\n");

        // uint32_t word = get_word(all_segments, 0, 1);
        // printf("This is a word before copying: %u \n", word);

        // word = get_word(all_segments, 5, 1);
        // printf("This is a word before copying: %u \n", word);

        // printf("expand\n");

        uint32_t new_capacity = seg_capacity * 2;
        seg_capacity = new_capacity;
        
        Array *new_segments = malloc(sizeof(Array) * new_capacity);

        for (uint32_t i = 0; i < seg_size; i++) {
            new_segments[i] = all_segments->segments[i]; //copying structs
        }

        free(all_segments->segments);

        all_segments->segments = new_segments;


        // uint32_t word_after = get_word(all_segments, 0, 1);
        // printf("This is a word after copying: %u \n", word_after);

        // word_after = get_word(all_segments, 5, 1);
        // printf("This is a word after copying: %u \n", word_after);

        all_segments->segments[seg_size] = *insert_seg;
        seg_size++;
        // printf("----expand complete\n");
    }
}

static inline uint32_t array_size(Array array1)
{
    return array1.size;
}

static inline uint32_t element_at(Array array1, uint32_t index)
{
    return array1.array[index];
}

static inline void replace_at(Array *array1, uint32_t insert_value, 
                                               uint32_t old_value_index)
{
    
    // printf("old_vale_index is %u insert_value is %u\n", old_value_index, insert_value);

    if (old_value_index == array1->size){
        push_at_back(array1, insert_value);
    }
    else{
        // int size = array1->size;
        // printf("size is %d\n", size);
        array1->array[old_value_index] = insert_value;
    }

    // printf("exit\n");
}

static inline void push_at_back(Array *array1, uint32_t insert_value)
{
    uint32_t size = array_size(*array1);
    array1->array[size] = insert_value;
    array1->size = size + 1;
}

uint32_t segmentlength(Memory *memory, uint32_t segment_index)
{
    uint32_t return_size = array_size(*(array_at(memory->segments, segment_index)));
    //printf("What is this return size? : %u\n", return_size);

    return return_size;
}

static inline Array *array_at(Array *arrays, uint32_t segment_index)
{
    return &(arrays[segment_index]);
}

void add_to_seg0(Memory *memory, uint32_t word)
{
    Array *segment_0 = (array_at(memory->segments, 0));

    push_at_back(segment_0, word);
}

uint32_t get_word(Memory *memory, uint32_t segment_index, 
                  uint32_t word_in_segment)
{
    // Seq_T find_segment = (Seq_T)(Seq_get(memory->segments, segment_index));

    // uint32_t *find_word = (uint32_t *)(Seq_get(find_segment, word_in_segment));

    Array target = *((array_at(memory->segments, segment_index)));

    uint32_t find_word = element_at(target, word_in_segment);

    return find_word;
}

void set_word(Memory *memory, uint32_t segment_index, 
              uint32_t word_index, uint32_t word)
{
    Array *target = array_at(memory->segments, segment_index);

    // printf("index is %u\n", segment_index);

    // printf("word is %u word_index is %u\n", word, word_index);
    replace_at(target, word, word_index);
}

uint32_t map_segment(Memory *memory, uint32_t num_words)
{
    // printf("^^^map seg\n");
    
    Array new_segment;
    new_segment.array = malloc((sizeof(uint32_t)) * num_words);
    new_segment.size = num_words;
    
    /*Initialize to ALL 0s*/
    for (uint32_t i = 0; i < num_words; i++) {
        replace_at(&new_segment, 0, i);
    }


    if (Seq_length(memory->map_queue) != 0) {
        uint32_t *seq_index = (uint32_t *)Seq_remlo(memory->map_queue);
        replace_array(memory, &new_segment, *seq_index);

        uint32_t segment_index = *seq_index;

        // printf("  reuse index %u\n", segment_index);
        free(seq_index);
        // printf("~ - ~map complete\n\n");
        return segment_index;
    }
    else {
        // printf("map %u\n", num_words);
        seg_push_at_back(memory, &new_segment);

        // printf("~ ~ ~map complete\n\n");
        return seg_size - 1;
    }
}

static inline void replace_array(Memory *all_segments, Array *insert_array, uint32_t old_array_index)
{
    // if (old_array_index != NULL){
    all_segments->segments[old_array_index] = *insert_array;
    // }
    // else {
        // seg_push_at_back(all_segments, insert_array);
    // }
    
    
    // if (old_array_index >= (seg_size - 1))
    // {
    //     printf("oh wow\n");
    //     seg_push_at_back(all_segments, insert_array);
    // }
    // else
    // {
    //     printf("lol\n");
    //     all_segments->segments[old_array_index] = *insert_array;
    // }
}

void unmap_segment(Memory *memory, uint32_t segment_index)
{
    // printf("-------unmap %u\n", segment_index);
    /* can't un-map segment 0 */
    assert(segment_index > 0);

    Array *seg_to_unmap = array_at(memory->segments, segment_index);
    /* can't un-map a segment that isn't mapped */

    free(seg_to_unmap->array);
    //free(seg_to_unmap);
    seg_to_unmap->array = NULL;


    uint32_t *num = malloc(sizeof(uint32_t));
    *num = segment_index;
    Seq_addhi(memory->map_queue, num);
    // printf("-------unmap complete\n");
}


void duplicate_segment(Memory *memory, uint32_t segment_to_copy)
{
    if (segment_to_copy != 0) {

        Array *seg_0 = array_at(memory->segments, 0);

        /*free all c array, that is  in segment 0*/
        free(seg_0->array);

        /*hard copy - duplicate array to create new segment 0*/
        Array *target = array_at(memory->segments, segment_to_copy);
        
        uint32_t target_seg_length = array_size(*target);

        seg_0->array = malloc(sizeof(uint32_t) * target_seg_length);
        seg_0->size = 0;

        /*Willl copy every single word onto the duplicate segment*/
        for (uint32_t i = 0; i < target_seg_length; i++) {
            uint32_t word = element_at(*target, i);
            push_at_back(seg_0, word);
        }

        /*replace segment0 with the duplicate*/
        //replace_array(memory, &duplicate, 0);                                     //this is problematic, duplicate goes out of scope
    } else {
        /*don't replace segment0 with itself --- do nothing*/
        return;
    }
}


void free_segments(Memory *memory)
{
    // printf("seg_size is %d\n", seg_size);
    for (uint32_t i = 0; i < seg_size; i++)
    {
        Array *target = array_at(memory->segments, i);

        if (target->array != NULL) {
            free(target->array);
        }
    }

    /* free map_queue Sequence that kept track of unmapped stuff*/
    uint32_t queue_length = Seq_length(memory->map_queue);
    
    /*free all the indexes words */
    for (uint32_t i = 0; i < queue_length; i++) {
        uint32_t *index = (uint32_t*)(Seq_get(memory->map_queue, i));
        free(index);
    }

    free(memory->segments);
    Seq_free(&(memory->map_queue));
}

struct Info *get_Info(uint32_t instruction)
{
    struct Info *info = malloc(sizeof(struct Info));
    assert(info != NULL); /*Check if heap allocation successful*/

    uint32_t op = instruction;
    op = op >> 28;

    info->op = op;

    if (op != 13){
        uint32_t registers_ins = instruction &= registers_mask;

        info->rA = registers_ins >> 6;
        info->rB = registers_ins << 26 >> 29;
        info->rC = registers_ins << 29 >> 29;
    }
    else {

        uint32_t registers_ins = instruction &= op13_mask;

        info->rA = registers_ins >> 25;
        info->value = registers_ins << 7 >> 7;
    }

    return info;
}


//change back for kcachegrind
static inline void instruction_executer(Info info, Memory *all_segments,
                           uint32_t *all_registers, uint32_t *counter)
{
    // printf("execute!\n");
    uint32_t code = info->op;

    /*We want a halt instruction to execute quicker*/
    if (code == HALT)
    {
        /////// halt /////
        free(info);
        
        free_segments(all_segments);

        exit(EXIT_SUCCESS);
    }

    /*Rest of instructions*/
    if (code == CMOV)
    {
        ///////// conditional_move /////
        uint32_t rA = info->rA;
        uint32_t rC = info->rC;

        uint32_t rC_val = get_register_value(all_registers, rC);

        if (rC_val != 0)
        {
            uint32_t rB = info->rB;

            uint32_t rB_val = get_register_value(all_registers, rB);

            //set_register_value(all_registers, rA, rB_val);
            uint32_t *word;
            word = &(all_registers[rA]);
            *word = rB_val;
        }
    }
    else if (code == SLOAD)
    {
        ////// segemented load /////
        /* Establishes which register indexes are being used */
        uint32_t rA = info->rA;
        uint32_t rB = info->rB;
        uint32_t rC = info->rC;

        /* Accesses the values at the register indexes*/
        uint32_t rB_val = get_register_value(all_registers, rB);
        uint32_t rC_val = get_register_value(all_registers, rC);

        uint32_t val = get_word(all_segments, rB_val, rC_val);

        //set_register_value(all_registers, rA, val);
        uint32_t *word;
        word = &(all_registers[rA]);
        *word = val;  
    }
    else if (code == SSTORE)
    {
        /////// segmented store /////
        /* Establishes which register indexes are being used */
        uint32_t rA = info->rA;
        uint32_t rB = info->rB;
        uint32_t rC = info->rC;

        /* Accesses the values at the register indexes*/
        uint32_t rA_val = get_register_value(all_registers, rA);
        uint32_t rB_val = get_register_value(all_registers, rB);
        uint32_t rC_val = get_register_value(all_registers, rC);

        set_word(all_segments, rA_val, rB_val, rC_val);
        
    }
    else if (code == ADD || code == MUL || code == DIV || code == NAND)
    {
        ////// arithetics /////
        uint32_t rA = info->rA;
        uint32_t rB = info->rB;
        uint32_t rC = info->rC;

        uint32_t rB_val = get_register_value(all_registers, rB);
        uint32_t rC_val = get_register_value(all_registers, rC);

        uint32_t value = 0;

        /*Determine which math operation to perform based on 4bit opcode*/
        if (code == ADD)
        {
            value = (rB_val + rC_val) % two_pow_32;
        }
        else if (code == MUL)
        {
            value = (rB_val * rC_val) % two_pow_32;
        }
        else if (code == DIV)
        {
            value = rB_val / rC_val;
        }
        else if (code == NAND)
        {
            value = ~(rB_val & rC_val);
        }
        else
        {
            exit(EXIT_FAILURE);
        }

        uint32_t *word;
        word = &(all_registers[rA]);
        *word = value;
    }
    else if (code == ACTIVATE)
    {
        ////////////// map_a_segment ////////
        uint32_t rB = info->rB;
        uint32_t rC_val = get_register_value(all_registers, info->rC);
        uint32_t mapped_index = map_segment(all_segments, rC_val);

        uint32_t *word;
        word = &(all_registers[rB]);
        *word = mapped_index;
    }
    else if (code == INACTIVATE)
    {
        ////////////// unmap_a_segment ////////
        uint32_t rC_val = get_register_value(all_registers, info->rC);
        unmap_segment(all_segments, rC_val);
    }
    else if (code == OUT)
    {
        ////////////// output ////////
        uint32_t rC = info->rC;
        uint32_t val = get_register_value(all_registers, rC);

        assert(val >= MIN);
        assert(val <= MAX);

        printf("%c", val);
    }
    else if (code == IN)
    {
        ////////////// input ////////
        uint32_t rC = info->rC;

        uint32_t input_value = (uint32_t)fgetc(stdin);
        uint32_t all_ones = ~0;

        /*Check if input value is EOF...aka: -1*/
        if (input_value == all_ones)
        {
            //set_register_value(all_registers, rC, all_ones);
            uint32_t *word;
            word = &(all_registers[rC]);
            *word = all_ones;
            
            return;
        }

        /* Check if the input value is in bounds */
        assert(input_value >= MIN);
        assert(input_value <= MAX);

        /* $r[C] gets loaded with input value */
        uint32_t *word;
        word = &(all_registers[rC]);
        *word = input_value;
    }
    else if (code == LOADP)
    {
        ////////////// load_program ////////
        uint32_t rB_val = all_registers[info->rB];

        uint32_t rC_val = get_register_value(all_registers, info->rC);
        duplicate_segment(all_segments, rB_val);
        *counter = rC_val;
    }
    else if (code == LV)
    {
        ////////////// load_value ////////
        uint32_t rA = info->rA;
        uint32_t val = info->value;

        uint32_t *word= &(all_registers[rA]);
        *word = val;
    }
    else
    {
        exit(EXIT_FAILURE);
    }

    free(info);
    return;
}

static inline uint32_t get_register_value(uint32_t *all_registers, uint32_t num_register)
{
    return all_registers[num_register]; // TODO ????
}
