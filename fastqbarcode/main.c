//
//  main.c
//  Filter_fastq_index
//
//  Created by Nicolas Rapin on 08/09/2015.
//  Copyright (c) 2015 Nicolas Rapin. All rights reserved.
//
#include "uthash.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_RESIDUE_PAIR 1000001
#define MAX_RECORD_LEN 300

typedef struct data {
    struct data *next_data;
    int index_data;
    char    barcode[6]; /*barcode*/
    char *record;   /*actual record, which size is dynamic (dynamic string)*/
} DATA;


typedef struct barcode_data {
    struct barcode_data *next_data;
    int index_data;
    int amount;
    char    barcode[6]; /*barcode*/
} BARCODE_DATA;


// Hash stuffs //
struct barcode_hash_map_struct {
    int id;                    /*! key , here index of the barcode */
    char name[20];             /*! here sequence of the barcode*/
    UT_hash_handle hh;         /*! makes this structure hashable */
};

struct barcode_hash_map_struct *hash_map_barcode;


void hash_add_barcode(char  seq[], int  index){
    struct barcode_hash_map_struct *s = NULL;
    s = malloc(sizeof(struct barcode_hash_map_struct));
    strcpy(s->name, seq);
    s->id = index;
    HASH_ADD_STR( hash_map_barcode, name, s );
}
/*!find one element*/
int  barcode_hash_find(char * user_id) {
    struct barcode_hash_map_struct *s;
    
    HASH_FIND_STR( hash_map_barcode, user_id, s );  /*! s: output pointer */
    if (s)
        return s->id;
    else
        return -1;
}
// Hash stuffs //

int  compare(const void *i,const void *j)
{
    //printf("data i is %5e",((DATA**)i)[0]->score_data);
    if ( ((BARCODE_DATA**)i)[0]->amount < ((BARCODE_DATA**)j)[0]->amount) return -1;
    else
        if ( ((BARCODE_DATA**)i)[0]->amount > ((BARCODE_DATA**)j)[0]->amount) return 1;
        else
            return 0;
}

BARCODE_DATA	*barcode_data_allocate( void ){
    BARCODE_DATA	*new_entry;
    if ( ( new_entry = ( BARCODE_DATA * ) malloc ( sizeof(BARCODE_DATA))) != NULL ){
        new_entry->next_data = NULL;
        new_entry->index_data = 0;
        new_entry->amount = 0;
    }
    return ( new_entry );
}


DATA	*data_allocate( void ){
    DATA	*new_entry;
    if ( ( new_entry = ( DATA * ) malloc ( sizeof(DATA))) != NULL ){
        new_entry->next_data = NULL;
        new_entry->index_data = 0;
        new_entry->record = (char*)malloc(sizeof(char)*MAX_RECORD_LEN);
    }
    return ( new_entry );
}



void write_file(char** barcodes, DATA *All_data, int MAX_BARCODE_OUT, char * prefix) {
    DATA *new;
    int i;
    for (i=0;i < MAX_BARCODE_OUT; i++){
        fprintf(stderr,"Exporting fastq file for %s\n",barcodes[i]);
        char outfilename[400]= "";
        strcat(outfilename, prefix);
        strcat(outfilename, barcodes[i]);
        strcat(outfilename,".fastq");
        FILE *f = fopen(outfilename, "w");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        for (new=All_data; new; new=new->next_data) {
            //printf("%s,%s\n",new->barcode , barcodes[i]);
            if (strcmp(new->barcode , barcodes[i] ) == 0){
                fprintf(f, "%s", new->record);
            }
        }
        fclose(f);
    }
}



int main (int argc, char *argv[]){
    int frompipe=0;
    char fastq_file[400]="seq.fastq";
    char prefix[400] = "";
    int disp=0;
    char line[200]="\0";
    int index;
    int MAX_REC_LEN = 1;
    int MAX_BARCODE_OUT = 10;
    /*---------------- parse command-line options --------------------*/
    int n;
    char *s;
    for(n=1; n<argc; n++){
        s = argv[n];
        if(*s++=='-')
            switch(*s++) {
                case 'f':
                    strcpy(fastq_file,argv[n+1]);
                    break;
                case 'p':
                    strcpy(prefix,argv[n+1]);
                    break;
                case 'm':
                    MAX_BARCODE_OUT = atoi(argv[n+1]);
                    break;
                case 'd':
                    disp=1;
                    break;
                    
                case '-': /*get data from pipe instead of a file*/
                    frompipe=1;
                    break;
                    
                case 'h':
                    printf("\n\nUsage: fileter_fastq_index -f FASTQ.\n   -- get data from pipe.\n   -p prefix for output files\n   -m number of barcode output\n   -h:Displays this help message.\n   -d Verbose\n");
                    exit(-1);
                    break;
                    
                default:
                    printf("\n%s: Illegal option [%c]\n\nUsage: %s -fdh\n  -f [fastq file]: custom  file, default: 'seq.fastq'.\n  -p prefix for output files\n  -m number of barcode output\n  -- get data from pipe.\n  -d: Turn on verbose mode.\n  -h:Displays this help message.\n", argv[0],*--s,argv[0]);
                    exit(-1); break;
            }
    }
    
    /* -------------- End of cmd-line parser ------------------------- */
    
    
    /*read file*/
    FILE * seq_file = NULL;
    
    if(!frompipe){
        
        seq_file=fopen(fastq_file, "r");
        
        if(seq_file==NULL) {
            printf("Error: Please provide an input file.\n");
            exit (EXIT_FAILURE);}
    }
    
    /* idline, seq, spacer,qual and so on... */
    DATA *new = NULL;
    DATA *last_entry=NULL;
    DATA *All_data = NULL;
    BARCODE_DATA *barcode_new = NULL;
    BARCODE_DATA *barcode_last_entry=NULL;
    BARCODE_DATA *barcode_All_data = NULL;
    BARCODE_DATA **barcode_dataptr = NULL;
    int i = 0;
    size_t ln;
    char id_record[200] = "";
    char seq[200] = "";
    char spacer[200] = "";
    char qual[200] = "";
    int barcode_index = 0;
    if (!frompipe) {
        while ( fgets (line, 200, seq_file)!=NULL) {
            //line[strcspn(line, "\n")] = 0;
            if (disp) 	fprintf(stderr,"Read data file: ##%s##\n",line );
            switch(i){
                case 0 :
                    MAX_REC_LEN = 0;
                    new = data_allocate();
                    if ( ! new ){
                        printf( "Error. Cannot allocate new DATA element. Exit\n" );
                        exit( 1 );
                    }
                    /*add a first entry to the list*/
                    if ( ! All_data  )
                        All_data = new;
                    else
                        last_entry->next_data = new;
                    last_entry = new;
                    new->index_data = (int)index++;
                    strcpy(id_record , line);
                    ln = strlen(line);
                    int j=0;
                    while (line[ln-7] != '\n') {
                        new->barcode[j++] = line[ln-7];
                        ln++;
                    }
                    new->barcode[j] = '\0';
                    // update barcode linked list here.
                    
                    if (barcode_hash_find(new->barcode) == -1){  // barcode is not in the list...
                        hash_add_barcode(/*&*/new->barcode,/*&*/barcode_index++);
                        //add to barcode linked list.
                        barcode_new = barcode_data_allocate();
                        /*add a first entry to the list if needed*/
                        if ( ! barcode_All_data ){
                            barcode_All_data = barcode_new;
                            barcode_dataptr = malloc( sizeof(BARCODE_DATA *));
                            if (barcode_dataptr == NULL){
                                puts("\nFailure to allocate room for barcode_data *pointers");
                                exit(EXIT_FAILURE);
                            }
                            barcode_dataptr[0] = barcode_new;
                        }
                        else
                            barcode_last_entry->next_data = barcode_new;
                        barcode_dataptr = realloc(barcode_dataptr,(barcode_index) * sizeof(BARCODE_DATA *));
                        barcode_dataptr[barcode_index-1] = barcode_new;
                        barcode_last_entry = barcode_new;
                        barcode_new->index_data  = barcode_index-1;
                        barcode_new->amount = 1;
                        strcpy(barcode_new->barcode , new->barcode );
                        
                    }else{
                        // just add +1 to the amount of the barcode in the list.
                        barcode_dataptr[barcode_hash_find(new->barcode)]->amount +=1;
                    }
                    i+=1;
                    break;
                case 1 :
                    strcpy(seq, line);
                    i+=1;
                    break;
                case 2 :
                    strcpy(spacer , line);
                    i+=1;
                    break;
                case 3 :
                    strcpy(qual , line);
                    //now put all index, quls etc... in in string and free the rest,
                    sprintf( new->record, "%s%s%s%s", id_record,seq,spacer,qual);
                    i=0;
                    break;
            }
            
        }
        fclose(seq_file);
    }else {
        while ( fgets (line, 200, stdin)!=NULL) {
            //line[strcspn(line, "\n")] = 0;
            if (disp) 	fprintf(stderr,"Read data file: ##%s##\n",line );
            switch(i){
                case 0 :
                    MAX_REC_LEN = 0;
                    new = data_allocate();
                    if ( ! new ){
                        printf( "Error. Cannot allocate new DATA element. Exit\n" );
                        exit( 1 );
                    }
                    /*add a first entry to the list*/
                    if ( ! All_data  )
                        All_data = new;
                    else
                        last_entry->next_data = new;
                    last_entry = new;
                    new->index_data = (int)index++;
                    strcpy(id_record , line);
                    ln = strlen(line);
                    int j=0;
                    while (line[ln-7] != '\n') {
                        new->barcode[j++] = line[ln-7];
                        ln++;
                    }
                    new->barcode[j] = '\0';
                    // update barcode linked list here.
                    
                    if (barcode_hash_find(new->barcode) == -1){  // barcode is not in the list...
                        hash_add_barcode(/*&*/new->barcode,/*&*/barcode_index++);
                        //add to barcode linked list.
                        barcode_new = barcode_data_allocate();
                        /*add a first entry to the list if needed*/
                        if ( ! barcode_All_data ){
                            barcode_All_data = barcode_new;
                            barcode_dataptr = malloc( sizeof(BARCODE_DATA *));
                            if (barcode_dataptr == NULL){
                                puts("\nFailure to allocate room for barcode_data *pointers");
                                exit(EXIT_FAILURE);
                            }
                            barcode_dataptr[0] = barcode_new;
                        }
                        else
                            barcode_last_entry->next_data = barcode_new;
                        barcode_dataptr = realloc(barcode_dataptr,(barcode_index) * sizeof(BARCODE_DATA *));
                        barcode_dataptr[barcode_index-1] = barcode_new;
                        barcode_last_entry = barcode_new;
                        barcode_new->index_data  = barcode_index-1;
                        barcode_new->amount = 1;
                        strcpy(barcode_new->barcode , new->barcode );
                        
                    }else{
                        // just add +1 to the amount of the barcode in the list.
                        barcode_dataptr[barcode_hash_find(new->barcode)]->amount +=1;
                    }
                    i+=1;
                    break;
                case 1 :
                    strcpy(seq, line);
                    i+=1;
                    break;
                case 2 :
                    strcpy(spacer , line);
                    i+=1;
                    break;
                case 3 :
                    strcpy(qual , line);
                    //now put all index, quls etc... in in string and free the rest,
                    sprintf( new->record, "%s%s%s%s", id_record,seq,spacer,qual);
                    i=0;
                    break;
            }
            
        }
    }
    /*data read, put address of entry in a table*/
    int entry_nb=0;
    if(last_entry != NULL){
        entry_nb=last_entry->index_data;
    }else {
        fprintf(stderr,"Last Entry pointer not valid\n");
        exit(EXIT_FAILURE);
    }

/*
    DATA **dataptr;
    
    dataptr = malloc(entry_nb * sizeof(DATA *));
    if (dataptr == NULL)
    {
        puts("\nFailure to allocate room for data *pointers");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr,"indexing records...\n");
    i=0;
    for (new=All_data; new; new=new->next_data) {
        //printf("%s\n%s\n%d\n\n",new->barcode,new->record,i++);
        dataptr[ i++ ] = new ;
        if (disp) {
            printf("%s %d\n",dataptr[i-1]->barcode,i);
        }
    }
*/
    fprintf(stderr,"sorting barcodes out...\n");
    qsort(barcode_dataptr , barcode_last_entry->index_data, sizeof(BARCODE_DATA *), compare);
    
    char *barcodes[MAX_BARCODE_OUT];
    for (i=0;i < MAX_BARCODE_OUT; i++){
        barcodes[i] = malloc(strlen("\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0") + 1);
    }
        
    int j ;
    fprintf(stderr,"Sorted barcodes: \n");
    for( i=barcode_index - MAX_BARCODE_OUT-1, j=0 ; i < barcode_index -1 ; i++,j++ ){
        strcpy( barcodes[j],barcode_dataptr[i]->barcode);
        fprintf(stderr,"%s %d\n", barcode_dataptr[i]->barcode,barcode_dataptr[i]->amount);
    }
    
    write_file(barcodes, All_data, MAX_BARCODE_OUT,prefix);
    
    
    //    for (barcode_new=barcode_All_data; barcode_new; barcode_new=barcode_new->next_data) {
    //       printf("%s\t%d\n",barcode_new->barcode,barcode_new->amount);
    //    }
    
    //    printf("%s",dataptr[ 999]->record);
    //now output the fastq  file of the n most frequent files.
    
    return 0;
    
}
