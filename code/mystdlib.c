#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"


void *xmalloc (size_t size){
	register void *value = malloc (size);
	if (value == 0){
		printf("virtual memory exhausted");
		exit(1);
	}
	return value;
}

void *xcalloc (size_t nelem,  size_t elsize){
    void *new_mem = (void*)calloc(nelem ? nelem : 1, elsize ? elsize : 1);
    if (new_mem == NULL) {
        fprintf(stderr,
                "xcalloc: request for %lu elements of size %lu failed.\n",
                (unsigned long)nelem, (unsigned long)elsize);
        exit(EXIT_FAILURE);
    }    
    return new_mem;
}

anytype loadparam(int argc, char flag[100], char *argv[], char defaultchar[100]){
	int j;
	char *temp;
	temp=defaultchar;
	for(j=1;j<argc;j++){
		if(argv[j][0]=='-'){
			if(strcmp(flag, argv[j])==0){
				temp=argv[j+1];
			}
		}
	}
	anytype tempstruct;
	tempstruct.s=temp;
	tempstruct.i=atoi(temp);
	tempstruct.ll=atoll(temp);
	tempstruct.f=atof(temp);
	return tempstruct;
}

void loadseries(int argc, char flag[100], char *argv[], char *chararray[100]){		// different from inelastic_hard_sphere
	int i, j, slot;
	for(i=0;i<argc;i++){
		if(argv[i][0]=='-'){
			if(strcmp(flag, argv[i])==0){
				slot=i+1;
				j=0;
				while(argv[slot][0]!='-'){
					chararray[j]=argv[slot];
					j++;
					slot++;
				}
			}
		}
	}
}

void my_exit(char *mychar){
	printf("%s\n", mychar);
	exit(1);
}

int mygetline(char s[], int lim, FILE *inp){
	int c=0, i;
	for(i=0;i<lim-1&&(c=getc(inp))!=EOF&&c!='\n';i++){
		s[i]=c;
	}
	if(c=='\n'){
		s[i]=c;
		++i;
	}
	s[i]='\0';
	return i;
}


