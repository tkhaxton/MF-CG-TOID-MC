# define declare_array(type, name, size) type *name; name=xcalloc(size, sizeof(type)); for(i=0;i<size;i++) name[i]=0;
# define declare_array_nozero(type, name, size) type *name; name=xcalloc(size, sizeof(type));
# define allocate_array(type, name, size) name=xcalloc(size, sizeof(type));

# define declare_matrix(type, name, sizex, sizey) type **name; name=xcalloc(sizex, sizeof(type *)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type)); for(j=0;j<sizey;j++){  name[i][j]=0;}}
# define declare_matrix_nozero(type, name, sizex, sizey) type **name; name=xcalloc(sizex, sizeof(type *)); for(i=0;i<sizex;i++) name[i]=xcalloc(sizey, sizeof(type));
# define allocate_matrix(type, name, sizex, sizey) name=xcalloc(sizex, sizeof(type *)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type)); for(j=0;j<sizey;j++){  name[i][j]=0;}}
# define allocate_matrix_nozero(type, name, sizex, sizey) name=xcalloc(sizex, sizeof(type *)+1); for(i=0;i<sizex;i++) name[i]=xcalloc(sizey, sizeof(type));

# define declare_3d_tensor(type, name, sizex, sizey, sizez) type ***name; name=xcalloc(sizex, sizeof(type **)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type *)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type)); for(k=0;k<sizez;k++) name[i][j][k]=0;}}
# define allocate_3d_tensor(type, name, sizex, sizey, sizez) name=xcalloc(sizex, sizeof(type **)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type *)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type)); for(k=0;k<sizez;k++) name[i][j][k]=0;}}
# define allocate_3d_tensor_nozero(type, name, sizex, sizey, sizez) name=xcalloc(sizex, sizeof(type **)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type *)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type));}}

# define declare_4d_tensor(type, name, sizex, sizey, sizez, sizea) type ****name; name=xcalloc(sizex, sizeof(type ***)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type **)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type *)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type)); for(l=0;l<sizea;l++) name[i][j][k][l]=0;}}}
# define allocate_4d_tensor_nozero(type, name, sizex, sizey, sizez, sizea) name=xcalloc(sizex, sizeof(type ***)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type **)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type *)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type));}}}
# define allocate_4d_tensor(type, name, sizex, sizey, sizez, sizea) name=xcalloc(sizex, sizeof(type ***)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type **)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type *)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type)); for(l=0;l<sizea;l++) name[i][j][k][l]=0;}}}

# define declare_5d_tensor(type, name, sizex, sizey, sizez, sizea, sizeb) type *****name; name=xcalloc(sizex, sizeof(type ****)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type ***)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type **)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type *)); for(l=0;l<sizea;l++){name[i][j][k][l]=xcalloc(sizeb, sizeof(type)); for(m=0;m<sizeb;m++) name[i][j][k][l][m]=0;}}}}
# define allocate_5d_tensor(type, name, sizex, sizey, sizez, sizea, sizeb) name=xcalloc(sizex, sizeof(type ****)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type ***)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type **)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type *)); for(l=0;l<sizea;l++){name[i][j][k][l]=xcalloc(sizeb, sizeof(type)); for(m=0;m<sizeb;m++) name[i][j][k][l][m]=0;}}}}
# define declare_6d_tensor(type, name, sizex, sizey, sizez, sizea, sizeb, sizec) type ******name; name=xcalloc(sizex, sizeof(type *****)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type ****)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type ***)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type **)); for(l=0;l<sizea;l++){name[i][j][k][l]=xcalloc(sizeb, sizeof(type *)); for(m=0;m<sizeb;m++){ name[i][j][k][l][m]=xcalloc(sizec, sizeof(type)); for(n=0;n<sizec;n++) name[i][j][k][l][m][n]=0;}}}}}
# define allocate_6d_tensor(type, name, sizex, sizey, sizez, sizea, sizeb, sizec) name=xcalloc(sizex, sizeof(type *****)); for(i=0;i<sizex;i++){ name[i]=xcalloc(sizey, sizeof(type ****)); for(j=0;j<sizey;j++){ name[i][j]=xcalloc(sizez, sizeof(type ***)); for(k=0;k<sizez;k++){ name[i][j][k]=xcalloc(sizea, sizeof(type **)); for(l=0;l<sizea;l++){name[i][j][k][l]=xcalloc(sizeb, sizeof(type *)); for(m=0;m<sizeb;m++){ name[i][j][k][l][m]=xcalloc(sizec, sizeof(type)); for(n=0;n<sizec;n++) name[i][j][k][l][m][n]=0;}}}}}

# define free_matrix(name, sizex) {for(i=0;i<sizex;i++){ free(name[i]);}; free(name);}
# define free_3d_tensor(name, sizex, sizey) {for(i=0;i<sizex;i++){ for(j=0;j<sizey;j++) free(name[i][j]); free(name[i]);}; free(name);}
# define free_4d_tensor(name, sizex, sizey, sizez) {for(i=0;i<sizex;i++){ for(j=0;j<sizey;j++){for(k=0;k<sizez;k++) free(name[i][j][k]); free(name[i][j]);}; free(name[i]);}; free(name);}

typedef struct anytype_tag{
	int i;
	long long ll;
	double f;
	char *s;
} anytype;

void *xmalloc (size_t size);
void *xcalloc (size_t nelem,  size_t elsize);
anytype loadparam(int argc, char flag[100], char *argv[], char defaultchar[100]);
void loadseries(int argc, char flag[100], char *argv[], char *chararray[100]);
void my_exit(char *mychar);
int mygetline(char s[], int lim, FILE *inp);

