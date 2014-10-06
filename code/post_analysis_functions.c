#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mymath.h"
#include "mystdlib.h"
#include "peptoid_functions.h"
#include "post_analysis_functions.h"

void assign_monomerid_and_chainlength(int Nchains, int uniformchainlength, monomernodes ***monomerid, int **chainlength){
	(*monomerid)=xmalloc(sizeof(monomernodes *)*Nchains+1);
    int i, j, nodecounter=0;
	(*chainlength)=xmalloc(sizeof(int)*Nchains+1);
    for(i=0;i<Nchains;i++){
		(*chainlength)[i]=uniformchainlength;
        (*monomerid)[i]=xmalloc(sizeof(monomernodes)*(*chainlength)[i]+1);
        for(j=0;j<(*chainlength)[i];j++){
            (*monomerid)[i][j].backbone=nodecounter;
			nodecounter++;
            (*monomerid)[i][j].sidechain=nodecounter;
			nodecounter++;
		}
	}
}

int input_trajectory(int Nchains, int chainlength, coord *coordarray, double_triple *pbox_dimension, FILE *inp, int *pcycle, int lim){
    int i, j, read, counter=0;
    declare_array(char, linechar, lim);
    for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength;j++){
            read=mygetline(linechar, lim, inp);
            sscanf(linechar, "%i %i %i %i %lf %lf %lf %lf %lf %lf", &(coordarray[counter].leafid), &(coordarray[counter].chainid), &(coordarray[counter].monomerid), &(coordarray[counter].nodetype), &(coordarray[counter].r.x), &(coordarray[counter].r.y), &(coordarray[counter].r.z), &(coordarray[counter].n.x), &(coordarray[counter].n.y), &(coordarray[counter].n.z));
            counter++;
            read=mygetline(linechar, lim, inp);
            sscanf(linechar, "%i %i %i %i %lf %lf %lf %lf %lf %lf", &(coordarray[counter].leafid), &(coordarray[counter].chainid), &(coordarray[counter].monomerid), &(coordarray[counter].nodetype), &(coordarray[counter].r.x), &(coordarray[counter].r.y), &(coordarray[counter].r.z), &(coordarray[counter].n.x), &(coordarray[counter].n.y), &(coordarray[counter].n.z));
            if(coordarray[counter].nodetype==1) coordarray[counter].orientationtype=1;
            else if(coordarray[counter].nodetype==2) coordarray[counter].orientationtype=0;
            else if(coordarray[counter].nodetype==3) coordarray[counter].orientationtype=0;
            else if(coordarray[counter].nodetype==-1) coordarray[counter].orientationtype=0;
			else{
				printf("sidechain nodetype=%i!\n", coordarray[counter].nodetype);
				exit(1);
			}
            counter++;
        }
    }
    read=mygetline(linechar, lim, inp);
    sscanf(linechar, "%i box %lf %lf %lf", &(*pcycle), &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
    read=mygetline(linechar, lim, inp);
	free(linechar);
    return read;
}

int input_trajectory_v0(int Nchains, int chainlength, coord *coordarray, double_triple *pbox_dimension, FILE *inp, int *pcycle, int lim, int leafflag){
    int i, j, read, counter=0;
    declare_array(char, linechar, lim);
    for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength;j++){
            read=mygetline(linechar, lim, inp);
			if(leafflag==1) (coordarray[counter].leafid)=0;
			else{
				if(i<Nchains/2) (coordarray[counter].leafid)=0;
				else (coordarray[counter].leafid)=1;
			}
            sscanf(linechar, "%i %i %i %lf %lf %lf %lf %lf %lf", &(coordarray[counter].chainid), &(coordarray[counter].monomerid), &(coordarray[counter].nodetype), &(coordarray[counter].r.x), &(coordarray[counter].r.y), &(coordarray[counter].r.z), &(coordarray[counter].n.x), &(coordarray[counter].n.y), &(coordarray[counter].n.z));
            counter++;
            read=mygetline(linechar, lim, inp);
			if(leafflag==1) (coordarray[counter].leafid)=0;
			else{
				if(i<Nchains/2) (coordarray[counter].leafid)=0;
				else (coordarray[counter].leafid)=1;
			}
            sscanf(linechar, "%i %i %i %lf %lf %lf %lf %lf %lf", &(coordarray[counter].chainid), &(coordarray[counter].monomerid), &(coordarray[counter].nodetype), &(coordarray[counter].r.x), &(coordarray[counter].r.y), &(coordarray[counter].r.z), &(coordarray[counter].n.x), &(coordarray[counter].n.y), &(coordarray[counter].n.z));
            if(coordarray[counter].nodetype==1) coordarray[counter].orientationtype=1;
            else if(coordarray[counter].nodetype==2) coordarray[counter].orientationtype=0;
            else if(coordarray[counter].nodetype==3) coordarray[counter].orientationtype=0;
			else{
				printf("sidechain nodetype=%i!\n", coordarray[counter].nodetype);
				exit(1);
			}
			counter++;
        }
    }
    read=mygetline(linechar, lim, inp);
    sscanf(linechar, "%i box %lf %lf %lf", &(*pcycle), &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
    read=mygetline(linechar, lim, inp);
	free(linechar);
    return read;
}

int input_trajectory_box(int Nchains, int chainlength, FILE *inp, int lim, int *pcycle, double_triple *pbox_dimension){
    int i, j, read;
    declare_array(char, linechar, lim);
    for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength;j++){
            read=mygetline(linechar, lim, inp);
            read=mygetline(linechar, lim, inp);
        }
    }
    read=mygetline(linechar, lim, inp);
    sscanf(linechar, "%i box %lf %lf %lf", &(*pcycle), &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
    read=mygetline(linechar, lim, inp);
    free(linechar);
    return read;
}

int input_trajectory_blank(int Nchains, int chainlength, FILE *inp, int lim){
    int i, j, read;
    declare_array(char, linechar, lim);
    for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength;j++){
            read=mygetline(linechar, lim, inp);
            read=mygetline(linechar, lim, inp);
        }
    }
    read=mygetline(linechar, lim, inp);
    read=mygetline(linechar, lim, inp);
    free(linechar);
    return read;
}

void poormansxrd_fromcoord_bin_leafid(int N, coord *coordarray, double_triple box_dimension, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int **edgeoncount, int *powdercount, int *radiallyaveragededgeoncount, double *****Iqedgeonhisto, double ****Iqpowderhisto, double ****Iqradiallyaveragededgeonhisto, int Nleaves){
    int i, j, k, l, x, y, z, inplaneradialbin, isotropicbin, sym;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, Nleaves);
	declare_matrix(double, Aqimag, Nnodetypes, Nleaves);
	double_triple resolution, q;
    resolution.x=2*M_PI/box_dimension.x;
    resolution.y=2*M_PI/box_dimension.y;
    resolution.z=2*M_PI/box_dimension.z;
	if(xrdqspace==0){
		xrdqspace=resolution.x;
		if((resolution.x!=resolution.y)||(resolution.x!=resolution.z)) my_exit("trying to use resolution-limited xrd, but box is not cubic");
	}
    int_triple range, bin;
    range.x=(int) floor((xrdqn+0.5)*xrdqspace/resolution.x);
    range.y=(int) floor((xrdqn+0.5)*xrdqspace/resolution.y);
    range.z=(int) floor((xrdqn+0.5)*xrdqspace/resolution.z);
    for(x=-range.x;x<=range.x;x++){
        q.x=resolution.x*x;
        bin.x=(int) floor(q.x/xrdqspace+0.5);
        for(y=-range.y;y<=range.y;y++){
            q.y=resolution.y*y;
            bin.y=(int) floor(q.y/xrdqspace+0.5);
            inplaneradialbin=(int) floor(sqrt(q.x*q.x+q.y*q.y)/xrdqspace+0.5);
            for(z=0;z<=xrdqn;z++){
                q.z=xrdqspace*z;                //  prescribed spacing, no periodic boundary conditions
				isotropicbin=(int) floor(sqrt(q.x*q.x+q.y*q.y+q.z*q.z)/xrdqspace+0.5);
                bin.z=(int) floor(q.z/xrdqspace+0.5);
                xrdcount[xrdqn+bin.x][xrdqn+bin.y][bin.z]++;
                if(inplaneradialbin<=xrdqn) edgeoncount[inplaneradialbin][bin.z]++;
                if(isotropicbin<=xrdqn){
                    powdercount[isotropicbin]++;
                    radiallyaveragededgeoncount[isotropicbin]++;
                }
				for(i=0;i<Nnodetypes;i++) for(j=0;j<Nleaves;j++) Aqreal[i][j]=Aqimag[i][j]=0;
				for(i=0;i<N;i++){
					qdotr=q.x*coordarray[i].r.x+q.y*coordarray[i].r.y+q.z*coordarray[i].r.z;
                    if(coordarray[i].nodetype>=0){
                        Aqreal[coordarray[i].nodetype][coordarray[i].leafid]+=cos(qdotr);
                        Aqimag[coordarray[i].nodetype][coordarray[i].leafid]+=sin(qdotr);
                    }
				}
				for(i=0;i<Nnodetypes;i++){
                    
                    for(k=0;k<Nleaves;k++){
                        for(l=0;l<Nleaves;l++){
                            if(k==l) sym=0;
                            else if(k/2==l/2) sym=1;
                            else sym=2;
                            val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                            Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][i][sym]+=val;
                            if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][i][sym]+=val;
                            if(isotropicbin<=xrdqn){
                                Iqpowderhisto[isotropicbin][i][i][sym]+=val;
                                if((x==0)&&(y==0)) Iqradiallyaveragededgeonhisto[isotropicbin][i][i][sym]+=val*(2./M_PI*sqrt((pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2)+q.z*q.z)/(pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2))));
                                else Iqradiallyaveragededgeonhisto[isotropicbin][i][i][sym]+=val*(2./M_PI*sqrt((q.x*q.x+q.y*q.y+q.z*q.z)/(q.x*q.x+q.y*q.y)));
                            }
                            for(j=i+1;j<Nnodetypes;j++){
                                val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                                Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][j][sym]+=val;
                                if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][j][sym]+=val;
                                if(isotropicbin<=xrdqn){
                                    Iqpowderhisto[isotropicbin][i][j][sym]+=val;
                                    if((x==0)&&(y==0)) Iqradiallyaveragededgeonhisto[isotropicbin][i][j][sym]+=val*(2./M_PI*sqrt((pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2)+q.z*q.z)/(pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2))));
                                    else Iqradiallyaveragededgeonhisto[isotropicbin][i][j][sym]+=val*(2./M_PI*sqrt((q.x*q.x+q.y*q.y+q.z*q.z)/(q.x*q.x+q.y*q.y)));
                                }
                            }
                        }
                    }
                    
				}
			}
		}
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
}

void sidechain_histo(int Nchains, int chainlength, coord *coordarray, double_triple box_dimension, sidechainhistoparams params, int **histo_rparallel, int **histo_rperp, int **histo_rdotn, int **histo_ndotn, int *count){
    int i, j, backindex, sideindex, type;
    double rparallel, r2, rperp, rdotn, ndotn;
    double_triple r;
    int rparallelbin, rperpbin, rdotnbin, ndotnbin;
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength;j++){
			backindex=2*(i*chainlength+j);
            sideindex=backindex+1;
            type=coordarray[sideindex].nodetype;
            if(type>-1){
                r=subtract_double_triple(coordarray[sideindex].r, coordarray[backindex].r);
                recenter_double_triple(&r, box_dimension);
                r2=dot_product(r, r);
                rparallel=dot_product(coordarray[backindex].n, r);
                rperp=sqrt(r2-rparallel*rparallel);
                rdotn=dot_product(coordarray[sideindex].n, r);
                ndotn=dot_product(coordarray[backindex].n, coordarray[sideindex].n);
                rparallelbin=(int) floor((rparallel-params.rparallelmin)/params.rparallelwidth);
                rperpbin=(int) floor(rperp/params.rperpwidth);
                rdotnbin=(int) floor((rdotn-params.rdotnmin)/params.rdotnwidth);
                ndotnbin=(int) floor((ndotn+1)/params.ndotnwidth);
                histo_rparallel[type][rparallelbin]++;
                histo_rperp[type][rperpbin]++;
                histo_rdotn[type][rdotnbin]++;
                histo_ndotn[type][ndotnbin]++;
                count[type]++;
            }
        }
    }
}

void backbone_histo(int Nchains, int chainlength, coord *coordarray, double_triple box_dimension, backbonehistoparams params, int *histo_r, int *histo_rleftdotn, int *histo_rrightdotn, int *histo_tripleproduct, int **histo_r_rleftdotn, int **histo_r_rrightdotn, int **histo_tripleproduct_rleftdotn, int **histo_tripleproduct_rrightdotn, int *pcount_r, int *pcount_rleftdotn, int *pcount_rrightdotn, int *pcount_tripleproduct, int *pcount_r_rleftdotn, int *pcount_r_rrightdotn, int *pcount_tripleproduct_rleftdotn, int *pcount_tripleproduct_rrightdotn){
	int i, j, backindex, rleftbin, rrightbin, rleftdotnbin, rrightdotnbin, tripleproductbin;
	double rrightnorm, rrightdotn, rleftdotn, tripleproduct;
	double_triple rright, rleft, n;
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength;j++){
			backindex=2*(i*chainlength+j);
			n=coordarray[backindex].n;
			if(j>0){
				rleft=rright;
				rleftbin=rrightbin;
				rleftdotn=dot_product(rleft, n);
				rleftdotnbin=(int) floor((rleftdotn+1)/params.dotwidth);
				if((rleftdotnbin<0)||(rleftdotnbin>=params.dotbins)){
					printf("rleftdotnbin=%i not in (0, %i)\n", rleftdotnbin, params.dotbins-1);
					exit(1);
				}
				histo_rleftdotn[rleftdotnbin]++;
				(*pcount_rleftdotn)++;
				histo_r_rleftdotn[rleftbin][rleftdotnbin]++;
				(*pcount_r_rleftdotn)++;
			}
			if(j<chainlength-1){
				rright=subtract_double_triple(coordarray[backindex+2].r, coordarray[backindex].r);
				recenter_double_triple(&rright, box_dimension);
				rrightnorm=norm(rright);
				rright=scalar_multiply_double_triple(rright, 1./rrightnorm);
				rrightbin=(int) floor((rrightnorm-params.rbottom)/params.rwidth);
				if((rrightbin<0)||(rrightbin>=params.rbins)){
					printf("rrightbin=%i not in (0, %i)\n", rrightbin, params.rbins-1);
					printf("rright=(%f, %f, %f)\n", rright.x, rright.y, rright.z);
					printf("rrightnorm=%f\n", rrightnorm);
					printf("box_dimension (%f, %f, %f)\n", box_dimension.x, box_dimension.y, box_dimension.z);
					exit(1);
				}
				histo_r[rrightbin]++;
				(*pcount_r)++;
				rrightdotn=dot_product(rright, n);
				rrightdotnbin=(int) floor((rrightdotn+1)/params.dotwidth);
				if((rrightdotnbin<0)||(rrightdotnbin>=params.dotbins)){
					printf("rrightdotnbin=%i not in (0, %i)\n", rrightdotnbin, params.dotbins-1);
					exit(1);
				}
				histo_rrightdotn[rrightdotnbin]++;
				(*pcount_rrightdotn)++;
				histo_r_rrightdotn[rrightbin][rrightdotnbin]++;
				(*pcount_r_rrightdotn)++;
				if(j>0){
					tripleproduct=dot_product(rleft, rright)-rleftdotn*rrightdotn;
					tripleproductbin=(int) floor((tripleproduct+1)/params.dotwidth);
					if((tripleproductbin<0)||(tripleproductbin>=params.dotbins)){
						printf("tripleproductbin=%i not in (0, %i)\n", tripleproductbin, params.dotbins-1);
						exit(1);
					}
					histo_tripleproduct[tripleproductbin]++;
					(*pcount_tripleproduct)++;
					histo_tripleproduct_rleftdotn[tripleproductbin][rleftdotnbin]++;
					(*pcount_tripleproduct_rleftdotn)++;
					histo_tripleproduct_rrightdotn[tripleproductbin][rrightdotnbin]++;
					(*pcount_tripleproduct_rrightdotn)++;
				}
                
			}
		}
	}
}

void gr_fromcoord(int N, coord *coordarray, double_triple box_dimension, double binwidth, double ****histo, int *pcount, int Nchains, int maxgrbin){
    int i, j, bin, type, leafi, leafj;
    double_triple sep;
    double r;
    for(i=0;i<N;i++){
        if(coordarray[i].nodetype>=0){
            leafi=coordarray[i].chainid/(Nchains/2);
            for(j=0;j<i;j++){
                if(coordarray[j].nodetype>=0){
                    leafj=coordarray[j].chainid/(Nchains/2);
                    sep=subtract_double_triple(coordarray[j].r, coordarray[i].r);
                    recenter_double_triple(&sep, box_dimension);
                    r=norm(sep);
                    bin=(int) floor(r/binwidth);
					if(bin<maxgrbin-1){
						if(coordarray[i].chainid==coordarray[j].chainid) type=0;        //  same poly
						else if(leafi==leafj) type=1;                                   //  same leaf
						else type=2;
						if(coordarray[i].nodetype<coordarray[j].nodetype) histo[bin][coordarray[i].nodetype][coordarray[j].nodetype][type]++;
						else histo[bin][coordarray[j].nodetype][coordarray[i].nodetype][type]++;
					}
                }
            }
        }
    }
    (*pcount)+=N;
}

void output_gr_3dnorm(char *filename, double ****histo, int Nsites, double binwidth, int Nbins, int count, double *fractions, int type){
    int i, j, k;
    FILE *outp;
    outp=fopen(filename, "w");
    for(i=0;i<Nbins;i++){
        fprintf(outp, "%.3f", binwidth*(0.5+i));
        for(j=0;j<Nsites;j++){
            fprintf(outp, " %.8g", 4*histo[i][j][j][type]/(fractions[j]*fractions[j]*4*M_PI*binwidth*pow((0.5+i)*count*binwidth, 2)));
            for(k=j+1;k<Nsites;k++){
                fprintf(outp, " %.8g", 4*histo[i][j][k][type]/(2*fractions[j]*fractions[k]*4*M_PI*binwidth*pow((0.5+i)*count*binwidth, 2)));
            }
        }
        fprintf(outp, "\n");
    }
    fclose(outp);
}

void output_gr_2dnorm(char *filename, double ****histo, int Nsites, double binwidth, int Nbins, int count, double *fractions, int type){
    int i, j, k;
    FILE *outp;
    outp=fopen(filename, "w");
    for(i=0;i<Nbins;i++){
        fprintf(outp, "%.3f", binwidth*(0.5+i));
        for(j=0;j<Nsites;j++){
            fprintf(outp, " %.8g", 4*histo[i][j][j][type]/(fractions[j]*fractions[j]*2*M_PI*binwidth*binwidth*(0.5+i)*count));        //  factor of 4 because I only did half of each sum
            for(k=j+1;k<Nsites;k++){
                fprintf(outp, " %.8g", 4*histo[i][j][k][type]/(2*fractions[j]*fractions[k]*2*M_PI*binwidth*binwidth*(0.5+i)*count));
            }
        }
        fprintf(outp, "\n");
    }
    fclose(outp);
}

void output_dotproduct(char *filename, int bins, int *histo, int count, double binwidth){
    FILE *outp;
    outp=fopen(filename, "w");
    int i;
    for(i=0;i<bins;i++){
        fprintf(outp, "%.3f %.8g\n", binwidth*(0.5+i)-1, histo[i]*2/(binwidth*count));
    }
    fclose(outp);
}

void output_unsigned_dotproduct_severaltypes(char *filename, int bins, int **histos, int *counts, double binwidth, int types, int countfactor){
    FILE *outp;
    outp=fopen(filename, "w");
    int i, j;
    for(i=0;i<bins;i++){
        fprintf(outp, "%.3f", binwidth*(0.5+i));
		for(j=0;j<types;j++) fprintf(outp, " %.8f", histos[j][i]/(binwidth*counts[j]*countfactor));
		fprintf(outp, "\n");
    }
    fclose(outp);
}

void dotproducthistosametype_nomixed(int N, coord *coordarray, double_triple box_dimension, int Nchains, double threshold, double binwidth, int *histo, int *pcount, int type){
    int i, j, leafi, leafj, bin;
    double_triple sep;
    double r, dotproduct;
    for(i=0;i<N;i++){
        if(coordarray[i].nodetype==type){
            leafi=coordarray[i].chainid/(Nchains/2);
            for(j=0;j<i;j++){
                if(coordarray[j].nodetype==type){
                    leafj=coordarray[j].chainid/(Nchains/2);
                    if(leafi==leafj){
                        if(coordarray[i].chainid!=coordarray[j].chainid){
                            sep=subtract_double_triple(coordarray[j].r, coordarray[i].r);
                            recenter_double_triple(&sep, box_dimension);
                            r=norm(sep);
                            if(r<threshold){
                                dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                histo[bin]++;
                                (*pcount)++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void dotproducthistosametype(int N, coord *coordarray, double_triple box_dimension, int Nchains, double threshold, double binwidth, int *histo, int *mixedhisto, int *pcount, int type){
    int i, j, leafi, leafj, bin;
    double_triple sep;
    double r, dotproduct, mixeddotproduct;
    for(i=0;i<N;i++){
        if(coordarray[i].nodetype==type){
            leafi=coordarray[i].chainid/(Nchains/2);
            for(j=0;j<i;j++){
                if(coordarray[j].nodetype==type){
                    leafj=coordarray[j].chainid/(Nchains/2);
                    if(leafi==leafj){
                        if(coordarray[i].chainid!=coordarray[j].chainid){
                            sep=subtract_double_triple(coordarray[j].r, coordarray[i].r);
                            recenter_double_triple(&sep, box_dimension);
                            r=norm(sep);
							sep=scalar_multiply_double_triple(sep, 1./r);
                            if(r<threshold){
                                dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                histo[bin]++;
								dotproduct=dot_product(coordarray[j].n, sep);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                mixedhisto[bin]++;
								dotproduct=-dot_product(coordarray[i].n, sep);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                mixedhisto[bin]++;
								(*pcount)++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void dotproducthistoalltypes(int N, coord *coordarray, double_triple box_dimension, int Nchains, double *thresholds, double binwidth, int **histo, int **mixedaddhisto, int **mixedsubtracthisto, int *count, int type, int issigned){
    int i, j, leafi, leafj, bin, mixedbinadd, mixedbinsubtract;
    double_triple sep;
    double r, dotproduct, mixeddotproduct, dotproduct2;
    for(i=0;i<N;i++){
        if(coordarray[i].nodetype==type){
            leafi=coordarray[i].chainid/(Nchains/2);
            for(j=0;j<i;j++){
                if(coordarray[j].nodetype==type){
                    leafj=coordarray[j].chainid/(Nchains/2);
					sep=subtract_double_triple(coordarray[j].r, coordarray[i].r);
					recenter_double_triple(&sep, box_dimension);
					r=norm(sep);
					if(leafi==leafj){
						if(coordarray[i].chainid==coordarray[j].chainid) type=0;
						else type=1;
					}
					else type=2;
					if(r<thresholds[type]){
						sep=scalar_multiply_double_triple(sep, 1./r);
						if(issigned==0){
							dotproduct=fabs(dot_product(coordarray[i].n, coordarray[j].n));
							bin=(int) floor((dotproduct)/binwidth);
							dotproduct=dot_product(coordarray[j].n, sep);
							dotproduct2=dot_product(coordarray[i].n, sep);
							mixedbinadd=(int) floor(fabs(dotproduct2+dotproduct)/binwidth);
							mixedbinsubtract=(int) floor(fabs(dotproduct2-dotproduct)/binwidth);
						}
						else{
							dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
							bin=(int) floor((dotproduct+1)/binwidth);
							dotproduct=dot_product(coordarray[j].n, sep);
							dotproduct2=dot_product(coordarray[i].n, sep);
							mixedbinadd=(int) floor((dotproduct2+dotproduct+1)/binwidth);
							mixedbinsubtract=(int) floor((dotproduct2-dotproduct+1)/binwidth);
						}
						histo[type][bin]++;
						mixedaddhisto[type][mixedbinadd]++;
						mixedsubtracthisto[type][mixedbinsubtract]++;
						count[type]++;
					}
                }
            }
        }
    }
}



void dotproducthistodifferenttypetwoshells(int N, coord *coordarray, double_triple box_dimension, int Nchains, double threshold, double outerthreshold, double binwidth, int *histo, int *pcount, int *outerhisto, int *poutercount, int typea, int typeb){
    int i, j, leafi, leafj, bin;
    double_triple sep;
    double r, dotproduct;
    for(i=0;i<N;i++){
        if(coordarray[i].nodetype==typea){
            leafi=coordarray[i].chainid/(Nchains/2);
            for(j=0;j<i;j++){
                if(coordarray[j].nodetype==typeb){
                    leafj=coordarray[j].chainid/(Nchains/2);
                    if(leafi==leafj){
                        if(coordarray[i].chainid!=coordarray[j].chainid){
                            sep=subtract_double_triple(coordarray[j].r, coordarray[i].r);
                            recenter_double_triple(&sep, box_dimension);
                            r=norm(sep);
                            if(r<threshold){
                                dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                histo[bin]++;
                                (*pcount)++;
                            }
                            else if(r<outerthreshold){
                                dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                outerhisto[bin]++;
                                (*poutercount)++;
                            }
                        }
                    }
                }
            }
        }
        else if(coordarray[i].nodetype==typeb){
            leafi=coordarray[i].chainid/(Nchains/2);
            for(j=0;j<i;j++){
                if(coordarray[j].nodetype==typea){
                    leafj=coordarray[j].chainid/(Nchains/2);
                    if(leafi==leafj){
                        if(coordarray[i].chainid!=coordarray[j].chainid){
                            sep=subtract_double_triple(coordarray[j].r, coordarray[i].r);
                            recenter_double_triple(&sep, box_dimension);
                            r=norm(sep);
                            if(r<threshold){
                                dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                histo[bin]++;
                                (*pcount)++;
                            }
                            else if(r<threshold){
                                dotproduct=dot_product(coordarray[i].n, coordarray[j].n);
                                bin=(int) floor((dotproduct+1)/binwidth);
                                outerhisto[bin]++;
                                (*poutercount)++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void output_backbone_histo(char *base, backbonehistoparams params, int *histo_r, int *histo_rleftdotn, int *histo_rrightdotn, int *histo_tripleproduct, int **histo_r_rleftdotn, int **histo_r_rrightdotn, int **histo_tripleproduct_rleftdotn, int **histo_tripleproduct_rrightdotn, int *pcount_r, int *pcount_rleftdotn, int *pcount_rrightdotn, int *pcount_tripleproduct, int *pcount_r_rleftdotn, int *pcount_r_rrightdotn, int *pcount_tripleproduct_rleftdotn, int *pcount_tripleproduct_rrightdotn){
	declare_array_nozero(char, filename, 400);
	int i, j;
	FILE *outp;
	sprintf(filename, "%s.normr", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.rbins;i++){
		fprintf(outp, "%.3f %.8g\n", params.rbottom+(i+0.5)*params.rwidth, (1.*histo_r[i])/(1.*(*pcount_r))/params.rwidth);
		histo_r[i]=0;
	}
	fclose(outp);
	sprintf(filename, "%s.rleftdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.dotbins;i++){
		fprintf(outp, "%.3f %.8g\n", -1+(i+0.5)*params.dotwidth, (1.*histo_rleftdotn[i])/(1.*(*pcount_rleftdotn))/params.dotwidth);
		histo_rleftdotn[i]=0;
	}
	fclose(outp);
	sprintf(filename, "%s.rrightdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.dotbins;i++){
		fprintf(outp, "%.3f %.8g\n", -1+(i+0.5)*params.dotwidth, (1.*histo_rrightdotn[i])/(1.*(*pcount_rrightdotn))/params.dotwidth);
		histo_rrightdotn[i]=0;
	}
	fclose(outp);
	sprintf(filename, "%s.tripleproduct", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.dotbins;i++){
		fprintf(outp, "%.3f %.8g\n", -1+(i+0.5)*params.dotwidth, (1.*histo_tripleproduct[i])/(1.*(*pcount_tripleproduct))/params.dotwidth);
		histo_tripleproduct[i]=0;
	}
	fclose(outp);
	sprintf(filename, "%s.normr.rleftdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.rbins;i++){
		for(j=0;j<params.dotbins;j++){
			fprintf(outp, "%.3f %.3f %.8g\n", params.rbottom+(i+0.5)*params.rwidth, -1+(j+0.5)*params.dotwidth, (1.*histo_r_rleftdotn[i][j])/(1.*(*pcount_r_rleftdotn))/params.rwidth/params.dotwidth);
			histo_r_rleftdotn[i][j]=0;
		}
	}
	fclose(outp);
	sprintf(filename, "%s.normr.rrightdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.rbins;i++){
		for(j=0;j<params.dotbins;j++){
			fprintf(outp, "%.3f %.3f %.8g\n", params.rbottom+(i+0.5)*params.rwidth, -1+(j+0.5)*params.dotwidth, (1.*histo_r_rrightdotn[i][j])/(1.*(*pcount_r_rrightdotn))/params.rwidth/params.dotwidth);
			histo_r_rrightdotn[i][j]=0;
		}
	}
	fclose(outp);
	sprintf(filename, "%s.tripleproduct.rleftdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.dotbins;i++){
		for(j=0;j<params.dotbins;j++){
			fprintf(outp, "%.3f %.3f %.8g\n", -1+(i+0.5)*params.dotwidth, -1+(j+0.5)*params.dotwidth, (1.*histo_tripleproduct_rleftdotn[i][j])/(1.*(*pcount_tripleproduct_rleftdotn))/params.dotwidth/params.dotwidth);
			histo_tripleproduct_rleftdotn[i][j]=0;
		}
	}
	fclose(outp);
	sprintf(filename, "%s.tripleproduct.rrightdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.dotbins;i++){
		for(j=0;j<params.dotbins;j++){
			fprintf(outp, "%.3f %.3f %.8g\n", -1+(i+0.5)*params.dotwidth, -1+(j+0.5)*params.dotwidth, (1.*histo_tripleproduct_rrightdotn[i][j])/(1.*(*pcount_tripleproduct_rrightdotn))/params.dotwidth/params.dotwidth);
			histo_tripleproduct_rrightdotn[i][j]=0;
		}
	}
	fclose(outp);
	(*pcount_r)=(*pcount_rleftdotn)=(*pcount_rrightdotn)=(*pcount_tripleproduct)=(*pcount_r_rrightdotn)=(*pcount_r_rleftdotn)=(*pcount_tripleproduct_rrightdotn)=(*pcount_tripleproduct_rleftdotn)=0;
}

void output_sidechain_histo(char *base, sidechainhistoparams params, int **histo_rparallel, int **histo_rperp, int **histo_rdotn, int **histo_ndotn, int *count){
	declare_array_nozero(char, filename, 400);
	int i, j, firsttype=1, lasttype=3;
	FILE *outp;
	sprintf(filename, "%s.rparallel", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.bins;i++){
        fprintf(outp, "%.3f", params.rparallelmin+(i+0.5)*params.rparallelwidth);
        for(j=firsttype;j<=lasttype;j++) fprintf(outp, " %.8g", (1.*histo_rparallel[j][i])/(1.*count[j]*params.rparallelwidth));
        fprintf(outp, "\n");
	}
	fclose(outp);
	sprintf(filename, "%s.rperp", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.bins;i++){
        fprintf(outp, "%.3f", (i+0.5)*params.rperpwidth);
        for(j=firsttype;j<=lasttype;j++) fprintf(outp, " %.8g", (1.*histo_rperp[j][i])/(1.*count[j]*params.rperpwidth));
        fprintf(outp, "\n");
	}
	fclose(outp);
	sprintf(filename, "%s.rdotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.bins;i++){
        fprintf(outp, "%.3f", params.rdotnmin+(i+0.5)*params.rdotnwidth);
        for(j=firsttype;j<=lasttype;j++) fprintf(outp, " %.8g", (1.*histo_rdotn[j][i])/(1.*count[j]*params.rdotnwidth));
        fprintf(outp, "\n");
	}
	fclose(outp);
	sprintf(filename, "%s.ndotn", base);
	outp=fopen(filename, "w");
	for(i=0;i<params.bins;i++){
        fprintf(outp, "%.3f", -1+(i+0.5)*params.ndotnwidth);
        for(j=firsttype;j<=lasttype;j++) fprintf(outp, " %.8g", (1.*histo_ndotn[j][i])/(1.*count[j]*params.ndotnwidth));
        fprintf(outp, "\n");
	}
	fclose(outp);
}

void output_3d_xrd_leaves(char *base, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension){
	if(xrdqspace==0) xrdqspace=2*M_PI/box_dimension.x;
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int nx, ny, nz, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(nx=-xrdqn;nx<=xrdqn;nx++){
        for(ny=-xrdqn;ny<=xrdqn;ny++){
            for(nz=0;nz<=xrdqn;nz++){
                fprintf(rawoutp, "%.4f\t%.4f\t%.4f", xrdqspace*nx, xrdqspace*ny, xrdqspace*nz);
                fprintf(weightedoutp, "%.4f\t%.4f\t%.4f", xrdqspace*nx, xrdqspace*ny, xrdqspace*nz);
                fprintf(weightedtotaloutp, "%.4f\t%.4f\t%.4f", xrdqspace*nx, xrdqspace*ny, xrdqspace*nz);
                total[0]=total[1]=0;
                for(i=0;i<Nnodetypes;i++){
                    for(j=i;j<Nnodetypes;j++){
                        for(l=0;l<symmetries;l++){
                            val=Iqhisto[xrdqn+nx][xrdqn+ny][nz][i][j][l]/(1.*xrdcount[xrdqn+nx][xrdqn+ny][nz]*Nresidues);
                            fprintf(rawoutp, "\t%.8g", val);
                            val=val*scatterfactor[i]*scatterfactor[j];
                            fprintf(weightedoutp, "\t%.8g", val);
                            total[l]+=val;
                        }
                    }
                }
                fprintf(rawoutp, "\n");
                fprintf(weightedoutp, "\n");
				for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
				fprintf(weightedtotaloutp, "\n");
            }
        }
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

void output_2d_xrd_from3d_z0_leaves(char *base, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension){
	if(xrdqspace==0) xrdqspace=2*M_PI/box_dimension.x;
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int nx, ny, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(nx=-xrdqn;nx<=xrdqn;nx++){
        for(ny=-xrdqn;ny<=xrdqn;ny++){
			fprintf(rawoutp, "%.4f\t%.4f", xrdqspace*nx, xrdqspace*ny);
			fprintf(weightedoutp, "%.4f\t%.4f", xrdqspace*nx, xrdqspace*ny);
			fprintf(weightedtotaloutp, "%.4f\t%.4f", xrdqspace*nx, xrdqspace*ny);
            total[0]=total[1]=0;
			for(i=0;i<Nnodetypes;i++){
				for(j=i;j<Nnodetypes;j++){
                    for(l=0;l<symmetries;l++){
                        val=Iqhisto[xrdqn+nx][xrdqn+ny][0][i][j][l]/(1.*xrdcount[xrdqn+nx][xrdqn+ny][0]*Nresidues);
                        fprintf(rawoutp, "\t%.8g", val);
                        val=val*scatterfactor[i]*scatterfactor[j];
                        fprintf(weightedoutp, "\t%.8g", val);
                        total[l]+=val;
                    }
				}
			}
			fprintf(rawoutp, "\n");
			fprintf(weightedoutp, "\n");
			for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
			fprintf(weightedtotaloutp, "\n");
        }
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

void output_2d_xrd_from3d_z0_leaves_resolutionlimited(char *base, double ******Iqhisto, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int ***xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple avbox, int boxcount){
	if(xrdqspace.x==0) xrdqspace.x=2*M_PI/avbox.x*boxcount;
	if(xrdqspace.y==0) xrdqspace.x=2*M_PI/avbox.y*boxcount;
	if(xrdqspace.z==0) xrdqspace.x=2*M_PI/avbox.z*boxcount;
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int nx, ny, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(nx=-xrdqn.x;nx<=xrdqn.x;nx++){
        for(ny=-xrdqn.y;ny<=xrdqn.y;ny++){
			fprintf(rawoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(weightedoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(weightedtotaloutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
            total[0]=total[1]=0;
			for(i=0;i<Nnodetypes;i++){
				for(j=i;j<Nnodetypes;j++){
                    for(l=0;l<symmetries;l++){
                        val=Iqhisto[xrdqn.x+nx][xrdqn.y+ny][0][i][j][l]/(1.*xrdcount[xrdqn.x+nx][xrdqn.y+ny][0]*Nresidues);
                        fprintf(rawoutp, "\t%.8g", val);
                        val=val*scatterfactor[i]*scatterfactor[j];
                        fprintf(weightedoutp, "\t%.8g", val);
                        total[l]+=val;
                    }
				}
			}
			fprintf(rawoutp, "\n");
			fprintf(weightedoutp, "\n");
			for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
			fprintf(weightedtotaloutp, "\n");
        }
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

void output_2d_xrd_leaves_resolutionlimited(char *base, double *****Iqhisto2d, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple avbox){

    printf("avbox=%f, %f, %f\n", avbox.x/(1.*xrdcount), avbox.y/(1.*xrdcount), avbox.z/(1.*xrdcount));
	if(xrdqspace.x==0) xrdqspace.x=2*M_PI/avbox.x*xrdcount;
	if(xrdqspace.y==0) xrdqspace.y=2*M_PI/avbox.y*xrdcount;

	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int nx, ny, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(nx=-xrdqn.x;nx<=xrdqn.x;nx++){
        for(ny=-xrdqn.y;ny<=xrdqn.y;ny++){
			fprintf(rawoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(weightedoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(weightedtotaloutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
            total[0]=total[1]=0;
			for(i=0;i<Nnodetypes;i++){
				for(j=i;j<Nnodetypes;j++){
                    for(l=0;l<symmetries;l++){
                        val=Iqhisto2d[xrdqn.x+nx][xrdqn.y+ny][i][j][l]/(1.*xrdcount*Nresidues);
                        fprintf(rawoutp, "\t%.8g", val);
                        val=val*scatterfactor[i]*scatterfactor[j];
                        fprintf(weightedoutp, "\t%.8g", val);
                        total[l]+=val;
                    }
				}
			}
			fprintf(rawoutp, "\n");
			fprintf(weightedoutp, "\n");
			for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
			fprintf(weightedtotaloutp, "\n");
        }
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

void output_2d_xrd_leaves_resolutionlimited_andradav(char *base, double *****Iqhisto2d, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, double *neutronscatterfactor, int symmetries, double_triple avbox, double binspace){
    
    printf("avbox=%f, %f, %f\n", avbox.x/(1.*xrdcount), avbox.y/(1.*xrdcount), avbox.z/(1.*xrdcount));
	if(xrdqspace.x==0) xrdqspace.x=2*M_PI/avbox.x*xrdcount;
	if(xrdqspace.y==0) xrdqspace.y=2*M_PI/avbox.y*xrdcount;
    
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
	declare_array_nozero(char, neutronweightedfile, charlength);
	declare_array_nozero(char, neutronweightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp, *neutronweightedoutp, *neutronweightedtotaloutp;
    int nx, ny, i, j, k, l, radialbin, maxradbin;
	maxradbin=(int) (xrdqspace.x*xrdqn.x/binspace);
	if(((int) (xrdqspace.y*xrdqn.y/binspace))<maxradbin) maxradbin=(int) (xrdqspace.y*xrdqn.y/binspace);
	declare_4d_tensor(double, radav, maxradbin, Nnodetypes, Nnodetypes, symmetries);
	declare_4d_tensor(int, radcount, maxradbin, Nnodetypes, Nnodetypes, symmetries);
    double val, Xrayval, neutronval;
    declare_array(double, total, symmetries);
    declare_array(double, neutrontotal, symmetries);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
    sprintf(neutronweightedfile, "%s.neutron.weighted", base);
    sprintf(neutronweightedtotalfile, "%s.neutron.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
	neutronweightedoutp=fopen(neutronweightedfile, "w");
	neutronweightedtotaloutp=fopen(neutronweightedtotalfile, "w");
    for(nx=-xrdqn.x;nx<=xrdqn.x;nx++){
        for(ny=-xrdqn.y;ny<=xrdqn.y;ny++){
            radialbin=(int) floor(sqrt(pow(xrdqspace.x*nx, 2)+pow(xrdqspace.y*ny, 2))/binspace+0.5);
			fprintf(rawoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(weightedoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(weightedtotaloutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(neutronweightedoutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			fprintf(neutronweightedtotaloutp, "%.4f\t%.4f", xrdqspace.x*nx, xrdqspace.y*ny);
			for(l=0;l<symmetries;l++){
				total[l]=0;
				neutrontotal[l]=0;
			}
			for(i=0;i<Nnodetypes;i++){
				for(j=i;j<Nnodetypes;j++){
                    for(l=0;l<symmetries;l++){
                        val=Iqhisto2d[xrdqn.x+nx][xrdqn.y+ny][i][j][l]/(1.*xrdcount*Nresidues);
						if(radialbin<maxradbin){
							radav[radialbin][i][j][l]+=val;
							radcount[radialbin][i][j][l]++;
						}
                        fprintf(rawoutp, "\t%.8g", val);
                        Xrayval=val*scatterfactor[i]*scatterfactor[j];
                        fprintf(weightedoutp, "\t%.8g", Xrayval);
                        total[l]+=Xrayval;
                        neutronval=val*neutronscatterfactor[i]*neutronscatterfactor[j];
                        fprintf(neutronweightedoutp, "\t%.8g", neutronval);
                        neutrontotal[l]+=neutronval;
                    }
				}
			}
			fprintf(rawoutp, "\n");
			fprintf(weightedoutp, "\n");
			for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
			fprintf(weightedtotaloutp, "\n");
			fprintf(neutronweightedoutp, "\n");
			for(l=0;l<symmetries;l++) fprintf(neutronweightedtotaloutp, "\t%.8g", neutrontotal[l]);
			fprintf(neutronweightedtotaloutp, "\n");
        }
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    fclose(neutronweightedoutp);
    fclose(neutronweightedtotaloutp);
    sprintf(rawfile, "%s.radav.raw", base);
    sprintf(weightedfile, "%s.radav.weighted", base);
    sprintf(weightedtotalfile, "%s.radav.total", base);
    sprintf(neutronweightedfile, "%s.radav.neutron.weighted", base);
    sprintf(neutronweightedtotalfile, "%s.radav.neutron.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
	neutronweightedoutp=fopen(neutronweightedfile, "w");
	neutronweightedtotaloutp=fopen(neutronweightedtotalfile, "w");
	for(radialbin=0;radialbin<maxradbin;radialbin++){
		fprintf(rawoutp, "%.4f", binspace*radialbin);
		fprintf(weightedoutp, "%.4f", binspace*radialbin);
		fprintf(weightedtotaloutp, "%.4f", binspace*radialbin);
		fprintf(neutronweightedoutp, "%.4f", binspace*radialbin);
		fprintf(neutronweightedtotaloutp, "%.4f", binspace*radialbin);
		for(l=0;l<symmetries;l++){
			total[l]=0;
			neutrontotal[l]=0;
		}
		for(i=0;i<Nnodetypes;i++){
			for(j=i;j<Nnodetypes;j++){
				for(l=0;l<symmetries;l++){
					val=radav[radialbin][i][j][l]/(1.*radcount[radialbin][i][j][l]);
					fprintf(rawoutp, "\t%.8g", val);
					Xrayval=val*scatterfactor[i]*scatterfactor[j];
					fprintf(weightedoutp, "\t%.8g", Xrayval);
					total[l]+=Xrayval;
					neutronval=val*neutronscatterfactor[i]*neutronscatterfactor[j];
					fprintf(neutronweightedoutp, "\t%.8g", neutronval);
					neutrontotal[l]+=neutronval;
				}
			}
		}
		fprintf(rawoutp, "\n");
		fprintf(weightedoutp, "\n");
		for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
		fprintf(weightedtotaloutp, "\n");
		fprintf(neutronweightedoutp, "\n");
		for(l=0;l<symmetries;l++) fprintf(neutronweightedtotaloutp, "\t%.8g", neutrontotal[l]);
		fprintf(neutronweightedtotaloutp, "\n");
	}	
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    fclose(neutronweightedoutp);
    fclose(neutronweightedtotaloutp);
	free(rawfile);
	free(weightedfile);
	free(weightedtotalfile);
	free(neutronweightedfile);
	free(neutronweightedtotalfile);
    free(total);
    free(neutrontotal);
	free_4d_tensor(radav, maxradbin, Nnodetypes, Nnodetypes);
	free_4d_tensor(radcount, maxradbin, Nnodetypes, Nnodetypes);
}

void output_1d_xrd_leaves_resolutionlimited(char *base, double ****Iqhisto1d, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, double *neutronscatterfactor, int symmetries, double_triple avbox){
    
    if(xrdqspace.z==0) xrdqspace.z=2*M_PI/avbox.z*xrdcount;
    
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
	declare_array_nozero(char, neutronweightedfile, charlength);
	declare_array_nozero(char, neutronweightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp, *neutronweightedoutp, *neutronweightedtotaloutp;
    int nz, i, j, l;
    double val, Xrayval, neutronval;
    declare_array(double, total, symmetries);
    declare_array(double, neutrontotal, symmetries);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
    sprintf(neutronweightedfile, "%s.neutron.weighted", base);
    sprintf(neutronweightedtotalfile, "%s.neutron.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
	neutronweightedoutp=fopen(neutronweightedfile, "w");
	neutronweightedtotaloutp=fopen(neutronweightedtotalfile, "w");
    for(nz=0;nz<=xrdqn.z;nz++){
        fprintf(rawoutp, "%.4f", xrdqspace.z*nz);
        fprintf(weightedoutp, "%.4f", xrdqspace.z*nz);
        fprintf(weightedtotaloutp, "%.4f", xrdqspace.z*nz);
        fprintf(neutronweightedoutp, "%.4f", xrdqspace.z*nz);
        fprintf(neutronweightedtotaloutp, "%.4f", xrdqspace.z*nz);
		for(l=0;l<symmetries;l++){
			total[l]=0;
			neutrontotal[l]=0;
		}
        for(i=0;i<Nnodetypes;i++){
            for(j=i;j<Nnodetypes;j++){
                for(l=0;l<symmetries;l++){
                    val=Iqhisto1d[nz][i][j][l]/(1.*xrdcount*Nresidues);
                    fprintf(rawoutp, "\t%.8g", val);
                    Xrayval=val*scatterfactor[i]*scatterfactor[j];
                    fprintf(weightedoutp, "\t%.8g", Xrayval);
                    total[l]+=Xrayval;
                    neutronval=val*neutronscatterfactor[i]*neutronscatterfactor[j];
                    fprintf(neutronweightedoutp, "\t%.8g", neutronval);
                    neutrontotal[l]+=neutronval;
                }
            }
        }
        fprintf(rawoutp, "\n");
        fprintf(weightedoutp, "\n");
        for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
        fprintf(weightedtotaloutp, "\n");
        fprintf(neutronweightedoutp, "\n");
        for(l=0;l<symmetries;l++) fprintf(neutronweightedtotaloutp, "\t%.8g", neutrontotal[l]);
        fprintf(neutronweightedtotaloutp, "\n");
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    fclose(neutronweightedoutp);
    fclose(neutronweightedtotaloutp);
	free(rawfile);
	free(weightedfile);
	free(weightedtotalfile);
	free(neutronweightedfile);
	free(neutronweightedtotalfile);
    free(total);
    free(neutrontotal);
}

void output_stack_xrd_leaves_resolutionlimited(char *base, double ****Iqhisto1d, int xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, double *neutronscatterfactor, int symmetries, double stackspace){
    
    double xrdqspace=2*M_PI/stackspace;
    
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
	declare_array_nozero(char, neutronweightedfile, charlength);
	declare_array_nozero(char, neutronweightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp, *neutronweightedoutp, *neutronweightedtotaloutp;
    int nz, i, j, l;
    double val, Xrayval, neutronval;
    declare_array(double, total, symmetries);
    declare_array(double, neutrontotal, symmetries);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
    sprintf(neutronweightedfile, "%s.neutron.weighted", base);
    sprintf(neutronweightedtotalfile, "%s.neutron.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
	neutronweightedoutp=fopen(neutronweightedfile, "w");
	neutronweightedtotaloutp=fopen(neutronweightedtotalfile, "w");
    for(nz=0;nz<=xrdqn;nz++){
        fprintf(rawoutp, "%.4f", xrdqspace*nz);
        fprintf(weightedoutp, "%.4f", xrdqspace*nz);
        fprintf(weightedtotaloutp, "%.4f", xrdqspace*nz);
        fprintf(neutronweightedoutp, "%.4f", xrdqspace*nz);
        fprintf(neutronweightedtotaloutp, "%.4f", xrdqspace*nz);
		for(l=0;l<symmetries;l++){
			total[l]=0;
			neutrontotal[l]=0;
		}
        for(i=0;i<Nnodetypes;i++){
            for(j=i;j<Nnodetypes;j++){
                for(l=0;l<symmetries;l++){
                    val=Iqhisto1d[nz][i][j][l]/(1.*xrdcount*Nresidues);
                    fprintf(rawoutp, "\t%.8g", val);
                    Xrayval=val*scatterfactor[i]*scatterfactor[j];
                    fprintf(weightedoutp, "\t%.8g", Xrayval);
                    total[l]+=Xrayval;
                    neutronval=val*neutronscatterfactor[i]*neutronscatterfactor[j];
                    fprintf(neutronweightedoutp, "\t%.8g", neutronval);
                    neutrontotal[l]+=neutronval;
                }
            }
        }
        fprintf(rawoutp, "\n");
        fprintf(weightedoutp, "\n");
        for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
        fprintf(weightedtotaloutp, "\n");
        fprintf(neutronweightedoutp, "\n");
        for(l=0;l<symmetries;l++) fprintf(neutronweightedtotaloutp, "\t%.8g", neutrontotal[l]);
        fprintf(neutronweightedtotaloutp, "\n");
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    fclose(neutronweightedoutp);
    fclose(neutronweightedtotaloutp);
	free(rawfile);
	free(weightedfile);
	free(weightedtotalfile);
	free(neutronweightedfile);
	free(neutronweightedtotalfile);
    free(total);
    free(neutrontotal);
}

void output_2d_xrd_leaves(char *base, double *****Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int **xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension){
	if(xrdqspace==0) xrdqspace=2*M_PI/box_dimension.x;
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int nr, nz, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(nr=0;nr<=xrdqn;nr++){
		for(nz=0;nz<=xrdqn;nz++){
			fprintf(rawoutp, "%.4f\t%.4f", xrdqspace*nr, xrdqspace*nz);
			fprintf(weightedoutp, "%.4f\t%.4f", xrdqspace*nr, xrdqspace*nz);
			fprintf(weightedtotaloutp, "%.4f\t%.4f", xrdqspace*nr, xrdqspace*nz);
            total[0]=total[1]=0;
			for(i=0;i<Nnodetypes;i++){
				for(j=i;j<Nnodetypes;j++){
                    for(l=0;l<symmetries;l++){
                        val=Iqhisto[nr][nz][i][j][l]/(1.*xrdcount[nr][nz]*Nresidues);
                        fprintf(rawoutp, "\t%.8g", val);
                        val=val*scatterfactor[i]*scatterfactor[j];
                        fprintf(weightedoutp, "\t%.8g", val);
                        total[l]+=val;
                    }
                }
			}
			fprintf(rawoutp, "\n");
			fprintf(weightedoutp, "\n");
			for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
			fprintf(weightedtotaloutp, "\n");
		}
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

void output_1d_xrd_from2d_z0_leaves(char *base, double *****Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int **xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension){
	if(xrdqspace==0) xrdqspace=2*M_PI/box_dimension.x;
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int nr, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(nr=0;nr<=xrdqn;nr++){
		fprintf(rawoutp, "%.4f", xrdqspace*nr);
		fprintf(weightedoutp, "%.4f", xrdqspace*nr);
		fprintf(weightedtotaloutp, "%.4f", xrdqspace*nr);
        total[0]=total[1]=0;
		for(i=0;i<Nnodetypes;i++){
			for(j=i;j<Nnodetypes;j++){
                for(l=0;l<symmetries;l++){
                    val=Iqhisto[nr][0][i][j][l]/(1.*xrdcount[nr][0]*Nresidues);
                    fprintf(rawoutp, "\t%.8g", val);
                    val=val*scatterfactor[i]*scatterfactor[j];
                    fprintf(weightedoutp, "\t%.8g", val);
                    total[l]+=val;
                }
            }
		}
		fprintf(rawoutp, "\n");
		fprintf(weightedoutp, "\n");
		for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
		fprintf(weightedtotaloutp, "\n");
    }
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

void output_1d_xrd_leaves(char *base, double ****Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int *xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension){
	if(xrdqspace==0) xrdqspace=2*M_PI/box_dimension.x;
	declare_array_nozero(char, rawfile, charlength);
	declare_array_nozero(char, weightedfile, charlength);
	declare_array_nozero(char, weightedtotalfile, charlength);
    FILE *rawoutp, *weightedoutp, *weightedtotaloutp;
    int n, i, j, l;
    double val;
    declare_array(double, total, 2);
    sprintf(rawfile, "%s.raw", base);
    sprintf(weightedfile, "%s.weighted", base);
    sprintf(weightedtotalfile, "%s.total", base);
	rawoutp=fopen(rawfile, "w");
	weightedoutp=fopen(weightedfile, "w");
	weightedtotaloutp=fopen(weightedtotalfile, "w");
    for(n=0;n<=xrdqn;n++){
		fprintf(rawoutp, "%.4f", xrdqspace*n);
		fprintf(weightedoutp, "%.4f", xrdqspace*n);
		fprintf(weightedtotaloutp, "%.4f", xrdqspace*n);
        total[0]=total[1]=0;
		for(i=0;i<Nnodetypes;i++){
			for(j=i;j<Nnodetypes;j++){
                for(l=0;l<symmetries;l++){
                    val=Iqhisto[n][i][j][l]/(1.*xrdcount[n]*Nresidues);
                    fprintf(rawoutp, "\t%.8g", val);
                    val=val*scatterfactor[i]*scatterfactor[j];
                    fprintf(weightedoutp, "\t%.8g", val);
                    total[l]+=val;
                }
			}
		}
		fprintf(rawoutp, "\n");
		fprintf(weightedoutp, "\n");
		for(l=0;l<symmetries;l++) fprintf(weightedtotaloutp, "\t%.8g", total[l]);
		fprintf(weightedtotaloutp, "\n");
	}
    fclose(rawoutp);
    fclose(weightedoutp);
    fclose(weightedtotaloutp);
    free(total);
}

int count_allatom_atoms_atomid(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, int *atomid, int *leaves){
	int i, j, total=0, sidechainindex, leaf;
	
	//  all atom groups
	
	//	0	back N
	//	1	back C
	//	2	back O
	//	3	back H
	//  4   first ethyl C
	//  5   first ethyl H
	//  6   phenyl second ethyl C
	//  7   phenyl second ethyl H
	//  8   phenyl aromatic C
	//  9   phenyl aromatic H
	//  10  amino second ethyl C
	//  11  amino second ethyl H
	//  12  amino N
	//  13  amino amino H
	//  14  carboxyl second ethyl C
	//  15  carboxyl second ethyl H
	//  16  carboxyl carboxyl C
	//  17  carboxyl O

	//	running through in same order as output_xyz_allatom and poormansxrd_fromcoord_bin_leafid_allatom
	
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength[i];j++){
			sidechainindex=monomerid[i][j].sidechain;
			leaf=coordarray[sidechainindex].leafid;
			atomid[total]=0;					//	backbone N
			leaves[total]=leaf;
			total++;
			if(coordarray[sidechainindex].nodetype>=0){
				atomid[total]=4;				//	Cbeta
				leaves[total]=leaf;
			}
			else{
				atomid[total]=3;				//	right terminal H (long terminus)
				leaves[total]=leaf;
			}
			total++;
            if(j==0){
				atomid[total]=3;				//	left terminal H
				leaves[total]=leaf;
				total++;
			}
			else{
				atomid[total]=1;				//	Cnyl
				atomid[total+1]=2;				//	O
				leaves[total]=leaves[total+1]=leaf;
				total+=2;
			}
			if(j==chainlength[i]-1){
				atomid[total]=3;				//	right terminal H
				leaves[total]=leaf;
				total++;
			}
			else{
				atomid[total]=1;				//	C
				atomid[total+1]=3;				//	H
				atomid[total+2]=3;				//	H
				leaves[total]=leaves[total+1]=leaves[total+2]=leaf;
				total+=3;
			}
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				if(coordarray[sidechainindex].nodetype==1){					//	phenyl
					atomid[total]=8;										//	aromatic C
					atomid[total+1]=8;										//	aromatic C
					atomid[total+2]=8;										//	aromatic C
					atomid[total+3]=8;										//	aromatic C
					atomid[total+4]=8;										//	aromatic C
					atomid[total+5]=8;										//	aromatic C
					atomid[total+6]=6;										//	phenyl Cgamma
					atomid[total+7]=9;										//	aromatic H
					atomid[total+8]=9;										//	aromatic H
					atomid[total+9]=9;										//	aromatic H
					atomid[total+10]=9;										//	aromatic H
					atomid[total+11]=9;										//	aromatic H
					leaves[total]=leaves[total+1]=leaves[total+2]=leaves[total+3]=leaves[total+4]=leaves[total+5]=leaves[total+6]=leaves[total+7]=leaves[total+8]=leaves[total+9]=leaves[total+10]=leaves[total+11]=leaf;
					total+=12;
				}
			}
			if(coordarray[sidechainindex].orientationtype==0){				//	parallel
				if(coordarray[sidechainindex].nodetype==2){					//	amino
					atomid[total]=12;										//	amino N
					atomid[total+1]=10;										//	amino Cgamma
					atomid[total+2]=13;										//	amino H
					atomid[total+3]=13;										//	amino H
					atomid[total+4]=13;										//	amino H
					leaves[total]=leaves[total+1]=leaves[total+2]=leaves[total+3]=leaves[total+4]=leaf;
					total+=5;				
				}
				if(coordarray[sidechainindex].nodetype==3){					//	carboxyl
					atomid[total]=14;										//	carboxyl second ethyl C
					atomid[total+1]=16;										//	carboxyl carboxyl C
					atomid[total+2]=17;										//	carboxyl O
					atomid[total+3]=17;										//	carboxyl O
					leaves[total]=leaves[total+1]=leaves[total+2]=leaves[total+3]=leaf;
					total+=4;
				}
				if(coordarray[sidechainindex].nodetype==4){ 
					total+=1;				//	methyl
					my_exit("need to write code for count_allatom_atoms_atomid methyl");
				}
			}
			if(coordarray[sidechainindex].nodetype>=0){
				atomid[total]=5;											//	Cbeta H
				atomid[total+1]=5;											//	Cbeta H
				leaves[total]=leaves[total+1]=leaf;
				total+=2;          
				if(coordarray[sidechainindex].nodetype==1){
					atomid[total]=7;										//	phenyl Cgamma H
					atomid[total+1]=7;										//	phenyl Cgamma H
				}
				else if(coordarray[sidechainindex].nodetype==2){
					atomid[total]=11;										//	amino Cgamma H
					atomid[total+1]=11;										//	amino Cgamma H
				}
				else if(coordarray[sidechainindex].nodetype==3){
					atomid[total]=15;										//	carboxyl Cgamma H
					atomid[total+1]=15;										//	carboxyl Cgamma H
				}
				leaves[total]=leaves[total+1]=leaf;
				total+=2;          
			}
		}
	}
	return total;
}

void poormansxrd_fromcoord_bin_leafid_allatom(int Nchains, coord *coordarray, double_triple box_dimension, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int **edgeoncount, int *powdercount, int *radiallyaveragededgeoncount, double *****Iqedgeonhisto, double ****Iqpowderhisto, double ****Iqradiallyaveragededgeonhisto, int allatoms, int *atomid, int *leaves, double_triple *position, int *chainlength, monomernodes **monomerid, allatomparams params, int Nleaves){
	int i, j, k, l, x, y, z, inplaneradialbin, isotropicbin, sym, backboneindex, sidechainindex, atomcount=0;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, Nleaves);
	declare_matrix(double, Aqimag, Nnodetypes, Nleaves);
	double_triple resolution, q;
    resolution.x=2*M_PI/box_dimension.x;
    resolution.y=2*M_PI/box_dimension.y;
    resolution.z=2*M_PI/box_dimension.z;
	if(xrdqspace==0){
		xrdqspace=resolution.x;
		if((resolution.x!=resolution.y)||(resolution.x!=resolution.z)) my_exit("trying to use resolution-limited xrd, but box is not cubic");
	}
    int_triple range, bin;
    range.x=(int) floor((xrdqn+0.5)*xrdqspace/resolution.x);
    range.y=(int) floor((xrdqn+0.5)*xrdqspace/resolution.y);
    range.z=(int) floor((xrdqn+0.5)*xrdqspace/resolution.z);
	double_triple sidevector, director, thirdvector, loc, backbonevector, backbonebonded, sidechainbonded, lastbackbonebonded, Cbetapos, Cgammapos, ethylvec, ethylplane, ethyloutofplane;
	
	//	assigning position, running through in same order as output_xyz_allatom and poormansxrd_fromcoord_bin_leafid_allatom
	
	for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength[i];j++){
			backboneindex=monomerid[i][j].backbone;
			sidechainindex=monomerid[i][j].sidechain;
			
			//	don't do center of mass calculation, but do use positions relative to bonded neighbors (don't worry about boundary conditions because XRD is independent of them)
			
			if(j==0) lastbackbonebonded=coordarray[backboneindex].r;
            backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			
            director=coordarray[backboneindex].n;
            if(j>0) backbonevector=subtract_double_triple(coordarray[backboneindex].r, coordarray[monomerid[i][j-1].backbone].r);
            else backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[backboneindex].r);
            
            recenter_double_triple(&backbonevector, box_dimension);
            normalize(&backbonevector);
            thirdvector=cross_product(backbonevector, director);
            normalize(&thirdvector);                                                            //  this is normal to the plane where I'll put all the atoms
            sidevector=cross_product(director, thirdvector);
            normalize(&sidevector);
            
			lastbackbonebonded=backbonebonded;
            
			position[atomcount]=backbonebonded;
			atomcount++;
            if(coordarray[sidechainindex].nodetype>=0){                                     
                loc=scalar_multiply_double_triple(director, params.Cbetaspacing);
                Cbetapos=loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            else{
                loc=scalar_multiply_double_triple(director, params.CterminuslongHdistance);
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            if(j==0){
                loc=add_double_triple(scalar_multiply_double_triple(director, params.NterminusHdepth), scalar_multiply_double_triple(sidevector, params.NterminusHwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            else{
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cnyldepth), scalar_multiply_double_triple(sidevector, params.Cnylwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Onyldepth), scalar_multiply_double_triple(sidevector, params.Onylwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            if(j==chainlength[i]-1){
                if(coordarray[sidechainindex].nodetype>=0){
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminusHdepth), scalar_multiply_double_triple(sidevector, params.CterminusHwidth));
                    loc=add_double_triple(backbonebonded, loc);
					position[atomcount]=loc;
					atomcount++;
                }
                else{
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminuslongHdepth), scalar_multiply_double_triple(sidevector, params.CterminuslongHwidth));
                    loc=add_double_triple(backbonebonded, loc);
                    //outputcoords_centered(outp, "H", loc, box_dimension);                       //  terminal H
					position[atomcount]=loc;
					atomcount++;
                }
            }
            else{                                                                               //  C
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cdepth), scalar_multiply_double_triple(sidevector, params.Cwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, -params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }            
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(sidevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==1){					//	phenyl
					loc=scalar_multiply_double_triple(sidevector, -params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(sidevector, params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
                    loc=scalar_multiply_double_triple(sidevector, -params.phenylCgammaspacing);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(sidevector, params.phenylHspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
                }
			}
			if(coordarray[sidechainindex].orientationtype==0){		//	parallel
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				if(j>0){
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
					else backbonevector=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
				}
				else{
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j].backbone].r);
					else{
						backbonevector.x=1;	//arbitrary
						backbonevector.y=backbonevector.z=0;
					}
				}
				recenter_double_triple(&backbonevector, box_dimension);
				normalize(&backbonevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(backbonevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==2){					//	amino
					loc=scalar_multiply_double_triple(director, params.aminoNdepth);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(director, params.aminoCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
					loc=add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
				}
				if(coordarray[sidechainindex].nodetype==3){					//	carboxyl
					loc=scalar_multiply_double_triple(director, params.carboxylbackCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(director, params.carboxylforwardCdepth);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, -params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
				}
			}
            if(coordarray[sidechainindex].nodetype>=0){
                
                //  C beta ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cbetapos, backbonebonded)), normed(subtract_double_triple(Cbetapos, Cgammapos))));
                ethylplane=subtract_double_triple(Cgammapos, backbonebonded);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                
                //  C gamma ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cgammapos, Cbetapos)), normed(subtract_double_triple(Cgammapos, sidechainbonded))));
                ethylplane=subtract_double_triple(sidechainbonded, Cbetapos);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
            }
		}
	}
	
	for(x=-range.x;x<=range.x;x++){
        q.x=resolution.x*x;
        bin.x=(int) floor(q.x/xrdqspace+0.5);
        for(y=-range.y;y<=range.y;y++){
            q.y=resolution.y*y;
            bin.y=(int) floor(q.y/xrdqspace+0.5);
            inplaneradialbin=(int) floor(sqrt(q.x*q.x+q.y*q.y)/xrdqspace+0.5);
            for(z=0;z<=xrdqn;z++){
                q.z=xrdqspace*z;                //  prescribed spacing, no periodic boundary conditions
				isotropicbin=(int) floor(sqrt(q.x*q.x+q.y*q.y+q.z*q.z)/xrdqspace+0.5);
                bin.z=(int) floor(q.z/xrdqspace+0.5);
                xrdcount[xrdqn+bin.x][xrdqn+bin.y][bin.z]++;
                if(inplaneradialbin<=xrdqn) edgeoncount[inplaneradialbin][bin.z]++;
                if(isotropicbin<=xrdqn){
                    powdercount[isotropicbin]++;
                    radiallyaveragededgeoncount[isotropicbin]++;
                }
				
				for(i=0;i<Nnodetypes;i++) for(j=0;j<Nleaves;j++) Aqreal[i][j]=Aqimag[i][j]=0;
				
				for(i=0;i<allatoms;i++){
					qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
					Aqreal[atomid[i]][leaves[i]]+=cos(qdotr);
					Aqimag[atomid[i]][leaves[i]]+=sin(qdotr);
				}
				for(i=0;i<Nnodetypes;i++){
                    
                    for(k=0;k<Nleaves;k++){
                        for(l=0;l<Nleaves;l++){
                            if(k==l) sym=0;
                            else if(k/2==l/2) sym=1;
                            else sym=2;
                            val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                            Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][i][sym]+=val;
                            if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][i][sym]+=val;
                            if(isotropicbin<=xrdqn){
                                Iqpowderhisto[isotropicbin][i][i][sym]+=val;
                                if((x==0)&&(y==0)) Iqradiallyaveragededgeonhisto[isotropicbin][i][i][sym]+=val*(2./M_PI*sqrt((pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2)+q.z*q.z)/(pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2))));
                                else Iqradiallyaveragededgeonhisto[isotropicbin][i][i][sym]+=val*(2./M_PI*sqrt((q.x*q.x+q.y*q.y+q.z*q.z)/(q.x*q.x+q.y*q.y)));
                            }
                            for(j=i+1;j<Nnodetypes;j++){
                                val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                                Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][j][sym]+=val;
                                if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][j][sym]+=val;
                                if(isotropicbin<=xrdqn){
                                    Iqpowderhisto[isotropicbin][i][j][sym]+=val;
                                    if((x==0)&&(y==0)) Iqradiallyaveragededgeonhisto[isotropicbin][i][j][sym]+=val*(2./M_PI*sqrt((pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2)+q.z*q.z)/(pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2))));
                                    else Iqradiallyaveragededgeonhisto[isotropicbin][i][j][sym]+=val*(2./M_PI*sqrt((q.x*q.x+q.y*q.y+q.z*q.z)/(q.x*q.x+q.y*q.y)));
                                }
                            }
                        }
                    }
                    
				}
			}
		}
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
}

void poormansxrd_fromcoord_bin_leafid_allatom_no3d(int Nchains, coord *coordarray, double_triple box_dimension, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int **edgeoncount, double *****Iqedgeonhisto, int allatoms, int *atomid, int *leaves, double_triple *position, int *chainlength, monomernodes **monomerid, allatomparams params, int Nleaves){

	int i, j, k, l, x, y, z, inplaneradialbin, sym, backboneindex, sidechainindex, atomcount=0;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, Nleaves);
	declare_matrix(double, Aqimag, Nnodetypes, Nleaves);
	double_triple resolution, q;
    resolution.x=2*M_PI/box_dimension.x;
    resolution.y=2*M_PI/box_dimension.y;
    resolution.z=2*M_PI/box_dimension.z;
	if(xrdqspace==0){
		xrdqspace=resolution.x;
		if((resolution.x!=resolution.y)||(resolution.x!=resolution.z)) my_exit("trying to use resolution-limited xrd, but box is not cubic");
	}
    int_triple range, bin;
    range.x=(int) floor((xrdqn+0.5)*xrdqspace/resolution.x);
    range.y=(int) floor((xrdqn+0.5)*xrdqspace/resolution.y);
    range.z=(int) floor((xrdqn+0.5)*xrdqspace/resolution.z);
	double_triple sidevector, director, thirdvector, loc, backbonevector, backbonebonded, sidechainbonded, lastbackbonebonded, Cbetapos, Cgammapos, ethylvec, ethylplane, ethyloutofplane;
	
	for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength[i];j++){
			backboneindex=monomerid[i][j].backbone;
			sidechainindex=monomerid[i][j].sidechain;
			
			//	don't do center of mass calculation, but do use positions relative to bonded neighbors (don't worry about boundary conditions because XRD is independent of them)
			
			if(j==0) lastbackbonebonded=coordarray[backboneindex].r;
            backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			
            director=coordarray[backboneindex].n;
            if(j>0) backbonevector=subtract_double_triple(coordarray[backboneindex].r, coordarray[monomerid[i][j-1].backbone].r);
            else backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[backboneindex].r);
            
            recenter_double_triple(&backbonevector, box_dimension);
            normalize(&backbonevector);
            thirdvector=cross_product(backbonevector, director);
            normalize(&thirdvector);                                                            //  this is normal to the plane where I'll put all the atoms
            sidevector=cross_product(director, thirdvector);
            normalize(&sidevector);
            
			lastbackbonebonded=backbonebonded;
            
			position[atomcount]=backbonebonded;
			atomcount++;
            if(coordarray[sidechainindex].nodetype>=0){
                loc=scalar_multiply_double_triple(director, params.Cbetaspacing);
                Cbetapos=loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            else{
                loc=scalar_multiply_double_triple(director, params.CterminuslongHdistance);
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            if(j==0){
                loc=add_double_triple(scalar_multiply_double_triple(director, params.NterminusHdepth), scalar_multiply_double_triple(sidevector, params.NterminusHwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            else{
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cnyldepth), scalar_multiply_double_triple(sidevector, params.Cnylwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Onyldepth), scalar_multiply_double_triple(sidevector, params.Onylwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            if(j==chainlength[i]-1){
                if(coordarray[sidechainindex].nodetype>=0){
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminusHdepth), scalar_multiply_double_triple(sidevector, params.CterminusHwidth));
                    loc=add_double_triple(backbonebonded, loc);
					position[atomcount]=loc;
					atomcount++;
                }
                else{
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminuslongHdepth), scalar_multiply_double_triple(sidevector, params.CterminuslongHwidth));
                    loc=add_double_triple(backbonebonded, loc);
					position[atomcount]=loc;
					atomcount++;
                }
            }
            else{                                                                               //  C
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cdepth), scalar_multiply_double_triple(sidevector, params.Cwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, -params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(sidevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==1){					//	phenyl
					loc=scalar_multiply_double_triple(sidevector, -params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(sidevector, params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
                    loc=scalar_multiply_double_triple(sidevector, -params.phenylCgammaspacing);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(sidevector, params.phenylHspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
                }
			}
			if(coordarray[sidechainindex].orientationtype==0){		//	parallel
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				if(j>0){
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
					else backbonevector=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
				}
				else{
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j].backbone].r);
					else{
						backbonevector.x=1;	//arbitrary
						backbonevector.y=backbonevector.z=0;
					}
				}
				recenter_double_triple(&backbonevector, box_dimension);
				normalize(&backbonevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(backbonevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==2){					//	amino
					loc=scalar_multiply_double_triple(director, params.aminoNdepth);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(director, params.aminoCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
					loc=add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
				}
				if(coordarray[sidechainindex].nodetype==3){					//	carboxyl
					loc=scalar_multiply_double_triple(director, params.carboxylbackCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(director, params.carboxylforwardCdepth);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, -params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
				}
			}
            if(coordarray[sidechainindex].nodetype>=0){
                
                //  C beta ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cbetapos, backbonebonded)), normed(subtract_double_triple(Cbetapos, Cgammapos))));
                ethylplane=subtract_double_triple(Cgammapos, backbonebonded);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                
                //  C gamma ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cgammapos, Cbetapos)), normed(subtract_double_triple(Cgammapos, sidechainbonded))));
                ethylplane=subtract_double_triple(sidechainbonded, Cbetapos);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
            }
		}
	}
	
    //  xy part
    
	for(x=-range.x;x<=range.x;x++){
        q.x=resolution.x*x;
        bin.x=(int) floor(q.x/xrdqspace+0.5);
        for(y=-range.y;y<=range.y;y++){
            q.y=resolution.y*y;
            bin.y=(int) floor(q.y/xrdqspace+0.5);
            inplaneradialbin=(int) floor(sqrt(q.x*q.x+q.y*q.y)/xrdqspace+0.5);
            z=0;
            q.z=xrdqspace*z;                //  prescribed spacing, no periodic boundary conditions
            bin.z=(int) floor(q.z/xrdqspace+0.5);
            xrdcount[xrdqn+bin.x][xrdqn+bin.y][bin.z]++;
            if(inplaneradialbin<=xrdqn) edgeoncount[inplaneradialbin][bin.z]++;
            
            for(i=0;i<Nnodetypes;i++) for(j=0;j<Nleaves;j++) Aqreal[i][j]=Aqimag[i][j]=0;
            
            for(i=0;i<allatoms;i++){
                qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
                Aqreal[atomid[i]][leaves[i]]+=cos(qdotr);
                Aqimag[atomid[i]][leaves[i]]+=sin(qdotr);
            }
            for(i=0;i<Nnodetypes;i++){
                
                for(k=0;k<Nleaves;k++){
                    for(l=0;l<Nleaves;l++){
                        if(k==l) sym=0;
                        else if(k/2==l/2) sym=1;
                        else sym=2;
                        val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                        Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][i][sym]+=val;
                        if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][i][sym]+=val;
                        for(j=i+1;j<Nnodetypes;j++){
                            val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                            Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][j][sym]+=val;
                            if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][j][sym]+=val;
                        }
                    }
                }
			}
		}
	}
    
    //  z
    
    x=0;
    y=0;
    q.x=resolution.x*x;
    bin.x=(int) floor(q.x/xrdqspace+0.5);
    q.y=resolution.y*y;
    bin.y=(int) floor(q.y/xrdqspace+0.5);
    inplaneradialbin=(int) floor(sqrt(q.x*q.x+q.y*q.y)/xrdqspace+0.5);
    for(z=0;z<=xrdqn;z++){
        q.z=xrdqspace*z;                //  prescribed spacing, no periodic boundary conditions
        bin.z=(int) floor(q.z/xrdqspace+0.5);
        xrdcount[xrdqn+bin.x][xrdqn+bin.y][bin.z]++;
        if(inplaneradialbin<=xrdqn) edgeoncount[inplaneradialbin][bin.z]++;
        
        for(i=0;i<Nnodetypes;i++) for(j=0;j<Nleaves;j++) Aqreal[i][j]=Aqimag[i][j]=0;
        
        for(i=0;i<allatoms;i++){
            qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
            Aqreal[atomid[i]][leaves[i]]+=cos(qdotr);
            Aqimag[atomid[i]][leaves[i]]+=sin(qdotr);
        }
        for(i=0;i<Nnodetypes;i++){
            
            for(k=0;k<Nleaves;k++){
                for(l=0;l<Nleaves;l++){
                    if(k==l) sym=0;
                    else if(k/2==l/2) sym=1;
                    else sym=2;
                    val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                    Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][i][sym]+=val;
                    if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][i][sym]+=val;
                    for(j=i+1;j<Nnodetypes;j++){
                        val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                        Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][j][sym]+=val;
                        if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][j][sym]+=val;
                    }
                }
            }
            
        }
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
    
}

void poormansxrd_fromposition_bin_leaves_no3d_resolutionlimited(int N, double_triple *position, double_triple box_dimension, double *****Iqhistoxy, double ****Iqhistoz, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int *pxrdcount, int *id, int *leaves, double_triple *pavbox){
    
    (*pxrdcount)++;
	(*pavbox)=add_double_triple((*pavbox), box_dimension);

    int i, j, x, y, z, inplaneradialbin, isotropicbin, k, l, sym;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, 2);
	declare_matrix(double, Aqimag, Nnodetypes, 2);
	double_triple resolution, q;

    if(xrdqspace.x==0) resolution.x=2*M_PI/box_dimension.x;
	else resolution.x=xrdqspace.x;
    if(xrdqspace.y==0) resolution.y=2*M_PI/box_dimension.y;
	else resolution.y=xrdqspace.y;
    if(xrdqspace.z==0) resolution.z=2*M_PI/box_dimension.z;
	else resolution.z=xrdqspace.z;
    
    //  xy part
    
	for(x=-xrdqn.x;x<=xrdqn.x;x++){
        q.x=resolution.x*x;
        for(y=-xrdqn.y;y<=xrdqn.y;y++){
            q.y=resolution.y*y;
            z=0;
            q.z=resolution.z*z;                //  prescribed spacing, no periodic boundary conditions
            
            for(i=0;i<Nnodetypes;i++) for(j=0;j<2;j++) Aqreal[i][j]=Aqimag[i][j]=0;
            for(i=0;i<N;i++){
                qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
                Aqreal[id[i]][leaves[i]]+=cos(qdotr);
                Aqimag[id[i]][leaves[i]]+=sin(qdotr);
            }
            for(i=0;i<Nnodetypes;i++){
                
                for(k=0;k<2;k++){
                    for(l=0;l<2;l++){
                        if(k==l) sym=0;
                        else sym=1;
                        val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                        Iqhistoxy[xrdqn.x+x][xrdqn.y+y][i][i][sym]+=val;
                        for(j=i+1;j<Nnodetypes;j++){
                            val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                            Iqhistoxy[xrdqn.x+x][xrdqn.y+y][i][j][sym]+=val;
                        }
                    }
                }
                
            }
        }
    }

    //  z part
    
    x=0;
    y=0;
    q.x=resolution.x*x;
    q.y=resolution.y*y;
    for(z=0;z<=xrdqn.z;z++){
        q.z=resolution.z*z;                //  prescribed spacing, no periodic boundary conditions
        for(i=0;i<Nnodetypes;i++) for(j=0;j<2;j++) Aqreal[i][j]=Aqimag[i][j]=0;
        for(i=0;i<N;i++){
            qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
            Aqreal[id[i]][leaves[i]]+=cos(qdotr);
            Aqimag[id[i]][leaves[i]]+=sin(qdotr);
        }
        for(i=0;i<Nnodetypes;i++){
            for(k=0;k<2;k++){
                for(l=0;l<2;l++){
                    if(k==l) sym=0;
                    else if(k/2==l/2) sym=1;
                    else sym=2;
                    val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                    Iqhistoz[z][i][i][sym]+=val;
                    for(j=i+1;j<Nnodetypes;j++){
                        val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                        Iqhistoz[z][i][j][sym]+=val;
                    }
                }
            }
        }
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
}

void poormansxrd_fromposition_bin_leaves_stack_resolutionlimited(int N, double_triple *position, double_triple box_dimension, double ****Iqhistostack, int xrdqnstack, int Nnodetypes, int *pxrdcount, int *id, int *leaves, double stackspace){
    
    (*pxrdcount)++;
    
    int i, j, z, inplaneradialbin, isotropicbin, k, l, sym;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, 2);
	declare_matrix(double, Aqimag, Nnodetypes, 2);
	double resolution, q;
    
    resolution=2*M_PI/stackspace;
    
    for(z=0;z<=xrdqnstack;z++){
        q=resolution*z;                //  prescribed spacing, no periodic boundary conditions
        for(i=0;i<Nnodetypes;i++) for(j=0;j<2;j++) Aqreal[i][j]=Aqimag[i][j]=0;
        for(i=0;i<N;i++){
            qdotr=q*position[i].z;
            Aqreal[id[i]][leaves[i]]+=cos(qdotr);
            Aqimag[id[i]][leaves[i]]+=sin(qdotr);
        }
        for(i=0;i<Nnodetypes;i++){
            for(k=0;k<2;k++){
                for(l=0;l<2;l++){
                    if(k==l) sym=0;
                    else if(k/2==l/2) sym=1;
                    else sym=2;
                    val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                    Iqhistostack[z][i][i][sym]+=val;
                    for(j=i+1;j<Nnodetypes;j++){
                        val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                        Iqhistostack[z][i][j][sym]+=val;
                    }
                }
            }
        }
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
}

void poormansxrd_fromcoord_bin_leafid_allatom_no3d_resolutionlimited(int Nchains, coord *coordarray, double_triple box_dimension, double *****Iqhistoxy, double ****Iqhistoz, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int *pxrdcount, int allatoms, int *atomid, int *leaves, double_triple *position, int *chainlength, monomernodes **monomerid, allatomparams params, int Nleaves, double_triple *pavbox){

	(*pxrdcount)++;
	(*pavbox)=add_double_triple((*pavbox), box_dimension);
	
	int i, j, k, l, x, y, z, sym, backboneindex, sidechainindex, atomcount=0;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, Nleaves);
	declare_matrix(double, Aqimag, Nnodetypes, Nleaves);
	double_triple resolution, q;
	
    if(xrdqspace.x==0) resolution.x=2*M_PI/box_dimension.x;
	else resolution.x=xrdqspace.x;
    if(xrdqspace.y==0) resolution.y=2*M_PI/box_dimension.y;
	else resolution.y=xrdqspace.y;
    if(xrdqspace.z==0) resolution.z=2*M_PI/box_dimension.z;
	else resolution.z=xrdqspace.z;

	double_triple sidevector, director, thirdvector, loc, backbonevector, backbonebonded, sidechainbonded, lastbackbonebonded, Cbetapos, Cgammapos, ethylvec, ethylplane, ethyloutofplane;
	
	for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength[i];j++){
			backboneindex=monomerid[i][j].backbone;
			sidechainindex=monomerid[i][j].sidechain;
			
			//	don't do center of mass calculation, but do use positions relative to bonded neighbors (don't worry about boundary conditions because XRD is independent of them)
			
			if(j==0) lastbackbonebonded=coordarray[backboneindex].r;
            backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			
            director=coordarray[backboneindex].n;
            if(j>0) backbonevector=subtract_double_triple(coordarray[backboneindex].r, coordarray[monomerid[i][j-1].backbone].r);
            else backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[backboneindex].r);
            
            recenter_double_triple(&backbonevector, box_dimension);
            normalize(&backbonevector);
            thirdvector=cross_product(backbonevector, director);
            normalize(&thirdvector);                                                            //  this is normal to the plane where I'll put all the atoms
            sidevector=cross_product(director, thirdvector);
            normalize(&sidevector);
            
			lastbackbonebonded=backbonebonded;
            
			position[atomcount]=backbonebonded;
			atomcount++;
            if(coordarray[sidechainindex].nodetype>=0){
                loc=scalar_multiply_double_triple(director, params.Cbetaspacing);
                Cbetapos=loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            else{
                loc=scalar_multiply_double_triple(director, params.CterminuslongHdistance);
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            if(j==0){
                loc=add_double_triple(scalar_multiply_double_triple(director, params.NterminusHdepth), scalar_multiply_double_triple(sidevector, params.NterminusHwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            else{
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cnyldepth), scalar_multiply_double_triple(sidevector, params.Cnylwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Onyldepth), scalar_multiply_double_triple(sidevector, params.Onylwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
            if(j==chainlength[i]-1){
                if(coordarray[sidechainindex].nodetype>=0){
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminusHdepth), scalar_multiply_double_triple(sidevector, params.CterminusHwidth));
                    loc=add_double_triple(backbonebonded, loc);
					position[atomcount]=loc;
					atomcount++;
                }
                else{
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminuslongHdepth), scalar_multiply_double_triple(sidevector, params.CterminuslongHwidth));
                    loc=add_double_triple(backbonebonded, loc);
					position[atomcount]=loc;
					atomcount++;
                }
            }
            else{                                                                               //  C
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cdepth), scalar_multiply_double_triple(sidevector, params.Cwidth));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, -params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
				position[atomcount]=loc;
				atomcount++;
            }
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(sidevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==1){					//	phenyl
					loc=scalar_multiply_double_triple(sidevector, -params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(sidevector, params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
                    loc=scalar_multiply_double_triple(sidevector, -params.phenylCgammaspacing);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(sidevector, params.phenylHspacing);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
                }
			}
			if(coordarray[sidechainindex].orientationtype==0){		//	parallel
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				if(j>0){
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
					else backbonevector=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
				}
				else{
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j].backbone].r);
					else{
						backbonevector.x=1;	//arbitrary
						backbonevector.y=backbonevector.z=0;
					}
				}
				recenter_double_triple(&backbonevector, box_dimension);
				normalize(&backbonevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(backbonevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==2){					//	amino
					loc=scalar_multiply_double_triple(director, params.aminoNdepth);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(director, params.aminoCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
                    
					loc=add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
				}
				if(coordarray[sidechainindex].nodetype==3){					//	carboxyl
					loc=scalar_multiply_double_triple(director, params.carboxylbackCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=scalar_multiply_double_triple(director, params.carboxylforwardCdepth);
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, -params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					position[atomcount]=loc;
					atomcount++;
				}
			}
            if(coordarray[sidechainindex].nodetype>=0){
                
                //  C beta ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cbetapos, backbonebonded)), normed(subtract_double_triple(Cbetapos, Cgammapos))));
                ethylplane=subtract_double_triple(Cgammapos, backbonebonded);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                
                //  C gamma ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cgammapos, Cbetapos)), normed(subtract_double_triple(Cgammapos, sidechainbonded))));
                ethylplane=subtract_double_triple(sidechainbonded, Cbetapos);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
				position[atomcount]=loc;
				atomcount++;
            }
		}
	}
	
    //  xy part
    
	for(x=-xrdqn.x;x<=xrdqn.x;x++){
        q.x=resolution.x*x;
        for(y=-xrdqn.y;y<=xrdqn.y;y++){
            q.y=resolution.y*y;
            z=0;
            q.z=resolution.z*z;                //  prescribed spacing, no periodic boundary conditions
            for(i=0;i<Nnodetypes;i++) for(j=0;j<Nleaves;j++) Aqreal[i][j]=Aqimag[i][j]=0;
            for(i=0;i<allatoms;i++){
                qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
                Aqreal[atomid[i]][leaves[i]]+=cos(qdotr);
                Aqimag[atomid[i]][leaves[i]]+=sin(qdotr);
            }
            for(i=0;i<Nnodetypes;i++){                
                for(k=0;k<Nleaves;k++){
                    for(l=0;l<Nleaves;l++){
                        if(k==l) sym=0;
                        else if(k/2==l/2) sym=1;
                        else sym=2;
                        val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                        Iqhistoxy[xrdqn.x+x][xrdqn.y+y][i][i][sym]+=val;
                        for(j=i+1;j<Nnodetypes;j++){
                            val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                            Iqhistoxy[xrdqn.x+x][xrdqn.y+y][i][j][sym]+=val;
                        }
                    }
                }
			}
		}
	}
    
    //  z
    
    x=0;
    y=0;
    q.x=resolution.x*x;
    q.y=resolution.y*y;
    for(z=0;z<=xrdqn.z;z++){
        q.z=resolution.z*z;                //  prescribed spacing, no periodic boundary conditions
        for(i=0;i<Nnodetypes;i++) for(j=0;j<Nleaves;j++) Aqreal[i][j]=Aqimag[i][j]=0;
        for(i=0;i<allatoms;i++){
            qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
            Aqreal[atomid[i]][leaves[i]]+=cos(qdotr);
            Aqimag[atomid[i]][leaves[i]]+=sin(qdotr);
        }
        for(i=0;i<Nnodetypes;i++){            
            for(k=0;k<Nleaves;k++){
                for(l=0;l<Nleaves;l++){
                    if(k==l) sym=0;
                    else if(k/2==l/2) sym=1;
                    else sym=2;
                    val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                    Iqhistoz[z][i][i][sym]+=val;
                    for(j=i+1;j<Nnodetypes;j++){
                        val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                        Iqhistoz[z][i][j][sym]+=val;
                    }
                }
            }            
        }
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
}

void poormansxrd_fromposition_bin_leaves(int N, double_triple *position, double_triple box_dimension, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int *id, int **edgeoncount, int *powdercount, int *radiallyaveragededgeoncount, double *****Iqedgeonhisto, double ****Iqpowderhisto, double ****Iqradiallyaveragededgeonhisto, int *leaves){
    int i, j, x, y, z, inplaneradialbin, isotropicbin, k, l, sym;
    double qdotr, val;
	declare_matrix(double, Aqreal, Nnodetypes, 2);
	declare_matrix(double, Aqimag, Nnodetypes, 2);
	double_triple resolution, q;
    resolution.x=2*M_PI/box_dimension.x;
    resolution.y=2*M_PI/box_dimension.y;
    resolution.z=2*M_PI/box_dimension.z;
    int_triple range, bin;
    range.x=(int) floor((xrdqn+0.5)*xrdqspace/resolution.x);
    range.y=(int) floor((xrdqn+0.5)*xrdqspace/resolution.y);
    range.z=(int) floor((xrdqn+0.5)*xrdqspace/resolution.z);
    for(x=-range.x;x<=range.x;x++){
        q.x=resolution.x*x;
        bin.x=(int) floor(q.x/xrdqspace+0.5);
        for(y=-range.y;y<=range.y;y++){
            q.y=resolution.y*y;
            bin.y=(int) floor(q.y/xrdqspace+0.5);
            inplaneradialbin=(int) floor(sqrt(q.x*q.x+q.y*q.y)/xrdqspace+0.5);
            for(z=0;z<=xrdqn;z++){
                q.z=xrdqspace*z;                //  prescribed spacing, no periodic boundary conditions
				isotropicbin=(int) floor(sqrt(q.x*q.x+q.y*q.y+q.z*q.z)/xrdqspace+0.5);
                bin.z=(int) floor(q.z/xrdqspace+0.5);
                xrdcount[xrdqn+bin.x][xrdqn+bin.y][bin.z]++;
                if(inplaneradialbin<=xrdqn) edgeoncount[inplaneradialbin][bin.z]++;
                if(isotropicbin<=xrdqn){
                    powdercount[isotropicbin]++;
                    radiallyaveragededgeoncount[isotropicbin]++;
                }
				for(i=0;i<Nnodetypes;i++) for(j=0;j<2;j++) Aqreal[i][j]=Aqimag[i][j]=0;
				for(i=0;i<N;i++){
					qdotr=q.x*position[i].x+q.y*position[i].y+q.z*position[i].z;
					Aqreal[id[i]][leaves[i]]+=cos(qdotr);
					Aqimag[id[i]][leaves[i]]+=sin(qdotr);
				}
				for(i=0;i<Nnodetypes;i++){
                    
                    for(k=0;k<2;k++){
                        for(l=0;l<2;l++){
                            if(k==l) sym=0;
                            else sym=1;
                            val=(Aqreal[i][k]*Aqreal[i][l]+Aqimag[i][k]*Aqimag[i][l]);
                            Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][i][sym]+=val;
                            if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][i][sym]+=val;
                            if(isotropicbin<=xrdqn){
                                Iqpowderhisto[isotropicbin][i][i][sym]+=val;
                                if((x==0)&&(y==0)) Iqradiallyaveragededgeonhisto[isotropicbin][i][i][sym]+=val*(2./M_PI*sqrt((pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2)+q.z*q.z)/(pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2))));
                                else Iqradiallyaveragededgeonhisto[isotropicbin][i][i][sym]+=val*(2./M_PI*sqrt((q.x*q.x+q.y*q.y+q.z*q.z)/(q.x*q.x+q.y*q.y)));
                            }
                            for(j=i+1;j<Nnodetypes;j++){
                                val=2*(Aqreal[i][k]*Aqreal[j][l]+Aqimag[i][k]*Aqimag[j][l]);
                                Iqhisto[xrdqn+bin.x][xrdqn+bin.y][bin.z][i][j][sym]+=val;
                                if(inplaneradialbin<=xrdqn) Iqedgeonhisto[inplaneradialbin][bin.z][i][j][sym]+=val;
                                if(isotropicbin<=xrdqn){
                                    Iqpowderhisto[isotropicbin][i][j][sym]+=val;
                                    if((x==0)&&(y==0)) Iqradiallyaveragededgeonhisto[isotropicbin][i][j][sym]+=val*(2./M_PI*sqrt((pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2)+q.z*q.z)/(pow(0.25*resolution.x, 2)+pow(0.25*resolution.y, 2))));
                                    else Iqradiallyaveragededgeonhisto[isotropicbin][i][j][sym]+=val*(2./M_PI*sqrt((q.x*q.x+q.y*q.y+q.z*q.z)/(q.x*q.x+q.y*q.y)));
                                }
                            }
                        }
                    }
                    
				}
			}
		}
	}
	free_matrix(Aqreal, Nnodetypes);
	free_matrix(Aqimag, Nnodetypes);
}

