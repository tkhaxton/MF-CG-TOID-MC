#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mymath.h"
#include "mystdlib.h"
#include "peptoid_functions.h"
#include "search_functions.h"

void reflect(double *pvar){
	if(*pvar>1) (*pvar)=2-(*pvar);
	if(*pvar<-1) (*pvar)=-2-(*pvar);
}

void output_monolayer(char *configname, monolayer config, double energy, int frame){
    FILE *outp;
    if(frame==0) outp=fopen(configname, "w");
    else outp=fopen(configname, "a");
    fprintf(outp, "%i", frame);
    fprintf(outp, "\t%.8f", energy);
    fprintf(outp, "\t%.8f", config.monomerspacing);
    fprintf(outp, "\t%.8f", config.monomerscaleoffsetfraction);
    fprintf(outp, "\t%.8f", config.terminusspacing);
    fprintf(outp, "\t%.8f", config.interchainspacing);
    
    int charged, chargetype, last, self;
    for(charged=0;charged<2;charged++){
        for(chargetype=2;chargetype<4;chargetype++){
            if(charged==1){
                last=0;
                self=chargetype;
            }
            else{
                last=chargetype;
                self=1;
            }
            fprintf(outp, "\t%.8f", config.backboneoffset[last][self].x);
            fprintf(outp, "\t%.8f", config.backboneoffset[last][self].y);
            fprintf(outp, "\t%.8f", config.backboneoffset[last][self].z);
            fprintf(outp, "\t%.8f", config.backbonenz[last][self]);
            fprintf(outp, "\t%.8f", config.backbonephi[last][self]);
            fprintf(outp, "\t%.8f", config.sidechainrelativepos[last][self].x);
            fprintf(outp, "\t%.8f", config.sidechainrelativepos[last][self].y);
            fprintf(outp, "\t%.8f", config.sidechainrelativepos[last][self].z);
            fprintf(outp, "\t%.8f", config.sidechainnpar[last][self]);
            fprintf(outp, "\t%.8f", config.sidechainnphi[last][self]);
        }
    }
    fprintf(outp, "\n");
    fclose(outp);
}

void copymonolayer(monolayer new, monolayer *pchanging, int Nnodetypes){
	(*pchanging).monomerspacing=new.monomerspacing;
	(*pchanging).interchainspacing=new.interchainspacing;
	(*pchanging).terminusspacing=new.terminusspacing;
	(*pchanging).offset=new.offset;
	(*pchanging).boxheight=new.boxheight;
	(*pchanging).phenylcode=new.phenylcode;
	(*pchanging).offsetfraction=new.offsetfraction;
	(*pchanging).monomerscaleoffsetfraction=new.monomerscaleoffsetfraction;
	(*pchanging).length=new.length;
	int i, j;
	for(i=0;i<Nnodetypes;i++){
		for(j=0;j<Nnodetypes;j++){
			(*pchanging).backboneoffset[i][j]=new.backboneoffset[i][j];
			(*pchanging).backbonenz[i][j]=new.backbonenz[i][j];
			(*pchanging).backbonephi[i][j]=new.backbonephi[i][j];
			(*pchanging).sidechainrelativepos[i][j]=new.sidechainrelativepos[i][j];
			(*pchanging).sidechainnpar[i][j]=new.sidechainnpar[i][j];
			(*pchanging).sidechainnphi[i][j]=new.sidechainnphi[i][j];
		}
	}
}

void initialize_periodic_monolayer(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config){
    int i, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	(*pbox_dimension).x=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).y=Nchains*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
	for(i=0;i<Nchains;i++){
		position.x=0.5*(*pbox_dimension).x+(i%2)*config.offset;
		position.y=(i+0.5)*config.interchainspacing;
		position.z=(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness;
		for(k=0;k<chainlength[chaincounter];k++){
			backboneid=monomerid[chaincounter][k].backbone;
			sidechainid=monomerid[chaincounter][k].sidechain;
            if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                lastpolartype=0;
                selftype=lastpolar;
            }
			else if(coordarray[sidechainid].nodetype>1){
                lastpolar=coordarray[sidechainid].nodetype;
                lastpolartype=0;                                //  self is polar
                selftype=coordarray[sidechainid].nodetype;
            }
            else{
                lastpolartype=lastpolar;
                selftype=coordarray[sidechainid].nodetype;
            }
			coordarray[backboneid].r=position;
			coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
            coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[lastpolartype][selftype]);
            fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
            nz=config.backbonenz[lastpolartype][selftype];
 			coordarray[backboneid].n.z=((k%2)*2-1)*nz;          //  amphiphilic
            phi=config.backbonephi[lastpolartype][selftype];
			coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
			coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
            if(coordarray[sidechainid].nodetype>=0){
                nback=coordarray[backboneid].n;
                nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
                normalize(&nalong);
                nside=cross_product(nback, nalong);
                npar=config.sidechainnpar[lastpolartype][coordarray[sidechainid].nodetype];
                nphi=config.sidechainnphi[lastpolartype][coordarray[sidechainid].nodetype];
                sidepos=config.sidechainrelativepos[lastpolartype][coordarray[sidechainid].nodetype];
                coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
                fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
                nx=sqrt(1-pow(npar, 2))*cos(nphi);
                ny=sqrt(1-pow(npar, 2))*sin(nphi);
                 coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
            }
		}
		chaincounter++;
		if(chaincounter==Nchains) return;
	}
}

void copybilayer(bilayer new, bilayer *pchanging, int Nnodetypes){
	(*pchanging).monomerspacing=new.monomerspacing;
	(*pchanging).interchainspacing=new.interchainspacing;
	(*pchanging).terminusspacing=new.terminusspacing;
	(*pchanging).offset=new.offset;
	(*pchanging).boxheight=new.boxheight;
	(*pchanging).phenylcode=new.phenylcode;
	(*pchanging).offsetfraction=new.offsetfraction;
	(*pchanging).monomerscaleoffsetfraction=new.monomerscaleoffsetfraction;
	(*pchanging).length=new.length;
	int k, i, j;
	for(k=0;k<2;k++){
		for(i=0;i<Nnodetypes;i++){
			for(j=0;j<Nnodetypes;j++){
				(*pchanging).backboneoffset[k][i][j]=new.backboneoffset[k][i][j];
				(*pchanging).backbonenz[k][i][j]=new.backbonenz[k][i][j];
				(*pchanging).backbonephi[k][i][j]=new.backbonephi[k][i][j];
				(*pchanging).sidechainrelativepos[k][i][j]=new.sidechainrelativepos[k][i][j];
				(*pchanging).sidechainnpar[k][i][j]=new.sidechainnpar[k][i][j];
				(*pchanging).sidechainnphi[k][i][j]=new.sidechainnphi[k][i][j];
			}
		}
	}
}

void initialize_periodic_bilayer(int Nmonomers, int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry){
    int j, i, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	(*pbox_dimension).x=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).y=(Nchains/2)*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
    if(Nchains%2!=0) my_exit("Nchains not divisible by 2!");
    for(j=0;j<2;j++){
        for(i=0;i<Nchains/2;i++){
            position.x=0.5*(*pbox_dimension).x+(i%2)*config.offset;
            position.y=(i+0.5)*config.interchainspacing;
            position.z=0.5*(*pbox_dimension).z;
            for(k=0;k<chainlength[chaincounter];k++){
                backboneid=monomerid[chaincounter][k].backbone;
                
                sidechainid=monomerid[chaincounter][k].sidechain;
				
				if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
					lastpolartype=0;
					selftype=lastpolar;
				}
				else if(coordarray[sidechainid].nodetype>1){
					lastpolar=coordarray[sidechainid].nodetype;
					lastpolartype=0;                                //  self is polar
					selftype=coordarray[sidechainid].nodetype;
				}
				else{
					lastpolartype=lastpolar;
					selftype=coordarray[sidechainid].nodetype;
				}
				
                coordarray[backboneid].r=position;
				if((j==1)&&(leafsymmetry==-1)) coordarray[backboneid].r.x-=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
				else coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
                coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[j][lastpolartype][selftype]);
                fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
				nz=config.backbonenz[j][lastpolartype][selftype];
                coordarray[backboneid].n.z=(1-2*j)*((k%2)*2-1)*nz;          //  leaves facing and amphiphilic
                phi=config.backbonephi[j][lastpolartype][selftype];
                coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
                coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
                nback=coordarray[backboneid].n;
                nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
                normalize(&nalong);
                nside=scalar_multiply_double_triple(cross_product(nback, nalong), (1-2*j));
                npar=config.sidechainnpar[j][lastpolartype][selftype];
                nphi=config.sidechainnphi[j][lastpolartype][selftype];
				sidepos=config.sidechainrelativepos[j][lastpolartype][selftype];
                coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
                fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
                nx=sqrt(1-pow(npar, 2))*cos(nphi);
                ny=sqrt(1-pow(npar, 2))*sin(nphi);
                coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
            }
			chaincounter++;
			if(chaincounter==Nchains) return;
        }
    }
}

void output_bilayer(char *configname, bilayer config, double energy, int frame){
    FILE *outp;
    if(frame==0) outp=fopen(configname, "w");
    else outp=fopen(configname, "a");
    fprintf(outp, "%i", frame);
    fprintf(outp, "\t%.8f", energy);
    fprintf(outp, "\t%.8f", config.monomerspacing);
    fprintf(outp, "\t%.8f", config.monomerscaleoffsetfraction);
    fprintf(outp, "\t%.8f", config.terminusspacing);
    fprintf(outp, "\t%.8f", config.interchainspacing);
    
    int charged, chargetype, last, self, leaf;
	for(leaf=0;leaf<2;leaf++){
		for(charged=0;charged<2;charged++){
			for(chargetype=2;chargetype<4;chargetype++){
				if(charged==1){
					last=0;
					self=chargetype;
				}
				else{
					last=chargetype;
					self=1;
				}
				fprintf(outp, "\t%.8f", config.backboneoffset[leaf][last][self].x);
				fprintf(outp, "\t%.8f", config.backboneoffset[leaf][last][self].y);
				fprintf(outp, "\t%.8f", config.backboneoffset[leaf][last][self].z);
				fprintf(outp, "\t%.8f", config.backbonenz[leaf][last][self]);
				fprintf(outp, "\t%.8f", config.backbonephi[leaf][last][self]);
				fprintf(outp, "\t%.8f", config.sidechainrelativepos[leaf][last][self].x);
				fprintf(outp, "\t%.8f", config.sidechainrelativepos[leaf][last][self].y);
				fprintf(outp, "\t%.8f", config.sidechainrelativepos[leaf][last][self].z);
				fprintf(outp, "\t%.8f", config.sidechainnpar[leaf][last][self]);
				fprintf(outp, "\t%.8f", config.sidechainnphi[leaf][last][self]);
			}
		}
	}
	fprintf(outp, "\n");
    fclose(outp);
}

