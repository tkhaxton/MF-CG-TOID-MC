#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mymath.h"
#include "mystdlib.h"
#include "peptoid_functions.h"

void parse_chemistry(char **chemistry_char, int *pchaintypes, int **Nchainsoftype, int **chainlengthoftype, int ***chainchemistry, int *pNmonomers, int *pNchains, monomertypes *pmonomercount, int sidechaintypesmax, int *pNnodetypes, int terminuscode){
	(*pNnodetypes)=0;
	(*pNmonomers)=0;
	(*pNchains)=0;
	(*pchaintypes)=atoi(chemistry_char[0]);
	(*Nchainsoftype)=xcalloc((*pchaintypes), sizeof(int));
	(*chainlengthoftype)=xcalloc((*pchaintypes), sizeof(int));
	int i, j, k, patterntype, counter=1, nonpolartype, firstpolartype, secondpolartype, Nblocks, blockcount, maxblocks=10, *polartype, *blocklength, count, extrachainlength;
	polartype=xcalloc(maxblocks, sizeof(int));
	blocklength=xcalloc(maxblocks, sizeof(int));
	(*chainchemistry)=xcalloc((*pchaintypes), sizeof(int *));
	(*pmonomercount).charged=(*pmonomercount).nonpolar=(*pmonomercount).monomers=0;
	for(i=0;i<sidechaintypesmax;i++) (*pmonomercount).type[i]=0;
	for(i=0;i<(*pchaintypes);i++){
		patterntype=atoi(chemistry_char[counter]);
		counter++;
		(*Nchainsoftype)[i]=atoi(chemistry_char[counter]);
		(*pNchains)+=(*Nchainsoftype)[i];
		counter++;
		if(patterntype==0){										//	alternating
			(*chainlengthoftype)[i]=atoi(chemistry_char[counter]);
			counter++;
			nonpolartype=atoi(chemistry_char[counter]);
			if(nonpolartype>=sidechaintypesmax) my_exit("type too big!");
			if(nonpolartype+1>(*pNnodetypes)) (*pNnodetypes)=nonpolartype+1;
			counter++;
			firstpolartype=atoi(chemistry_char[counter]);
			if(firstpolartype>=sidechaintypesmax) my_exit("type too big!");
			if(firstpolartype+1>(*pNnodetypes)) (*pNnodetypes)=firstpolartype+1;
			counter++;
			secondpolartype=atoi(chemistry_char[counter]);
			if(secondpolartype+1>(*pNnodetypes)) (*pNnodetypes)=secondpolartype+1;
			if(secondpolartype>=sidechaintypesmax) my_exit("type too big!");
			counter++;
            if(terminuscode==1) extrachainlength=1;
            else extrachainlength=0;
            (*chainlengthoftype)[i]+=extrachainlength;
			(*chainchemistry)[i]=xcalloc((*chainlengthoftype)[i], sizeof(int));
			for(j=0;j<(*chainlengthoftype)[i]-extrachainlength;j++){
				if(j%4==0){
					(*chainchemistry)[i][j]=firstpolartype;				//	start with polar
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];						//	backbone site
					(*pmonomercount).type[firstpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];							//	backbone
				}
				if(j%2==1){
					(*chainchemistry)[i][j]=nonpolartype;
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=2*(*Nchainsoftype)[i];						//	one for sidechain, one for backbone
					(*pmonomercount).type[nonpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];							//	backbone
				}
				if(j%4==3){
					(*chainchemistry)[i][j]=secondpolartype;
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];						//	backbone site
					(*pmonomercount).type[secondpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];							//	backbone
				}
			}
			(*pNmonomers)+=((*chainlengthoftype)[i])*((*Nchainsoftype)[i]);
            if(terminuscode==1){
                (*chainchemistry)[i][j]=-1;                              //  chainchemistry[i][j]=-1 denotes backbone site with no sidechain site
            }
		}
		if(patterntype==1){																//	block
			Nblocks=atoi(chemistry_char[counter]);
			if(Nblocks>maxblocks) my_exit("Nblocks too big!");
			counter++;
			nonpolartype=atoi(chemistry_char[counter]);
			if(nonpolartype>=sidechaintypesmax) my_exit("type too big!");
			if(nonpolartype+1>(*pNnodetypes)) (*pNnodetypes)=nonpolartype+1;
			counter++;
			blockcount=0;
            (*chainlengthoftype)[i]=0;
			while(blockcount<Nblocks){
				polartype[blockcount]=atoi(chemistry_char[counter]);
				if(polartype[blockcount]>=sidechaintypesmax) my_exit("type too big!");
				if(polartype[blockcount]+1>(*pNnodetypes)) (*pNnodetypes)=polartype[blockcount]+1;
				counter++;
				blocklength[blockcount]=atoi(chemistry_char[counter]);			//	 in dimers
				(*chainlengthoftype)[i]+=2*blocklength[blockcount];
				counter++;
				blockcount++;
			}
			count=0;
            
            if(terminuscode==1) extrachainlength=1;
            else extrachainlength=0;
            (*chainlengthoftype)[i]+=extrachainlength;
            
			(*chainchemistry)[i]=xcalloc((*chainlengthoftype)[i], sizeof(int));
			for(j=0;j<blockcount;j++){
				for(k=0;k<blocklength[j];k++){
					(*chainchemistry)[i][count]=polartype[j];				//	start with polar
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];							//	backbone site
					(*pmonomercount).type[polartype[j]]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];								//	backbone
					count++;
					(*chainchemistry)[i][count]=nonpolartype;
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=2*(*Nchainsoftype)[i];							//	one for sidechain, one for backbone
					(*pmonomercount).type[nonpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];								//	backbone
					count++;
				}
			}
			(*pNmonomers)+=((*chainlengthoftype)[i])*((*Nchainsoftype)[i]);
            if(terminuscode==1){
                (*chainchemistry)[i][count]=-1;                              //  chainchemistry[i][j]=-1 denotes backbone site with no sidechain site
            }
		}
		if(patterntype==2){																//	block, not necessarily alternating hydrophobic-hydrophilic
			Nblocks=atoi(chemistry_char[counter]);
			if(Nblocks>maxblocks) my_exit("Nblocks too big!");
			counter++;
			blockcount=0;
			(*chainlengthoftype)[i]=0;
			while(blockcount<Nblocks){
				polartype[blockcount]=atoi(chemistry_char[counter]);
				if(polartype[blockcount]>=sidechaintypesmax) my_exit("type too big!");
				if(polartype[blockcount]+1>(*pNnodetypes)) (*pNnodetypes)=polartype[blockcount]+1;
				counter++;
				blocklength[blockcount]=atoi(chemistry_char[counter]);			//	 in MONOMERS
				counter++;
				(*chainlengthoftype)[i]+=blocklength[blockcount];
				blockcount++;
			}
			count=0;
            
            if(terminuscode==1) extrachainlength=1;
            else extrachainlength=0;
            (*chainlengthoftype)[i]+=extrachainlength;
            
			(*chainchemistry)[i]=xcalloc((*chainlengthoftype)[i], sizeof(int));
			for(j=0;j<blockcount;j++){
				for(k=0;k<blocklength[j];k++){
					(*chainchemistry)[i][count]=polartype[j];				//	start with polar
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];							//	backbone site
					(*pmonomercount).type[polartype[j]]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];								//	backbone
					count++;
				}
			}
			(*pNmonomers)+=((*chainlengthoftype)[i])*((*Nchainsoftype)[i]);
            if(terminuscode==1){
                (*chainchemistry)[i][count]=-1;                              //  chainchemistry[i][j]=-1 denotes backbone site with no sidechain site
            }
		}
	}
}

void associate_chains_and_nodes(int Nmonomers, int Nchains, int chaintypes, int *chainsoftype, int *chainlengthoftype, int **chainchemistry, int *pNnodes, coord **coordarray, monomernodes ***monomerid, int **chainlength, bonded_params my_bonded_params){
	int i, j,k,  nodespermonomer=2, chaincounter=0, nodecounter=0;
	(*pNnodes)=Nmonomers*nodespermonomer;
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*monomerid)=xcalloc(Nchains, sizeof(monomernodes *));
	(*chainlength)=xcalloc(Nchains, sizeof(int));
	for(i=0;i<chaintypes;i++){
		for(j=0;j<chainsoftype[i];j++){
			(*monomerid)[chaincounter]=xcalloc(chainlengthoftype[i], sizeof(monomernodes));
			(*chainlength)[chaincounter]=chainlengthoftype[i];
			for(k=0;k<chainlengthoftype[i];k++){
				(*monomerid)[chaincounter][k].backbone=nodecounter;
				(*coordarray)[nodecounter].chainid=chaincounter;
				(*coordarray)[nodecounter].monomerid=k;
				(*coordarray)[nodecounter].nodetype=0;							//	backbone
                (*coordarray)[nodecounter].leafid=0;                            //  default
				nodecounter++;
				(*monomerid)[chaincounter][k].sidechain=nodecounter;
				(*coordarray)[nodecounter].chainid=chaincounter;
				(*coordarray)[nodecounter].monomerid=k;
				(*coordarray)[nodecounter].nodetype=chainchemistry[i][k];		//	sidechain
                if((*coordarray)[nodecounter].nodetype==-1) (*coordarray)[nodecounter].orientationtype=-1;
				else (*coordarray)[nodecounter].orientationtype=my_bonded_params.sidechainorientationtype[(*coordarray)[nodecounter].nodetype];
                (*coordarray)[nodecounter].leafid=0;                            //  default
				nodecounter++;
			}
			chaincounter++;
		}
	}
}

void initialize_polymer_config(polymer *pconfig, double backbonezzigzag, double backbonenz, double backbonephi, int Nnodetypes, bonded_params my_bonded_params, double phenylheight){
    int i, j, last, self;
    allocate_matrix_nozero(double_triple, (*pconfig).backboneoffset, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonenz, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonephi, Nnodetypes, Nnodetypes);
    allocate_matrix_nozero(double_triple, (*pconfig).sidechainrelativepos, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnpar, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnphi, Nnodetypes, Nnodetypes);
    (*pconfig).backbonenz[0][2]=backbonenz;
    (*pconfig).backbonenz[0][3]=backbonenz;
    (*pconfig).backbonenz[2][1]=backbonenz;
    (*pconfig).backbonenz[3][1]=backbonenz;                                      //  perpendicular, keeping phi at 0
    (*pconfig).backbonephi[0][2]=backbonephi;
    (*pconfig).backbonephi[0][3]=backbonephi;
    (*pconfig).backbonephi[2][1]=backbonephi;
    (*pconfig).backbonephi[3][1]=backbonephi;
    (*pconfig).sidechainrelativepos[2][1]=(*pconfig).sidechainrelativepos[3][1]=rsol(my_bonded_params.sidechain[1]);
    (*pconfig).sidechainrelativepos[0][2]=rsol(my_bonded_params.sidechain[2]);
    (*pconfig).sidechainrelativepos[0][3]=rsol(my_bonded_params.sidechain[3]);
    (*pconfig).sidechainnpar[2][1]=(*pconfig).sidechainnpar[3][1]=nparhardsol(my_bonded_params.sidechainorientationtype[1], &((*pconfig).sidechainrelativepos[2][1]), my_bonded_params.sidechain[1]);
    (*pconfig).sidechainnpar[0][2]=nparhardsol(my_bonded_params.sidechainorientationtype[2], &((*pconfig).sidechainrelativepos[0][2]), my_bonded_params.sidechain[2]);
    (*pconfig).sidechainnpar[0][3]=nparhardsol(my_bonded_params.sidechainorientationtype[3], &((*pconfig).sidechainrelativepos[0][3]), my_bonded_params.sidechain[3]);
    (*pconfig).sidechainnphi[2][1]=(*pconfig).sidechainnphi[3][1]=nphihardsol(my_bonded_params.sidechainorientationtype[1], (*pconfig).sidechainnpar[2][1], (*pconfig).sidechainrelativepos[2][1], my_bonded_params.sidechain[1]);
    (*pconfig).sidechainnphi[0][2]=nphihardsol(my_bonded_params.sidechainorientationtype[2], (*pconfig).sidechainnpar[0][2], (*pconfig).sidechainrelativepos[0][2], my_bonded_params.sidechain[2]);
    (*pconfig).sidechainnphi[0][3]=nphihardsol(my_bonded_params.sidechainorientationtype[3], (*pconfig).sidechainnpar[0][3], (*pconfig).sidechainrelativepos[0][3], my_bonded_params.sidechain[3]);
    for(i=0;i<2;i++){
        for(j=2;j<4;j++){
            if(i==0){
                last=0;
                self=j;
            }
            else{
                last=j;
                self=1;
            }
            (*pconfig).backboneoffset[last][self].z=phenylheight-(*pconfig).sidechainrelativepos[2][1].z;
            if(last==0) (*pconfig).backboneoffset[last][self].z-=backbonezzigzag;
            (*pconfig).backboneoffset[last][self].x=(*pconfig).backboneoffset[last][self].y=0;
        }
    }

}

void initialize_2dlattice_fromconfig(int Nmonomers, int Nchains, double density, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, polymer config){
	int chainsonaside, i, j, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	(*pbox_dimension).x=pow(Nmonomers/density, 1.0/3.0);
	(*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
	chainsonaside=(int) floor(sqrt(1.*Nchains))+1;
	double chainspacing=(*pbox_dimension).x/chainsonaside;
	double chainspacingz, offsetz;
	if(my_nonbonded_params.solvationparams.interface==1){
		chainspacingz=((*pbox_dimension).z-my_nonbonded_params.vacuumthickness)/chainsonaside;
		offsetz=0.5*my_nonbonded_params.vacuumthickness;
	}
	else{
		chainspacingz=chainspacing;
		offsetz=0;
	}
	for(i=0;i<chainsonaside;i++){
		for(j=0;j<chainsonaside;j++){
			position.x=0.5*(*pbox_dimension).x;
			position.y=(i+0.5)*chainspacing;
			position.z=offsetz+(j+0.5)*chainspacingz;
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
					sidepos=config.sidechainrelativepos[lastpolartype][coordarray[sidechainid].nodetype];
					nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));  //  projection onto r_{back-side} but perp to nhat, assuming r_{back-side} has no y component
					normalize(&nalong);
					nside=cross_product(nback, nalong);
					npar=config.sidechainnpar[lastpolartype][coordarray[sidechainid].nodetype];
					nphi=config.sidechainnphi[lastpolartype][coordarray[sidechainid].nodetype];
					
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
}

void input_and_copy_single_polymer(char *filename, int *pNnodes, int Nchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, double molarity, double cutoff2, int *prunningtime){
	int i, j, inputchains, inputnodes, success, k, l;
	FILE *inp;
	inp=fopen(filename, "r");
	double_triple input_box_dimension, com, sep, try;
    sep.x=sep.y=sep.z=0;
	double maxr=0, packingfraction;
	
	fscanf(inp, "%lf %lf %lf", &((input_box_dimension).x), &((input_box_dimension).y), &((input_box_dimension).z));

	fscanf(inp, "%i", &(*pinterface));
    if((*pinterface)==1){
        fscanf(inp, " %lf", pvacuumthickness);
        (*pbox_dimension).x=(*pbox_dimension).y=pow((1.*Nchains)/(molarity), 1.0/2.0);
        (*pbox_dimension).z=(input_box_dimension).z;
    }
    else{
        (*pbox_dimension).x=pow((1.*Nchains)/(molarity*molardensity), 1.0/3.0);
        (*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
    }
	fscanf(inp, "%i %i %i %i %i %i", &inputnodes, &inputchains, &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	
	((*pmonomercount).charged)*=Nchains;
	((*pmonomercount).nonpolar)*=Nchains;
	((*pmonomercount).monomers)*=Nchains;
	
	if(inputchains!=1){
		printf("inputchains=%i! (%f %f %f %i)\n", inputchains, input_box_dimension.x, input_box_dimension.y, input_box_dimension.z, inputnodes);
		exit(1);
	}
	(*pNnodes)=Nchains*inputnodes;
	
	for(i=0;i<(*pNnodetypes);i++){
		fscanf(inp, " %i", &((*pmonomercount).type[i]));
		((*pmonomercount).type[i])*=Nchains;
	}
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((Nchains), sizeof(int));
	(*monomerid)=xcalloc((Nchains), sizeof(monomernodes *));
	fscanf(inp, "%i", &((*chainlength)[0]));								//	scan in one polymer, assign its attributes to all
	((*monomerid)[0])=xcalloc(((*chainlength)[0]), sizeof(monomernodes));
	for(j=0;j<(*chainlength)[0];j++){
		fscanf(inp, "%i %i", &((*monomerid)[0][j].backbone), &((*monomerid)[0][j].sidechain));
	}
	for(i=1;i<Nchains;i++){
		(*chainlength)[i]=(*chainlength)[0];
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[0];j++){
			(*monomerid)[i][j].backbone=(*monomerid)[0][j].backbone+i*inputnodes;
			(*monomerid)[i][j].sidechain=(*monomerid)[0][j].sidechain+i*inputnodes;
		}
	}
	com.x=com.y=com.z=0;
    int realnodes=0;
	for(i=0;i<inputnodes;i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
        if((*coordarray)[i].nodetype>=0){
            realnodes++;
            if(i==0){
                com=(*coordarray)[i].r;
            }
            else{
                sep=subtract_double_triple((*coordarray)[i].r, (*coordarray)[i-1].r);
                recenter_double_triple(&sep, input_box_dimension);
                (*coordarray)[i].r=add_double_triple((*coordarray)[i-1].r, sep);
                com=add_double_triple(com, (*coordarray)[i].r);
            }
        }
	}
	com=scalar_multiply_double_triple(com, (1./realnodes));
    if((*pinterface)==1) com.z=0;
    double norm2;
	for(i=0;i<inputnodes;i++){
        if((*coordarray)[i].nodetype>=0){
            (*coordarray)[i].r=subtract_double_triple((*coordarray)[i].r, com);
            recenter_double_triple(&((*coordarray)[i].r), input_box_dimension);
            norm2=(*coordarray)[i].r.x*(*coordarray)[i].r.x+(*coordarray)[i].r.y*(*coordarray)[i].r.y;
            if((*pinterface)==0) norm2+=(*coordarray)[i].r.z*(*coordarray)[i].r.z;
            if(norm2>maxr*maxr) maxr=sqrt(norm2);
        }
	}
	maxr=maxr+sqrt(cutoff2);
    if((*pinterface)==1) packingfraction=Nchains*M_PI*pow(maxr, 2)/((*pbox_dimension).x*(*pbox_dimension).y);
	else packingfraction=4.*Nchains*M_PI/3.*pow(maxr, 3)/((*pbox_dimension).x*(*pbox_dimension).y*(*pbox_dimension).z);
	if(packingfraction>0.5){
        printf("packing fraction=%f!\n", packingfraction);
        exit(1);
    }
	double_triple *centers;
	centers=xcalloc(Nchains, sizeof(double_triple));
	centers[0].x=centers[0].y=centers[0].z=0;
    for(i=1;i<Nchains;i++){
		success=0;
		while(success==0){
            if((*pinterface)==1){
                try.x=(*pbox_dimension).x*rand_double;
                try.y=(*pbox_dimension).y*rand_double;
                try.z=0;
            }
			else{
                try=rand_unit_cube();
                try=scalar_multiply_double_triple(try, (*pbox_dimension).x);
			}
            success=1;
			for(j=0;j<i;j++){
                sep=subtract_double_triple(centers[j], try);
				recenter_double_triple(&sep, (*pbox_dimension));
				if(norm(sep)<2*maxr) success=0;
			}
		}
		centers[i]=try;
	}
    double_triple axis_vector;
	declare_matrix(double, one_matrix, 3, 3);
	declare_matrix(double, new_matrix, 3, 3);
	declare_matrix(double, total_matrix, 3, 3);
    double angle;
    for(i=0;i<Nchains;i++){
        if((*pinterface)==0){
            for(j=0;j<3;j++){
                axis_vector.x=axis_vector.y=axis_vector.z=0;
                if(j==0) axis_vector.x=1;
                else if(j==1) axis_vector.y=1;
                else axis_vector.z=1;
                angle=2*M_PI*rand_double;
                forward_matrix(angle, axis_vector, one_matrix);
                if(j==0){
                    for(k=0;k<3;k++){
                        for(l=0;l<3;l++){
                            total_matrix[k][l]=one_matrix[k][l];
                        }
                    }
                }
                else{
                    matrix_multiply(total_matrix, one_matrix, new_matrix);
                    for(k=0;k<3;k++){
                        for(l=0;l<3;l++){
                            total_matrix[k][l]=new_matrix[k][l];
                        }
                    }
                }
            }
        }
        else{
            angle=2*M_PI*rand_double;
            axis_vector.x=axis_vector.y=0;
            axis_vector.z=1;
            forward_matrix(angle, axis_vector, total_matrix);
        }
		for(j=0;j<inputnodes;j++){
			(*coordarray)[i*inputnodes+j]=(*coordarray)[j];
			(*coordarray)[i*inputnodes+j].r=add_double_triple(centers[i], rotate_by_matrix((*coordarray)[j].r, total_matrix));
			(*coordarray)[i*inputnodes+j].n=rotate_by_matrix((*coordarray)[j].n, total_matrix);
			(*coordarray)[i*inputnodes+j].chainid=i;
		}
	}
    for(i=0;i<(*pNnodes);i++){
        fmod_double_triple(&((*coordarray)[i].r), (*pbox_dimension));
    }
    free(centers);
    free_matrix(one_matrix, 3);
    free_matrix(new_matrix, 3);
    free_matrix(total_matrix, 3);
	fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
}

void input_bilayer(char *configname, bilayer *pconfig, int chainlength, int Nnodetypes, double residueshift){
 	FILE *fp;
	fp=fopen(configname, "r");
	int size=1000, i, j, k;
	char *buff;
	buff=xcalloc(size, sizeof(char));
	while(fgets(buff, size, fp));
	fclose(fp);
	int offset;
	sscanf(buff, "%*i%n", &offset);
	buff+=offset;
	sscanf(buff, "%*f%n", &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerscaleoffsetfraction), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).terminusspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).interchainspacing), &offset);
	buff+=offset;
    
    (*pconfig).length=chainlength*(*pconfig).monomerspacing+(*pconfig).terminusspacing;
    (*pconfig).offset=(*pconfig).length*(*pconfig).offsetfraction+((*pconfig).monomerscaleoffsetfraction+residueshift)*(*pconfig).monomerspacing;
    
    allocate_3d_tensor_nozero(double_triple, (*pconfig).backboneoffset, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonenz, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonephi, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor_nozero(double_triple, (*pconfig).sidechainrelativepos, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnpar, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnphi, 2, Nnodetypes, Nnodetypes);
    
    int leaf, charged, chargetype, last, self;
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
                sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[leaf][last][self].x), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[leaf][last][self].y), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[leaf][last][self].z), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backbonenz[leaf][last][self]), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backbonephi[leaf][last][self]), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[leaf][last][self].x), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[leaf][last][self].y), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[leaf][last][self].z), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainnpar[leaf][last][self]), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainnphi[leaf][last][self]), &offset);
                buff+=offset;
            }
        }
    }
}

void initialize_bilayer_fragment(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry, int xchains, double molarity){
    int j, l, i, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback, periodicbox;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	periodicbox.x=spaceperchain*xchains;
	periodicbox.y=(Nchains/2/xchains)*config.interchainspacing;
	periodicbox.z=config.boxheight;
    
	(*pbox_dimension).x=pow((1.*Nchains)/(molarity*molardensity), 1.0/3.0);
	(*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
	if((periodicbox.x>(*pbox_dimension).x)||(periodicbox.y>(*pbox_dimension).y)) my_exit("molarity too high; sheet doesn't fit\n");
	
	if(Nchains%2!=0) my_exit("Nchains not divisible by 2!");
	int ychains=Nchains/2/xchains;
    for(j=0;j<2;j++){
		for(i=0;i<ychains;i++){
			for(l=0;l<xchains;l++){
				position.x=l*spaceperchain+(i%2)*config.offset-0.5*periodicbox.x+0.5*(*pbox_dimension).x;
				position.y=(i+0.5)*config.interchainspacing-0.5*periodicbox.y+0.5*(*pbox_dimension).y;
				position.z=0.5*(*pbox_dimension).z;
				for(k=0;k<chainlength[chaincounter];k++){
					backboneid=monomerid[chaincounter][k].backbone;
					sidechainid=monomerid[chaincounter][k].sidechain;
                    coordarray[backboneid].leafid=j;
                    coordarray[sidechainid].leafid=j;
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
					if(coordarray[sidechainid].nodetype>=0){
						nback=coordarray[backboneid].n;
						nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
						normalize(&nalong);
						nside=scalar_multiply_double_triple(cross_product(nback, nalong), (1-2*j));
						npar=config.sidechainnpar[j][lastpolartype][coordarray[sidechainid].nodetype];
						nphi=config.sidechainnphi[j][lastpolartype][coordarray[sidechainid].nodetype];
						sidepos=config.sidechainrelativepos[j][lastpolartype][coordarray[sidechainid].nodetype];
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
    }
}

void create_ordered_bilayer(bilayer *pconfig, int chainlength, int Nnodetypes, double residueshift, double chargedbackbonezoffset, double backbonezzigzag, double backbonenz, double backbonephi, bonded_params my_bonded_params, double leafspacing, double leafxoffsetfrac, double leafyoffsetfrac, int leafsymmetry){
	int i, j, k;
    
    (*pconfig).length=chainlength*(*pconfig).monomerspacing+(*pconfig).terminusspacing;
    (*pconfig).offset=(*pconfig).length*(*pconfig).offsetfraction+((*pconfig).monomerscaleoffsetfraction+residueshift)*(*pconfig).monomerspacing;
    
    allocate_3d_tensor_nozero(double_triple, (*pconfig).backboneoffset, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonenz, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonephi, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor_nozero(double_triple, (*pconfig).sidechainrelativepos, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnpar, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnphi, 2, Nnodetypes, Nnodetypes);
    
    int leaf, charged, chargetype, last, self;
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
				(*pconfig).backboneoffset[leaf][last][self].x=0+leaf*leafxoffsetfrac*(*pconfig).monomerspacing;
				(*pconfig).backboneoffset[leaf][last][self].y=0+leaf*leafyoffsetfrac*(*pconfig).interchainspacing;
				(*pconfig).backboneoffset[leaf][last][self].z=(1-leaf*2)*chargedbackbonezoffset+(leaf-0.5)*leafspacing;
				if(self==1) (*pconfig).backboneoffset[leaf][last][self].z+=(1-leaf*2)*backbonezzigzag;
				(*pconfig).backbonenz[leaf][last][self]=backbonenz;
				if((leaf==1)&&(leafsymmetry==-1)){
					(*pconfig).backbonephi[leaf][last][self]=M_PI-backbonephi;
				}
				else{
					(*pconfig).backbonephi[leaf][last][self]=backbonephi;
				}
                (*pconfig).sidechainrelativepos[leaf][last][self]=rsol(my_bonded_params.sidechain[self]);
                if((leaf==1)&&(leafsymmetry==-1)) (*pconfig).sidechainrelativepos[leaf][last][self].x=-(*pconfig).sidechainrelativepos[leaf][last][self].x;
				(*pconfig).sidechainnpar[leaf][last][self]=nparhardsol(my_bonded_params.sidechainorientationtype[self], &((*pconfig).sidechainrelativepos[leaf][last][self]), my_bonded_params.sidechain[self]);
				if((i==1)&&(leafsymmetry==-1)){
					(*pconfig).sidechainnphi[leaf][last][self]=M_PI-nphihardsol(my_bonded_params.sidechainorientationtype[self], (*pconfig).sidechainnpar[leaf][last][self], (*pconfig).sidechainrelativepos[leaf][last][self], my_bonded_params.sidechain[self]);
				}
				else{
					(*pconfig).sidechainnphi[leaf][last][self]=nphihardsol(my_bonded_params.sidechainorientationtype[self], (*pconfig).sidechainnpar[leaf][last][self], (*pconfig).sidechainrelativepos[leaf][last][self], my_bonded_params.sidechain[self]);
				}
			}
        }
    }
}

void initialize_periodic_bilayer_xchains(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry, int xchains){
    int j, l, i, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).x=spaceperchain*xchains;
	(*pbox_dimension).y=(Nchains/2/xchains)*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
    if(Nchains%2!=0) my_exit("Nchains not divisible by 2!");
	int ychains=Nchains/2/xchains;
    for(j=0;j<2;j++){
		for(i=0;i<ychains;i++){
			for(l=0;l<xchains;l++){
				position.x=l*spaceperchain+(i%2)*config.offset;
				position.y=(i+0.5)*config.interchainspacing;
				position.z=0.5*(*pbox_dimension).z;
				for(k=0;k<chainlength[chaincounter];k++){
					backboneid=monomerid[chaincounter][k].backbone;
					sidechainid=monomerid[chaincounter][k].sidechain;
                    coordarray[backboneid].leafid=j;
                    coordarray[sidechainid].leafid=j;
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
					if(coordarray[sidechainid].nodetype>=0){
						nback=coordarray[backboneid].n;
						nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
						normalize(&nalong);
						nside=scalar_multiply_double_triple(cross_product(nback, nalong), (1-2*j));
						npar=config.sidechainnpar[j][lastpolartype][coordarray[sidechainid].nodetype];
						nphi=config.sidechainnphi[j][lastpolartype][coordarray[sidechainid].nodetype];
						sidepos=config.sidechainrelativepos[j][lastpolartype][coordarray[sidechainid].nodetype];
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
    }
}

void input_monolayer(char *configname, monolayer *pconfig, int chainlength, int Nnodetypes, double residueshift){
 	FILE *fp;
	fp=fopen(configname, "r");
	int size=1000, i, j;
	char *buff;
	buff=xcalloc(size, sizeof(char));
	while(fgets(buff, size, fp));
	fclose(fp);
	int offset;
	sscanf(buff, "%*i%n", &offset);
	buff+=offset;
	sscanf(buff, "%*f%n", &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerscaleoffsetfraction), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).terminusspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).interchainspacing), &offset);
	buff+=offset;
    
    (*pconfig).length=chainlength*(*pconfig).monomerspacing+(*pconfig).terminusspacing;
    (*pconfig).offset=(*pconfig).length*(*pconfig).offsetfraction+((*pconfig).monomerscaleoffsetfraction+residueshift)*(*pconfig).monomerspacing;
    
	if(Nnodetypes<4) Nnodetypes=4;
    
    allocate_matrix_nozero(double_triple, ((*pconfig).backboneoffset), Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonenz, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonephi, Nnodetypes, Nnodetypes);
    allocate_matrix_nozero(double_triple, (*pconfig).sidechainrelativepos, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnpar, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnphi, Nnodetypes, Nnodetypes);
    
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
            sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[last][self].x), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[last][self].y), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[last][self].z), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backbonenz[last][self]), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backbonephi[last][self]), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[last][self].x), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[last][self].y), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[last][self].z), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainnpar[last][self]), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainnphi[last][self]), &offset);
            buff+=offset;
        }
    }
}

double_triple rsol(sidechain_params sidechain){
    double arctan=quartic_global_minimum(sidechain.J14, sidechain.J13, sidechain.J12, sidechain.J11);
    double rparoverrperp=tan(arctan);
    
    //  solve quadratic
    
    double a=pow(rparoverrperp, 2)+1;
    double b=-2*(rparoverrperp*sidechain.rpar0+sidechain.rperp0);
    double c=pow(sidechain.rperp0, 2)+pow(sidechain.rpar0, 2)-pow(sidechain.r0, 2);
    double rperp=(-b+sqrt(b*b-4*a*c))/(2*a);
    
    double rpar=rparoverrperp*rperp;
    double_triple result;
    result.z=rpar;
    result.x=rperp;
    result.y=0;
    return result;
}

double nparhardsol(int orientationtype, double_triple *pr, sidechain_params sidechain){
    if(orientationtype==1) return 0;
    double result=sidechain.r20+(*pr).z*sidechain.r21+pow((*pr).z, 2)*sidechain.r22, rperpdif, rpardif;
    if(result>1){
        printf("changing sidechainrelativepos.z from %f to", (*pr).z);
        (*pr).z=(-sidechain.r21+sqrt(pow(sidechain.r21, 2)-4*sidechain.r22*(sidechain.r20-1)))/(2*sidechain.r22);
        rpardif=(*pr).z-sidechain.rpar0;
        rperpdif=sqrt(pow(sidechain.r0, 2)-pow(rpardif, 2));
        (*pr).x=sidechain.rperp0+rperpdif;
        printf(" %f\n", (*pr).z);
        return 1;
    }
    else if(result<-1){
        printf("changing sidechainrelativepos.z from %f to", (*pr).z);
        (*pr).z=(-sidechain.r21+sqrt(pow(sidechain.r21, 2)-4*sidechain.r22*(sidechain.r20-1)))/(2*sidechain.r22);
        rpardif=(*pr).z-sidechain.rpar0;
        rperpdif=sqrt(pow(sidechain.r0, 2)-pow(rpardif, 2));
        (*pr).x=sidechain.rperp0+rperpdif;
        printf(" %f\n", (*pr).z);
        return -1;
    }
    return result;
}

double nphihardsol(int orientationtype, double npar, double_triple relativepos, sidechain_params sidechain){
    double rpar=relativepos.z;
    double rperp=relativepos.x;
    double rdotn, rdotnperphat;//, rmag=norm(relativepos);
    double nparmag=fabs(npar);
    if(orientationtype==0){
        if(nparmag==1){
            rdotn=0;
            return 0;           //  phi is unconstrained; choose 0 arbitrarily
        }
        rdotn=sidechain.r30+npar*sidechain.r31+npar*npar*sidechain.r32;
    }
    else if(orientationtype==1){
        double p4, p3, p2, p1;
        if(nparmag==1){
            rdotn=0;
            return 0;           //  phi is unconstrained; choose 0 arbitrarily
        }
        else{
            
            //  find minimum of quartic
            
            p4=sidechain.J340;
            p3=sidechain.J331*nparmag;
            p2=(sidechain.J320+sidechain.J322*nparmag*nparmag);
            p1=sidechain.J311*nparmag+sidechain.J313*nparmag*nparmag*nparmag;
            rdotn=quartic_global_minimum(p4, p3, p2, p1);
            if(npar<0) rdotn=-rdotn;
        }
    }
    rdotnperphat=(rdotn-rpar*npar)/(rperp*sqrt(1-npar*npar));
    if(rdotnperphat>1){
        rdotnperphat=1;
        return 0;
        
        //  can't adjust relativepos.x to fix this; would have to make it bigger
        
    }
    else if(rdotnperphat<-1){
        rdotnperphat=-1;
        
        //  can't adjust relativepos.x to fix this; would have to make it bigger
        
        return M_PI;
    }
    return acos(rdotnperphat);
}

void initialize_periodic_monolayer_xchains(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config, int xchains){
    int i, j, k, lastpolar=2, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).x=spaceperchain*xchains;
	(*pbox_dimension).y=(Nchains/xchains)*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
	
	int ychains=Nchains/xchains;
	for(i=0;i<ychains;i++){
		for(j=0;j<xchains;j++){
			position.x=j*spaceperchain+(i%2)*config.offset;
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
}

void initialize_periodic_monolayer_xchains_bothinterfaces(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config, int xchains){
    int i, j, k, interface, lastpolar=2, lastpolartype, chaincounter=0, backboneid, sidechainid, firstdirection, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).x=spaceperchain*xchains;
	int ychains=Nchains/xchains/2;
	(*pbox_dimension).y=ychains*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
	
    for(interface=0;interface<2;interface++){
        for(i=0;i<ychains;i++){
            for(j=0;j<xchains;j++){
                position.x=j*spaceperchain+(i%2)*config.offset;
                position.y=(i+0.5)*config.interchainspacing;
                if(interface==0) position.z=(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness;
                else position.z=0.5*my_nonbonded_params.vacuumthickness;
                for(k=0;k<chainlength[chaincounter];k++){
                    backboneid=monomerid[chaincounter][k].backbone;
                    sidechainid=monomerid[chaincounter][k].sidechain;
                    coordarray[backboneid].leafid=interface;
                    coordarray[sidechainid].leafid=interface;
                    if(coordarray[sidechainid].nodetype==1) firstdirection=1;           //  phenyl toward air
                    else firstdirection=-1;
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
                    nz=(2*interface-1)*firstdirection*config.backbonenz[lastpolartype][selftype];
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
    }
}

void initialize_2dxylattice_fromconfig(int Nmonomers, int Nchains, double density, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, polymer config){
	int chainsonaside, i, j, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	(*pbox_dimension).x=pow(Nmonomers/density, 1.0/3.0);
	(*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
	chainsonaside=(int) floor(sqrt(1.*Nchains))+1;
	double chainspacing=(*pbox_dimension).x/chainsonaside;
	position.z=(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness;
    for(i=0;i<chainsonaside;i++){
		for(j=0;j<chainsonaside;j++){
			position.x=(i+0.5)*chainspacing;
			position.y=(j+0.5)*chainspacing;
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
}

void input_coords(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime){
	int i, j;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &(*pinterface));
	if((*pinterface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &(*pNnodes), &(*pNchains), &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	for(i=0;i<(*pNnodetypes);i++) fscanf(inp, " %i", &((*pmonomercount).type[i]));
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<*pNchains;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	for(i=0;i<(*pNnodes);i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
	}
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
}

void input_coords_and_replicate_stacks(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, int replicates, double_triple displacement, double minsep){
	int i, j, Nchainsinput, Nnodesinput, r, replicateindex, count=0;
    double dif, av=0, max=0;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &(*pinterface));
	if((*pinterface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &Nnodesinput, &Nchainsinput, &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
    (*pNnodes)=replicates*Nnodesinput;
    (*pNchains)=replicates*Nchainsinput;
    ((*pmonomercount).charged)*=replicates;
    ((*pmonomercount).nonpolar)*=replicates;
    ((*pmonomercount).monomers)*=replicates;
	for(i=0;i<(*pNnodetypes);i++){
        fscanf(inp, " %i", &((*pmonomercount).type[i]));
        ((*pmonomercount).type[i])*=replicates;
    }
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<Nchainsinput;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	for(i=0;i<Nnodesinput;i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
        if((*coordarray)[i].nodetype>=0){                           //  not empty
            dif=(*coordarray)[i].r.z-(*coordarray)[0].r.z;
            recenter(dif, ((*pbox_dimension).z));
            av+=dif;
            count++;
        }
	}
    av=av/(1.*count)+(*coordarray)[0].r.z;
    recenter(av, ((*pbox_dimension).z));
    for(i=0;i<Nnodesinput;i++){
        if((*coordarray)[i].nodetype>=0){                           //  not empty
            dif=(*coordarray)[i].r.z-av;
            recenter(dif, ((*pbox_dimension).z));
            if(dif*dif>max) max=dif*dif;
        }
    }
    max=sqrt(max);
	printf("max displacement 2*%f+%f=%f\n", max, minsep, 2*max+minsep);
    if(2*max+minsep>displacement.z){
        printf("intersheet z displacement %f not bigger than 2*%f+%f!\n", displacement.z, max, minsep);
        exit(1);
    }
    
    ((*pbox_dimension).z)=replicates*displacement.z;
    
    for(i=0;i<Nnodesinput;i++){
        fmod_double_triple(&((*coordarray)[i].r), (*pbox_dimension));
    }
    
    replicateindex=Nnodesinput;
    for(r=1;r<replicates;r++){
        for(i=0;i<Nchainsinput;i++){
            (*chainlength)[r*Nchainsinput+i]=(*chainlength)[i];
           ((*monomerid)[r*Nchainsinput+i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
            for(j=0;j<(*chainlength)[i];j++){
                (*coordarray)[replicateindex]=(*coordarray)[(*monomerid)[i][j].backbone];
                if(i!=(*coordarray)[replicateindex].chainid){
                    printf("input chainids %i and %i don't agree!\n", i, (*coordarray)[replicateindex].chainid);
                    exit(1);
                }
                (*coordarray)[replicateindex].chainid=i+r*Nchainsinput;
				(*coordarray)[replicateindex].leafid+=2*r;
                (*monomerid)[i+r*Nchainsinput][j].backbone=replicateindex;
                (*coordarray)[replicateindex].r=add_double_triple((*coordarray)[replicateindex].r, scalar_multiply_double_triple(displacement, r));
                fmod_double_triple(&((*coordarray)[replicateindex].r), (*pbox_dimension));
                replicateindex++;
                (*coordarray)[replicateindex]=(*coordarray)[(*monomerid)[i][j].sidechain];
                if(i!=(*coordarray)[replicateindex].chainid){
                    printf("input chainids %i and %i don't agree!\n", i, (*coordarray)[replicateindex].chainid);
                    exit(1);
                }
                (*coordarray)[replicateindex].chainid=i+r*Nchainsinput;
 				(*coordarray)[replicateindex].leafid+=2*r;
                (*monomerid)[i+r*Nchainsinput][j].sidechain=replicateindex;
                (*coordarray)[replicateindex].r=add_double_triple((*coordarray)[replicateindex].r, scalar_multiply_double_triple(displacement, r));
                fmod_double_triple(&((*coordarray)[replicateindex].r), (*pbox_dimension));
                replicateindex++;
            }
        }
    }
    
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));

}

void input_coords_v0_bilayer(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime){
	int i, j;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &(*pinterface));
	if((*pinterface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &(*pNnodes), &(*pNchains), &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	for(i=0;i<(*pNnodetypes);i++) fscanf(inp, " %i", &((*pmonomercount).type[i]));
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<*pNchains;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	for(i=0;i<(*pNnodes);i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype));
	}
    for(i=0;i<(*pNnodes)/2;i++) (*coordarray)[i].leafid=0;
    for(i=(*pNnodes)/2;i<(*pNnodes);i++) (*coordarray)[i].leafid=1;
    fscanf(inp, "%lli %i", &(*pt), &(*pframe));
    (*prunningtime)=0;
}

void output_coords(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, int interface, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime){
	int i, j;
	FILE *outp;
	outp=fopen(filename, "w");
	fprintf(outp, "%.8f %.8f %.8f\n", box_dimension.x, box_dimension.y, box_dimension.z);
	fprintf(outp, "%i", interface);
	if(interface==1)fprintf(outp, " %.8f\n", vacuumthickness);
	else fprintf(outp, "\n");
	fprintf(outp, "%i %i %i %i %i %i", Nnodes, Nchains, monomercount.charged, monomercount.nonpolar, monomercount.monomers, Nnodetypes);
	for(i=0;i<Nnodetypes;i++) fprintf(outp, " %i", monomercount.type[i]);
	fprintf(outp, "\n");
	for(i=0;i<Nchains;i++){
		fprintf(outp, "%i", chainlength[i]);
		for(j=0;j<chainlength[i];j++){
			fprintf(outp, " %i %i", monomerid[i][j].backbone, monomerid[i][j].sidechain);
		}
		fprintf(outp, "\n");
	}
	for(i=0;i<Nnodes;i++){
		fprintf(outp, "%.8f %.8f %.8f %.8f %.8f %.8f %i %i %i %i %i\n", coordarray[i].r.x, coordarray[i].r.y, coordarray[i].r.z, coordarray[i].n.x, coordarray[i].n.y, coordarray[i].n.z, coordarray[i].nodetype, coordarray[i].chainid, coordarray[i].monomerid, coordarray[i].orientationtype, coordarray[i].leafid);
	}
    fprintf(outp, "%lli %i %i\n", t, frame, runningtime);
	fclose(outp);
}

void configure_cells_struct(linkedlistfull *plink, double_triple box_dimension){
	(*plink).core.cellsperside.x=floor(box_dimension.x/(*plink).mincellwidth)+1;
	(*plink).core.cellsperside.y=floor(box_dimension.y/(*plink).mincellwidth)+1;
	(*plink).core.cellsperside.z=floor(box_dimension.z/(*plink).mincellwidth)+1;
	if(((*plink).core.cellsperside.x<3)||((*plink).core.cellsperside.y<3)||((*plink).core.cellsperside.y<3)){
		printf("code not written for small cellsperside (%i x %i x %i)\n", (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
		printf("box_dimension=%f x %f x %f\n", box_dimension.x, box_dimension.y, box_dimension.z);
		exit(1);
	}
	(*plink).cellwidth.x=box_dimension.x/(*plink).core.cellsperside.x;
	(*plink).cellwidth.y=box_dimension.y/(*plink).core.cellsperside.y;
	(*plink).cellwidth.z=box_dimension.z/(*plink).core.cellsperside.z;
}

void allocate_linklist(linkedlistfull *plink, int Nnodes){
	int i, j, k;
	allocate_3d_tensor(int, (*plink).core.head, (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
	allocate_array(int, (*plink).core.list, Nnodes);
	allocate_array(int, (*plink).reverselist, Nnodes);
	allocate_array(int_triple, (*plink).core.cell, Nnodes);
	
	allocate_array(int, (*plink).core.number_neighbors, Nnodes);
	allocate_matrix(int, (*plink).core.neighbor, Nnodes, (*plink).maxneighbors);
}

void allocate_linklist_pairenergies(linkedlistfull *plink, int Nnodes){
	int i, j, k;
	allocate_3d_tensor(int, (*plink).core.head, (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
	allocate_array(int, (*plink).core.list, Nnodes);
	allocate_array(int, (*plink).reverselist, Nnodes);
	allocate_array(int_triple, (*plink).core.cell, Nnodes);
	
	allocate_array(int, (*plink).core.number_neighbors, Nnodes);
	allocate_matrix(int, (*plink).core.neighbor, Nnodes, (*plink).maxneighbors);
	allocate_matrix(double, (*plink).core.pair_energy, Nnodes, (*plink).maxneighbors);
}

void cellconstructdoublylinked(coord *coordarray, int N, linkedlistfull *plink, int lownodetype, int highnodetype){
	int i, j, k, n;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				(*plink).core.head[i][j][k]=-1;						//	empty
			}
		}
	}
	for(n=0;n<N;n++) (*plink).reverselist[n]=-1;
	for(n=0;n<N;n++){
		if((coordarray[n].nodetype>=lownodetype)&&(coordarray[n].nodetype<=highnodetype)){
			i=(int) floor(coordarray[n].r.x/(*plink).cellwidth.x);
			j=(int) floor(coordarray[n].r.y/(*plink).cellwidth.y);
			k=(int) floor(coordarray[n].r.z/(*plink).cellwidth.z);
			i=mod(i, (*plink).core.cellsperside.x);
			j=mod(j, (*plink).core.cellsperside.y);
			k=mod(k, (*plink).core.cellsperside.z);
			(*plink).core.cell[n].x=i;
			(*plink).core.cell[n].y=j;
			(*plink).core.cell[n].z=k;
			(*plink).core.list[n]=(*plink).core.head[i][j][k];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
			(*plink).core.head[i][j][k]=n;
		}
	}
}

void constructneighborlist(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=coordarray[self].r;
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
                    
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void constructneighborlist_pairenergies_phenyl(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double energy;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=coordarray[self].r;
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
							
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;

									energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
									
									(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
									(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;

								energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
								
								(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
								(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
					
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
							
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void constructneighborlist_pairenergies_charged(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double energy;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=coordarray[self].r;
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
									if(coordarray[self].nodetype==2){
										if(coordarray[neighbor].nodetype==2){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
										}
										else if(coordarray[neighbor].nodetype==3){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
										}
									}
									else if(coordarray[self].nodetype==3){
										if(coordarray[neighbor].nodetype==2){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
										}
										else if(coordarray[neighbor].nodetype==3){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
										}
									}
									(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
									(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
								if(coordarray[self].nodetype==2){
									if(coordarray[neighbor].nodetype==2){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
									}
									else if(coordarray[neighbor].nodetype==3){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
									}
								}
								else if(coordarray[self].nodetype==3){
									if(coordarray[neighbor].nodetype==2){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
									}
									else if(coordarray[neighbor].nodetype==3){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
									}
								}
								(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
								(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
					
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void update_neighbor_list(int Nchains, int *chainlength, coord *coordarray, nonbonded_params my_nonbonded_params, linkedlistset linkset, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_neighbors, int *assigned){
	int i, j, index, target, targetpolymer, firstnode, lastnode;
	double rsq;
	double_triple r;
	for(index=0;index<Nchains;index++) number_neighbors[index]=0;
	for(index=0;index<Nchains;index++){
		for(i=index+1;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
				for(j=0;j<linkset.phenylphenyl.core.number_neighbors[i];j++){
					target=linkset.phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(targetpolymer>index){							//	to avoid double counting
						if(assigned[targetpolymer]==0){
							r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
							rsq=dot_product(r, r);
							if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
								assigned[targetpolymer]=1;
								neighbor[index][number_neighbors[index]]=targetpolymer;
								number_neighbors[index]++;
								if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough");
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough");
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
				for(j=0;j<linkset.chargedcharged.core.number_neighbors[i];j++){
					target=linkset.chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(targetpolymer>index){							//	to avoid double counting
						if(assigned[targetpolymer]==0){
							r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
							rsq=dot_product(r, r);
							if(rsq<my_nonbonded_params.cutoff2){				//	interacting
								assigned[targetpolymer]=1;
								neighbor[index][number_neighbors[index]]=targetpolymer;
								number_neighbors[index]++;
								if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough");
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
 								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough");
							}
						}
					}
				}
			}
		}
	}
}

void deletefile(char *filename){
	FILE *outp;
	outp=fopen(filename, "w");
	fclose(outp);
}

void initialize_energycomponents_nofiles(energycomponents *pmy_energycomponents, int Nnodetypes){
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
    
	(*pmy_energycomponents).nncross=0;
	(*pmy_energycomponents).ccunlikecross=0;
	(*pmy_energycomponents).cclikecross=0;
	(*pmy_energycomponents).cncross=0;

	(*pmy_energycomponents).nndifferent=0;
	(*pmy_energycomponents).ccunlikedifferent=0;
	(*pmy_energycomponents).cclikedifferent=0;
	(*pmy_energycomponents).cndifferent=0;

	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
	int i;
	for(i=0;i<Nnodetypes;i++){
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
	}
}

void output_xyz_cg(int Nchains, int *chainlength, monomernodes **monomerid, int Natoms, coord *coordarray, double_triple box_dimension, char *xyzname, char *sourcename, cgparams params, int *pframe){
	int i, j, jmid, backboneindex, sidechainindex, atomcount=0, lastbackbonecount, lastbackboneindex, firstbackbonecount, newcgsourcefile=0;
	double_triple com, backbonebonded, sidevector, sidechainbonded, lastbackbonebonded, firstbackbonebonded, sep, lastbackbone;
	FILE *outp, *sourceoutp;
	
	if((*pframe)==0){
        newcgsourcefile=1;
		sourceoutp=fopen(sourcename, "w");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fprintf(sourceoutp, "topo clearbonds\n");
		fprintf(sourceoutp, "set sel [atomselect top \" name C\"]\n$sel set radius %.6f\n", params.backrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name H\"]\n$sel set radius %.6f\n", params.phenrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name N\"]\n$sel set radius %.6f\n", params.aminrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name O\"]\n$sel set radius %.6f\n", params.carbrad);
		fprintf(sourceoutp, "color Name C white\n");
		fprintf(sourceoutp, "color Name H yellow\n");
		fprintf(sourceoutp, "color Name N blue\n");
		fprintf(sourceoutp, "color Name O red\n");
		outp=fopen(xyzname, "w");
	}
	else{

		if ((sourceoutp = fopen(sourcename, "r")) == NULL) {
            newcgsourcefile=1;
			
			sourceoutp=fopen(sourcename, "w");
			fprintf(sourceoutp, "topo clearbonds\n");
			fprintf(sourceoutp, "set sel [atomselect top \" name C\"]\n$sel set radius %.6f\n", params.backrad);
			fprintf(sourceoutp, "set sel [atomselect top \" name H\"]\n$sel set radius %.6f\n", params.phenrad);
			fprintf(sourceoutp, "set sel [atomselect top \" name N\"]\n$sel set radius %.6f\n", params.aminrad);
			fprintf(sourceoutp, "set sel [atomselect top \" name O\"]\n$sel set radius %.6f\n", params.carbrad);
			
		} else {
			fclose(sourceoutp);
			sourceoutp=fopen(sourcename, "a");
		}
		
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		if(newcgsourcefile==0) fclose(sourceoutp);
		outp=fopen(xyzname, "a");
	}
	fprintf(outp, "%i\n%i\t%.6f\t%.6f\t%.6f\n", Natoms, Nchains, box_dimension.x, box_dimension.y, box_dimension.z);
	for(i=0;i<Nchains;i++){
		com.x=com.y=com.z=0;
		for(j=0;j<chainlength[i];j++){
			if(j==0){
				com=coordarray[monomerid[i][j].backbone].r;
				lastbackbone=coordarray[monomerid[i][j].backbone].r;
			}
			else{
				sep=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, lastbackbone);
				recenter_double_triple(&sep, box_dimension);
				lastbackbone=add_double_triple(lastbackbone, sep);
				com=add_double_triple(com, lastbackbone);
			}
		}
		com=scalar_multiply_double_triple(com, (1./chainlength[i]));
		fmod_double_triple(&com, box_dimension);
		for(jmid=0;jmid<chainlength[i];jmid++){
			if(jmid+chainlength[i]/2<chainlength[i]) j=chainlength[i]/2+jmid;
			else j=chainlength[i]/2-(1+jmid-(chainlength[i]-chainlength[i]/2));
			backboneindex=monomerid[i][j].backbone;
			if(j==chainlength[i]/2) firstbackbonebonded=backbonebonded=nearest_image(coordarray[backboneindex].r, com, box_dimension);
			else if(j==chainlength[i]/2-1) backbonebonded=nearest_image(coordarray[backboneindex].r, firstbackbonebonded, box_dimension);
			else backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			lastbackbonebonded=backbonebonded;
			lastbackboneindex=backboneindex;
			outputcoords_centered(outp, "C", backbonebonded, box_dimension);				//	backbone N
            if(newcgsourcefile==1){
				if(jmid==0) firstbackbonecount=atomcount;
				else if(jmid==chainlength[i]-chainlength[i]/2){
					fprintf(sourceoutp, "topo addbond %i %i\n", firstbackbonecount, atomcount);
				}
				else{
					fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount, atomcount);
				}
				lastbackbonecount=atomcount;
				atomcount++;
			}
			sidechainindex=monomerid[i][j].sidechain;
            sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
            recenter_double_triple(&sidevector, box_dimension);
            sidechainbonded=add_double_triple(backbonebonded, sidevector);
            if(coordarray[sidechainindex].nodetype>0){
                if(coordarray[sidechainindex].nodetype==1) outputcoords_centered(outp, "H", sidechainbonded, box_dimension);
                else if(coordarray[sidechainindex].nodetype==2) outputcoords_centered(outp, "N", sidechainbonded, box_dimension);
                else if(coordarray[sidechainindex].nodetype==3) outputcoords_centered(outp, "O", sidechainbonded, box_dimension);
                else my_exit("haven't written output_xyz_cg for other nodetypes");
                if(newcgsourcefile==1){
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount-1, atomcount);
                    atomcount+=1;
                }
			}
		}
	}
    if(newcgsourcefile==1){
		fclose(sourceoutp);
	}
	fclose(outp);
	(*pframe)++;
}

int rand_int(int max){
	int max_allowed=RAND_MAX-RAND_MAX%max;
	int temp=max_allowed;
	while(temp>=max_allowed){
		temp=rand()%max;
	}
	return temp;
}

void mc_translate_single_hard(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob=0, difference=0;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	difference=calc_energy_difference_changer_hard(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink).core);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		fmod_double_triple(&(coord_new.r), box_dimension);
		updatecell(coord_new.r, &(*plink), mover);
		coordarray[mover]=coord_new;
		updateneighborlist(&(*plink), mover, coordarray, box_dimension);
	}
}

void mc_translate_single_hard_pull(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, int pulledchain, double force){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	difference=calc_energy_difference_changer_hard(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink).core);
    if(coord_old.chainid==pulledchain){
        if(coord_old.monomerid==0) difference+=shift.x*force;
        else if(coord_old.monomerid==chainlength[pulledchain]-1) difference-=shift.x*force;
    }
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		fmod_double_triple(&(coord_new.r), box_dimension);
		updatecell(coord_new.r, &(*plink), mover);
		coordarray[mover]=coord_new;
		updateneighborlist(&(*plink), mover, coordarray, box_dimension);
	}
}

void mc_translate_single_phenyl(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink1, linkedlistfull *plink2, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	difference=calc_energy_difference_changer_phenyl(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink1).core, (*plink2).core);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		fmod_double_triple(&(coord_new.r), box_dimension);
		updatecell(coord_new.r, &(*plink1), mover);
		updatecell(coord_new.r, &(*plink2), mover);
		coordarray[mover]=coord_new;
		updateneighborlist(&(*plink1), mover, coordarray, box_dimension);
		updateneighborlist_pairenergies_phenyl(&(*plink2), mover, coordarray, box_dimension, my_nonbonded_params);
	}
}

void mc_translate_single_charged(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink1, linkedlistfull *plink2, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	difference=calc_energy_difference_changer_charged(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink1).core, (*plink2).core);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		fmod_double_triple(&(coord_new.r), box_dimension);
		updatecell(coord_new.r, &(*plink1), mover);
		updatecell(coord_new.r, &(*plink2), mover);
		coordarray[mover]=coord_new;
		updateneighborlist(&(*plink1), mover, coordarray, box_dimension);
		updateneighborlist_pairenergies_charged(&(*plink2), mover, coordarray, box_dimension, my_nonbonded_params);
	}
}

void mc_rotate_single_hard(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple newn=scalar_multiply_double_triple(rand_unit_cube(), max_rotate);
	newn=add_double_triple(coord_new.n, newn);
	normalize(&newn);
	coord_new.n=newn;
	difference=calc_energy_difference_changen_hard(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, link);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		coordarray[mover]=coord_new;
	}
}

void mc_rotate_single_phenyl(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple newn=scalar_multiply_double_triple(rand_unit_cube(), max_rotate);
	newn=add_double_triple(coord_new.n, newn);
	normalize(&newn);
	coord_new.n=newn;
	difference=calc_energy_difference_changen_phenyl(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, link);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		coordarray[mover]=coord_new;
	}
}

void mc_rotate_single_charged(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple newn=scalar_multiply_double_triple(rand_unit_cube(), max_rotate);
	newn=add_double_triple(coord_new.n, newn);
	normalize(&newn);
	coord_new.n=newn;
	difference=calc_energy_difference_changen_charged(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, link);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		coordarray[mover]=coord_new;
	}
}

int mc_aspect_ratios_cellstruct(int Nnodes, int Nchains, coord *coordarray, coord *newcoordarray, double_triple *pbox_dimension, double maxlogaspect, bonded_params my_bonded_params, nonbonded_params *pmy_nonbonded_params, int *chainlength, monomernodes **monomerid, linkedlistset linkset, pressureparam my_pressureparam, double temperature){
	double shiftlog, factor, complement, xfactor, yfactor, zfactor, newenergy, oldenergy, acceptance_prob, surfacepressureterm;
	double_triple new_box_dimension;
	int result=0, i, dim=rand_int(3);
    shiftlog=(rand_double*2-1)*maxlogaspect;
	factor=exp(shiftlog);
	complement=exp(-0.5*shiftlog);
	if(my_pressureparam.pressuretype==0){			//	isotropic
		if(dim==0){
			xfactor=factor;
			yfactor=zfactor=complement;
		}
		else if(dim==1){
			yfactor=factor;
			xfactor=zfactor=complement;
		}
		else{
			zfactor=factor;
			xfactor=yfactor=complement;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
	}
	else if(my_pressureparam.pressuretype==1){		//	planar (still isotropic, but only changing two sides per move)
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
        surfacepressureterm+=0.5*(xfactor+1.)*0.5*(yfactor+1.)*(zfactor-1)*(*pbox_dimension).x*(*pbox_dimension).y*(*pbox_dimension).z*my_pressureparam.normalforceperunitarea;
	}
	else if(my_pressureparam.pressuretype==2){		//	uniaxial in x direction
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(xfactor-1.)*0.5*(yfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
	}
	else if(my_pressureparam.pressuretype==3){		//	uniaxial in y direction
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(yfactor-1.)*0.5*(xfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
	}
	else if(my_pressureparam.pressuretype==4){		//	fixed z, only swapping x and y
		xfactor=factor;
		yfactor=1./factor;
		zfactor=1.;
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
	}
	else if(my_pressureparam.pressuretype==5){		//	fixed z, x and y fluctuating independently
		dim=rand_int(2);
		if(dim==0){
			xfactor=factor;
			yfactor=1.;
			zfactor=1.;
		}
		else if(dim==1){
			xfactor=1.;
			yfactor=factor;
			zfactor=1.;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
	}
    else if(my_pressureparam.pressuretype==6){		//	both x and y, different pressures
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(xfactor-1.)*0.5*(yfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
		surfacepressureterm+=((yfactor-1.)*0.5*(xfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.ypressure);
	}
	else if(my_pressureparam.pressuretype==7){			//	not conserving volume; isotropic surface pressure; separate vertical pressure; changing only one dimension at a time
		if(dim==0){
			xfactor=factor;
			yfactor=zfactor=1;
		}
		else if(dim==1){
			yfactor=factor;
			xfactor=zfactor=1;
		}
		else{
			zfactor=factor;
			xfactor=yfactor=1;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*my_pressureparam.surfacepressure;
        surfacepressureterm+=0.5*(xfactor+1.)*0.5*(yfactor+1.)*(zfactor-1)*(*pbox_dimension).x*(*pbox_dimension).y*(*pbox_dimension).z*my_pressureparam.normalforceperunitarea;
	}
    
	new_box_dimension.x=(*pbox_dimension).x*xfactor;
	new_box_dimension.y=(*pbox_dimension).y*yfactor;
	if((*pmy_nonbonded_params).solvationparams.interface==1) new_box_dimension.z=((*pbox_dimension).z-(*pmy_nonbonded_params).vacuumthickness)*zfactor+(*pmy_nonbonded_params).vacuumthickness;		//	dilate/contract only water part									//	dilate/contract only water part
	else new_box_dimension.z=(*pbox_dimension).z*zfactor;
	for(i=0;i<Nnodes;i++){
		newcoordarray[i]=coordarray[i];
		newcoordarray[i].r.x*=xfactor;
		newcoordarray[i].r.y*=yfactor;
		if((*pmy_nonbonded_params).solvationparams.interface==1) newcoordarray[i].r.z=(coordarray[i].r.z-0.5*(*pbox_dimension).z)*zfactor+0.5*new_box_dimension.z;				//	dilate/contract only water part by dilating/contracting relative to center of water slab
		else newcoordarray[i].r.z*=zfactor;
		fmod_double_triple(&(newcoordarray[i].r), new_box_dimension);																											//	in case went out of bounds due to numerical drift
	}
	nonbonded_params new_nonbonded_params=(*pmy_nonbonded_params);
	if((*pmy_nonbonded_params).solvationparams.interface==1) new_nonbonded_params.solvationparams.interfaceheight2=new_box_dimension.z-0.5*(*pmy_nonbonded_params).vacuumthickness;														//	 interface1 stays put
    
	newenergy=total_energy_cellstruct(newcoordarray, my_bonded_params, new_nonbonded_params, Nchains, chainlength, monomerid, new_box_dimension, linkset);			//	shouldn't need to redo cell until accepted, as long as maxlogaspect isn't too big
    if(newenergy<1.0/0.0){
		oldenergy=total_energy_cellstruct(coordarray, my_bonded_params, (*pmy_nonbonded_params), Nchains, chainlength, monomerid, *pbox_dimension, linkset);		//	shouldn't need to redo cell until accepted, as long as maxlogaspect isn't too big
		acceptance_prob=exp((oldenergy-newenergy-surfacepressureterm)/temperature);
		if((acceptance_prob>=1)||(rand_double<acceptance_prob)){
			result=1;
			(*pbox_dimension)=new_box_dimension;
			for(i=0;i<Nnodes;i++){
				coordarray[i]=newcoordarray[i];
			}
			(*pmy_nonbonded_params)=new_nonbonded_params;
		}
	}
	return result;
}

int change_cells_doublylinked(linkedlistfull *plink, double_triple box_dimension){
	int flag=0, i, j;
	int_triple newcellsperside;
	newcellsperside.x=floor(box_dimension.x/(*plink).mincellwidth)+1;
	newcellsperside.y=floor(box_dimension.y/(*plink).mincellwidth)+1;
	newcellsperside.z=floor(box_dimension.z/(*plink).mincellwidth)+1;
	if(newcellsperside.x!=(*plink).core.cellsperside.x) flag=1;
	if(newcellsperside.y!=(*plink).core.cellsperside.y) flag=1;
	if(newcellsperside.z!=(*plink).core.cellsperside.z) flag=1;
	if((newcellsperside.x<3)||(newcellsperside.y<3)||(newcellsperside.y<3)) my_exit("code not written for small cellsperside.");
	if(flag==1){
		for(i=0;i<(*plink).core.cellsperside.x;i++){
			for(j=0;j<(*plink).core.cellsperside.y;j++){
				free((*plink).core.head[i][j]);
			}
			free((*plink).core.head[i]);
		}
		free((*plink).core.head);
		(*plink).core.head=xcalloc(newcellsperside.x, sizeof(int **));
		for(i=0;i<newcellsperside.x;i++){
			(*plink).core.head[i]=xcalloc(newcellsperside.y, sizeof(int *));
			for(j=0;j<newcellsperside.y;j++){
				(*plink).core.head[i][j]=xcalloc(newcellsperside.z, sizeof(int));
			}
		}
	}
	(*plink).core.cellsperside=newcellsperside;
	(*plink).cellwidth.x=box_dimension.x/(*plink).core.cellsperside.x;
	(*plink).cellwidth.y=box_dimension.y/(*plink).core.cellsperside.y;
	(*plink).cellwidth.z=box_dimension.z/(*plink).core.cellsperside.z;
	return flag;
}

void vmmc_translate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, linkedlistset linkset, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned, int move3d){
	double individual_energy, collective_energy, probability, reverse_individual_energy, reverse_probability;
	int seed=rand_int(Nchains), cluster_complete=0, number_interior_members=0, number_exterior_members=1, cluster_number=1, i, j, k, l, m, candidate, index, firstnode, lastnode, exterior_member_rank=0, number_unchecked_neighbors, neighbor_rank, frustrated_link_index=0, temp, firstnodecandidate, lastnodecandidate;
	
	for(i=0;i<Nchains;i++) in_cluster[i]=0;
	in_cluster[seed]=1;
    double_triple shift;
    if(move3d==0){
        shift.x=(rand_int(2)*2-1)*max_translate*rand_double;
        shift.y=(rand_int(2)*2-1)*max_translate*rand_double;
        shift.z=0;
    }
	else shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
    
	index=seed;
    
	exterior_member[0]=index;
	firstnode=monomerid[index][0].backbone;
	lastnode=monomerid[index][chainlength[index]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;						//	in case long terminus
	for(i=firstnode;i<=lastnode;i++){
 		newcoordarray[i]=coordarray[i];
		reversecoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		reversecoordarray[i].r=subtract_double_triple(coordarray[i].r, shift);											//	don't need to apply b.c. because calc_pair computes minimum image
	}
	while(cluster_complete==0){
		number_unchecked_neighbors=number_neighbors[index];
		while((number_unchecked_neighbors>0)&&(cluster_complete==0)){
			neighbor_rank=rand_int(number_unchecked_neighbors);						//	choose randomly from neighbors
			candidate=neighbor[index][neighbor_rank];
			for(i=neighbor_rank;i<number_unchecked_neighbors-1;i++){				//	shift neighbor list to close gap
				neighbor[index][i]=neighbor[index][i+1];
			}
			neighbor[index][i]=candidate;
			number_unchecked_neighbors--;
			if(in_cluster[candidate]==0){											//	not already in cluster
				individual_energy=collective_energy=0;

				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnodecandidate].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;

				for(i=firstnode;i<=lastnode;i++){
					individual_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					collective_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);

					if(newcoordarray[i].nodetype==1){
						individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
					else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
						individual_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
				}
				probability=1.0-exp(-(individual_energy-collective_energy)/temperature);
				if(probability>0){
					if((rand_double<probability)||(probability==1)){						//	pre-link
						reverse_individual_energy=0;
						for(i=firstnode;i<=lastnode;i++){
							reverse_individual_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
							if(newcoordarray[i].nodetype==1){
								reverse_individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
							}
							else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
								reverse_individual_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
							}
						}
						reverse_probability=1.0-exp(-(reverse_individual_energy-collective_energy)/temperature);
						
						if((reverse_probability>probability)||(rand_double<reverse_probability/probability)){		//	form link
							cluster_number++;
							exterior_member[number_exterior_members]=candidate;
							number_exterior_members++;
							in_cluster[candidate]=1;
							for(i=firstnodecandidate;i<=lastnodecandidate;i++){
								newcoordarray[i]=coordarray[i];
								reversecoordarray[i]=coordarray[i];
								newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
								reversecoordarray[i].r=subtract_double_triple(coordarray[i].r, shift);					//	don't need to apply b.c. because calc_pair computes minimum image
							}
						}
						else{
							frustrated_link[frustrated_link_index].x=index;				//	mark LINK as frustrated
							frustrated_link[frustrated_link_index].y=candidate;
							frustrated_link_index++;
							if(frustrated_link_index==max_frustrated_links){
								printf("max_frustrated_links=%i not big enough.\n", max_frustrated_links);
								exit(1);
							}
						}
					}
				}
			}
		}
		interior_member[number_interior_members]=index;								//	searched all neighbors; shift index from exterior to interior
		number_interior_members++;
		for(i=exterior_member_rank;i<number_exterior_members-1;i++){
			exterior_member[i]=exterior_member[i+1];
		}
		number_exterior_members--;
		if(number_exterior_members==0) cluster_complete=1;
		else{																		//	assign new cluster member to check neighbors of
			exterior_member_rank=rand_int(number_exterior_members);
			index=exterior_member[exterior_member_rank];
			firstnode=monomerid[index][0].backbone;
			lastnode=monomerid[index][chainlength[index]-1].sidechain;
			if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		}
	}
	for(i=0;i<frustrated_link_index;i++){
		if((in_cluster[frustrated_link[i].x]==0)||(in_cluster[frustrated_link[i].y]==0)){
			return;																	//	if there are any frustrated LINKS external to the cluster, reject
		}
	}
	
	//	Acceptance prob is min(1, x), where x is a product over pairs {i, j}
	//	where i is in cluster and j is outside
	//	energy either
	//		initially overlapping (positive, neighbor) and finally noninteracting (0)
	//		or initially noninteracting (0, not a neighbor but check in case it is a neighbor) and finally overlapping (positive)
	//	Note that if I use D(C) not equal to 1, then that becomes the acceptance probability. But for now I am putting the size dependence in the distribution of max_number.
    
	double sum_energy=0, initial_energy, final_energy;
	declare_array(int, neighborofindex, Nchains);
	declare_array(double, initiallynoninteractingenergiesbypolymer, Nchains);
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(m=0;m<Nchains;m++){
			neighborofindex[m]=0;
			initiallynoninteractingenergiesbypolymer[m]=0;
		}
		neighborofindex[index]=1;		//	self
		for(m=0;m<number_neighbors[index];m++){
			candidate=neighbor[index][m];
			neighborofindex[candidate]=1;
			if(in_cluster[candidate]==0){
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnode].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				initial_energy=0;

				for(i=firstnode;i<=lastnode;i++){
					for(j=firstnodecandidate;j<=lastnodecandidate;j++){
						initial_energy+=calc_energy_pair(coordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
					}
				}
				
				if(initial_energy>0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					
					if(final_energy==0) sum_energy-=initial_energy;
				}
				else if(initial_energy==0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					if(final_energy>0) sum_energy+=final_energy;
				}
			}
		}
		//	next, add up energies for any polymers that were not initially neighbors
		//	bin the energies by polymers, because only positive pairwise polymer-polymer energies contribute to acceptance criterion
		
		for(i=firstnode;i<lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);                   //  need to apply boundary conditions because I look up cells, not just neighbors, here
            add_hardenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).shortrange, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			if(coordarray[i].nodetype==1) add_phenylphenylenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).phenylphenyl, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
                add_chargedchargedenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).chargedcharged, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
            }
		}
		for(m=0;m<Nchains;m++){
			if(initiallynoninteractingenergiesbypolymer[m]>0) sum_energy+=initiallynoninteractingenergiesbypolymer[m];
		}
	}
	free(neighborofindex);
	free(initiallynoninteractingenergiesbypolymer);
	
	//	external potential
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			sum_energy+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
		}
	}
	probability=exp(-sum_energy/temperature);
	if((probability<1.0)&&(rand_double>probability)){			//	reject
		return;
	}
    
	//	otherwise accept
	
	//	update coordinates and cells
    
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);
			updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
			if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
			coordarray[i]=newcoordarray[i];
		}
	}
	
	//	update site-site neighbor list
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
			if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_phenyl(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
		}
	}
	
	//	Update polymer-polymer neighbor lists. First, wipe out neighbor lists starting and ending in the cluster
	
	for(i=0;i<number_interior_members;i++){
		index=interior_member[i];
		for(j=0;j<number_neighbors[index];j++){
			temp=neighbor[index][j];
			if(in_cluster[temp]==0){
				for(k=0;k<number_neighbors[temp];k++){
					if(neighbor[temp][k]==index) break;								//	find the reverse neighbor map; can I do this more efficiently?
				}
				if(k==number_neighbors[temp]){
					printf("Error: did not find reverse map.\n");
					exit(1);
				}
				number_neighbors[temp]--;											//	remove reverse neighbor map
				for(l=k;l<number_neighbors[temp];l++){
					neighbor[temp][l]=neighbor[temp][l+1];
				}
			}
		}
		number_neighbors[index]=0;
	}
    
	//	Now update neighbor lists
    
	int target, targetpolymer, found;
	double rsq;
	double_triple r;
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=0;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
        
        for(i=0;i<number_neighbors[index];i++) assigned[neighbor[index][i]]=1;
		
        for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
				for(j=0;j<linkset.phenylphenyl.core.number_neighbors[i];j++){
					target=linkset.phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (a)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (b)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (c)");
								}
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
				for(j=0;j<linkset.chargedcharged.core.number_neighbors[i];j++){
					target=linkset.chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.cutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (d)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
 								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (e)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
 									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (f)");
								}
							}
						}
					}
				}
			}
		}
	}
}		

void vmmc_translate_wholepolymer_slow(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, linkedlistset linkset, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned){
	double individual_energy, collective_energy, probability, reverse_individual_energy, reverse_probability;
	int seed=rand_int(Nchains), cluster_complete=0, number_interior_members=0, number_exterior_members=1, cluster_number=1, i, j, k, l, m, candidate, index, firstnode, lastnode, firstnodecandidate, lastnodecandidate, exterior_member_rank=0, number_unchecked_neighbors, neighbor_rank, frustrated_link_index=0, temp;
	
	for(i=0;i<Nchains;i++) in_cluster[i]=0;
	in_cluster[seed]=1;
    double_triple shift;
    if(my_nonbonded_params.solvationparams.interface==1){
        shift.x=(rand_int(2)*2-1)*max_translate*rand_double;
        shift.y=(rand_int(2)*2-1)*max_translate*rand_double;
        shift.z=0;
    }
	else shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
    
	index=seed;
	exterior_member[0]=index;
	firstnode=monomerid[index][0].backbone;
	lastnode=monomerid[index][chainlength[index]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;						//	in case long terminus
	for(i=firstnode;i<=lastnode;i++){
 		newcoordarray[i]=coordarray[i];
		reversecoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		reversecoordarray[i].r=subtract_double_triple(coordarray[i].r, shift);											//	don't need to apply b.c. because calc_pair computes minimum image
	}
	while(cluster_complete==0){
		number_unchecked_neighbors=number_neighbors[index];
		while((number_unchecked_neighbors>0)&&(cluster_complete==0)){
			neighbor_rank=rand_int(number_unchecked_neighbors);						//	choose randomly from neighbors
			candidate=neighbor[index][neighbor_rank];
			for(i=neighbor_rank;i<number_unchecked_neighbors-1;i++){				//	shift neighbor list to close gap
				neighbor[index][i]=neighbor[index][i+1];
			}
			neighbor[index][i]=candidate;
			number_unchecked_neighbors--;
			if(in_cluster[candidate]==0){											//	not already in cluster
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnodecandidate].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				individual_energy=collective_energy=0;
				for(i=firstnode;i<=lastnode;i++){
					for(j=firstnodecandidate;j<=lastnodecandidate;j++){
						individual_energy+=calc_energy_pair(newcoordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
						collective_energy+=calc_energy_pair(coordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
					}
				}
				probability=1.0-exp(-(individual_energy-collective_energy)/temperature);
				if(probability>0){
					if((rand_double<probability)||(probability==1)){						//	pre-link
						reverse_individual_energy=0;
						for(i=firstnode;i<=lastnode;i++){
							for(j=firstnodecandidate;j<=lastnodecandidate;j++){
								reverse_individual_energy+=calc_energy_pair(reversecoordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
								reverse_probability=1.0-exp(-(reverse_individual_energy-collective_energy)/temperature);
							}
						}
						if((reverse_probability>probability)||(rand_double<reverse_probability/probability)){		//	form link
							cluster_number++;
							exterior_member[number_exterior_members]=candidate;
							number_exterior_members++;
							in_cluster[candidate]=1;
							for(i=firstnodecandidate;i<=lastnodecandidate;i++){
								newcoordarray[i]=coordarray[i];
								reversecoordarray[i]=coordarray[i];
								newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
								reversecoordarray[i].r=subtract_double_triple(coordarray[i].r, shift);					//	don't need to apply b.c. because calc_pair computes minimum image
							}
						}
						else{
							frustrated_link[frustrated_link_index].x=index;				//	mark LINK as frustrated
							frustrated_link[frustrated_link_index].y=candidate;
							frustrated_link_index++;
							if(frustrated_link_index==max_frustrated_links){
								printf("max_frustrated_links=%i not big enough.\n", max_frustrated_links);
								exit(1);
							}
						}
					}
				}
			}
		}
		interior_member[number_interior_members]=index;								//	searched all neighbors; shift index from exterior to interior
		number_interior_members++;
		for(i=exterior_member_rank;i<number_exterior_members-1;i++){
			exterior_member[i]=exterior_member[i+1];
		}
		number_exterior_members--;
		if(number_exterior_members==0) cluster_complete=1;
		else{																		//	assign new cluster member to check neighbors of
			exterior_member_rank=rand_int(number_exterior_members);
			index=exterior_member[exterior_member_rank];
			firstnode=monomerid[index][0].backbone;
			lastnode=monomerid[index][chainlength[index]-1].sidechain;
			if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		}
	}
	for(i=0;i<frustrated_link_index;i++){
		if((in_cluster[frustrated_link[i].x]==0)||(in_cluster[frustrated_link[i].y]==0)){
			return;																	//	if there are any frustrated LINKS external to the cluster, reject
		}
	}
	
	//	Acceptance prob is min(1, x), where x is a product over pairs {i, j}
	//	where i is in cluster and j is outside
	//	energy either
	//		initially overlapping (positive) and finally noninteracting (0)
	//		or initially noninteracting (0) and finally overlapping (positive)
	//	Note that if I use D(C) not equal to 1, then that becomes the acceptance probability. But for now I am putting the size dependence in the distribution of max_number.
    
	double sum_energy=0, initial_energy, final_energy;
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(m=0;m<number_neighbors[index];m++){
			candidate=neighbor[index][m];
			if(in_cluster[candidate]==0){
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnode].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				initial_energy=0;
				
				for(i=firstnode;i<=lastnode;i++){
					for(j=firstnodecandidate;j<=lastnodecandidate;j++){
						initial_energy+=calc_energy_pair(coordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
					}
				}
				
				if(initial_energy>0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						for(j=firstnodecandidate;j<=lastnodecandidate;j++){
							final_energy+=calc_energy_pair(newcoordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
						}
					}
					if(final_energy==0) sum_energy-=initial_energy;
				}
				else if(initial_energy==0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						for(j=firstnodecandidate;j<=lastnodecandidate;j++){
							final_energy+=calc_energy_pair(newcoordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
						}
					}
					if(final_energy>0) sum_energy+=final_energy;
				}
			}
		}
	}
	
	//	external potential
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			sum_energy+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
		}
	}
	probability=exp(-sum_energy/temperature);
	if((probability<1.0)&&(rand_double>probability)){			//	reject
		return;
	}
    
	//	otherwise accept
	
	//	update coordinates and cells
    
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);
			updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
			if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
			coordarray[i]=newcoordarray[i];
		}
	}
	
	//	update site-site neighbor list
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
			if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_phenyl(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
		}
	}
	
	//	Update polymer-polymer neighbor lists. First, wipe out neighbor lists starting and ending in the cluster
	
	for(i=0;i<number_interior_members;i++){
		index=interior_member[i];
		for(j=0;j<number_neighbors[index];j++){
			temp=neighbor[index][j];
			if(in_cluster[temp]==0){
				for(k=0;k<number_neighbors[temp];k++){
					if(neighbor[temp][k]==index) break;								//	find the reverse neighbor map; can I do this more efficiently?
				}
				if(k==number_neighbors[temp]){
					printf("Error: did not find reverse map.\n");
					exit(1);
				}
				number_neighbors[temp]--;											//	remove reverse neighbor map
				for(l=k;l<number_neighbors[temp];l++){
					neighbor[temp][l]=neighbor[temp][l+1];
				}
			}
		}
		number_neighbors[index]=0;
	}
    
	//	Now update neighbor lists
    
	int target, targetpolymer, found;
	double rsq;
	double_triple r;
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=0;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
        
        for(i=0;i<number_neighbors[index];i++) assigned[neighbor[index][i]]=1;
		
        for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
				for(j=0;j<linkset.phenylphenyl.core.number_neighbors[i];j++){
					target=linkset.phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (a)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (b)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (c)");
								}
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
				for(j=0;j<linkset.chargedcharged.core.number_neighbors[i];j++){
					target=linkset.chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.cutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (d)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
 								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (e)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
 									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (f)");
								}
							}
						}
					}
				}
			}
		}
	}
}		

double_triple com(int firstnode, int lastnode, coord *coordarray, double_triple box_dimension){
	int realnodes=0, i;
	double_triple result, sep, pos;
	for(i=firstnode;i<=lastnode;i++){
        if(coordarray[i].nodetype>=0){
            realnodes++;
            if(i==firstnode){
                pos=coordarray[i].r;
                result=pos;
            }
            else{
                sep=subtract_double_triple(coordarray[i].r, coordarray[i-1].r);
                recenter_double_triple(&sep, box_dimension);
                pos=add_double_triple(pos, sep);
                result=add_double_triple(result, pos);
            }
        }
	}
	result=scalar_multiply_double_triple(result, (1./realnodes));
    fmod_double_triple(&result, box_dimension);
	return result;
}

int calccomcheckr(double_triple *pcom, int firstnode, int lastnode, coord *coordarray, double_triple box_dimension, double maxr){
	int realnodes=0, i;
	double_triple sep, pos;
	for(i=firstnode;i<=lastnode;i++){
        if(coordarray[i].nodetype>=0){
            realnodes++;
            if(i==firstnode){
                pos=coordarray[i].r;
                (*pcom)=pos;
            }
            else{
                sep=subtract_double_triple(coordarray[i].r, coordarray[i-1].r);
                recenter_double_triple(&sep, box_dimension);
                pos=add_double_triple(pos, sep);
                (*pcom)=add_double_triple((*pcom), pos);
            }
        }
	}
	(*pcom)=scalar_multiply_double_triple((*pcom), (1./realnodes));
    fmod_double_triple(&(*pcom), box_dimension);
	for(i=firstnode;i<=lastnode;i++){
        if(coordarray[i].nodetype>=0){
			sep=subtract_double_triple(coordarray[i].r, (*pcom));
			recenter_double_triple(&sep, box_dimension);
			if(norm(sep)>maxr){
				//printf("site distance from com=%f!\n", norm(sep));
				return 1;
			}
		}
	}
	return 0;
}

int_triple cellofpos(double_triple pos, double_triple cellwidth){
	int_triple new;
	new.x=(int) floor(pos.x/cellwidth.x);
	new.y=(int) floor(pos.y/cellwidth.y);
	new.z=(int) floor(pos.z/cellwidth.z);
	return new;
}
	
void vmmc_rotate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, linkedlistset linkset, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned, int move3d, double maxr){
	double individual_energy, collective_energy, probability, reverse_individual_energy, reverse_probability;
	int checkrflag, seed=rand_int(Nchains), cluster_complete=0, number_interior_members=0, number_exterior_members=1, cluster_number=1, i, j, k, l, m, candidate, index, firstnode, lastnode, exterior_member_rank=0, number_unchecked_neighbors, neighbor_rank, frustrated_link_index=0, temp, firstnodecandidate, lastnodecandidate, checkselfoverlap=0;
	for(i=0;i<Nchains;i++) in_cluster[i]=0;
	in_cluster[seed]=1;
    double_triple axis_vector, sep;
	double angle=(2.0*rand_double-1.0)*max_rotate;
	declare_matrix(double, forward_matrix, 3, 3);
	declare_matrix(double, backward_matrix, 3, 3);
	declare_array_nozero(double_triple, mycom, Nchains);
    if(move3d==0){
		axis_vector.z=1;
		axis_vector.x=0;
		axis_vector.y=0;
    }
	else axis_vector=rand_unit_sphere();
	forward_and_backward_matrix(angle, axis_vector, forward_matrix, backward_matrix);

	double minbox=box_dimension.x;
	if(box_dimension.y<minbox) minbox=box_dimension.y;
	if(box_dimension.z<minbox) minbox=box_dimension.z;
	
	index=seed;

	exterior_member[0]=index;
	firstnode=monomerid[index][0].backbone;
	lastnode=monomerid[index][chainlength[index]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;						//	in case long terminus
	
	checkrflag=calccomcheckr(&(mycom[index]), firstnode, lastnode, coordarray, box_dimension, maxr);										//	self-consistently, make sure no site is too far from com or could defeat selfoverlap check
	if(checkrflag==1){
		free_matrix(forward_matrix, 3);
		free_matrix(backward_matrix, 3);
		free(mycom);
        //printf("\treturn due to maxr\n");
		return;																										
	}
	
	for(i=firstnode;i<=lastnode;i++){
 		newcoordarray[i]=coordarray[i];
		reversecoordarray[i]=coordarray[i];
		sep=subtract_double_triple(coordarray[i].r, mycom[index]);
		recenter_double_triple(&sep, box_dimension);
		newcoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, forward_matrix));
		reversecoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, backward_matrix));
		newcoordarray[i].n=rotate_by_matrix(coordarray[i].n, forward_matrix);
		reversecoordarray[i].n=rotate_by_matrix(coordarray[i].n, backward_matrix);
	}
	while(cluster_complete==0){
		number_unchecked_neighbors=number_neighbors[index];
		while((number_unchecked_neighbors>0)&&(cluster_complete==0)){
			neighbor_rank=rand_int(number_unchecked_neighbors);						//	choose randomly from neighbors
			candidate=neighbor[index][neighbor_rank];
			for(i=neighbor_rank;i<number_unchecked_neighbors-1;i++){				//	shift neighbor list to close gap
				neighbor[index][i]=neighbor[index][i+1];
			}
			neighbor[index][i]=candidate;
			number_unchecked_neighbors--;
			if(in_cluster[candidate]==0){											//	not already in cluster
				individual_energy=collective_energy=0;				
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnodecandidate].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				
				for(i=firstnode;i<=lastnode;i++){
					individual_energy+=calc_hard_energy_cellstruct_givenpolymer_recreatecell((*plinkset).shortrange, newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					collective_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					if(newcoordarray[i].nodetype==1){
						individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer_recreatecell((*plinkset).phenylphenyl, newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
					else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
						individual_energy+=calc_charged_energy_cellstruct_givenpolymer_recreatecell((*plinkset).chargedcharged, newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
				}
				probability=1.0-exp(-(individual_energy-collective_energy)/temperature);
				if(probability>0){
					if((rand_double<probability)||(probability==1)){						//	pre-link
						reverse_individual_energy=0;
						for(i=firstnode;i<=lastnode;i++){
							reverse_individual_energy+=calc_hard_energy_cellstruct_givenpolymer_recreatecell((*plinkset).shortrange, reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
							if(newcoordarray[i].nodetype==1){
								reverse_individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer_recreatecell((*plinkset).phenylphenyl, reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
							}
							else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
								reverse_individual_energy+=calc_charged_energy_cellstruct_givenpolymer_recreatecell((*plinkset).chargedcharged, reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
							}
						}
						reverse_probability=1.0-exp(-(reverse_individual_energy-collective_energy)/temperature);
						
						if((reverse_probability>probability)||(rand_double<reverse_probability/probability)){		//	form link
							cluster_number++;
							exterior_member[number_exterior_members]=candidate;
							number_exterior_members++;
							in_cluster[candidate]=1;
														
							checkrflag=calccomcheckr(&(mycom[candidate]), firstnodecandidate, lastnodecandidate, coordarray, box_dimension, maxr);
							if(checkrflag==1){
								free_matrix(forward_matrix, 3);
								free_matrix(backward_matrix, 3);
								free(mycom);
								return;
							}							
							mycom[candidate]=subtract_double_triple(mycom[candidate], mycom[index]);
							recenter_double_triple(&(mycom[candidate]), box_dimension);
							mycom[candidate]=add_double_triple(mycom[candidate], mycom[index]);
														
							for(i=firstnodecandidate;i<=lastnodecandidate;i++){
								newcoordarray[i]=coordarray[i];
								reversecoordarray[i]=coordarray[i];
								sep=subtract_double_triple(coordarray[i].r, mycom[candidate]);
								recenter_double_triple(&sep, box_dimension);
								sep=subtract_double_triple(add_double_triple(sep, mycom[candidate]), mycom[seed]);
								newcoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, forward_matrix));
								reversecoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, backward_matrix));
								newcoordarray[i].n=rotate_by_matrix(coordarray[i].n, forward_matrix);
								reversecoordarray[i].n=rotate_by_matrix(coordarray[i].n, backward_matrix);
							}
							if(checkselfoverlap==0){
								if(norm(mycom[candidate])+2*maxr+sqrt(my_nonbonded_params.cutoff2)>minbox) checkselfoverlap=1;
							}
						}
						else{
							frustrated_link[frustrated_link_index].x=index;				//	mark LINK as frustrated
							frustrated_link[frustrated_link_index].y=candidate;
							frustrated_link_index++;
							if(frustrated_link_index==max_frustrated_links){
								printf("max_frustrated_links=%i not big enough.\n", max_frustrated_links);
								exit(1);
							}
						}
					}
				}
			}
		}
		interior_member[number_interior_members]=index;								//	searched all neighbors; shift index from exterior to interior
		number_interior_members++;
		for(i=exterior_member_rank;i<number_exterior_members-1;i++){
			exterior_member[i]=exterior_member[i+1];
		}
		number_exterior_members--;
		if(number_exterior_members==0) cluster_complete=1;
		else{																		//	assign new cluster member to check neighbors of
			exterior_member_rank=rand_int(number_exterior_members);
			index=exterior_member[exterior_member_rank];
			firstnode=monomerid[index][0].backbone;
			lastnode=monomerid[index][chainlength[index]-1].sidechain;
			if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		}
	}
	for(i=0;i<frustrated_link_index;i++){
		if((in_cluster[frustrated_link[i].x]==0)||(in_cluster[frustrated_link[i].y]==0)){
			free_matrix(forward_matrix, 3);
			free_matrix(backward_matrix, 3);
			free(mycom);
			return;																	//	if there are any frustrated LINKS external to the cluster, reject
		}
	}
	
	//	Acceptance prob is min(1, x), where x is a product over pairs {i, j}
	//	where i is in cluster and j is outside
	//	energy either
	//		initially overlapping (positive) and finally noninteracting (0)
	//		or initially noninteracting (0) and finally overlapping (positive)
	//	Note that if I use D(C) not equal to 1, then that becomes the acceptance probability. But for now I am putting the size dependence in the distribution of max_number.
    
	//	BECAUSE OF POSSIBLE LARGE DISPLACEMENTS IN ROTATIONS, NEED TO ACCOUNT FOR POSSIBILITY THAT FINALLY POSITIVE INTERACTING PAIRS WERE NOT INITIALLY NEIGHBORS
		
	double sum_energy=0, initial_energy, final_energy;
	declare_array(int, neighborofindex, Nchains);
	declare_array(double, initiallynoninteractingenergiesbypolymer, Nchains);
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
 		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;

		//	first, do regular sum over neighbors of polymer #index
		
		for(m=0;m<Nchains;m++){
			neighborofindex[m]=0;
			initiallynoninteractingenergiesbypolymer[m]=0;
		}
		neighborofindex[index]=1;		//	self
		for(m=0;m<number_neighbors[index];m++){
			candidate=neighbor[index][m];
			neighborofindex[candidate]=1;
			
			//	calculating contribution from initial neighbors index and candidate
			
			if(in_cluster[candidate]==0){
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnode].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				initial_energy=0;
				
				for(i=firstnode;i<=lastnode;i++){
					for(j=firstnodecandidate;j<=lastnodecandidate;j++){
						initial_energy+=calc_energy_pair(coordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
					}
				}
				
				if(initial_energy>0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					
					if(final_energy==0) sum_energy-=initial_energy;
				}
				else if(initial_energy==0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					if(final_energy>0) sum_energy+=final_energy;
				}
			}
		}
		
		//	next, add up energies for any polymers that were not initially neighbors
		//	bin the energies by polymers, because only positive pairwise polymer-polymer energies contribute to acceptance criterion
		
		for(i=firstnode;i<lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);                   //  need to apply boundary conditions because I look up cells, not just neighbors, here
            add_hardenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).shortrange, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			if(coordarray[i].nodetype==1) add_phenylphenylenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).phenylphenyl, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
                add_chargedchargedenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).chargedcharged, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
            }
		}
		for(m=0;m<Nchains;m++){
			if(initiallynoninteractingenergiesbypolymer[m]>0) sum_energy+=initiallynoninteractingenergiesbypolymer[m];
		}
	}
	free(neighborofindex);
	free(initiallynoninteractingenergiesbypolymer);
	
	//	external potential
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			sum_energy+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
		}
	}
	probability=exp(-sum_energy/temperature);
	if((probability<1.0)&&(rand_double>probability)){			//	reject
		free_matrix(forward_matrix, 3);
		free_matrix(backward_matrix, 3);
		free(mycom);
		return;
	}
    
	//	otherwise accept
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);
			updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
			if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
			coordarray[i]=newcoordarray[i];
		}
	}
	
	//	update site-site neighbor list
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
			if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_phenyl(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
		}
	}
	
	//	Update polymer-polymer neighbor lists. First, wipe out neighbor lists starting and ending in the cluster
	
	for(i=0;i<number_interior_members;i++){
		index=interior_member[i];
		for(j=0;j<number_neighbors[index];j++){
			temp=neighbor[index][j];
			if(in_cluster[temp]==0){
				for(k=0;k<number_neighbors[temp];k++){
					if(neighbor[temp][k]==index) break;								//	find the reverse neighbor map; can I do this more efficiently?
				}
				if(k==number_neighbors[temp]){
					printf("Error: did not find reverse map.\n");
					exit(1);
				}
				number_neighbors[temp]--;											//	remove reverse neighbor map
				for(l=k;l<number_neighbors[temp];l++){
					neighbor[temp][l]=neighbor[temp][l+1];
				}
			}
		}
		number_neighbors[index]=0;                                                  //  remove all forward maps
	}
    
	//	Now update neighbor lists
    
	int target, targetpolymer, found;
	double rsq;
	double_triple r;
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=0;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
        
        for(i=0;i<number_neighbors[index];i++) assigned[neighbor[index][i]]=1;
		
        for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
				for(j=0;j<linkset.phenylphenyl.core.number_neighbors[i];j++){
					target=linkset.phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (a)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
                                
                                //  okay because in_cluster equivalent to interior_member, since there should be no exterior_members
                                
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (b)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (c)");
								}
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
				for(j=0;j<linkset.chargedcharged.core.number_neighbors[i];j++){
					target=linkset.chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.cutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (d)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
 								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (e)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
 									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (f)");
								}
							}
						}
					}
				}
			}
		}
	}
	free_matrix(forward_matrix, 3);
	free_matrix(backward_matrix, 3);
	free(mycom);
}

double calc_hard_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
    fmod_double_triple(&(coord_new.r), box_dimension);
	int_triple newcell=cellofpos(coord_new.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, a, b, c;
    double result=0;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					if(coordarray[trialneighbor].chainid==polymer){
                         
						if(coord_new.nodetype==0){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, 2.*my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(coord_new.nodetype==1){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(coord_new.nodetype==2){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						else if(coord_new.nodetype==3){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
                    
				}
			}
		}
	}
    return result;
}

double calc_phenyl_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
    fmod_double_triple(&(coord_new.r), box_dimension);
	int_triple newcell=cellofpos(coord_new.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, a, b, c;
    double result=0;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					if(coordarray[trialneighbor].chainid==polymer){
						result+=phenylphenyl_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
					
				}
			}
		}
	}
    return result;
}

double calc_charged_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
    fmod_double_triple(&(coord_new.r), box_dimension);
	int_triple newcell=cellofpos(coord_new.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, a, b, c;
    double result=0;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
 					if(coordarray[trialneighbor].chainid==polymer){
						if(coord_new.nodetype==2){
							if(coordarray[trialneighbor].nodetype==2){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
						else if(coord_new.nodetype==3){
							if(coordarray[trialneighbor].nodetype==2){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
					}
                    trialneighbor=mylinkedlistfull.core.list[trialneighbor];
				}
			}
		}
	}
    return result;
}

void add_hardenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	int_triple newcell=cellofpos(newcoord.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, trialpolymer, a, b, c;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					trialpolymer=coordarray[trialneighbor].chainid;
					if(list[trialpolymer]==0){
							
						if(newcoord.nodetype==0){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, 2.*my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(newcoord.nodetype==1){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(newcoord.nodetype==2){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						else if(newcoord.nodetype==3){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];

				}
			}
		}
	}
}

void add_phenylphenylenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	int_triple newcell=cellofpos(newcoord.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, trialpolymer, a, b, c;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					trialpolymer=coordarray[trialneighbor].chainid;
					if(list[trialpolymer]==0){						
						energy[trialpolymer]+=phenylphenyl_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
					
				}
			}
		}
	}
}

void add_chargedchargedenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	int_triple newcell=cellofpos(newcoord.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, trialpolymer, a, b, c;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					trialpolymer=coordarray[trialneighbor].chainid;
					if(list[trialpolymer]==0){		
						if(newcoord.nodetype==2){
							if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
						else if(newcoord.nodetype==3){
							if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
					}
                    trialneighbor=mylinkedlistfull.core.list[trialneighbor];					
				}
			}
		}
	}
}

void mc_translate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid){
	int movingpolymer=rand_int(Nchains), firstnode, lastnode, i;
	firstnode=monomerid[movingpolymer][0].backbone;
	lastnode=monomerid[movingpolymer][chainlength[movingpolymer]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[movingpolymer][chainlength[movingpolymer]-1].backbone;
	double acceptance_prob, difference=0;
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	for(i=firstnode;i<=lastnode;i++){
		newcoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		difference+=calc_hard_energy_cellstruct_otherpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
		difference-=calc_hard_energy_cellstruct_otherpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_otherpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_phenyl_energy_cellstruct_otherpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_otherpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_charged_energy_cellstruct_otherpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
		}
		difference+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	}
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        for(i=firstnode;i<=lastnode;i++){
            fmod_double_triple(&(newcoordarray[i].r), box_dimension);
            updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
            if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
            coordarray[i]=newcoordarray[i];
        }
        for(i=firstnode;i<=lastnode;i++){
            updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
            if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_phenyl(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
        }
	}
}

void mc_translate_wholebilayer(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid){
	int movingsheet=rand_int(Nsheets), firstnode, lastnode, i, foundfirst=0;
    for(i=0;i<Nnodes;i++){
        if(coordarray[i].leafid/2==movingsheet){
            if(foundfirst==0){
                firstnode=i;
                foundfirst=1;
            }
            lastnode=i;
        }
    }
	if(coordarray[lastnode].nodetype<0) lastnode--;                         //  accounting for possibility of long right terminus
	double acceptance_prob, difference=0;
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	for(i=firstnode;i<=lastnode;i++){
		newcoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		difference+=calc_hard_energy_cellstruct_otherbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
		difference-=calc_hard_energy_cellstruct_otherbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_otherbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_phenyl_energy_cellstruct_otherbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_otherbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_charged_energy_cellstruct_otherbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
		}
		difference+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	}
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        for(i=firstnode;i<=lastnode;i++){
            fmod_double_triple(&(newcoordarray[i].r), box_dimension);
            updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
            if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
            coordarray[i]=newcoordarray[i];
        }
        for(i=firstnode;i<=lastnode;i++){
            updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
            if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_phenyl(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
        }
	}
}

int mc_shiftbilayergap(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple *pbox_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, double normalforceperunitarea, int sitespersheet){
	int changinginterval=rand_int(Nsheets), result=0, lowermover, uppermover, i, r;
    lowermover=changinginterval;
    uppermover=(changinginterval+1);                                            //  don't modulate this for purposes of shifting; only calculating energy
	double sep;
    int countlower=0, countupper=0;
	declare_array(int, count, Nsheets);
	declare_array(double, com, Nsheets);
	declare_array(double, first, Nsheets);

    for(r=0;r<Nsheets;r++){
        first[r]=coordarray[r*sitespersheet].r.z;
        for(i=0;i<sitespersheet;i++){
            if(coordarray[i].nodetype>=0){										//  accounting for possibility of long right terminus
                sep=coordarray[r*sitespersheet+i].r.z-first[r];
                recenter((sep), ((*pbox_dimension).z));
                com[r]+=sep;
                count[r]++;
            }
        }
        com[r]=com[r]/(1.*count[r])+first[r];
        fmod((com[r]), ((*pbox_dimension).z));
    }
    for(r=0;r<lowermover;r++){
        sep=com[lowermover]-com[r];
        fmod((sep), ((*pbox_dimension).z));
        com[r]=com[lowermover]-sep;                                         //  no boundary should separate sheets 0, ..., lowermover-1 from lowermover
    }
    for(r=uppermover;r<Nsheets;r++){
        sep=com[r]-com[lowermover];
        fmod((sep), ((*pbox_dimension).z));
        com[r]=com[lowermover]+sep;                                         //  no boundary should separate sheets uppermover, ..., Nsheets-1 from lowermover
    }
    for(r=0;r<Nsheets;r++){
        for(i=0;i<sitespersheet;i++){
            newcoordarray[r*sitespersheet+i]=coordarray[r*sitespersheet+i];
            sep=coordarray[r*sitespersheet+i].r.z-com[r];
            recenter(sep, (*pbox_dimension).z);
            newcoordarray[r*sitespersheet+i].r.z=com[r]+sep;
        }
    }                                                                       //  boundary is between sheet Nsheets-1 and 0
           	
    double acceptance_prob, difference=0;
    double shift=rand_double*2-1;
    double_triple new_box_dimension=(*pbox_dimension);
    new_box_dimension.z+=shift;
    for(i=uppermover*sitespersheet;i<Nnodes;i++){
		newcoordarray[i].r.z+=shift;    //  move upper sheets and sheets above it by shift
    }
    
    uppermover=uppermover%Nsheets;                                          //  now shift uppermover from Nsheets to 0 if equal to Nsheets
	    
    for(i=lowermover*sitespersheet;i<(lowermover+1)*sitespersheet;i++){
		difference+=calc_hard_energy_cellstruct_specificbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params, uppermover);
		difference-=calc_hard_energy_cellstruct_specificbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params, uppermover);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_specificbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params, uppermover);
			difference-=calc_phenyl_energy_cellstruct_specificbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params, uppermover);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_specificbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params, uppermover);
			difference-=calc_charged_energy_cellstruct_specificbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params, uppermover);
		}        
        if(my_nonbonded_params.solvationparams.interface!=0) my_exit("mc_translate_wholebilyaer_changez not written for interface!\n");
	}
	acceptance_prob=exp(-(difference+normalforceperunitarea*(new_box_dimension.z-(*pbox_dimension).z)*(*pbox_dimension).x*(*pbox_dimension).y)/temperature);
	
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        result=1;
		printf("box from %f to ", (*pbox_dimension).z);
        (*pbox_dimension)=new_box_dimension;
		printf("%f\n", (*pbox_dimension).z);
		for(i=0;i<Nnodes;i++){
            fmod_double_triple(&(newcoordarray[i].r), (*pbox_dimension));
            coordarray[i]=newcoordarray[i];
        }
	}
	free(count);
	free(com);
	free(first);
    return result;
}

int mc_translate_wholebilayer_changez_old(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple *pbox_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, double normalforceperunitarea){
	int movingsheet=rand_int(Nsheets), firstnode, lastnode, i, foundfirst=0, result=0, sheet;
	double sep;
	declare_array(double, first, Nsheets);
	declare_array(int, found, Nsheets);
	declare_array(int, count, Nsheets);
	declare_array(double, com, Nsheets);
    for(i=0;i<Nnodes;i++){
		if(coordarray[i].nodetype>=0){										//  accounting for possibility of long right terminus
			newcoordarray[i]=coordarray[i];
			sheet=coordarray[i].leafid/2;
			if(found[sheet]==0){
				found[sheet]=1;
				first[sheet]=coordarray[i].r.z;
			}
			sep=coordarray[i].r.z-first[sheet];
			recenter(sep, (*pbox_dimension).z);
			com[sheet]+=sep;
			count[sheet]++;
			if(coordarray[i].leafid/2==movingsheet){
				if(foundfirst==0){
					firstnode=i;
					foundfirst=1;
				}
				lastnode=i;
			}
		}
    }
	for(sheet=0;sheet<Nsheets;sheet++){
		com[sheet]=com[sheet]/count[sheet]+first[sheet];
		fmod(com[sheet], (*pbox_dimension).z);
	}
	for(i=0;i<Nnodes;i++){
		if(coordarray[i].nodetype>=0){										//  accounting for possibility of long right terminus
			sheet=coordarray[i].leafid/2;
			sep=coordarray[i].r.z-com[sheet];
			recenter(sep, (*pbox_dimension).z);
			newcoordarray[i].r.z=com[sheet]+sep;							//	in newcoordarray, making sure no boundary separates z coordinates of same sheet			
		}
	}
	
	double acceptance_prob, difference=0;
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
    
    int boxchangedirection=2*rand_int(2)-1;
    double_triple new_box_dimension=(*pbox_dimension);
    new_box_dimension.z+=boxchangedirection*shift.z;
	
    for(i=firstnode;i<=lastnode;i++){
		newcoordarray[i].r=add_double_triple(newcoordarray[i].r, shift);
		difference+=calc_hard_energy_cellstruct_otherbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params);
		difference-=calc_hard_energy_cellstruct_otherbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_otherbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params);
			difference-=calc_phenyl_energy_cellstruct_otherbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_otherbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params);
			difference-=calc_charged_energy_cellstruct_otherbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params);
		}
        
        if(my_nonbonded_params.solvationparams.interface!=0) my_exit("mc_translate_wholebilyaer_changez not written for interface!\n");
	}
    
	acceptance_prob=exp(-(difference+normalforceperunitarea*(new_box_dimension.z-(*pbox_dimension).z)*(*pbox_dimension).x*(*pbox_dimension).y)/temperature);
	
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        result=1;
        (*pbox_dimension)=new_box_dimension;
		for(i=0;i<Nnodes;i++){            
            fmod_double_triple(&(newcoordarray[i].r), (*pbox_dimension));
            coordarray[i]=newcoordarray[i];
        }
	}
    return result;
}

void mc_reptate(int Nchains, int *chainlength, coord *coordarray, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, reptationparams leftreptation, reptationparams rightreptation, sidereptationparams *sidereptation){
    int debug=0;
	int movingpolymer=rand_int(Nchains), i, reversesecondsidenodetype, reverseendsidenodetype;
    int direction=2*rand_int(2)-1;
    int last=chainlength[movingpolymer]-1;
    if(coordarray[monomerid[movingpolymer][last].sidechain].nodetype>0) my_exit("mc_reptate not written to account for empty sidechains!");
	reptationparams reptation, reversereptation;
	int reverseendindex, reversesecondindex, reversethirdindex, reversefourthindex, endindex, secondindex;
    if(direction==1){		//	right
		reverseendindex=0;
		reversesecondindex=1;
		reversethirdindex=2;
		reversefourthindex=3;
		endindex=last;
		secondindex=last-1;
		reptation=rightreptation;
		reversereptation=leftreptation;
	}
	else{					//	left
		reverseendindex=last;
		reversesecondindex=last-1;
		reversethirdindex=last-2;
		reversefourthindex=last-3;
		endindex=0;
		secondindex=1;
		reptation=leftreptation;
		reversereptation=rightreptation;
	}
    
	//  make sure that the reverse move would be proposed
    
	int reversefourthsite=monomerid[movingpolymer][reversefourthindex].backbone;
	int reversethirdsite=monomerid[movingpolymer][reversethirdindex].backbone;
	int reversesecondsite=monomerid[movingpolymer][reversesecondindex].backbone;
	if(in_neighborhood(coordarray[reversefourthsite], coordarray[reversethirdsite], coordarray[reversesecondsite], reversereptation, box_dimension)==0){
		return;
	}
	int reverseendsite=monomerid[movingpolymer][reverseendindex].backbone;
	if(in_neighborhood(coordarray[reversethirdsite], coordarray[reversesecondsite], coordarray[reverseendsite], reversereptation, box_dimension)==0){
		return;
	}
	int reversesecondsidesite=monomerid[movingpolymer][reversesecondindex].sidechain;
	reversesecondsidenodetype=coordarray[reversesecondsidesite].nodetype;
	if(in_side_neighborhood(coordarray[reversesecondsite], coordarray[monomerid[movingpolymer][reversesecondindex].sidechain], sidereptation[reversesecondsidenodetype], box_dimension)==0){
		return;
	}
	int reverseendsidesite=monomerid[movingpolymer][reverseendindex].sidechain;
	reverseendsidenodetype=coordarray[monomerid[movingpolymer][reverseendindex].sidechain].nodetype;
	if(in_side_neighborhood(coordarray[reverseendsite], coordarray[reverseendsidesite], sidereptation[reverseendsidenodetype], box_dimension)==0){
		return;
	}
    if(debug==1) printf("passed neighborhood criteria\n");
    
	int second_number_neighbors, end_number_neighbors, secondside_shortrange_number_neighbors, secondside_longrange_number_neighbors, endside_shortrange_number_neighbors, endside_longrange_number_neighbors;
    double difference=0;
	
	//  reptate: calculate new nonbonded interactions
	
	int endsite=monomerid[movingpolymer][endindex].backbone;
	int secondsite=monomerid[movingpolymer][secondindex].backbone;
	
	//	new second-to-end backbone
	
    if(debug==1) printf("new second...\n");
	declare_array(int, second_neighborlist, (*plinkset).shortrange.maxneighbors);
	coord newsecondcoord=reptate(coordarray[monomerid[movingpolymer][secondindex].backbone], coordarray[monomerid[movingpolymer][endindex].backbone], rightreptation, box_dimension);
	int_triple newsecondcell=cell(newsecondcoord.r, (*plinkset).shortrange.cellwidth);
	create_temporary_neighborlist((*plinkset).shortrange, coordarray, box_dimension, &second_number_neighbors, second_neighborlist, newsecondcell, newsecondcoord);
	difference+=calc_hard_energy_cellstruct(endsite, -1, -1, second_number_neighbors, second_neighborlist, newsecondcoord, coordarray, box_dimension, my_nonbonded_params);
    free(second_neighborlist);
	
	//	new end backbone
	
    if(debug==1) printf("\tnew end (difference %f)...\n", difference);
    declare_array(int, end_neighborlist, (*plinkset).shortrange.maxneighbors);
	coord newendcoord=reptate(coordarray[monomerid[movingpolymer][endindex].backbone], newsecondcoord, rightreptation, box_dimension);
	int_triple newendcell=cell(newendcoord.r, (*plinkset).shortrange.cellwidth);
	create_temporary_neighborlist((*plinkset).shortrange, coordarray, box_dimension, &end_number_neighbors, end_neighborlist, newendcell, newendcoord);
	difference+=calc_hard_energy_cellstruct(-1, -1, -1, end_number_neighbors, end_neighborlist, newendcoord, coordarray, box_dimension, my_nonbonded_params);
    free(end_neighborlist);
    
	//	new second-to-end sidechain
	
    if(debug==1) printf("\tnew second sidechain (difference %f)...\n", difference);
    int secondsidesite=monomerid[movingpolymer][secondindex].sidechain;
	int secondsidenodetype=coordarray[secondsidesite].nodetype;
	coord newsecondsidecoord=side_reptate(newsecondcoord, coordarray[secondsidesite], sidereptation[secondsidenodetype], box_dimension);
    declare_array(int, secondside_shortrange_neighborlist, (*plinkset).shortrange.maxneighbors);
	int *secondside_longrange_neighborlist;
    if(secondsidenodetype==1){
		allocate_array(int, secondside_longrange_neighborlist, (*plinkset).phenylphenyl.maxneighbors);
	}
	else if((secondsidenodetype==2)||(secondsidenodetype==3)){
		allocate_array(int, secondside_longrange_neighborlist, (*plinkset).chargedcharged.maxneighbors);
	}
	else my_exit("code not written for other nodetypes");
	int_triple newsecondsideshortrangecell=cell(newsecondsidecoord.r, (*plinkset).shortrange.cellwidth);
	int_triple newsecondsidelongrangecell;
	if(secondsidenodetype==1) newsecondsidelongrangecell=cell(newsecondsidecoord.r, (*plinkset).phenylphenyl.cellwidth);
	else newsecondsidelongrangecell=cell(newsecondsidecoord.r, (*plinkset).chargedcharged.cellwidth);
	create_temporary_neighborlist((*plinkset).shortrange, coordarray, box_dimension, &secondside_shortrange_number_neighbors, secondside_shortrange_neighborlist, newsecondsideshortrangecell, newsecondsidecoord);
	if(secondsidenodetype==1){
		create_temporary_neighborlist((*plinkset).phenylphenyl, coordarray, box_dimension, &secondside_longrange_number_neighbors, secondside_longrange_neighborlist, newsecondsidelongrangecell, newsecondsidecoord);
	}
	else{
		create_temporary_neighborlist((*plinkset).chargedcharged, coordarray, box_dimension, &secondside_longrange_number_neighbors, secondside_longrange_neighborlist, newsecondsidelongrangecell, newsecondsidecoord);
	}
	difference+=calc_hard_energy_cellstruct(secondsite, -1, -1, secondside_shortrange_number_neighbors, secondside_shortrange_neighborlist, newsecondsidecoord, coordarray, box_dimension, my_nonbonded_params);
	if(secondsidenodetype==1) difference+=calc_phenyl_energy_cellstruct(secondsite, -1, -1, secondside_longrange_number_neighbors, secondside_longrange_neighborlist, newsecondsidecoord, coordarray, box_dimension, my_nonbonded_params);
	else difference+=calc_charged_energy_cellstruct(secondsite, -1, -1, secondside_longrange_number_neighbors, secondside_longrange_neighborlist, newsecondsidecoord, coordarray, box_dimension, my_nonbonded_params);
    free(secondside_longrange_neighborlist);
	
	//	new end sidechain
	
    if(debug==1) printf("\tnew end sidechain (difference %f)...\n", difference);
    int endsidesite=monomerid[movingpolymer][endindex].sidechain;
	int endsidenodetype=coordarray[endsidesite].nodetype;
	coord newendsidecoord=side_reptate(newendcoord, coordarray[endsidesite], sidereptation[endsidenodetype], box_dimension);
    declare_array(int, endside_shortrange_neighborlist, (*plinkset).shortrange.maxneighbors);
	int *endside_longrange_neighborlist;
    if(endsidenodetype==1){
		allocate_array(int, endside_longrange_neighborlist, (*plinkset).phenylphenyl.maxneighbors);
	}
	else if((endsidenodetype==2)||(endsidenodetype==3)){
		allocate_array(int, endside_longrange_neighborlist, (*plinkset).chargedcharged.maxneighbors);
	}
	else my_exit("code not written for other nodetypes");
	int_triple newendsideshortrangecell=cell(newendsidecoord.r, (*plinkset).shortrange.cellwidth);
	int_triple newendsidelongrangecell;
	if(endsidenodetype==1) newendsidelongrangecell=cell(newendsidecoord.r, (*plinkset).phenylphenyl.cellwidth);
	else newendsidelongrangecell=cell(newendsidecoord.r, (*plinkset).chargedcharged.cellwidth);
	create_temporary_neighborlist((*plinkset).shortrange, coordarray, box_dimension, &endside_shortrange_number_neighbors, endside_shortrange_neighborlist, newendsideshortrangecell, newendsidecoord);
	if(endsidenodetype==1){
		create_temporary_neighborlist((*plinkset).phenylphenyl, coordarray, box_dimension, &endside_longrange_number_neighbors, endside_longrange_neighborlist, newendsidelongrangecell, newendsidecoord);
	}
	else{
		create_temporary_neighborlist((*plinkset).chargedcharged, coordarray, box_dimension, &endside_longrange_number_neighbors, endside_longrange_neighborlist, newendsidelongrangecell, newendsidecoord);
	}
	difference+=calc_hard_energy_cellstruct(endsite, -1, -1, endside_shortrange_number_neighbors, endside_shortrange_neighborlist, newendsidecoord, coordarray, box_dimension, my_nonbonded_params);
	if(endsidenodetype==1) difference+=calc_phenyl_energy_cellstruct(endsite, -1, -1, endside_longrange_number_neighbors, endside_longrange_neighborlist, newendsidecoord, coordarray, box_dimension, my_nonbonded_params);
	else difference+=calc_charged_energy_cellstruct(endsite, -1, -1, endside_longrange_number_neighbors, endside_longrange_neighborlist, newendsidecoord, coordarray, box_dimension, my_nonbonded_params);
    free(endside_longrange_neighborlist);
    if(debug==1) printf("\tdifference %f\n", difference);
    
	//	add nonbonded interactions among new sites
    
	if(newsecondsidecoord.nodetype==1){
		difference+=hard_energy_sumradii(newendcoord.r, newsecondsidecoord.r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
	}
	else if(newsecondsidecoord.nodetype==2){
		difference+=hard_energy_sumradii(newendcoord.r, newsecondsidecoord.r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
	}
	else if(newsecondsidecoord.nodetype==3){
		difference+=hard_energy_sumradii(newendcoord.r, newsecondsidecoord.r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
	}
	if(newendsidecoord.nodetype==1){
		difference+=hard_energy_sumradii(newendsidecoord.r, newsecondcoord.r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
		if(newsecondsidecoord.nodetype==1){
			//difference+=phenylphenyl_energy(newendsidecoord, newsecondcoord, my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2);
			
			difference+=phenylphenyl_energy(newendsidecoord, newsecondcoord, my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		
		}
		else if(newsecondsidecoord.nodetype==2){
			difference+=hard_energy_sumradii(newendsidecoord.r, newsecondsidecoord.r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
		}
		else if(newsecondsidecoord.nodetype==3){
			difference+=hard_energy_sumradii(newendsidecoord.r, newsecondsidecoord.r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
		}
	}
	else if(newendsidecoord.nodetype==2){
		difference+=hard_energy_sumradii(newendsidecoord.r, newsecondcoord.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
		if(newsecondsidecoord.nodetype==1){
			difference+=hard_energy_sumradii(newendsidecoord.r, newsecondcoord.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(newsecondsidecoord.nodetype==2){
			difference+=electrostatic_energy(newendsidecoord, newsecondcoord, my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
		}
		else if(newsecondsidecoord.nodetype==3){
			difference+=electrostatic_energy(newendsidecoord, newsecondcoord, my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
		}
	}
	else if(newendsidecoord.nodetype==3){
		difference+=hard_energy_sumradii(newendsidecoord.r, newsecondcoord.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
		if(newsecondsidecoord.nodetype==1){
			difference+=hard_energy_sumradii(newendsidecoord.r, newsecondcoord.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(newsecondsidecoord.nodetype==2){
			difference+=electrostatic_energy(newendsidecoord, newsecondcoord, my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
		}
		else if(newsecondsidecoord.nodetype==3){
			difference+=electrostatic_energy(newendsidecoord, newsecondcoord, my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
		}
	}
    if(debug==1) printf("\tafter new site nonbonded, difference %f\n", difference);
	
	//	calculate proposed deleted nonbonded interactions; using that backbone sites and cross terms are hard-core
	
	if(reverseendsidenodetype==1) difference-=calc_phenyl_energy_cellstruct(reverseendsite, -1, -1, (*plinkset).phenylphenyl.core.number_neighbors[reverseendsidesite], (*plinkset).phenylphenyl.core.neighbor[reverseendsidesite], coordarray[reverseendsidesite], coordarray, box_dimension, my_nonbonded_params);
	else if((reverseendsidenodetype==2)||(reverseendsidenodetype==3)) difference-=calc_charged_energy_cellstruct(reverseendsite, -1, -1, (*plinkset).chargedcharged.core.number_neighbors[reverseendsidesite], (*plinkset).chargedcharged.core.neighbor[reverseendsidesite], coordarray[reverseendsidesite], coordarray, box_dimension, my_nonbonded_params);
	else my_exit("code not written for other nodetypes");
    
	if(reversesecondsidenodetype==1) difference-=calc_phenyl_energy_cellstruct(reversesecondsite, -1, -1, (*plinkset).phenylphenyl.core.number_neighbors[reversesecondsidesite], (*plinkset).phenylphenyl.core.neighbor[reversesecondsidesite], coordarray[reversesecondsidesite], coordarray, box_dimension, my_nonbonded_params);
	else if((reversesecondsidenodetype==2)||(reversesecondsidenodetype==3)) difference-=calc_charged_energy_cellstruct(reversesecondsite, -1, -1, (*plinkset).chargedcharged.core.number_neighbors[reversesecondsidesite], (*plinkset).chargedcharged.core.neighbor[reversesecondsidesite], coordarray[reversesecondsidesite], coordarray, box_dimension, my_nonbonded_params);
	else my_exit("code not written for other nodetypes");
    if(debug==1) printf("\tafter delete old site nonbonded, difference %f\n", difference);
	
	//	change in onebody interactions
	
	difference+=onebody_energy(newsecondcoord, chooseonebody(newsecondcoord.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	difference+=onebody_energy(newendcoord, chooseonebody(newendcoord.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	difference+=onebody_energy(newsecondsidecoord, chooseonebody(newsecondsidecoord.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	difference+=onebody_energy(newendsidecoord, chooseonebody(newendsidecoord.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
    
	difference-=onebody_energy(coordarray[reversesecondsite], chooseonebody(coordarray[reversesecondsite].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	difference-=onebody_energy(coordarray[reverseendsite], chooseonebody(coordarray[reverseendsite].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	difference-=onebody_energy(coordarray[reversesecondsidesite], chooseonebody(coordarray[reversesecondsidesite].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	difference-=onebody_energy(coordarray[reverseendsidesite], chooseonebody(coordarray[reverseendsidesite].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
    if(debug==1) printf("\tafter onebody, difference %f\n", difference);
	
	//	change in bonded interactions
	
	if(direction==1){
		difference+=calc_backbone_bonded_energy(my_bonded_params, newsecondcoord, coordarray[endsite], newendcoord, box_dimension);
		difference+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[endsite], coordarray[secondsite], newsecondcoord, box_dimension);
		difference-=calc_backbone_bonded_energy(my_bonded_params, coordarray[reversesecondsite], coordarray[reverseendsite], coordarray[reversethirdsite], box_dimension);
		difference-=calc_backbone_bonded_energy_noriip1(my_bonded_params, coordarray[reversethirdsite], coordarray[reversesecondsite], coordarray[reversefourthsite], box_dimension);
	}
	else{
		difference+=calc_backbone_bonded_energy(my_bonded_params, newsecondcoord, newendcoord, coordarray[endsite], box_dimension);
		difference+=calc_backbone_bonded_energy_noriip1(my_bonded_params, coordarray[endsite], newsecondcoord, coordarray[secondsite], box_dimension);
		difference-=calc_backbone_bonded_energy(my_bonded_params, coordarray[reversesecondsite], coordarray[reversethirdsite], coordarray[reverseendsite], box_dimension);
		difference-=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[reversethirdsite], coordarray[reversefourthsite], coordarray[reversesecondsite], box_dimension);
	}
	if(newendsidecoord.nodetype>=0) difference+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[newendsidecoord.nodetype], newendcoord, newendsidecoord, box_dimension);
	if(newsecondsidecoord.nodetype>=0) difference+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[newsecondsidecoord.nodetype], newsecondcoord, newsecondsidecoord, box_dimension);
	if(coordarray[reverseendsidesite].nodetype>=0) difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[reverseendsidesite].nodetype], coordarray[reverseendsite], coordarray[reverseendsidesite], box_dimension);
	if(coordarray[reversesecondsidesite].nodetype>=0) difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[reversesecondsidesite].nodetype], coordarray[reversesecondsite], coordarray[reversesecondsidesite], box_dimension);
    if(debug==1) printf("\tafter bonded, difference %f\n", difference);
    
	//	change in interactions for middle sidechain that switches identities (assuming single identity for backbone sites)
	
    int oldtype, newtype, oldindex, startsearch, endsearch, newindex;
    coord changingcoord, changingbackbonecoord;
    if(direction==1){
        startsearch=2;
        endsearch=chainlength[movingpolymer]-1;
    }
    else{
        startsearch=0;
        endsearch=chainlength[movingpolymer]-3;
    }
    for(i=startsearch;i<=endsearch;i++){          //  check all of them, not assuming it is in block configuration
        oldindex=monomerid[movingpolymer][i].sidechain;
        oldtype=coordarray[oldindex].nodetype;
        newtype=coordarray[monomerid[movingpolymer][i-2*direction].sidechain].nodetype;
        if(oldtype!=newtype){
            changingcoord=coordarray[oldindex];
            changingbackbonecoord=coordarray[monomerid[movingpolymer][i].backbone];
            difference+=calc_charged_energy_difference_swapidentity_cellstruct(monomerid[movingpolymer][i].backbone, -1, -1, (*plinkset).chargedcharged.core.number_neighbors[oldindex], (*plinkset).chargedcharged.core.neighbor[oldindex], changingcoord, coordarray, box_dimension, my_nonbonded_params, newtype, oldtype);
            difference+=onebody_energy(changingcoord, chooseonebody(newtype, my_nonbonded_params), my_nonbonded_params.solvationparams);
            difference-=onebody_energy(changingcoord, chooseonebody(oldtype, my_nonbonded_params), my_nonbonded_params.solvationparams);
            if(newtype>=0) difference+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[newtype], changingbackbonecoord, changingcoord, box_dimension);
            if(oldtype>=0) difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[oldtype], changingbackbonecoord, changingcoord, box_dimension);
		}
    }
    if(debug==1) printf("\tafter middle identity change, difference %f\n", difference);
 	
	double acceptance_prob=exp(-difference/temperature);
    if(debug==1) printf("prop=%f=exp(-%f/%f)\n", acceptance_prob, difference, temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        
        printf("accepted\n");
        
        //  shift interior coord, cell, and neighbor data
        
        if(direction==1){
            for(i=0;i<chainlength[movingpolymer]-2;i++){
                oldindex=monomerid[movingpolymer][i+2].backbone;
                newindex=monomerid[movingpolymer][i].backbone;
                coordarray[newindex]=coordarray[oldindex];
                shiftcell(&((*plinkset).shortrange), newindex, oldindex);
                
                oldindex=monomerid[movingpolymer][i+2].sidechain;
                newindex=monomerid[movingpolymer][i].sidechain;
                newtype=coordarray[newindex].nodetype;
                coordarray[newindex]=coordarray[oldindex];
                coordarray[newindex].nodetype=newtype;
                shiftcell(&((*plinkset).shortrange), newindex, oldindex);
                if(newtype==1) shiftcell(&((*plinkset).phenylphenyl), newindex, oldindex);
                else shiftcell(&((*plinkset).chargedcharged), newindex, oldindex);
            }
        }
        else{
            for(i=chainlength[movingpolymer]-3;i>=0;i--){
                oldindex=monomerid[movingpolymer][i].backbone;
                newindex=monomerid[movingpolymer][i+2].backbone;
                coordarray[newindex]=coordarray[oldindex];
                shiftcell(&((*plinkset).shortrange), newindex, oldindex);
                
                oldindex=monomerid[movingpolymer][i].sidechain;
                newindex=monomerid[movingpolymer][i+2].sidechain;
                newtype=coordarray[newindex].nodetype;
                coordarray[newindex]=coordarray[oldindex];
                coordarray[newindex].nodetype=newtype;
                shiftcell(&((*plinkset).shortrange), newindex, oldindex);
                if(newtype==1) shiftcell(&((*plinkset).phenylphenyl), newindex, oldindex);
                else shiftcell(&((*plinkset).chargedcharged), newindex, oldindex);
            }
        }
        
        //  delete cell and neighbor data for reverse end sites
        
        deletecell(&((*plinkset).shortrange), reverseendsite);
        deletecell(&((*plinkset).shortrange), reversesecondsite);
        deletecell(&((*plinkset).shortrange), reverseendsidesite);
        if(coordarray[reverseendsidesite].nodetype==1) deletecell(&((*plinkset).phenylphenyl), reverseendsidesite);
        else if((coordarray[reverseendsidesite].nodetype==2)||(coordarray[reverseendsidesite].nodetype==3)) deletecell(&((*plinkset).chargedcharged), reverseendsidesite);
        deletecell(&((*plinkset).shortrange), reversesecondsidesite);
        if(coordarray[reversesecondsidesite].nodetype==1) deletecell(&((*plinkset).phenylphenyl), reversesecondsidesite);
        else if((coordarray[reversesecondsidesite].nodetype==2)||(coordarray[reversesecondsidesite].nodetype==3)) deletecell(&((*plinkset).chargedcharged), reversesecondsidesite);
        
        //  update growing end coord, cell, and neighbor data
        
        coordarray[endsite]=newendcoord;
        coordarray[secondsite]=newsecondcoord;
        coordarray[endsidesite]=newendsidecoord;
        coordarray[secondsidesite]=newsecondsidecoord;
        
        updatecellnocheck(newendcoord.r, &((*plinkset).shortrange), endsite);
        updatecellnocheck(newsecondcoord.r, &((*plinkset).shortrange), secondsite);
        updatecellnocheck(newendsidecoord.r, &((*plinkset).shortrange), endsidesite);
        if(newendsidecoord.nodetype==1) updatecellnocheck(newendsidecoord.r, &((*plinkset).phenylphenyl), endsidesite);
        else if((newendsidecoord.nodetype==2)||(newendsidecoord.nodetype==3)) updatecellnocheck(newendsidecoord.r, &((*plinkset).chargedcharged), endsidesite);
        updatecellnocheck(newsecondsidecoord.r, &((*plinkset).shortrange), secondsidesite);
        if(newsecondsidecoord.nodetype==1) updatecellnocheck(newsecondsidecoord.r, &((*plinkset).phenylphenyl), secondsidesite);
        else if((newsecondsidecoord.nodetype==2)||(newsecondsidecoord.nodetype==3)) updatecellnocheck(newsecondsidecoord.r, &((*plinkset).chargedcharged), secondsidesite);
        
        updateneighborlist_nowipeout(&((*plinkset).shortrange), endsite, coordarray, box_dimension);
        updateneighborlist_nowipeout(&((*plinkset).shortrange), secondsite, coordarray, box_dimension);
        updateneighborlist_nowipeout(&((*plinkset).shortrange), endsidesite, coordarray, box_dimension);
        if(newendsidecoord.nodetype==1) updateneighborlist_pairenergies_phenyl_nowipeout(&((*plinkset).shortrange), endsidesite, coordarray, box_dimension, my_nonbonded_params);
        else if((newendsidecoord.nodetype==2)||(newendsidecoord.nodetype==3)) updateneighborlist_pairenergies_charged_nowipeout(&((*plinkset).shortrange), endsidesite, coordarray, box_dimension, my_nonbonded_params);
        updateneighborlist_nowipeout(&((*plinkset).shortrange), secondsidesite, coordarray, box_dimension);
        if(newsecondsidecoord.nodetype==1) updateneighborlist_pairenergies_phenyl_nowipeout(&((*plinkset).shortrange), secondsidesite, coordarray, box_dimension, my_nonbonded_params);
        else if((newsecondsidecoord.nodetype==2)||(newsecondsidecoord.nodetype==3)) updateneighborlist_pairenergies_charged_nowipeout(&((*plinkset).shortrange), secondsidesite, coordarray, box_dimension, my_nonbonded_params);
	}
}

void mc_row_surgery(int Nchains, int *chainlength, coord *coordarray, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, reptationparams leftreptation, reptationparams rightreptation, int maxsurgery){
    int initialpolymer=rand_int(Nchains);
    int direction=2*rand_int(2)-1;
    reptationparams reptation, reversereptation;
    int wrapped=0, i=0, foundneighbor, j, debug=0;
    declare_array(int, polymer, maxsurgery+1);
    declare_array(int, endsite, maxsurgery+1);
    declare_array(int, secondsite, maxsurgery+1);
    declare_array(int, reverseendsite, maxsurgery+1);
    declare_array(int, reversesecondsite, maxsurgery+1);
    declare_array(int, reversethirdsite, maxsurgery+1);
    declare_array(int, reversefourthsite, maxsurgery+1);
	int last, reverseendindex, reversesecondindex, reversethirdindex, reversefourthindex, endindex, secondindex, index;
    if(direction==1){		//	right
        reptation=rightreptation;
        reversereptation=leftreptation;
    }
    else{					//	left
        reptation=leftreptation;
        reversereptation=rightreptation;
    }
    polymer[0]=initialpolymer;
    i=0;
    if(debug==1){
        if(direction==1) printf("direction %i, initial polymer %i: %i (%f) %i (%f) ... %i (%f) %i (%f)\n", direction, polymer[0], monomerid[polymer[0]][0].backbone, coordarray[monomerid[polymer[0]][0].backbone].r.x, monomerid[polymer[0]][1].backbone, coordarray[monomerid[polymer[0]][1].backbone].r.x, monomerid[polymer[0]][chainlength[polymer[0]]-2].backbone, coordarray[monomerid[polymer[0]][chainlength[polymer[0]]-2].backbone].r.x, monomerid[polymer[0]][chainlength[polymer[0]]-1].backbone, coordarray[monomerid[polymer[0]][chainlength[polymer[0]]-1].backbone].r.x);
        else printf("direction %i, initial polymer %i: %i (%f) %i (%f) ... %i (%f) %i (%f)\n", direction, polymer[0], monomerid[polymer[0]][chainlength[polymer[0]]-1].backbone, coordarray[monomerid[polymer[0]][chainlength[polymer[0]]-1].backbone].r.x, monomerid[polymer[0]][chainlength[polymer[0]]-2].backbone, coordarray[monomerid[polymer[0]][chainlength[polymer[0]]-2].backbone].r.x, monomerid[polymer[0]][1].backbone, coordarray[monomerid[polymer[0]][1].backbone].r.x, monomerid[polymer[0]][0].backbone, coordarray[monomerid[polymer[0]][0].backbone].r.x);
    }
    while((wrapped==0)&&(i<maxsurgery)){
        
        //  recalculate indices each time, in principle allowing for different lengths of polymers
        
        last=chainlength[polymer[i]]-1;
        if(coordarray[monomerid[polymer[i]][last].sidechain].nodetype>0) my_exit("mc_row_surgery not written to account for empty sidechains!");
        if(direction==1){		//	right
            reverseendindex=0;
            reversesecondindex=1;
            reversethirdindex=2;
            reversefourthindex=3;
            endindex=last;
            secondindex=last-1;
        }
        else{					//	left
            reverseendindex=last;
            reversesecondindex=last-1;
            reversethirdindex=last-2;
            reversefourthindex=last-3;
            endindex=0;
            secondindex=1;
            reptation=leftreptation;
            reversereptation=rightreptation;
        }
        
        endsite[i]=monomerid[polymer[i]][endindex].backbone;
        secondsite[i]=monomerid[polymer[i]][secondindex].backbone;
        reverseendsite[i]=monomerid[polymer[i]][reverseendindex].backbone;
        reversesecondsite[i]=monomerid[polymer[i]][reversesecondindex].backbone;
        reversethirdsite[i]=monomerid[polymer[i]][reversethirdindex].backbone;
        reversefourthsite[i]=monomerid[polymer[i]][reversefourthindex].backbone;
        
        //  reverse neighborhood criterion
        
        if(in_neighborhood(coordarray[reversefourthsite[i]], coordarray[reversethirdsite[i]], coordarray[reversesecondsite[i]], reversereptation, box_dimension)==0) {
            if(debug==1) printf("\treverse move not in neighborhood\n");
            return;
        }
        if(debug==1) printf("\treverse move in neighborhood\n");
        
        //  forward neighborhood criterion
        
        foundneighbor=0;
        if(debug==1) printf("\tsearching through %i neighbors of %i\n", (*plinkset).shortrange.core.number_neighbors[endsite[i]], endsite[i]);
        for(j=0;(j<(*plinkset).shortrange.core.number_neighbors[endsite[i]])&&(foundneighbor==0);j++){
            index=(*plinkset).shortrange.core.neighbor[endsite[i]][j];
            if(debug==1) printf("\t\tneighbor %i (type %i) (%i %i)\n", index, coordarray[index].nodetype, coordarray[index].chainid, coordarray[index].monomerid);
            if((coordarray[index].monomerid==reverseendindex)&&(coordarray[index].nodetype==0)){          //  backbone and end site
                if(in_neighborhood(coordarray[secondsite[i]], coordarray[endsite[i]], coordarray[index], reptation, box_dimension)==1) foundneighbor=1;
            }
        }
        if(foundneighbor==0){
            if(debug==1) printf("\tno neighbor found\n");
            if(debug==1) printf("\tcoordarray[%i].monomerid=%i, coordarray[%i].nodetype=%i\n", 1796, coordarray[1796].monomerid, 1796, coordarray[1796].nodetype);
            return;
        }
        polymer[i+1]=coordarray[index].chainid;
        
        if(polymer[i+1]==initialpolymer) wrapped=1;
        else{
            if(debug==1){
                if(direction==1) printf("\tpolymer %i: %i (%f) %i (%f) ... %i (%f) %i (%f)\n", polymer[i+1], monomerid[polymer[i+1]][0].backbone, coordarray[monomerid[polymer[i+1]][0].backbone].r.x, monomerid[polymer[i+1]][1].backbone, coordarray[monomerid[polymer[i+1]][1].backbone].r.x, monomerid[polymer[i+1]][chainlength[polymer[i+1]]-2].backbone, coordarray[monomerid[polymer[i+1]][chainlength[polymer[i+1]]-2].backbone].r.x, monomerid[polymer[i+1]][chainlength[polymer[i+1]]-1].backbone, coordarray[monomerid[polymer[i+1]][chainlength[polymer[i+1]]-1].backbone].r.x);
                else printf("\tpolymer %i: %i (%f) %i (%f) ... %i (%f) %i (%f)\n", polymer[i+1], monomerid[polymer[i+1]][chainlength[polymer[i+1]]-1].backbone, coordarray[monomerid[polymer[i+1]][chainlength[polymer[i+1]]-1].backbone].r.x, monomerid[polymer[i+1]][chainlength[polymer[i+1]]-2].backbone, coordarray[monomerid[polymer[i+1]][chainlength[polymer[i+1]]-2].backbone].r.x, monomerid[polymer[i+1]][1].backbone, coordarray[monomerid[polymer[i+1]][1].backbone].r.x, monomerid[polymer[i+1]][0].backbone, coordarray[monomerid[polymer[i+1]][0].backbone].r.x);
            }
            i++;
        }
    }
    double difference=0, acceptance_prob;
    int oldtype, newtype, oldindex, next, looplength=i+1, start, end, newmonomerid, oldbackindex;
    coord changingcoord, changingbackbonecoord;
    if(wrapped==1){
        if(debug==1) printf("\twrapped, looplength=%i\n", looplength);
        for(j=0;j<looplength;j++){
            if(j==looplength-1) next=0;
            else next=j+1;
            if(direction==1){
                start=chainlength[polymer[j]]-1;
                end=-1;
            }
            else{
                start=0;
                end=chainlength[polymer[j]];
            }
            
            //  change in energy (bonded and nonbonded interactions that turn on and off due to change in bond topology)
            
            if(direction==1){
                difference+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[endsite[j]], coordarray[secondsite[j]], coordarray[reverseendsite[next]], box_dimension);
                if(debug==1) printf("\t\t%i: after new bond energy, difference %f\n", j, difference);
                difference+=calc_backbone_bonded_energy_onlythreebody(my_bonded_params, coordarray[reverseendsite[next]], coordarray[endsite[j]], coordarray[reversesecondsite[next]], box_dimension);
                if(debug==1) printf("\t\t%i: after new bond energy, difference %f\n", j, difference);
                difference-=calc_backbone_bonded_energy_noriip1(my_bonded_params, coordarray[reversethirdsite[j]], coordarray[reversesecondsite[j]], coordarray[reversefourthsite[j]], box_dimension);
                if(debug==1) printf("\t\t%i: after remove old bond energy, difference %f\n", j, difference);
                difference-=calc_backbone_bonded_energy_onlythreebody(my_bonded_params, coordarray[reversesecondsite[j]], coordarray[reverseendsite[j]], coordarray[reversethirdsite[j]], box_dimension);
                if(debug==1) printf("\t\t%i: after remove old bond energy, difference %f\n", j, difference);
                difference-=calc_nonbonded_energy_pair(coordarray[endsite[j]], coordarray[reverseendsite[next]], box_dimension, my_nonbonded_params);
                if(debug==1) printf("\t\t%i: after remove old nonbond energy, difference %f\n", j, difference);
                difference+=calc_nonbonded_energy_pair(coordarray[reversesecondsite[j]], coordarray[reversethirdsite[j]], box_dimension, my_nonbonded_params);
                if(debug==1) printf("\t\t%i: after add new nonbond energy, difference %f\n", j, difference);
            }
            else{
                difference+=calc_backbone_bonded_energy_noriip1(my_bonded_params, coordarray[endsite[j]], coordarray[reverseendsite[next]], coordarray[secondsite[j]], box_dimension);
                
                if(debug==1) printf("\t\t%i: after new bond energy, difference %f\n", j, difference);
                difference+=calc_backbone_bonded_energy_onlythreebody(my_bonded_params, coordarray[reverseendsite[next]], coordarray[reversesecondsite[next]], coordarray[endsite[j]], box_dimension);
                if(debug==1) printf("\t\t%i: after new bond energy, difference %f\n", j, difference);
                difference-=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[reversethirdsite[j]], coordarray[reversefourthsite[j]], coordarray[reversesecondsite[j]], box_dimension);
                if(debug==1) printf("\t\t%i: after remove old bond energy, difference %f\n", j, difference);
                difference-=calc_backbone_bonded_energy_onlythreebody(my_bonded_params, coordarray[reversesecondsite[j]], coordarray[reversethirdsite[j]], coordarray[reverseendsite[j]], box_dimension);
                if(debug==1) printf("\t\t%i: after remove old bond energy, difference %f\n", j, difference);
                difference-=calc_nonbonded_energy_pair(coordarray[endsite[j]], coordarray[reverseendsite[next]], box_dimension, my_nonbonded_params);
                if(debug==1) printf("\t\t%i: after remove old nonbond energy, difference %f\n", j, difference);
                difference+=calc_nonbonded_energy_pair(coordarray[reversesecondsite[j]], coordarray[reversethirdsite[j]], box_dimension, my_nonbonded_params);
                if(debug==1) printf("\t\t%i: after add new nonbond energy, difference %f\n", j, difference);
            }
            if(debug==1)  printf("\t\t%i: after change bond energy, difference %f\n", j, difference);
            
            //  change in energy due to change in nodetype
            
            for(i=start;direction*i>direction*end;i-=direction){
                oldindex=monomerid[polymer[j]][i].sidechain;
                oldtype=coordarray[oldindex].nodetype;
                if(direction*i>direction*end+2){
                    newmonomerid=coordarray[oldindex].monomerid-direction*2;
                    newtype=coordarray[monomerid[polymer[j]][newmonomerid].sidechain].nodetype;
                }
                else{
                    if(direction==1) newmonomerid=chainlength[polymer[next]]-2+i;
                    else newmonomerid=i-(chainlength[polymer[j]]-2);
                    newtype=coordarray[monomerid[polymer[next]][newmonomerid].sidechain].nodetype;
                }
                if(oldtype!=newtype){
                    changingcoord=coordarray[oldindex];
                    oldbackindex=monomerid[polymer[j]][i].backbone;
                    changingbackbonecoord=coordarray[oldbackindex];
                    difference+=calc_charged_energy_difference_swapidentity_cellstruct(oldbackindex, -1, -1, (*plinkset).chargedcharged.core.number_neighbors[oldindex], (*plinkset).chargedcharged.core.neighbor[oldindex], changingcoord, coordarray, box_dimension, my_nonbonded_params, newtype, oldtype);
                    if(debug==1) printf("\t\tafter nonbonded difference %f\n", difference);
                    difference+=onebody_energy(changingcoord, chooseonebody(newtype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    difference-=onebody_energy(changingcoord, chooseonebody(oldtype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    if(debug==1) printf("\t\tafter onebody difference %f\n", difference);
                    if(newtype>=0) difference+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[newtype], changingbackbonecoord, changingcoord, box_dimension);
                    if(oldtype>=0) difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[oldtype], changingbackbonecoord, changingcoord, box_dimension);
                    if(debug==1) printf("\t\tafter sidechain difference %f\n", difference);
                }
            }
            if(debug==1) printf("\t\tafter change nodetype energy, difference %f\n", difference);
        }
        acceptance_prob=exp(-difference/temperature);
        monomernodes firstmonomerid, secondmonomerid;
        if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
            printf("\taccepted\n");
            
            if(direction==1){
                firstmonomerid=monomerid[polymer[0]][0];
                secondmonomerid=monomerid[polymer[0]][1];
            }
            else{
                firstmonomerid=monomerid[polymer[0]][chainlength[polymer[0]]-1];
                secondmonomerid=monomerid[polymer[0]][chainlength[polymer[0]]-2];
            }
            for(j=0;j<looplength;j++){
                if(j==looplength-1) next=0;
                else next=j+1;
                if(direction==1){
                    for(i=0;i<chainlength[polymer[j]]-2;i++){
                        coordarray[monomerid[polymer[j]][i+2].sidechain].nodetype=coordarray[monomerid[polymer[j]][i].sidechain].nodetype;
                        monomerid[polymer[j]][i]=monomerid[polymer[j]][i+2];
                    }
                    if(j<looplength-1){
                        for(i=chainlength[polymer[j]]-2;i<chainlength[polymer[j]];i++){
                            coordarray[monomerid[polymer[next]][i-(chainlength[polymer[j]]-2)].sidechain].nodetype=coordarray[monomerid[polymer[j]][i].sidechain].nodetype;
                            monomerid[polymer[j]][i]=monomerid[polymer[next]][i-(chainlength[polymer[j]]-2)];
                        }
                    }
                    else{
                        coordarray[firstmonomerid.sidechain].nodetype=coordarray[monomerid[polymer[j]][chainlength[polymer[j]]-2].sidechain].nodetype;
                        coordarray[secondmonomerid.sidechain].nodetype=coordarray[monomerid[polymer[j]][chainlength[polymer[j]]-1].sidechain].nodetype;
                        monomerid[polymer[j]][chainlength[polymer[j]]-2]=firstmonomerid;
                        monomerid[polymer[j]][chainlength[polymer[j]]-1]=secondmonomerid;
                    }
                }
                else{
                    for(i=chainlength[polymer[j]]-1;i>1;i--){
                        coordarray[monomerid[polymer[j]][i-2].sidechain].nodetype=coordarray[monomerid[polymer[j]][i].sidechain].nodetype;
                        monomerid[polymer[j]][i]=monomerid[polymer[j]][i-2];
                    }
                    if(j<looplength-1){
                        for(i=1;i>=0;i--){
                            coordarray[monomerid[polymer[next]][(chainlength[polymer[j]]-2)+i].sidechain].nodetype=coordarray[monomerid[polymer[j]][i].sidechain].nodetype;
                            monomerid[polymer[j]][i]=monomerid[polymer[next]][(chainlength[polymer[j]]-2)+i];
                        }
                    }
                    else{
                        coordarray[firstmonomerid.sidechain].nodetype=coordarray[monomerid[polymer[j]][1].sidechain].nodetype;
                        coordarray[secondmonomerid.sidechain].nodetype=coordarray[monomerid[polymer[j]][0].sidechain].nodetype;
                        monomerid[polymer[j]][1]=firstmonomerid;
                        monomerid[polymer[j]][0]=secondmonomerid;
                    }
                }
            }
            for(j=0;j<looplength;j++){
                for(i=0;i<chainlength[polymer[j]];i++){
                    coordarray[monomerid[polymer[j]][i].backbone].monomerid=i;
                    coordarray[monomerid[polymer[j]][i].sidechain].monomerid=i;
                    coordarray[monomerid[polymer[j]][i].backbone].chainid=polymer[j];
                    coordarray[monomerid[polymer[j]][i].sidechain].chainid=polymer[j];
                }
            }
            if(debug==1){
                printf("\tshifted polymers are now:\n");
                for(j=0;j<looplength;j++){
                    printf("\t\t%i: %i (%f) %i (%f) ... %i (%f) %i (%f)\n", polymer[j], monomerid[polymer[j]][0].backbone, coordarray[monomerid[polymer[j]][0].backbone].r.x, monomerid[polymer[j]][1].backbone, coordarray[monomerid[polymer[j]][1].backbone].r.x, monomerid[polymer[j]][chainlength[polymer[j]]-2].backbone, coordarray[monomerid[polymer[j]][chainlength[polymer[j]]-2].backbone].r.x, monomerid[polymer[j]][chainlength[polymer[j]]-1].backbone, coordarray[monomerid[polymer[j]][chainlength[polymer[j]]-1].backbone].r.x);
                }
                printf("\tcoordarray[%i].monomerid=%i, coordarray[%i].nodetype=%i\n", 1796, coordarray[1796].monomerid, 1796, coordarray[1796].nodetype);
            }
        }
        else if(debug==1) printf("\trejected (prob %f)\n", acceptance_prob);
    }
    else if(debug==1) printf("no wrap\n");
    free(polymer);
    free(endsite);
    free(secondsite);
    free(reverseendsite);
    free(reversesecondsite);
    free(reversethirdsite);
    free(reversefourthsite);
}

void calc_cellstruct_leafid(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double *avheight, double_triple *pavboxdimension){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	double height;
	(*pavboxdimension).x+=box_dimension.x;
	(*pavboxdimension).y+=box_dimension.y;
	(*pavboxdimension).z+=box_dimension.z;
	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){
				
                backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				
				//	using old sidechainid:
				//	backbone only has short-ranged interactions (zero)
                
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1){
                        height=coordarray[backboneid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[backboneid].nodetype]+=height;
                        height=coordarray[sidechainid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[sidechainid].nodetype]+=height;
                    }
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1){
						calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    }
					if(coordarray[sidechainid].nodetype>1){
						calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
					}
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				else if(j>=3){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions (zero)
            
		}
		else if(chainlength[i]==2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1){
                        height=coordarray[backboneid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[backboneid].nodetype]+=height;
                        height=coordarray[sidechainid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[sidechainid].nodetype]+=height;
                    }
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==1){
					
					//	backbone only has short-ranged interactions (zero)
                    
					(*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_onlyleft(my_bonded_params, coordarray[backboneid], coordarray[left], box_dimension);
				}
				left=backboneid;
			}
			
			//	backbone only has short-ranged interactions (zero)
            
		}
		else{
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1){
                        height=coordarray[backboneid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[backboneid].nodetype]+=height;
                        height=coordarray[sidechainid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[sidechainid].nodetype]+=height;
                    }
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
			}
        }
	}
}

void calc_nonbonded_energy_components_cellstruct_leafid(int n, int neighbor1, int neighbor2, int neighbor3, linkedlist link, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, energycomponents *pmy_energycomponents){
	int i, j, k, i_image, j_image, k_image, index;
	coord coord_image=coordarray[n];
	for(i=-1;i<=1;i++){																	//	loop through neighboring cells
		if((link.cell[n].x==0)&&(i==-1)){
			i_image=link.cellsperside.x-1;
			coord_image.r.x=coordarray[n].r.x+box_dimension.x;
		}
		else if((link.cell[n].x==link.cellsperside.x-1)&&(i==1)){
			i_image=0;
			coord_image.r.x=coordarray[n].r.x-box_dimension.x;
		}
		else{
            coord_image.r.x=coordarray[n].r.x;
            i_image=link.cell[n].x+i;
        }
		for(j=-1;j<=1;j++){
			if((link.cell[n].y==0)&&(j==-1)){
				j_image=link.cellsperside.y-1;
				coord_image.r.y=coordarray[n].r.y+box_dimension.y;
			}
			else if((link.cell[n].y==link.cellsperside.y-1)&&(j==1)){
				j_image=0;
				coord_image.r.y=coordarray[n].r.y-box_dimension.y;
			}
			else{
                coord_image.r.y=coordarray[n].r.y;
                j_image=link.cell[n].y+j;
            }
			for(k=-1;k<=1;k++){
				if((link.cell[n].z==0)&&(k==-1)){
					k_image=link.cellsperside.z-1;
					coord_image.r.z=coordarray[n].r.z+box_dimension.z;
				}
				else if((link.cell[n].z==link.cellsperside.z-1)&&(k==1)){
					k_image=0;
					coord_image.r.z=coordarray[n].r.z-box_dimension.z;
				}
				else{
                    coord_image.r.z=coordarray[n].r.z;
                    k_image=link.cell[n].z+k;
                }
				index=link.head[i_image][j_image][k_image];
				while(index>=0){
					if((index!=n)&&(index!=neighbor1)&&(index!=neighbor2)&&(index!=neighbor3)){					//	don't calculate nonbonded energy for bonded neighbors
                        
						calc_nonbonded_energy_pair_components_leafid(coord_image, coordarray[index], box_dimension, my_nonbonded_params, pmy_energycomponents);                                   
					}
					index=link.list[index];
				}
			}
		}
	}
}

void calc_nonbonded_energy_pair_components_leafid(coord a, coord b, double_triple box_dimension, nonbonded_params params, energycomponents *pmy_energycomponents){
    if(a.leafid==b.leafid){
		if(a.chainid==b.chainid){
			if(a.nodetype==0){
				if(b.nodetype==0){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==1){
				if(b.nodetype==0){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					
					(*pmy_energycomponents).nnsamepoly+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
					
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==2){
				if(b.nodetype==0){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cclikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).ccunlikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==3){
				if(b.nodetype==0){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).ccunlikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cclikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else my_exit("don't have energy functions for these node types!");
		}
		else{
			if(a.nodetype==0){
				if(b.nodetype==0){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==1){
				if(b.nodetype==0){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nn+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
					
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==2){
				if(b.nodetype==0){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cclike+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).ccunlike+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==3){
				if(b.nodetype==0){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).ccunlike+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cclike+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else my_exit("don't have energy functions for these node types!");
		}
    }
    else if((a.leafid)/2==(b.leafid)/2){
        if(a.nodetype==0){
            if(b.nodetype==0){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==1){
            if(b.nodetype==0){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                
				(*pmy_energycomponents).nncross+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
                
				return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==2){
            if(b.nodetype==0){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cclikecross+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).ccunlikecross+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==3){
            if(b.nodetype==0){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).ccunlikecross+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cclikecross+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else my_exit("don't have energy functions for these node types!");
    }
    else{
        if(a.nodetype==0){
            if(b.nodetype==0){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==1){
            if(b.nodetype==0){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                
				(*pmy_energycomponents).nndifferent+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
                
				return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==2){
            if(b.nodetype==0){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cclikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).ccunlikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==3){
            if(b.nodetype==0){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).ccunlikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cclikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else my_exit("don't have energy functions for these node types!");
    }
}

void output_timeseries(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double *avheight, double_triple *pavboxdimension, int runningtime){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers);
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
	for(i=0;i<Nnodetypes;i++){
		fprintf(outp, "\t%.8f", avheight[i]*factor/monomercount.type[i]);
        avheight[i]=0;
	}
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
}

void output_timeseries_leaves(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double *avheight, double_triple *pavboxdimension, int runningtime){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikecross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikecross*factor/monomercount.monomers);
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
	for(i=0;i<Nnodetypes;i++){
		fprintf(outp, "\t%.8f", avheight[i]*factor/monomercount.type[i]);
        avheight[i]=0;
	}
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
    
	(*pmy_energycomponents).nncross=0;
	(*pmy_energycomponents).ccunlikecross=0;
	(*pmy_energycomponents).cclikecross=0;
	(*pmy_energycomponents).cncross=0;

	(*pmy_energycomponents).nndifferent=0;
	(*pmy_energycomponents).ccunlikedifferent=0;
	(*pmy_energycomponents).cclikedifferent=0;
	(*pmy_energycomponents).cndifferent=0;
	
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
}

void output_timeseries_sheets(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double *avheight, double_triple *pavboxdimension, int runningtime){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikecross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikecross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nndifferent*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cndifferent*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikedifferent*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikedifferent*factor/monomercount.monomers);
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
	for(i=0;i<Nnodetypes;i++){
		fprintf(outp, "\t%.8f", avheight[i]*factor/monomercount.type[i]);
        avheight[i]=0;
	}
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
    
	(*pmy_energycomponents).nncross=0;
	(*pmy_energycomponents).ccunlikecross=0;
	(*pmy_energycomponents).cclikecross=0;
	(*pmy_energycomponents).cncross=0;

	(*pmy_energycomponents).nndifferent=0;
	(*pmy_energycomponents).ccunlikedifferent=0;
	(*pmy_energycomponents).cclikedifferent=0;
	(*pmy_energycomponents).cndifferent=0;
	
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;	
}

void output_trajectory(int Nchains, int *chainlength, coord *coordarray, monomernodes **monomerid, double_triple box_dimension, char *filename, int *pframe, int cycle){
	int i, j;
    FILE *outp;
    if((*pframe)==0) outp=fopen(filename, "w");
    else outp=fopen(filename, "a");
    for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength[i];j++){
            fprintf(outp, "%i %i %i %i %.6f %.6f %.6f %.6f %.6f %.6f\n", coordarray[monomerid[i][j].backbone].leafid, i, j, coordarray[monomerid[i][j].backbone].nodetype, coordarray[monomerid[i][j].backbone].r.x, coordarray[monomerid[i][j].backbone].r.y, coordarray[monomerid[i][j].backbone].r.z, coordarray[monomerid[i][j].backbone].n.x, coordarray[monomerid[i][j].backbone].n.y, coordarray[monomerid[i][j].backbone].n.z);
            if(coordarray[monomerid[i][j].sidechain].nodetype>=0) fprintf(outp, "%i %i %i %i %.6f %.6f %.6f %.6f %.6f %.6f\n", coordarray[monomerid[i][j].sidechain].leafid, i, j, coordarray[monomerid[i][j].sidechain].nodetype, coordarray[monomerid[i][j].sidechain].r.x, coordarray[monomerid[i][j].sidechain].r.y, coordarray[monomerid[i][j].sidechain].r.z, coordarray[monomerid[i][j].sidechain].n.x, coordarray[monomerid[i][j].sidechain].n.y, coordarray[monomerid[i][j].sidechain].n.z);
            else fprintf(outp, "%i %i %i %i 0 0 0 0 0 0\n", coordarray[monomerid[i][j].sidechain].leafid, i, j, coordarray[monomerid[i][j].sidechain].nodetype);
        }
    }
    fprintf(outp, "%i box %.6f %.6f %.6f\n\n", cycle, box_dimension.x, box_dimension.y, box_dimension.z);
    fclose(outp);
    (*pframe)++;
}

void calc_displacement(coord left, coord right, displacementparams *pdisplacement, double boxwidth){
    double sep=right.r.x-left.r.x;
	double dif=sep-(*pdisplacement).lastsep;
	recenter(dif, boxwidth);
	(*pdisplacement).lastsep+=dif;
	sep=(*pdisplacement).lastsep;
    (*pdisplacement).av+=sep;
    (*pdisplacement).count++;
}

void output_displacement_timeseries(char *filename, int cycle, displacementparams *pdisplacement){
	(*pdisplacement).av/=(1.*(*pdisplacement).count);
	FILE *outp;
	outp=fopen(filename, "a");
	fprintf(outp, "%i\t%.8f\n", cycle, (*pdisplacement).av);
	(*pdisplacement).av=(*pdisplacement).count=0;
	fclose(outp);
}

double calc_sidechain_bonded_energy(sidechain_params sidechain, coord backbonecoord, coord sidechaincoord, double_triple box_dimension){
 	double_triple sep=subtract_double_triple(sidechaincoord.r, backbonecoord.r);
	recenter_double_triple(&sep, box_dimension);
	double norm2=dot_product(sep, sep);
	double rparallel=dot_product(sep, backbonecoord.n);
	double rperp=sqrt(norm2-rparallel*rparallel);
    double npar=dot_product(sidechaincoord.n, backbonecoord.n);
	double rdotn=dot_product(sep, sidechaincoord.n);
    
	double energy=sidechain.k1*pow(sqrt(pow(rperp-sidechain.rperp0, 2)+pow(rparallel-sidechain.rpar0, 2))-sidechain.r0, 2);
	double arctan=atan(rparallel/rperp);
	energy+=(sidechain.J10+arctan*(sidechain.J11+arctan*(sidechain.J12+arctan*(sidechain.J13+arctan*sidechain.J14))));
	if(sidechaincoord.orientationtype==0){              //  parallel
		energy+=(sidechain.k2*pow(npar-(sidechain.r20+rparallel*(sidechain.r21+rparallel*sidechain.r22)), 2));
		energy+=(sidechain.k3*pow(rdotn-(sidechain.r30+npar*(sidechain.r31+npar*sidechain.r32)), 2));
	}
	if(sidechaincoord.orientationtype==1){              //  perpendicular
        energy+=(rparallel*(sidechain.J201+rparallel*(sidechain.J202+rparallel*(sidechain.J203+rparallel*(sidechain.J204+rparallel*(sidechain.J205+rparallel*sidechain.J206))))));
		double npar2=npar*npar;
		double npar4=npar2*npar2;
		energy+=(npar2*(sidechain.J220+rparallel*(sidechain.J221+rparallel*sidechain.J222)));
		energy+=(npar4*sidechain.J240);
		
		//	J30x are redundant, 300, 302, 304 added to 200, 220, and 240, respectively
		
		energy+=(rdotn*npar*(sidechain.J311+npar2*sidechain.J313));
		double rdotnpower=rdotn*rdotn;
		energy+=(rdotnpower*(sidechain.J320+npar2*sidechain.J322));
		rdotnpower*=rdotn;
		energy+=(rdotnpower*npar*sidechain.J331);
		rdotnpower*=rdotn;
		energy+=(rdotnpower*sidechain.J340);
	}
	return energy;
}

double scalar_triple_product(double_triple a, double_triple b, double_triple c){
	return a.x*b.y*c.z+a.y*b.z*c.x+a.z*b.x*c.y-a.x*b.z*c.y-a.z*b.y*c.x-a.y*b.x*c.z;
}

double calc_backbone_bonded_energy(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension){
	double_triple riip1=subtract_double_triple(rightcoord.r, centercoord.r);
	recenter_double_triple(&riip1, box_dimension);
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double riip1norm=norm(riip1);
	double rim1inorm=norm(rim1i);
	double_triple rhatiip1=scalar_multiply_double_triple(riip1, 1./riip1norm);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+riip1norm*(my_bonded_params.K11+riip1norm*(my_bonded_params.K12+riip1norm*(my_bonded_params.K13+riip1norm*my_bonded_params.K14)));	//	onlyright
    
    energy+=(my_bonded_params.K10+rim1inorm*(my_bonded_params.K11+rim1inorm*(my_bonded_params.K12+rim1inorm*(my_bonded_params.K13+rim1inorm*my_bonded_params.K14))));	//	onlyleft
	
	double riip1dotni=dot_product(rhatiip1, centercoord.n);
	double riip1dotnip1=dot_product(rhatiip1, rightcoord.n);
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	double rim1idotnim1=dot_product(rhatim1i, leftcoord.n);

	energy+=(my_bonded_params.kr*pow(riip1dotni-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));			//	onlyright, right interaction
	energy+=(my_bonded_params.kl*pow(riip1dotnip1-(riip1norm-my_bonded_params.r0l)/my_bonded_params.sl, 2));		//	onlyright, left interaction
	energy+=(my_bonded_params.kl*pow(rim1idotni-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));			//	onlyright, left interaction
	energy+=(my_bonded_params.kr*pow(rim1idotnim1-(rim1inorm-my_bonded_params.r0r)/my_bonded_params.sr, 2));		//	onlyright, right interaction
	
	double dotcrossproduct=dot_product(rhatim1i, rhatiip1)-rim1idotni*riip1dotni;
	
	energy+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));

	return my_bonded_params.factor*energy;
}

double calc_backbone_bonded_energy_onlythreebody(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension){
	double_triple riip1=subtract_double_triple(rightcoord.r, centercoord.r);
	recenter_double_triple(&riip1, box_dimension);
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double riip1norm=norm(riip1);
	double rim1inorm=norm(rim1i);
	double_triple rhatiip1=scalar_multiply_double_triple(riip1, 1./riip1norm);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	double riip1dotni=dot_product(rhatiip1, centercoord.n);
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	double dotcrossproduct=dot_product(rhatim1i, rhatiip1)-rim1idotni*riip1dotni;
	
	double energy=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
	
	return energy;
}

double calc_backbone_bonded_energy_norim1i(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension){
	double_triple riip1=subtract_double_triple(rightcoord.r, centercoord.r);
	recenter_double_triple(&riip1, box_dimension);
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double riip1norm=norm(riip1);
	double rim1inorm=norm(rim1i);
	double_triple rhatiip1=scalar_multiply_double_triple(riip1, 1./riip1norm);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+riip1norm*(my_bonded_params.K11+riip1norm*(my_bonded_params.K12+riip1norm*(my_bonded_params.K13+riip1norm*my_bonded_params.K14)));	//	onlyright
    	
	double riip1dotni=dot_product(rhatiip1, centercoord.n);
	double riip1dotnip1=dot_product(rhatiip1, rightcoord.n);
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	
	energy+=(my_bonded_params.kr*pow(riip1dotni-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));			//	onlyright, right interaction
	
    energy+=(my_bonded_params.kl*pow(riip1dotnip1-(riip1norm-my_bonded_params.r0l)/my_bonded_params.sl, 2));		//	onlyright, left interaction
		
	double dotcrossproduct=dot_product(rhatim1i, rhatiip1)-rim1idotni*riip1dotni;
	
	energy+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
	
	return my_bonded_params.factor*energy;
}

double calc_backbone_bonded_energy_noriip1(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension){
	double_triple riip1=subtract_double_triple(rightcoord.r, centercoord.r);
	recenter_double_triple(&riip1, box_dimension);
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double riip1norm=norm(riip1);
	double rim1inorm=norm(rim1i);
	double_triple rhatiip1=scalar_multiply_double_triple(riip1, 1./riip1norm);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+rim1inorm*(my_bonded_params.K11+rim1inorm*(my_bonded_params.K12+rim1inorm*(my_bonded_params.K13+rim1inorm*my_bonded_params.K14)));	//	onlyleft
	
 	double riip1dotni=dot_product(rhatiip1, centercoord.n);
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	double rim1idotnim1=dot_product(rhatim1i, leftcoord.n);
	
	energy+=(my_bonded_params.kl*pow(rim1idotni-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));			//	onlyright, left interaction
	energy+=(my_bonded_params.kr*pow(rim1idotnim1-(rim1inorm-my_bonded_params.r0r)/my_bonded_params.sr, 2));		//	onlyright, right interaction
	
	double dotcrossproduct=dot_product(rhatim1i, rhatiip1)-rim1idotni*riip1dotni;
	
	energy+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
	
	return my_bonded_params.factor*energy;
}

double calc_backbone_bonded_energy_onlyleft(bonded_params my_bonded_params, coord centercoord, coord leftcoord, double_triple box_dimension){
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double rim1inorm=norm(rim1i);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+rim1inorm*(my_bonded_params.K11+rim1inorm*(my_bonded_params.K12+rim1inorm*(my_bonded_params.K13+rim1inorm*my_bonded_params.K14)));	//	onlyleft
	
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	double rim1idotnim1=dot_product(rhatim1i, leftcoord.n);

	energy+=(my_bonded_params.kl*pow(rim1idotni-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));			//	onlyleft, left interaction
	energy+=(my_bonded_params.kr*pow(rim1idotnim1-(rim1inorm-my_bonded_params.r0r)/my_bonded_params.sr, 2));		//	onlyleft, right interaction

	return my_bonded_params.factor*energy;
}

double calc_energy_difference_changer_hard(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, twoleftneighbor=-1, leftneighbor=-1, rightneighbor=-1, tworightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
			difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
			difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0){
			leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
			if(movermonomer>1) twoleftneighbor=monomerid[moverchain][movermonomer-2].backbone;
		}
		if(movermonomer<chainlength[moverchain]-1){
			rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
			if(movermonomer<chainlength[moverchain]-2) tworightneighbor=monomerid[moverchain][movermonomer+2].backbone;
		}
		difference+=calc_backbone_bonded_difference_fourneighbors_changer(twoleftneighbor, leftneighbor, rightneighbor, tworightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_hard_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
		difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
		difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_hard_energy_cellstruct(neighbor, -1, -1, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changer_hard)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return my_bonded_params.factor*difference;
}

double calc_energy_difference_changer_phenyl(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link1, linkedlist link2){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, twoleftneighbor=-1, leftneighbor=-1, rightneighbor=-1, tworightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0){
			leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
			if(movermonomer>1) twoleftneighbor=monomerid[moverchain][movermonomer-2].backbone;
		}
		if(movermonomer<chainlength[moverchain]-1){
			rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
			if(movermonomer<chainlength[moverchain]-2) tworightneighbor=monomerid[moverchain][movermonomer+2].backbone;
		}
		difference+=calc_backbone_bonded_difference_fourneighbors_changer(twoleftneighbor, leftneighbor, rightneighbor, tworightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_hard_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_phenyl_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_hard_energy_cellstruct(neighbor, -1, -1, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
 		difference+=calc_phenyl_energy_cellstruct(neighbor, -1, -1, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changer_phenyl)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changer_charged(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link1, linkedlist link2){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, twoleftneighbor=-1, leftneighbor=-1, rightneighbor=-1, tworightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0){
			leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
			if(movermonomer>1) twoleftneighbor=monomerid[moverchain][movermonomer-2].backbone;
		}
		if(movermonomer<chainlength[moverchain]-1){
			rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
			if(movermonomer<chainlength[moverchain]-2) tworightneighbor=monomerid[moverchain][movermonomer+2].backbone;
		}
		difference+=calc_backbone_bonded_difference_fourneighbors_changer(twoleftneighbor, leftneighbor, rightneighbor, tworightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_hard_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_charged_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_hard_energy_cellstruct(neighbor, -1, -1, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_charged_energy_cellstruct(neighbor, -1, -1, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changer_charged)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changen_hard(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, leftneighbor=-1, rightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0) leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
		if(movermonomer<chainlength[moverchain]-1) rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
		difference+=calc_backbone_bonded_difference_twoneighbors_changen(leftneighbor, rightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changen_hard)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changen_phenyl(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, leftneighbor=-1, rightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
		difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
		difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0) leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
		if(movermonomer<chainlength[moverchain]-1) rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
		difference+=calc_backbone_bonded_difference_twoneighbors_changen(leftneighbor, rightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_phenyl_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
		difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
		difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
		difference+=calc_phenyl_energy_cellstruct(neighbor, -1, -1, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changen_phenyl)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changen_charged(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, leftneighbor=-1, rightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0) leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
		if(movermonomer<chainlength[moverchain]-1) rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
		difference+=calc_backbone_bonded_difference_twoneighbors_changen(leftneighbor, rightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_charged_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_charged_energy_cellstruct(neighbor, -1, -1, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changen_charged)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_nonbonded_energy_pair(coord a, coord b, double_triple box_dimension, nonbonded_params params){
	if(a.nodetype==0){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			return hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			return hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else if(a.nodetype==1){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			
			return phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
	
		}
		else if(b.nodetype==2){
			return hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			return hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else if(a.nodetype==2){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			return electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
		}
		else if(b.nodetype==3){
			return electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else if(a.nodetype==3){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			return electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
		}
		else if(b.nodetype==3){
			return electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else my_exit("don't have energy functions for these node types!");
	return 0;
}

double hard_energy(double_triple ar, double_triple br, double arhard, double brhard, double_triple box_dimension){
	double_triple sep=subtract_double_triple(br, ar);
	recenter_double_triple(&sep, box_dimension);
	double rmag=norm(sep);
	if(rmag>(arhard+brhard)){
		return 0;
	}
	else{
		return 1.0/0.0;
	}
}

double hard_energy_sumradii(double_triple ar, double_triple br, double sum, double_triple box_dimension){
	double_triple sep=subtract_double_triple(br, ar);
	recenter_double_triple(&sep, box_dimension);
	double rmag=norm(sep);
	if(rmag>(sum)){
		return 0;
	}
	else{
		return 1.0/0.0;
	}
}

double calc_hard_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
    int i, index;
    double result=0;
    for(i=0;i<numberneighbors;i++){
        index=neighborlist[i];
        if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
            if(coord_new.nodetype==0){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==1){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
                }
                else if(coordarray[index].nodetype==2){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
                }
                else if(coordarray[index].nodetype==3){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
                }
            }
            else if(coord_new.nodetype==1){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==2){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
                }
                else if(coordarray[index].nodetype==3){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
                }
            }
            else if(coord_new.nodetype==2){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==1){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
                }
            }
            else if(coord_new.nodetype==3){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==1){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
                }
            }
        }
    }
    return result;
}

double calc_hard_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid!=coord_new.chainid){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double calc_hard_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid==polymer){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double calc_hard_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index, sheet=coord_new.leafid/2;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2!=sheet){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double calc_hard_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2==othersheet){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double electrostatic_energy(coord a, coord b, solvation_parameters solvationparams, electrostaticparam pa, electrostaticparam pb, double_triple box_dimension, double cutoff2){

	// assuming a is amino
	
    double_triple r=subtract_double_triple(b.r, a.r);
	recenter_double_triple(&r, box_dimension);
	double rsq=dot_product(r, r);
	if(rsq>cutoff2) return 0;
	double rmag=sqrt(rsq);
    double rminusrhard=rmag-(pa.rhard+pb.rhard);
	if(rminusrhard<0) return 1.0/0.0;
	double eps;
	r=subtract_double_triple(add_double_triple(b.r, scalar_multiply_double_triple(b.n, pb.shift)), add_double_triple(a.r, scalar_multiply_double_triple(a.n, pa.shift)));

	recenter_double_triple(&r, box_dimension);
	
    rsq=dot_product(r, r);
	rmag=sqrt(rsq);
    double rmagminushard=rmag-(pa.rhard+pb.rhard-pa.shift-pb.shift);
	if(rmagminushard>solvationparams.saturationrange) eps=waterpermittivity;
	else eps=solvationparams.shortrangepermittivity+(waterpermittivity-solvationparams.shortrangepermittivity)*rmagminushard/solvationparams.saturationrange;
	double energy=pa.charge*pb.charge/rmag*(coulombconstant/eps);
    if(rminusrhard<solvationparams.waterradius){
        energy+=(pa.solvationenergy+pb.solvationenergy)*solvationparams.firstshellfraction*(1-4*pow((rminusrhard-0.5*solvationparams.waterradius)/solvationparams.waterradius, 2));
    }
	return pa.factor*pb.factor*energy;
}

int count_cg_atoms(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray){
	int i, j, total=0, sidechainindex;
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength[i];j++){
			total+=1;		//	backbone
			sidechainindex=monomerid[i][j].sidechain;
			if(coordarray[sidechainindex].nodetype>=0) total++;
		}
	}
	return total;
}

double calc_onebody_difference(coord coord_old, coord coord_new, onebodyparam p, solvation_parameters solvp){
	if(coord_old.r.z==coord_new.r.z) return 0;
	double affinityheight, interfaceheight;
	if(coord_new.r.z>0.5*(solvp.interfaceheight1+solvp.interfaceheight2)){
		affinityheight=coord_new.r.z-(solvp.interfaceheight2+p.z0);
		interfaceheight=coord_new.r.z-(solvp.interfaceheight2+p.zinterface);
	}
	else{
		affinityheight=solvp.interfaceheight1-p.z0-coord_new.r.z;
		interfaceheight=solvp.interfaceheight1-p.zinterface-coord_new.r.z;
	}
	double newaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double newinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	if(coord_old.r.z>0.5*(solvp.interfaceheight1+solvp.interfaceheight2)){
		affinityheight=coord_old.r.z-(solvp.interfaceheight2+p.z0);
		interfaceheight=coord_old.r.z-(solvp.interfaceheight2+p.zinterface);
	}
	else{
		affinityheight=solvp.interfaceheight1-p.z0-coord_old.r.z;
		interfaceheight=solvp.interfaceheight1-p.zinterface-coord_old.r.z;
	}
	double oldaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double oldinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	return (newaffinity+newinterface-oldaffinity-oldinterface);
}

void calc_onebody_energy_components(coord mycoord, onebodyparam p, solvation_parameters solvp, energycomponents *pmy_energycomponents){
	double affinityheight, interfaceheight;
	if(mycoord.r.z>0.5*(solvp.interfaceheight1+solvp.interfaceheight2)){
		affinityheight=mycoord.r.z-(solvp.interfaceheight2+p.z0);
		interfaceheight=mycoord.r.z-(solvp.interfaceheight2+p.zinterface);
	}
	else{
		affinityheight=solvp.interfaceheight1-p.z0-mycoord.r.z;
		interfaceheight=solvp.interfaceheight1-p.zinterface-mycoord.r.z;
	}
	double affinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double interface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	(*pmy_energycomponents).solvation[mycoord.nodetype]+=affinity;
	(*pmy_energycomponents).interface[mycoord.nodetype]+=interface;
}

onebodyparam chooseonebody(int type, nonbonded_params p){
	if(type==0) return p.one0;
	if(type==1) return p.one1;
	if(type==2) return p.one2;
	if(type==3) return p.one3;
	printf("don't have onebody param for type %i!", type);
    exit(1);
    return p.one0;
}

double onebody_energy(coord mycoord, onebodyparam p, solvation_parameters solvp){
    
	double affinityheight, interfaceheight;
	if(mycoord.r.z>0.5*(solvp.interfaceheight1+solvp.interfaceheight2)){
		affinityheight=mycoord.r.z-(solvp.interfaceheight2+p.z0);
		interfaceheight=mycoord.r.z-(solvp.interfaceheight2+p.zinterface);
	}
	else{
		affinityheight=solvp.interfaceheight1-p.z0-mycoord.r.z;
		interfaceheight=solvp.interfaceheight1-p.zinterface-mycoord.r.z;
	}
	double affinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double interface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	return affinity+interface;
}

void outputcoords_centered(FILE *outp, char *name, double_triple pos, double_triple box_dimension){
	fprintf(outp, "%s\t%.6f\t%.6f\t%.6f\n", name, pos.x-0.5*box_dimension.x, pos.y-0.5*box_dimension.y, pos.z-0.5*box_dimension.z);
}

double_triple nearest_image(double_triple a, double_triple b, double_triple box){
	double_triple result=subtract_double_triple(a, b);
	recenter_double_triple(&result, box);
	return add_double_triple(b, result);
}

void updatecell(double_triple pos, linkedlistfull *plink, int n){
	int_triple new;
	new.x=(int) floor(pos.x/(*plink).cellwidth.x);
	new.y=(int) floor(pos.y/(*plink).cellwidth.y);
	new.z=(int) floor(pos.z/(*plink).cellwidth.z);
	if((new.x!=(*plink).core.cell[n].x)||(new.y!=(*plink).core.cell[n].y)||(new.z!=(*plink).core.cell[n].z)){
		if((*plink).reverselist[n]==-1){					//	n is head
			
			(*plink).core.head[(*plink).core.cell[n].x][(*plink).core.cell[n].y][(*plink).core.cell[n].z]=(*plink).core.list[n];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=-1;
		}
		else{
			(*plink).core.list[(*plink).reverselist[n]]=(*plink).core.list[n];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=(*plink).reverselist[n];
		}
		(*plink).core.list[n]=(*plink).core.head[new.x][new.y][new.z];
		if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
		
		(*plink).core.head[new.x][new.y][new.z]=n;
		(*plink).reverselist[n]=-1;
		(*plink).core.cell[n]=new;
	}
}

double total_energy_cellstruct(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, linkedlistset linkset){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	double energy=0;
	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				if(j>0){
                    
					//	backbone only has short-ranged interactions
					
                    energy+=0.5*calc_hard_energy_cellstruct(sidechainid, twoleft, backboneid, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
                }
				sidechainid=monomerid[i][j].sidechain;
				if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    energy+=0.5*calc_hard_energy_cellstruct(backboneid, -1, -1, linkset.shortrange.core.number_neighbors[sidechainid], linkset.shortrange.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype==1) energy+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) energy+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype>=0) energy+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2) energy+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
				else if(j>=3) energy+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions
            
            energy+=calc_hard_energy_cellstruct(sidechainid, twoleft, -1, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
            
            //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
		}
		else if(chainlength[i]==2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				if(j==1){
					//	backbone only has short-ranged interactions

                    energy+=calc_hard_energy_cellstruct(sidechainid, twoleft, backboneid, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
				}
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    
                    energy+=0.5*calc_hard_energy_cellstruct(backboneid, -1, -1, linkset.shortrange.core.number_neighbors[sidechainid], linkset.shortrange.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype==1) energy+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) energy+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    
                    if(coordarray[sidechainid].nodetype>=0) energy+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==1){
                    
                    //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
                    
					energy+=calc_backbone_bonded_energy_onlyleft(my_bonded_params, coordarray[backboneid], coordarray[left], box_dimension);
				}
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions
            
            energy+=calc_hard_energy_cellstruct(sidechainid, twoleft, -1, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
            
            //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
		}
        else{
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    
                    energy+=0.5*calc_hard_energy_cellstruct(backboneid, -1, -1, linkset.shortrange.core.number_neighbors[sidechainid], linkset.shortrange.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype==1) energy+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) energy+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                }
				if(coordarray[sidechainid].nodetype>=0) energy+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                energy+=calc_hard_energy_cellstruct(sidechainid, -1, -1, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
            }
        }
	}
 	return energy;
}

double total_energy_cellstruct_debug(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, linkedlistset linkset){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	energycomponents my_energycomponents;
    
    my_energycomponents.backbone=0;
	my_energycomponents.sidechain=0;
    
	my_energycomponents.nnsamepoly=0;
	my_energycomponents.ccunlikesamepoly=0;
	my_energycomponents.cclikesamepoly=0;
	my_energycomponents.cnsamepoly=0;
    
	my_energycomponents.nn=0;
	my_energycomponents.ccunlike=0;
	my_energycomponents.cclike=0;
	my_energycomponents.cn=0;
	
	my_energycomponents.nncross=0;
	my_energycomponents.ccunlikecross=0;
	my_energycomponents.cclikecross=0;
	my_energycomponents.cncross=0;
    
	my_energycomponents.nndifferent=0;
	my_energycomponents.ccunlikedifferent=0;
	my_energycomponents.cclikedifferent=0;
	my_energycomponents.cndifferent=0;
    
    my_energycomponents.solvation=xmalloc(sizeof(double)*4+1);
    my_energycomponents.interface=xmalloc(sizeof(double)*4+1);
    for(i=0;i<4;i++){
        my_energycomponents.solvation[i]=0;
        my_energycomponents.interface[i]=0;
    }

	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) my_energycomponents.solvation[coordarray[backboneid].nodetype]+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				if(j>0){
                    
					//	backbone only has short-ranged interactions
					
                }
				sidechainid=monomerid[i][j].sidechain;
				if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) my_energycomponents.solvation[coordarray[sidechainid].nodetype]+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    //printf("%i %i after sidechain onebody %f\n", i, j, energy);
                    if(coordarray[sidechainid].nodetype==1) my_energycomponents.nn+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) my_energycomponents.ccunlike+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    
                    //  factor of 0.5 to account for double counting
                    
                    if(coordarray[sidechainid].nodetype>=0) my_energycomponents.sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2) my_energycomponents.backbone+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
				else if(j>=3) my_energycomponents.backbone+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions
                        
            //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
		}
	}
    
    double totalenergy=0;
    printf("total_energy_debug: \n");
    printf("\t%f\n", my_energycomponents.backbone);
	printf("\t%f\n", my_energycomponents.sidechain);
	printf("\t%f\n", my_energycomponents.nn);
	printf("\t%f\n", my_energycomponents.ccunlike);
    for(i=0;i<4;i++){
        printf("\t%f\n", my_energycomponents.solvation[i]+my_energycomponents.interface[i]);
    }
    
    totalenergy+=my_energycomponents.backbone;
	totalenergy+=my_energycomponents.sidechain;
    
	totalenergy+=my_energycomponents.nnsamepoly;
	totalenergy+=my_energycomponents.ccunlikesamepoly;
	totalenergy+=my_energycomponents.cclikesamepoly;
	totalenergy+=my_energycomponents.cnsamepoly;
    
	totalenergy+=my_energycomponents.nn;
	totalenergy+=my_energycomponents.ccunlike;
	totalenergy+=my_energycomponents.cclike;
	totalenergy+=my_energycomponents.cn;
	
	totalenergy+=my_energycomponents.nncross;
	totalenergy+=my_energycomponents.ccunlikecross;
	totalenergy+=my_energycomponents.cclikecross;
	totalenergy+=my_energycomponents.cncross;
    
	totalenergy+=my_energycomponents.nndifferent;
	totalenergy+=my_energycomponents.ccunlikedifferent;
	totalenergy+=my_energycomponents.cclikedifferent;
	totalenergy+=my_energycomponents.cndifferent;
    for(i=0;i<4;i++){
        totalenergy+=my_energycomponents.solvation[i];
        totalenergy+=my_energycomponents.interface[i];
    }
    
    printf("\ttotal %f\n", totalenergy);

	return totalenergy;
}

double phenylphenyl_energy(coord a, coord b, GBQvac pvac, double_triple box_dimension, double cutoff2, solvation_parameters solvp){
	double epsprime, strength, GBeps, sigmaGB, LJGBarg, sixthpower, twelfthpower, energy=0, normr5, rhatd1sq, rhatd2sq;
	double_triple r=subtract_double_triple(b.r, a.r);
	recenter_double_triple(&r, box_dimension);
	double normr2=dot_product(r, r);
	if(normr2>cutoff2) return 0;
	double normr=sqrt(normr2);
	double_triple rhat=scalar_multiply_double_triple(r, 1./normr);
	double rhatd1=dot_product(rhat, a.n);
	double rhatd2=dot_product(rhat, b.n);
	double d1d2=dot_product(a.n, b.n);
	double rhatd1plusrhatd2sq=(rhatd1+rhatd2)*(rhatd1+rhatd2);
	double rhatd1minusrhatd2sq=(rhatd1-rhatd2)*(rhatd1-rhatd2);
    double solvated, heighta, heightb, sigma, epsprimesolvent;
	if((solvp.interface==0)||(pvac.code==1)) solvated=1;
	else{
        if(pvac.code==0) solvated=0;
        else{
			fmod((a.r.z), (box_dimension.z));
			fmod((b.r.z), (box_dimension.z));
            if(a.r.z>0.5*(solvp.interfaceheight1+solvp.interfaceheight2)) heighta=a.r.z-solvp.interfaceheight2;
            else heighta=solvp.interfaceheight1-a.r.z;
            if(b.r.z>0.5*(solvp.interfaceheight1+solvp.interfaceheight2)) heightb=b.r.z-solvp.interfaceheight2;
            else heightb=solvp.interfaceheight1-b.r.z;
            heighta-=pvac.p11solv.interpolatemiddle;
            heightb-=pvac.p11solv.interpolatemiddle;
            solvated=0.5*(1./(1+exp(4*heighta/pvac.p11solv.interpolatewidth))+1./(1+exp(4*heightb/pvac.p11solv.interpolatewidth)));
        }
	}
	if(solvated<1){
		epsprime=1-0.5*pvac.chiprime*(rhatd1plusrhatd2sq/(1+pvac.chiprime*d1d2)+rhatd1minusrhatd2sq/(1-pvac.chiprime*d1d2));
		strength=1./sqrt(1-pvac.chi*pvac.chi*d1d2*d1d2);
		GBeps=epsprime/(strength*strength);			//	strength^nu*epsprime^mu, fixed mu=1, nu=-2
		sigmaGB=1./sqrt(1.-0.5*pvac.chi*(rhatd1plusrhatd2sq/(1+pvac.chi*d1d2)+rhatd1minusrhatd2sq/(1-pvac.chi*d1d2)));
		LJGBarg=pvac.sigma0*pvac.xi/(normr+pvac.sigma0*(pvac.xi-sigmaGB));
		sixthpower=LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg;
		twelfthpower=sixthpower*sixthpower;
		energy+=(1.-solvated)*4*pvac.eps0*GBeps*(twelfthpower-sixthpower);
	}
	if(solvated>0){
		epsprime=1-0.5*pvac.p11solv.chiprime*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chiprime*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chiprime*d1d2));
		strength=1./sqrt(1-pvac.p11solv.chi*pvac.p11solv.chi*d1d2*d1d2);
		GBeps=epsprime/(strength*strength);			//	strength^nu*epsprime^mu, fixed mu=1, nu=-2
		sigmaGB=1./sqrt(1.-0.5*pvac.p11solv.chi*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chi*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chi*d1d2)));
		LJGBarg=pvac.p11solv.sigma0*pvac.p11solv.xi/(normr+pvac.p11solv.sigma0*(pvac.p11solv.xi-sigmaGB));
		sixthpower=LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg;
		twelfthpower=sixthpower*sixthpower;
		energy+=solvated*4*pvac.p11solv.eps0*GBeps*(twelfthpower-sixthpower);
	}
	if(solvated>0){
		sigma=pvac.p11solv.sigma0*sigmaGB;
		if(normr<sigma+2*pvac.p11solv.w){
			if(normr<0.5*sigma) energy+=1./0.;
			else if(normr<sigma+pvac.p11solv.w){
				epsprimesolvent=1-0.5*pvac.p11solv.chirep*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chirep*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chirep*d1d2));
				energy+=solvated*pvac.p11solv.amprep*epsprimesolvent*(1-4/(pvac.p11solv.w*pvac.p11solv.w)*pow(normr-(sigma+0.5*pvac.p11solv.w), 2));
			}
			else{
				epsprimesolvent=1-0.5*pvac.p11solv.chiattr*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chiattr*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chiattr*d1d2));
				energy-=solvated*pvac.p11solv.ampattr*epsprimesolvent*(1-4/(pvac.p11solv.w*pvac.p11solv.w)*pow(normr-(sigma+1.5*pvac.p11solv.w), 2));
			}
		}
	}
    
	if((solvated<1)&&(normr>0.5)){
		normr5=normr2*normr2*normr;
		rhatd1sq=rhatd1*rhatd1;
		rhatd2sq=rhatd2*rhatd2;
		energy+=((1.-solvated)*coulombconstant*0.75*pvac.Q*pvac.Q/normr5*(1+2*d1d2*d1d2-5*(rhatd1sq+rhatd2sq)-20*rhatd1*rhatd2*d1d2+35*rhatd1sq*rhatd2sq));
	}

    return pvac.factor*energy;
}

void updateneighborlist(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension){
	int i, j, k, i_image, j_image, k_image, neighbor;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_pairenergies_phenyl(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, i_image, j_image, k_image, neighbor;
	double energy;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			(*plink).core.pair_energy[neighbor][j-1]=(*plink).core.pair_energy[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
							
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_pairenergies_charged(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, i_image, j_image, k_image, neighbor;
	double energy;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			(*plink).core.pair_energy[neighbor][j-1]=(*plink).core.pair_energy[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}		
		
	}
}

void updateneighborlist_nowipeout(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension){
	int i, j, k, i_image, j_image, k_image, neighbor;
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_pairenergies_phenyl_nowipeout(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, i_image, j_image, k_image, neighbor;
	double energy;
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
							
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_pairenergies_charged_nowipeout(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, i_image, j_image, k_image, neighbor;
	double energy;
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

double calc_phenyl_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
		if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_charged_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
		if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid!=coord_new.chainid){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid==polymer){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid!=coord_new.chainid){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid==polymer){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index, sheet=coord_new.leafid/2;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2!=sheet){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2==othersheet){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index, sheet=coord_new.leafid/2;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2!=sheet){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2==othersheet){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

int_triple cell(double_triple pos, double_triple cellwidth){
    int_triple result;
    result.x=(int) floor(pos.x/cellwidth.x);
    result.y=(int) floor(pos.y/cellwidth.y);
    result.z=(int) floor(pos.z/cellwidth.z);
    return result;
}

void create_temporary_neighborlist(linkedlistfull link, coord *coordarray, double_triple box_dimension, int *pnumber_neighbors, int *neighborlist, int_triple selfcell, coord selfcoord){
	int i, j, k, i_image, j_image, k_image, neighbor;
	(*pnumber_neighbors)=0;
	double_triple selfpos=selfcoord.r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=link.core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==link.core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=link.core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==link.core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=link.core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==link.core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=link.core.head[i_image][j_image][k_image];
				while(neighbor>=0){
                    if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<link.mincellwidth){
                        neighborlist[*pnumber_neighbors]=neighbor;
                        (*pnumber_neighbors)++;
                    }
					neighbor=link.core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==link.core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==link.core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==link.core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

int in_neighborhood(coord second, coord last, coord subject, reptationparams reptation, double_triple box_dimension){
	double_triple sidechaindirector=last.n;
	double_triple backbonedirector=subtract_double_triple(last.r, second.r);
	recenter_double_triple(&backbonedirector, box_dimension);
	backbonedirector=subtract_double_triple(backbonedirector, scalar_multiply_double_triple(sidechaindirector, dot_product(sidechaindirector, backbonedirector)));
	normalize(&backbonedirector);
	double_triple center=add_double_triple(last.r, add_double_triple(scalar_multiply_double_triple(backbonedirector, reptation.rside0), scalar_multiply_double_triple(sidechaindirector, reptation.rforward0)));
	double_triple sep=subtract_double_triple(subject.r, center);
	recenter_double_triple(&sep, box_dimension);
	if(sep.x*sep.x+sep.y*sep.y+sep.z*sep.z<reptation.range2){
		double dotproduct=dot_product(last.n, subject.n);
		if(dotproduct<reptation.maxconsecutivedirectordotproduct) return 1;
        else return 0;
    }
	else{
        return 0;
    }
}

int in_side_neighborhood(coord backbone, coord subject, sidereptationparams reptation, double_triple box_dimension){
	double_triple center=add_double_triple(backbone.r, scalar_multiply_double_triple(backbone.n, reptation.r0));
	double_triple sep=subtract_double_triple(subject.r, center);
	recenter_double_triple(&sep, box_dimension);
 	if(sep.x*sep.x+sep.y*sep.y+sep.z*sep.z<reptation.range2){
        double dotproduct=dot_product(backbone.n, subject.n);
		if(dotproduct>reptation.mindotproduct) return 1;
        else return 0;
    }
	else return 0;
}

coord reptate(coord second, coord last, reptationparams reptation, double_triple box_dimension){
	coord result=second;
	double_triple sidechaindirector=last.n;
	double_triple backbonedirector=subtract_double_triple(last.r, second.r);
	recenter_double_triple(&backbonedirector, box_dimension);
	backbonedirector=subtract_double_triple(backbonedirector, scalar_multiply_double_triple(sidechaindirector, dot_product(sidechaindirector, backbonedirector)));
	normalize(&backbonedirector);
	double_triple center=add_double_triple(scalar_multiply_double_triple(backbonedirector, reptation.rside0), scalar_multiply_double_triple(sidechaindirector, reptation.rforward0));
	result.r=add_double_triple(center, scalar_multiply_double_triple(rand_unit_ball(), reptation.range));
	result.n=rand_unit_sphere();
	if(reptation.maxconsecutivedirectordotproduct<1){
		double dotproduct=dot_product(result.n, last.n);
		while(dotproduct>reptation.maxconsecutivedirectordotproduct){
			result.n=rand_unit_sphere();
			dotproduct=dot_product(result.n, last.n);
		}
	}
    fmod_double_triple(&(result.r), box_dimension);
	return result;
}

coord side_reptate(coord backbone, coord image, sidereptationparams reptation, double_triple box_dimension){
	coord result=image;
	double_triple center=add_double_triple(backbone.r, scalar_multiply_double_triple(backbone.n, reptation.r0));
	result.r=add_double_triple(center, scalar_multiply_double_triple(rand_unit_ball(), reptation.range));
	result.n=rand_unit_sphere();
	if(reptation.mindotproduct>-1){
		double dotproduct=dot_product(result.n, backbone.n);
		while(dotproduct<reptation.mindotproduct){
			result.n=rand_unit_sphere();
			dotproduct=dot_product(result.n, backbone.n);
		}
	}
    fmod_double_triple(&(result.r), box_dimension);
	return result;
}

double calc_charged_energy_difference_swapidentity_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int newtype, int oldtype){
    int i, index;
    double result=0;
    for(i=0;i<numberneighbors;i++){
        index=neighborlist[i];
        if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
            if((newtype==2)&&(oldtype==3)){
                if(coordarray[index].nodetype==2){
                    result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
                    result-=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
                }
                else if(coordarray[index].nodetype==3){
                    result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
                    result-=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
                }
            }
            else if((newtype==3)&&(oldtype==2)){
                if(coordarray[index].nodetype==2){
                    result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
                    result-=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
                }
                else if(coordarray[index].nodetype==3){
                    result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
                    result-=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
                }
            }
            else{
                printf("swapping unknown types %i and %i", oldtype, newtype);
                exit(1);
            }
        }
    }
    return result;
}

void deletecell(linkedlistfull *plink, int n){
    if((*plink).reverselist[n]==-1){					//	n is head
        (*plink).core.head[(*plink).core.cell[n].x][(*plink).core.cell[n].y][(*plink).core.cell[n].z]=(*plink).core.list[n];
        if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=-1;
    }
    else{
        (*plink).core.list[(*plink).reverselist[n]]=(*plink).core.list[n];
        if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=(*plink).reverselist[n];
    }
    int i, j, neighbor;
	for(i=0;i<(*plink).core.number_neighbors[n];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[n][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=n) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			(*plink).core.pair_energy[neighbor][j-1]=(*plink).core.pair_energy[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[n]=0;
}

void shiftcell(linkedlistfull *plink, int newindex, int oldindex){
    if((*plink).reverselist[oldindex]==-1){					//	n is head
        (*plink).core.head[(*plink).core.cell[oldindex].x][(*plink).core.cell[oldindex].y][(*plink).core.cell[oldindex].z]=newindex;
        (*plink).reverselist[newindex]=-1;
    }
    else{
        (*plink).core.list[(*plink).reverselist[oldindex]]=newindex;
        (*plink).reverselist[newindex]=(*plink).reverselist[oldindex];
    }
    (*plink).core.list[newindex]=(*plink).core.list[oldindex];
    if((*plink).core.list[newindex]!=-1) (*plink).reverselist[(*plink).core.list[newindex]]=newindex;
    (*plink).core.cell[newindex]=(*plink).core.cell[oldindex];
    
    (*plink).core.number_neighbors[newindex]=(*plink).core.number_neighbors[oldindex];
    int i;
    for(i=0;i<(*plink).core.number_neighbors[newindex];i++){
        (*plink).core.neighbor[newindex][i]=(*plink).core.neighbor[oldindex][i];
        (*plink).core.pair_energy[newindex][i]=(*plink).core.pair_energy[oldindex][i];
    }
}

void updatecellnocheck(double_triple pos, linkedlistfull *plink, int n){
	int_triple new;
	new.x=(int) floor(pos.x/(*plink).cellwidth.x);
	new.y=(int) floor(pos.y/(*plink).cellwidth.y);
	new.z=(int) floor(pos.z/(*plink).cellwidth.z);
    if((*plink).reverselist[n]==-1){					//	n is head
        (*plink).core.head[(*plink).core.cell[n].x][(*plink).core.cell[n].y][(*plink).core.cell[n].z]=(*plink).core.list[n];
        if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=-1;
    }
    else{
        (*plink).core.list[(*plink).reverselist[n]]=(*plink).core.list[n];
        if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=(*plink).reverselist[n];
    }
    (*plink).core.list[n]=(*plink).core.head[new.x][new.y][new.z];
    if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
    (*plink).core.head[new.x][new.y][new.z]=n;
    (*plink).reverselist[n]=-1;
    (*plink).core.cell[n]=new;
}

double calc_energy_pair(coord a, coord b, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	double result=0;
	if(a.nodetype==0){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, 2.*my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
		}
	}
	else if(a.nodetype==1){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			
			result+=phenylphenyl_energy(a, b, my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
		else if(b.nodetype==2){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
		}
	}
	else if(a.nodetype==2){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
		}
		else if(b.nodetype==3){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
		}
	}
	else if(a.nodetype==3){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
		}
		else if(b.nodetype==3){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
		}
	}
	return result;
}

double calc_backbone_bonded_difference_fourneighbors_changer(int twoleft, int left, int right, int tworight, bonded_params my_bonded_params, coord *coordarray, coord oldcoord, coord newcoord, double_triple box_dimension){
	double difference=0;
	double riip1normnew, riip1normold, riip1dotninew, riip1dotniold, riip1dotnip1new, riip1dotnip1old, rim1inormnew, rim1inormold, rim1idotninew, rim1idotniold, rim1idotnim1new, rim1idotnim1old, dotcrossproduct;
	double_triple riip1new, riip1old, rhatiip1new, rhatiip1old, rim1inew, rim1iold, rhatim1inew, rhatim1iold;
	double rip1ip2dotnip1, rip1ip2dotnip2, rim2im1dotnim1, rim2im1dotnim2;
	double_triple rip1ip2, rhatip1ip2, rim2im1, rhatim2im1;
	if(right>=0){
		riip1new=subtract_double_triple(coordarray[right].r, newcoord.r);
		riip1old=subtract_double_triple(coordarray[right].r, oldcoord.r);
		recenter_double_triple(&riip1new, box_dimension);
		recenter_double_triple(&riip1old, box_dimension);
		riip1normnew=norm(riip1new);
		riip1normold=norm(riip1old);
		rhatiip1new=scalar_multiply_double_triple(riip1new, 1./riip1normnew);
		rhatiip1old=scalar_multiply_double_triple(riip1old, 1./riip1normold);
		
        difference+=my_bonded_params.factor*(my_bonded_params.K10+riip1normnew*(my_bonded_params.K11+riip1normnew*(my_bonded_params.K12+riip1normnew*(my_bonded_params.K13+riip1normnew*my_bonded_params.K14))));	//	onlyright
        difference-=my_bonded_params.factor*(my_bonded_params.K10+riip1normold*(my_bonded_params.K11+riip1normold*(my_bonded_params.K12+riip1normold*(my_bonded_params.K13+riip1normold*my_bonded_params.K14))));	//	onlyright
		
		riip1dotninew=dot_product(rhatiip1new, newcoord.n);
		riip1dotniold=dot_product(rhatiip1old, newcoord.n);
		riip1dotnip1new=dot_product(rhatiip1new, coordarray[right].n);
		riip1dotnip1old=dot_product(rhatiip1old, coordarray[right].n);

        difference+=(my_bonded_params.kr*pow(riip1dotninew-(riip1normnew-my_bonded_params.r0r)/my_bonded_params.sr, 2));			
        difference-=(my_bonded_params.kr*pow(riip1dotniold-(riip1normold-my_bonded_params.r0r)/my_bonded_params.sr, 2));			
        difference+=(my_bonded_params.kl*pow(riip1dotnip1new-(riip1normnew-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference-=(my_bonded_params.kl*pow(riip1dotnip1old-(riip1normold-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        
		if(tworight>=0){
			rip1ip2=subtract_double_triple(coordarray[tworight].r, coordarray[right].r);
			recenter_double_triple(&rip1ip2, box_dimension);
			rhatip1ip2=scalar_multiply_double_triple(rip1ip2, 1./norm(rip1ip2));            
			rip1ip2dotnip1=dot_product(rhatip1ip2, coordarray[right].n);
			rip1ip2dotnip2=dot_product(rhatip1ip2, coordarray[tworight].n);
            
			dotcrossproduct=dot_product(rhatiip1new, rhatip1ip2)-riip1dotnip1new*rip1ip2dotnip1;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));

			dotcrossproduct=dot_product(rhatiip1old, rhatip1ip2)-riip1dotnip1old*rip1ip2dotnip1;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
		}
	}
	if(left>=0){
		rim1inew=subtract_double_triple(newcoord.r, coordarray[left].r);
		rim1iold=subtract_double_triple(oldcoord.r, coordarray[left].r);
		recenter_double_triple(&rim1inew, box_dimension);
		recenter_double_triple(&rim1iold, box_dimension);
		rim1inormnew=norm(rim1inew);
		rim1inormold=norm(rim1iold);
		rhatim1inew=scalar_multiply_double_triple(rim1inew, 1./rim1inormnew);
		rhatim1iold=scalar_multiply_double_triple(rim1iold, 1./rim1inormold);
		
        difference+=my_bonded_params.factor*(my_bonded_params.K10+rim1inormnew*(my_bonded_params.K11+rim1inormnew*(my_bonded_params.K12+rim1inormnew*(my_bonded_params.K13+rim1inormnew*my_bonded_params.K14))));
        difference-=my_bonded_params.factor*(my_bonded_params.K10+rim1inormold*(my_bonded_params.K11+rim1inormold*(my_bonded_params.K12+rim1inormold*(my_bonded_params.K13+rim1inormold*my_bonded_params.K14))));
        
		rim1idotninew=dot_product(rhatim1inew, newcoord.n);
		rim1idotniold=dot_product(rhatim1iold, oldcoord.n);
		rim1idotnim1new=dot_product(rhatim1inew, coordarray[left].n);
		rim1idotnim1old=dot_product(rhatim1iold, coordarray[left].n);
        
        difference+=(my_bonded_params.kl*pow(rim1idotninew-(rim1inormnew-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference-=(my_bonded_params.kl*pow(rim1idotniold-(rim1inormold-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference+=(my_bonded_params.kr*pow(rim1idotnim1new-(rim1inormnew-my_bonded_params.r0r)/my_bonded_params.sr, 2));
        difference-=(my_bonded_params.kr*pow(rim1idotnim1old-(rim1inormold-my_bonded_params.r0r)/my_bonded_params.sr, 2));
		
		if(right>=0){
			dotcrossproduct=dot_product(rhatim1inew, rhatiip1new)-rim1idotninew*riip1dotninew;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));

			dotcrossproduct=dot_product(rhatim1iold, rhatiip1old)-rim1idotniold*riip1dotniold;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
		}
		if(twoleft>=0){
			rim2im1=subtract_double_triple(coordarray[left].r, coordarray[twoleft].r);
			recenter_double_triple(&rim2im1, box_dimension);
			rhatim2im1=scalar_multiply_double_triple(rim2im1, 1./norm(rim2im1));
			
			rim2im1dotnim1=dot_product(rhatim2im1, coordarray[left].n);
			rim2im1dotnim2=dot_product(rhatim2im1, coordarray[twoleft].n);
            
			dotcrossproduct=dot_product(rhatim2im1, rhatim1inew)-rim2im1dotnim1*rim1idotnim1new;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));

			dotcrossproduct=dot_product(rhatim2im1, rhatim1iold)-rim2im1dotnim1*rim1idotnim1old;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
		}
	}
	return difference;
}

double calc_backbone_bonded_difference_twoneighbors_changen(int left, int right, bonded_params my_bonded_params, coord *coordarray, coord oldcoord, coord newcoord, double_triple box_dimension){
	double difference=0;
	double riip1dotninew, riip1dotniold, rim1idotninew, rim1idotniold, dotcrossproduct, mydot, riip1norm, rim1inorm;
	double_triple riip1, rhatiip1, rim1i, rhatim1i;
	if(right>=0){
		riip1=subtract_double_triple(coordarray[right].r, newcoord.r);
		recenter_double_triple(&riip1, box_dimension);
		riip1norm=norm(riip1);
		rhatiip1=scalar_multiply_double_triple(riip1, 1./norm(riip1));
        
		riip1dotninew=dot_product(rhatiip1, newcoord.n);
		riip1dotniold=dot_product(rhatiip1, oldcoord.n);

        difference+=(my_bonded_params.kr*pow(riip1dotninew-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));
        difference-=(my_bonded_params.kr*pow(riip1dotniold-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));
	}
	if(left>=0){
		rim1i=subtract_double_triple(newcoord.r, coordarray[left].r);
		recenter_double_triple(&rim1i, box_dimension);
		rim1inorm=norm(rim1i);
		rhatim1i=scalar_multiply_double_triple(rim1i, 1./norm(rim1i));
        
		rim1idotninew=dot_product(rhatim1i, newcoord.n);
		rim1idotniold=dot_product(rhatim1i, oldcoord.n);
        difference+=(my_bonded_params.kl*pow(rim1idotninew-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference-=(my_bonded_params.kl*pow(rim1idotniold-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));
		
		if(right>=0){
			mydot=dot_product(rhatim1i, rhatiip1);
			dotcrossproduct=mydot-rim1idotninew*riip1dotninew;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
			
            dotcrossproduct=mydot-rim1idotniold*riip1dotniold;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
        }
	}
	return difference;
}

int count_allatom_atoms(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray){
	int i, j, total=0, sidechainindex;
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength[i];j++){
            if(j==0) total+=6;
            else if(j==chainlength[i]-1){
                total+=5;
            }
			else total+=7;		//	backbone
			sidechainindex=monomerid[i][j].sidechain;
            if(coordarray[sidechainindex].nodetype>=0) total+=4;           //  ethyl H's
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				if(coordarray[sidechainindex].nodetype==1) total+=12;				//	phenyl
			}
			if(coordarray[sidechainindex].orientationtype==0){				//	parallel
				if(coordarray[sidechainindex].nodetype==2) total+=5;				//	amino
				if(coordarray[sidechainindex].nodetype==3) total+=4;				//	carboxyl
				if(coordarray[sidechainindex].nodetype==4) total+=1;				//	methyl
			}
		}
	}
	return total;
}

void output_xyz_allatom(int Nchains, int *chainlength, monomernodes **monomerid, int Natoms, coord *coordarray, double_triple box_dimension, char *xyzname, char *sourcename, allatomparams params, int *pframe){
	int i, j, k, backboneindex, sidechainindex, atomcount=0, lastbackbonecount, Cgammacount;
	double_triple sidevector, director, thirdvector, loc, backbonevector, com, backbonebonded, sidechainbonded, lastbackbonebonded, firstbackbonebonded, sep, lastbackbone, Cbetapos, Cgammapos, ethylvec, ethylplane, ethyloutofplane;
    
	FILE *outp, *sourceoutp;
    if((*pframe)==0){
		sourceoutp=fopen(sourcename, "w");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fprintf(sourceoutp, "topo clearbonds\n");
		fprintf(sourceoutp, "set sel [atomselect top \" name B\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name G\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name P\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name A\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name L\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name D\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name C\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name O\"]\n$sel set radius %.6f\n", params.Orad);
		fprintf(sourceoutp, "set sel [atomselect top \" name R\"]\n$sel set radius %.6f\n", params.Orad);
		fprintf(sourceoutp, "set sel [atomselect top \" name N\"]\n$sel set radius %.6f\n", params.Nrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name M\"]\n$sel set radius %.6f\n", params.Nrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name H\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name Q\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name I\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name J\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name K\"]\n$sel set radius %.6f\n", params.Hrad);
		fprintf(sourceoutp, "color Name N white\n");      // backbone N
		fprintf(sourceoutp, "color Name B white\n");     //  C beta
		fprintf(sourceoutp, "color Name G white\n");     //  phenyl C gamma
		fprintf(sourceoutp, "color Name H white\n");    //  backbone H
		fprintf(sourceoutp, "color Name C white\n");     //  backbone C
		fprintf(sourceoutp, "color Name O white\n");      //  backbone (carbonyl) O
		fprintf(sourceoutp, "color Name P yellow\n");     //  phenyl C
		fprintf(sourceoutp, "color Name Q yellow\n");     //  phenyl H
		fprintf(sourceoutp, "color Name A white\n");     //  amino C gamma
		fprintf(sourceoutp, "color Name M blue\n");     //  amino N
		fprintf(sourceoutp, "color Name I blue\n");     //  amino H
		fprintf(sourceoutp, "color Name L white\n");     //  carboxyl C gamma
		fprintf(sourceoutp, "color Name D red\n");     //  carboxyl C
		fprintf(sourceoutp, "color Name R red\n");     //  carboxyl O
		fprintf(sourceoutp, "color Name J white\n");     //  beta H
		fprintf(sourceoutp, "color Name K white\n");     //  gamma H
		outp=fopen(xyzname, "w");
	}
	else{
		sourceoutp=fopen(sourcename, "a");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fclose(sourceoutp);
		outp=fopen(xyzname, "a");
	}
	fprintf(outp, "%i\n%i\t%.6f\t%.6f\t%.6f\n", Natoms, Nchains, box_dimension.x, box_dimension.y, box_dimension.z);
	for(i=0;i<Nchains;i++){
		com.x=com.y=com.z=0;
		for(j=0;j<chainlength[i];j++){
			if(j==0){
				com=coordarray[monomerid[i][j].backbone].r;
				lastbackbone=coordarray[monomerid[i][j].backbone].r;
			}
			else{
				sep=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, lastbackbone);
				recenter_double_triple(&sep, box_dimension);
				lastbackbone=add_double_triple(lastbackbone, sep);
				com=add_double_triple(com, lastbackbone);
			}
		}
		com=scalar_multiply_double_triple(com, (1./chainlength[i]));
		fmod_double_triple(&com, box_dimension);
        
        //  work left from middle to keep entire polymer in same box as center of mass
        
        for(j=chainlength[i]/2;j>=0;j--){
			backboneindex=monomerid[i][j].backbone;
			if(j==chainlength[i]/2) firstbackbonebonded=backbonebonded=nearest_image(coordarray[backboneindex].r, com, box_dimension);
			else backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			lastbackbonebonded=backbonebonded;
        }
        
        //  now go forward, keeping polymer in same unit cell
        
        for(j=0;j<chainlength[i];j++){
			backboneindex=monomerid[i][j].backbone;
			sidechainindex=monomerid[i][j].sidechain;
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
            
			outputcoords_centered(outp, "N", backbonebonded, box_dimension);				//	backbone N
            if(coordarray[sidechainindex].nodetype>=0){                                     
                loc=scalar_multiply_double_triple(director, params.Cbetaspacing);
                Cbetapos=loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "B", loc, box_dimension);                       //  Cbeta
            }
            else{
                loc=scalar_multiply_double_triple(director, params.CterminuslongHdistance);
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);                       //  one of two terminal H's for short R terminus
            }
            if(j==0){
                loc=add_double_triple(scalar_multiply_double_triple(director, params.NterminusHdepth), scalar_multiply_double_triple(sidevector, params.NterminusHwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);                       //  left terminal H
            }
            else{
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cnyldepth), scalar_multiply_double_triple(sidevector, params.Cnylwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "C", loc, box_dimension);
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Onyldepth), scalar_multiply_double_triple(sidevector, params.Onylwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "O", loc, box_dimension);                       //  carbonyl
            }
            if(j==chainlength[i]-1){
                if(coordarray[sidechainindex].nodetype>=0){
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminusHdepth), scalar_multiply_double_triple(sidevector, params.CterminusHwidth));
                    loc=add_double_triple(backbonebonded, loc);
                    outputcoords_centered(outp, "H", loc, box_dimension);                       //  terminal H
                }
                else{
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminuslongHdepth), scalar_multiply_double_triple(sidevector, params.CterminuslongHwidth));
                    loc=add_double_triple(backbonebonded, loc);
                    outputcoords_centered(outp, "H", loc, box_dimension);                       //  terminal H
                }
            }
            else{                                                                               //  C
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cdepth), scalar_multiply_double_triple(sidevector, params.Cwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "C", loc, box_dimension);
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, -params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);
            }
            
            if((*pframe)==0){
                if(j==0){
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - terminal H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+3);                //  N - C
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+3, atomcount+4);              //  C - H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+3, atomcount+5);              //  C - H
                }
                else if(j==1){
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+3, atomcount+2);  //  previous C - current Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+2, atomcount+3);              //  Cnyl - O
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);                //  N - C
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+5);              //  C - H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+6);              //  C - H
                }
                else if(j==chainlength[i]-1){
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+4, atomcount+2);  //  previous C - current Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+2, atomcount+3);              //  Cnyl - O
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);                //  N - terminal H
                }
                else{
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+4, atomcount+2);  //  previous C - current Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+2, atomcount+3);              //  Cnyl - O
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);                //  N - C
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+5);              //  C - H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+6);              //  C - H
                }
            }
            
            lastbackbonecount=atomcount;
            if(j==0) atomcount+=6;
            else if(j==chainlength[i]-1) atomcount+=5;
			else atomcount+=7;		//	backbone
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
                    Cgammacount=atomcount+6;
                    if((*pframe)==0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount+6);          //  beta C to gamma C
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+6, atomcount);            //  gamma C to aromatic ring
                        for(k=0;k<5;k++){
                            fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+k, atomcount+k+1);    //  ring
                        }
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+5, atomcount);            //  complete ring
                        for(k=1;k<6;k++){
                            fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+k, atomcount+k+6);    //  hydrogens
                        }
                    }
                    atomcount+=12;
					loc=scalar_multiply_double_triple(sidevector, -params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=scalar_multiply_double_triple(sidevector, params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
                    
                    loc=scalar_multiply_double_triple(sidevector, -params.phenylCgammaspacing);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "G", loc, box_dimension);
                    
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=scalar_multiply_double_triple(sidevector, params.phenylHspacing);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
                    
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
                    Cgammacount=atomcount+1;
                    if((*pframe)==0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount+1);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+1, atomcount);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+3);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);
                    }
                    atomcount+=5;
					loc=scalar_multiply_double_triple(director, params.aminoNdepth);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "M", loc, box_dimension);
					loc=scalar_multiply_double_triple(director, params.aminoCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "A", loc, box_dimension);
                    
					loc=add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "I", loc, box_dimension);
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "I", loc, box_dimension);
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "I", loc, box_dimension);
				}
				if(coordarray[sidechainindex].nodetype==3){					//	carboxyl
                    Cgammacount=atomcount;
                    if((*pframe)==0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+1, atomcount+2);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+1, atomcount+3);
                    }
                    atomcount+=4;
					loc=scalar_multiply_double_triple(director, params.carboxylbackCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "L", loc, box_dimension);
					loc=scalar_multiply_double_triple(director, params.carboxylforwardCdepth);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "D", loc, box_dimension);
					//loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(thirdvector, -params.carboxylOwidth));
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, -params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "R", loc, box_dimension);
					//loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(thirdvector, params.carboxylOwidth));
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "R", loc, box_dimension);
				}
			}
            if(coordarray[sidechainindex].nodetype>=0){
                
                //  C beta ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cbetapos, backbonebonded)), normed(subtract_double_triple(Cbetapos, Cgammapos))));
                ethylplane=subtract_double_triple(Cgammapos, backbonebonded);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
                outputcoords_centered(outp, "J", loc, box_dimension);
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
                outputcoords_centered(outp, "J", loc, box_dimension);
                if((*pframe)==0){
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount);
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount+1);
                }
                
                //  C gamma ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cgammapos, Cbetapos)), normed(subtract_double_triple(Cgammapos, sidechainbonded))));
                ethylplane=subtract_double_triple(sidechainbonded, Cbetapos);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
                outputcoords_centered(outp, "K", loc, box_dimension);
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
                outputcoords_centered(outp, "K", loc, box_dimension);
                if((*pframe)==0){
                    fprintf(sourceoutp, "topo addbond %i %i\n", Cgammacount, atomcount+2);
                    fprintf(sourceoutp, "topo addbond %i %i\n", Cgammacount, atomcount+3);
                }
                atomcount+=4;
            }
		}
	}
    if((*pframe)==0) fclose(sourceoutp);
	fclose(outp);
    (*pframe)++;
}

void calc_com_trajectory(int Nnodes, coord *coordarray, double_triple box_dimension, double_triple *pavcom, int *pcomcount){
    int i, realnodes=0;
    double_triple com, sep, pos;
    for(i=0;i<Nnodes;i++){
        if(coordarray[i].nodetype>=0){
            realnodes++;
            if(i==0){
                pos=coordarray[i].r;
                com=pos;
            }
            else{
                sep=subtract_double_triple(coordarray[i].r, coordarray[i-1].r);
                recenter_double_triple(&sep, box_dimension);
                pos=add_double_triple(pos, sep);
                com=add_double_triple(com, pos);
            }
        }
	}
	com=scalar_multiply_double_triple(com, (1./realnodes));
    fmod_double_triple(&com, box_dimension);
    (*pavcom)=add_double_triple((*pavcom), com);
    (*pcomcount)++;
}

void output_com(char *filename, double_triple *pavcom, int *pcomcount, int cycle){
    FILE *outp;
    outp=fopen(filename, "a");
    fprintf(outp, "%i %f %f %f\n", cycle, (*pavcom).x/(*pcomcount), (*pavcom).y/(*pcomcount), (*pavcom).z/(*pcomcount));
    (*pavcom).x=(*pavcom).y=(*pavcom).z=0;
    (*pcomcount)=0;
    fclose(outp);
}

void expose_input_polymer_to_interface(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double newheight, double newwidth, double newdepth, double newvacuumthickness){
	int i, j;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &(*pinterface));
	if((*pinterface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &(*pNnodes), &(*pNchains), &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	if((*pNchains)!=1){
		printf("Nchains=%i in expose_input_polymer_to_interface!\n", *pNchains);
		exit(1);
	}
	for(i=0;i<(*pNnodetypes);i++) fscanf(inp, " %i", &((*pmonomercount).type[i]));
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<*pNchains;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	double_triple com, sep, newpos;
	int realnodes=0;
	for(i=0;i<(*pNnodes);i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
        if((*coordarray)[i].nodetype>=0){
            realnodes++;
            if(i==0){
                com=(*coordarray)[i].r;
            }
            else{
                sep=subtract_double_triple((*coordarray)[i].r, (*coordarray)[i-1].r);
                recenter_double_triple(&sep, *pbox_dimension);
                (*coordarray)[i].r=add_double_triple((*coordarray)[i-1].r, sep);
                com=add_double_triple(com, (*coordarray)[i].r);
            }
        }
	}
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
	com=scalar_multiply_double_triple(com, (1./realnodes));
	
	for(i=0;i<(*pNnodes);i++){
        if((*coordarray)[i].nodetype>=0){			
            (*coordarray)[i].r=subtract_double_triple((*coordarray)[i].r, com);
			recenter_double_triple(&((*coordarray)[i].r), (*pbox_dimension));				//	putting coordinates near 0 before changing box_dimension
        }
	}

	(*pbox_dimension).x=(*pbox_dimension).y=newwidth;
	(*pvacuumthickness)=newvacuumthickness;
	(*pbox_dimension).z=newheight+(*pvacuumthickness);
	newpos.x=0.5*(*pbox_dimension).x;
	newpos.y=0.5*(*pbox_dimension).y;
	newpos.z=(*pbox_dimension).z-0.5*(*pvacuumthickness)-newdepth;	
	for(i=0;i<(*pNnodes);i++){
        if((*coordarray)[i].nodetype>=0){
			(*coordarray)[i].r=add_double_triple((*coordarray)[i].r, newpos);
        }
	}
				
}
