#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"
#include "mymath.h"
#include "peptoid_functions.h"
#include "post_analysis_functions.h"

int main(int argc, char *argv[]){
	int charlength=400, i, j, k, l, m, n, Nchains, read, uniformchainlength, cycle, mincycle, foundfirstcycle=0, firstcycle, count_r=0, count_rleftdotn=0, count_rrightdotn=0, count_tripleproduct=0, count_r_rleftdotn=0, count_r_rrightdotn=0, count_tripleproduct_rleftdotn=0, count_tripleproduct_rrightdotn=0, dobackbone, dosidechain, dogr, dophendotproduct, dochargeddotproduct, leafflag, allatoms, allatomframe=0, doallatomxyz, *chainlength, doxrdallatom, ngrouptypes=18, symmetries, Nleaves=0, Nresidues, Nsites, terminuscode, nframes, totalframes=0, framecounter=0, xrdcode, xrdcount=0, nsitetypes=4;
    double xrdqspace, rtop, grbinwidth, chargedthreshold, chargedouterthreshold, dpbinwidth, maxq;
	declare_array(double, phenthresholds, 3);
    int_triple xrdqnarray;
	double_triple box_dimension, xrdqspacearray, avbox;
	avbox.x=avbox.y=avbox.z=0;
	backbonehistoparams histoparams;
	declare_array(double, scatterfactorallatom, ngrouptypes);
	declare_array(double, neutronscatterfactorallatom, ngrouptypes);
	FILE *inp;
	declare_array(char, trajectoryfile, charlength);
	declare_array(char, directory, charlength);
	declare_array(char, filename, charlength);
	declare_array(char, allatommoviefile, charlength);
	declare_array(char, allatommoviesourcefile, charlength);
    allatomparams myallatomparams;
	monomernodes **monomerid;
	
	strcpy(trajectoryfile, loadparam(argc, "-trajectoryfile", argv, "0").s);
	printf("trajectoryfile: %s\n", trajectoryfile);
	strcpy(directory, loadparam(argc, "-directory", argv, "0").s);
	uniformchainlength=loadparam(argc, "-chainlength", argv, "28").i;
	terminuscode=loadparam(argc, "-terminuscode", argv, "0").i;
	Nchains=loadparam(argc, "-Nchains", argv, "96").i;
	
	xrdcode=loadparam(argc, "-xrdcode", argv, "2").i;
	xrdqspace=loadparam(argc, "-xrdqspace", argv, "0.1").f;
	maxq=loadparam(argc, "-maxq", argv, "3").f;
	mincycle=loadparam(argc, "-mincycle", argv, "5000000").i;
	nframes=loadparam(argc, "-nframes", argv, "1").i;
	histoparams.dotwidth=loadparam(argc, "-dotwidth", argv, "0.1").f;
	histoparams.rbottom=loadparam(argc, "-rbottom", argv, "2").f;
	rtop=loadparam(argc, "-rtop", argv, "4.5").f;
	histoparams.rwidth=loadparam(argc, "-rwidth", argv, "0.1").f;
    grbinwidth=loadparam(argc, "-grbinwidth", argv, "0.05").f;
    phenthresholds[0]=loadparam(argc, "-phenthresholdsamepoly", argv, "10").f;
    phenthresholds[1]=loadparam(argc, "-phenthresholdsameleaf", argv, "6.5").f;
    phenthresholds[2]=loadparam(argc, "-phenthresholdoppositeleaf", argv, "6.5").f;
    chargedthreshold=loadparam(argc, "-chargedthreshold", argv, "4.5").f;
    chargedouterthreshold=loadparam(argc, "-chargedouterthreshold", argv, "6").f;
    dpbinwidth=loadparam(argc, "-dpbinwidth", argv, "0.05").f;

    doxrdallatom=loadparam(argc, "-doxrdallatom", argv, "1").i;
    dobackbone=loadparam(argc, "-dobackbone", argv, "1").i;
    dosidechain=loadparam(argc, "-dosidechain", argv, "1").i;
    dogr=loadparam(argc, "-dogr", argv, "1").i;
    doallatomxyz=loadparam(argc, "-doallatomxyz", argv, "0").i;
	leafflag=loadparam(argc, "-leafflag", argv, "0").i;
    dophendotproduct=loadparam(argc, "-dophendotproduct", argv, "1").i;					//	dot product of phenethyl sidechain sites within given range
    dochargeddotproduct=loadparam(argc, "-dochargeddotproduct", argv, "1").i;

	histoparams.rbins=(int) floor((rtop-histoparams.rbottom)/histoparams.rwidth);
	histoparams.dotbins=(int) floor(2/histoparams.dotwidth);
	declare_array(int, histo_r, histoparams.rbins);
	declare_array(int, histo_rleftdotn, histoparams.dotbins);
	declare_array(int, histo_rrightdotn, histoparams.dotbins);
	declare_array(int, histo_tripleproduct, histoparams.dotbins);
	declare_matrix(int, histo_r_rleftdotn, histoparams.rbins, histoparams.dotbins);
	declare_matrix(int, histo_r_rrightdotn, histoparams.rbins, histoparams.dotbins);
	declare_matrix(int, histo_tripleproduct_rleftdotn, histoparams.dotbins, histoparams.dotbins);
	declare_matrix(int, histo_tripleproduct_rrightdotn, histoparams.dotbins, histoparams.dotbins);
    
    sidechainhistoparams sidehistoparams;
    sidehistoparams.bins=100;
    sidehistoparams.rparallelmin=-2;
    sidehistoparams.rparallelwidth=(6-sidehistoparams.rparallelmin)/(1.*sidehistoparams.bins);
    sidehistoparams.rperpwidth=(5.)/(1.*sidehistoparams.bins);
    sidehistoparams.rdotnmin=-1;
    sidehistoparams.rdotnwidth=(6-sidehistoparams.rdotnmin)/(1.*sidehistoparams.bins);
    sidehistoparams.ndotnwidth=(2.)/(1.*sidehistoparams.bins);
    
    declare_matrix(int, histo_rparallel, 4, sidehistoparams.bins);
    declare_matrix(int, histo_rperp, 4, sidehistoparams.bins);
    declare_matrix(int, histo_rdotn, 4, sidehistoparams.bins);
    declare_matrix(int, histo_ndotn, 4, sidehistoparams.bins);
    declare_array(int, sidehistocount, 4);
    
    Nresidues=Nchains*uniformchainlength;                                               //  for purposes of normalizing XRD spectra
	if(terminuscode==1) uniformchainlength++;											//	for purposes of inputting trajectory file
    Nsites=2*Nchains*uniformchainlength;												//  for purposes of searching through coordarray
	
    declare_array_nozero(coord, coordarray, Nsites);
    
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
    
    //  18 atom groups, so 18*(17/2+1)*3=171*3=513 columns in .weighted (2 to 514 in 1d)
    //  1d (transverse) columns for same leaf:
    //  back-back   2, 5, 8, 11, 14, 17, 56, 59, 62, 65, 68, 107, 110, 113, 116, 155, 158, 161, 200, 203, 242,
    //  back-phen   20, 23, 26, 29, 71, 74, 77, 80, 119, 122, 125, 128, 164, 167, 170, 173, 206, 209, 212, 215, 245, 248, 251, 254,
    //  back-amin   32, 35, 38, 41, 83, 86, 89, 92, 131, 134, 137, 140, 176, 179, 182, 185, 218, 221, 224, 227, 257, 260, 263, 266,
    //  back-carb   44, 47, 50, 53, 95, 98, 101, 104, 143, 146, 149, 152, 188, 191, 194, 197, 230, 233, 236, 239, 269, 272, 275, 278,
    //  phen-phen   281, 284, 287, 290, 317, 320, 323, 350, 353, 380
    //  phen-amin   293, 296, 299, 302, 326, 329, 332, 335, 356, 359, 362, 365, 383, 386, 389, 392
    //  phen-carb   305, 308, 311, 314, 338, 341, 344, 347, 368, 371, 374, 377, 395, 398, 401, 404
    //  amin-amin   407, 410, 413, 416, 431, 434, 437, 452, 455, 470
    //  amin-carb   419, 422, 425, 428, 440, 443, 446, 449, 458, 461, 464, 467, 473, 476, 479, 482
    //  carb-carb   485, 488, 491, 494, 497, 500, 503, 506, 509, 512
	
	//	monolayer: 171 columns in .weighted (3 to 173 in 2d)
	//	phen-phen 1d	93, 94, 95, 96, 107, 108, 109, 118, 119, 128
	//	phen-phen 2d	94, 95, 96, 97, 108, 109, 110, 119, 120, 129
	
	//	looked up real part of the atomic scattering factor at 11111 eV (Nam et al powder xrd) at http://henke.lbl.gov/optical_constants/pert_form.html; nearly identical to atomic number; just use atomic number
	
	scatterfactorallatom[0]=scatterfactorallatom[12]=7; // 7.018;
	scatterfactorallatom[1]=scatterfactorallatom[4]=scatterfactorallatom[6]=scatterfactorallatom[8]=scatterfactorallatom[10]=scatterfactorallatom[14]=scatterfactorallatom[16]=6; // 6.010;
	scatterfactorallatom[2]=scatterfactorallatom[17]=8; // 8.030;
	scatterfactorallatom[3]=scatterfactorallatom[5]=scatterfactorallatom[7]=scatterfactorallatom[9]=scatterfactorallatom[11]=scatterfactorallatom[13]=scatterfactorallatom[15]=1; // 1.0000;

	neutronscatterfactorallatom[0]=neutronscatterfactorallatom[12]=11.51;
	neutronscatterfactorallatom[1]=neutronscatterfactorallatom[4]=neutronscatterfactorallatom[6]=neutronscatterfactorallatom[8]=neutronscatterfactorallatom[10]=neutronscatterfactorallatom[14]=neutronscatterfactorallatom[16]=5.551;
	neutronscatterfactorallatom[2]=neutronscatterfactorallatom[17]=4.232;
	neutronscatterfactorallatom[3]=neutronscatterfactorallatom[5]=neutronscatterfactorallatom[7]=neutronscatterfactorallatom[9]=neutronscatterfactorallatom[11]=neutronscatterfactorallatom[13]=neutronscatterfactorallatom[15]=82.02;

    myallatomparams.Cbetaspacing=1.43;
    myallatomparams.NterminusHdepth=-0.5*0.997;
    myallatomparams.NterminusHwidth=-0.5*sqrt(3.)*0.997;
    myallatomparams.Cnyldepth=-0.5*1.345;
    myallatomparams.Cnylwidth=-0.5*sqrt(3.)*1.345;
    myallatomparams.Onyldepth=-0.5*sqrt(3.)-1.23;                   //  assuming trans
    myallatomparams.Onylwidth=-0.5*sqrt(3.)*0.997;
    myallatomparams.CterminusHdepth=0.997*cos(117./180.*M_PI);
    myallatomparams.CterminusHwidth=0.997*sin(117./180.*M_PI);
    myallatomparams.CterminuslongHdistance=1.000;
    myallatomparams.CterminuslongHdepth=-1.000*0.5;
    myallatomparams.CterminuslongHwidth=1.000*0.5*sqrt(3.);
    myallatomparams.Cdepth=-0.5*1.43;
    myallatomparams.Cwidth=0.5*sqrt(3.)*1.43;
    myallatomparams.CHdepth=-0.5*1.43-0.5*1.08;
    myallatomparams.CHwidth=0.5*sqrt(3.)*1.43;
    myallatomparams.CHoutofplane=0.5*sqrt(2.)*1.08;                 //  tetrahedral, assuming trans
	myallatomparams.phenylspacing=1.375;
	myallatomparams.phenylHspacing=1.375+1.08;
	myallatomparams.phenylCgammaspacing=1.375+1.49;
	myallatomparams.aminoNdepth=0.696;
	myallatomparams.aminoCdepth=-0.784;
	myallatomparams.aminoHdepth=0.696+1.04/3.;
	myallatomparams.aminoHwidth=1.04*sqrt(2.)*2./3.;
	myallatomparams.carboxylbackCdepth=-1.371;
	myallatomparams.carboxylforwardCdepth=0.151;
	myallatomparams.carboxylOdepth=0.621;
	myallatomparams.carboxylOwidth=0.883;
    myallatomparams.ethylHdepth=0.5*1.111;
    myallatomparams.ethylHoutofplane=0.5*sqrt(3.)*1.111;
	myallatomparams.Crad=0.7;
	myallatomparams.Orad=0.6;
	myallatomparams.Nrad=0.65;
    myallatomparams.Hrad=0.25;      //	atomic radii

	assign_monomerid_and_chainlength(Nchains, uniformchainlength, &monomerid, &chainlength);	
	
    sprintf(allatommoviefile, "%s/vmd.allatom.fromtrajectory.xyz", directory);
	sprintf(allatommoviesourcefile, "%s/vmd.allatom.fromtrajectory.source", directory);

    //  calculate totalframes
    
 	inp=fopen(trajectoryfile, "r");
    read=input_trajectory_blank(Nchains, uniformchainlength, inp, charlength);
    while(read>0){
        totalframes++;
        read=input_trajectory_blank(Nchains, uniformchainlength, inp, charlength);
    }
    fclose(inp);
	
	printf("total frames %i\n", totalframes);

    //  calculate average box size during production interval
    
    inp=fopen(trajectoryfile, "r");
    framecounter=0;
    read=input_trajectory_box(Nchains, uniformchainlength, inp, charlength, &cycle, &box_dimension);
    while(read>0){
        if((cycle>=mincycle)&&(framecounter>=totalframes-nframes)){
            avbox=add_double_triple(avbox, box_dimension);
            xrdcount++;
        }
        read=input_trajectory_box(Nchains, uniformchainlength, inp, charlength, &cycle, &box_dimension);
        framecounter++;
    }
    fclose(inp);

	inp=fopen(trajectoryfile, "r");
    if(leafflag>0) read=input_trajectory_v0(Nchains, uniformchainlength, coordarray, &box_dimension, inp, &cycle, charlength, leafflag);
    else read=input_trajectory(Nchains, uniformchainlength, coordarray, &box_dimension, inp, &cycle, charlength);
    
    for(i=0;i<Nsites;i++){
        if(coordarray[i].nodetype>=0){
            if(coordarray[i].leafid>Nleaves) Nleaves=coordarray[i].leafid;
        }
    }
    Nleaves++;
    if(Nleaves<3) symmetries=Nleaves;
    else symmetries=3;
	
	double padfactor=1.2;
    printf("average box dimensions %f, %f, %f\n", avbox.x/(1.*xrdcount), avbox.y/(1.*xrdcount), avbox.z/(1.*xrdcount));
	if(xrdcode==2){
		xrdqspacearray.x=0;
		xrdqspacearray.y=0;
		xrdqspacearray.z=xrdqspace;
		xrdqnarray.x=(int) avbox.x/(1.*xrdcount)*maxq*padfactor/(2*M_PI);
		xrdqnarray.y=(int) avbox.y/(1.*xrdcount)*maxq*padfactor/(2*M_PI);
		xrdqnarray.z=(int) maxq/xrdqspace*padfactor;
	}
	else if(xrdcode==3){
		xrdqspacearray.x=0;
		xrdqspacearray.y=0;
		xrdqspacearray.z=0;
		xrdqnarray.x=(int) avbox.x/(1.*xrdcount)*maxq*padfactor/(2*M_PI);
		xrdqnarray.y=(int) avbox.y/(1.*xrdcount)*maxq*padfactor/(2*M_PI);
		xrdqnarray.z=(int) avbox.z/(1.*xrdcount)*maxq*padfactor/(2*M_PI);
	}
	else my_exit("haven't written code for that xrdcode");
	
    declare_5d_tensor(double, Iqhistoallatomxy, 2*xrdqnarray.x+1, 2*xrdqnarray.y+1, ngrouptypes, ngrouptypes, symmetries);
    declare_4d_tensor(double, Iqhistoallatomz, xrdqnarray.z+1, ngrouptypes, ngrouptypes, symmetries);
		
	allatoms=count_allatom_atoms(Nchains, chainlength, monomerid, coordarray);
	declare_array(int, atomid, allatoms);
	declare_array(int, leaves, allatoms);
	declare_array_nozero(double_triple, position, allatoms);
	if(doxrdallatom==1){
		if(count_allatom_atoms_atomid(Nchains, chainlength, monomerid, coordarray, atomid, leaves)!=allatoms) my_exit("allatom count doesn't agree");
	}
	printf("allatoms=%i\n", allatoms);
	
    double maxr=0.5*norm(scalar_multiply_double_triple(avbox, (1./(1.*xrdcount))))+2*grbinwidth;
    int maxgrbin=(int) floor(maxr/grbinwidth)+1;
    declare_4d_tensor(double, grhisto, maxgrbin, nsitetypes, nsitetypes, 3);
    int grcount=0;
    declare_array(double, fractions, nsitetypes);
    fractions[0]=1.;
    fractions[1]=0.5;
    fractions[2]=fractions[3]=0.25;
    
    int dpbins=(int) floor(2./dpbinwidth);
	int dpbinsunsigned=(int) floor(1./dpbinwidth);
	int dpbinsunsignedadd=(int) floor(2./dpbinwidth);
    declare_matrix(int, phenhistos, 3, dpbinsunsigned);
    declare_matrix(int, phenmixedaddhistos, 3, dpbinsunsignedadd);
    declare_matrix(int, phenmixedsubtracthistos, 3, dpbinsunsignedadd);
    declare_array(int, phencounts, 3);
    declare_array(int, chargedhisto, dpbins);
    int chargedcount=0;
    declare_array(int, chargedouterhisto, dpbins);
    int chargedoutercount=0;

    framecounter=0;
    xrdcount=0;
    avbox.x=avbox.y=avbox.z=0;      //  no need to recalculate, but doing it anyway
    while(read>0){
        if(doallatomxyz==1) output_xyz_allatom(Nchains, chainlength, monomerid, allatoms, coordarray, box_dimension, allatommoviefile, allatommoviesourcefile, myallatomparams, &allatomframe);    //  haven't written com version
        if((cycle>=mincycle)&&(framecounter>=totalframes-nframes)){
            if(foundfirstcycle==0) firstcycle=cycle;
            foundfirstcycle=1;
            printf("frame %i of %i\n", framecounter, totalframes);
			if(doxrdallatom==1){
                poormansxrd_fromcoord_bin_leafid_allatom_no3d_resolutionlimited(Nchains, coordarray, box_dimension, Iqhistoallatomxy, Iqhistoallatomz, xrdqspacearray, xrdqnarray, ngrouptypes, &xrdcount, allatoms, atomid, leaves, position, chainlength, monomerid, myallatomparams, Nleaves, &avbox);
			}
			if(dobackbone==1) backbone_histo(Nchains, uniformchainlength, coordarray, box_dimension, histoparams, histo_r, histo_rleftdotn, histo_rrightdotn, histo_tripleproduct, histo_r_rleftdotn, histo_r_rrightdotn, histo_tripleproduct_rleftdotn, histo_tripleproduct_rrightdotn, &count_r, &count_rleftdotn, &count_rrightdotn, &count_tripleproduct, &count_r_rleftdotn, &count_r_rrightdotn, &count_tripleproduct_rleftdotn, &count_tripleproduct_rrightdotn);
			if(dosidechain==1) sidechain_histo(Nchains, uniformchainlength, coordarray, box_dimension, sidehistoparams, histo_rparallel, histo_rperp, histo_rdotn, histo_ndotn, sidehistocount);
			if(dogr==1) gr_fromcoord(Nsites, coordarray, box_dimension, grbinwidth, grhisto, &grcount, Nchains, maxgrbin);
			if(dophendotproduct==1) dotproducthistoalltypes(Nsites, coordarray, box_dimension, Nchains, phenthresholds, dpbinwidth, phenhistos, phenmixedaddhistos, phenmixedsubtracthistos, phencounts, 1, 0);
			if(dochargeddotproduct==1) dotproducthistodifferenttypetwoshells(Nsites, coordarray, box_dimension, Nchains, chargedthreshold, chargedouterthreshold, dpbinwidth, chargedhisto, &chargedcount, chargedouterhisto, &chargedoutercount, 2, 3);
		}
        if(leafflag>0) read=input_trajectory_v0(Nchains, uniformchainlength, coordarray, &box_dimension, inp, &cycle, charlength, leafflag);
        else read=input_trajectory(Nchains, uniformchainlength, coordarray, &box_dimension, inp, &cycle, charlength);
        framecounter++;
    }
    fclose(inp);

    mincycle=firstcycle;
    if(dobackbone==1){
        sprintf(filename, "%s/histo.%i-%i", directory, mincycle, cycle);
        output_backbone_histo(filename, histoparams, histo_r, histo_rleftdotn, histo_rrightdotn, histo_tripleproduct, histo_r_rleftdotn, histo_r_rrightdotn, histo_tripleproduct_rleftdotn, histo_tripleproduct_rrightdotn, &count_r, &count_rleftdotn, &count_rrightdotn, &count_tripleproduct, &count_r_rleftdotn, &count_r_rrightdotn, &count_tripleproduct_rleftdotn, &count_tripleproduct_rrightdotn);
	}
    if(dosidechain==1){
        sprintf(filename, "%s/sidehisto.%i-%i", directory, mincycle, cycle);
        output_sidechain_histo(filename, sidehistoparams, histo_rparallel, histo_rperp, histo_rdotn, histo_ndotn, sidehistocount);
    }
    if(doxrdallatom==1){
        sprintf(filename, "%s/xrd.allatom.%i-%i.faceon", directory, mincycle, cycle);
        output_2d_xrd_leaves_resolutionlimited_andradav(filename, Iqhistoallatomxy, xrdqspacearray, xrdqnarray, ngrouptypes, xrdcount, Nchains*uniformchainlength, charlength, scatterfactorallatom, neutronscatterfactorallatom, symmetries, avbox, xrdqspace);
        sprintf(filename, "%s/xrd.allatom.%i-%i.transverse", directory, mincycle, cycle);
        output_1d_xrd_leaves_resolutionlimited(filename, Iqhistoallatomz, xrdqspacearray, xrdqnarray, ngrouptypes, xrdcount, Nchains*uniformchainlength, charlength, scatterfactorallatom, neutronscatterfactorallatom, symmetries, avbox);
	}
    if(dogr==1){
        sprintf(filename, "%s/gr.%i-%i.samepoly.3dnorm", directory, mincycle, cycle);
        output_gr_3dnorm(filename, grhisto, nsitetypes, grbinwidth, maxgrbin, grcount, fractions, 0);
        sprintf(filename, "%s/gr.%i-%i.samepoly.2dnorm", directory, mincycle, cycle);
        output_gr_2dnorm(filename, grhisto, nsitetypes, grbinwidth, maxgrbin, grcount, fractions, 0);
        sprintf(filename, "%s/gr.%i-%i.sameleaf.3dnorm", directory, mincycle, cycle);
        output_gr_3dnorm(filename, grhisto, nsitetypes, grbinwidth, maxgrbin, grcount, fractions, 1);
        sprintf(filename, "%s/gr.%i-%i.sameleaf.2dnorm", directory, mincycle, cycle);
        output_gr_2dnorm(filename, grhisto, nsitetypes, grbinwidth, maxgrbin, grcount, fractions, 1);
        sprintf(filename, "%s/gr.%i-%i.oppositeleaf.3dnorm", directory, mincycle, cycle);
        output_gr_3dnorm(filename, grhisto, nsitetypes, grbinwidth, maxgrbin, grcount, fractions, 2);
        sprintf(filename, "%s/gr.%i-%i.oppositeleaf.2dnorm", directory, mincycle, cycle);
        output_gr_2dnorm(filename, grhisto, nsitetypes, grbinwidth, maxgrbin, grcount, fractions, 2);
    }
    if(dophendotproduct==1){
        sprintf(filename, "%s/phenhisto.%i-%i", directory, mincycle, cycle);
        output_unsigned_dotproduct_severaltypes(filename, dpbinsunsigned, phenhistos, phencounts, dpbinwidth, 3, 1);
        sprintf(filename, "%s/phenmixedaddhisto.%i-%i", directory, mincycle, cycle);
        output_unsigned_dotproduct_severaltypes(filename, dpbinsunsignedadd, phenmixedaddhistos, phencounts, dpbinwidth, 3, 2);
        sprintf(filename, "%s/phenmixedsubtracthisto.%i-%i", directory, mincycle, cycle);
        output_unsigned_dotproduct_severaltypes(filename, dpbinsunsignedadd, phenmixedsubtracthistos, phencounts, dpbinwidth, 3, 2);
    }
    if(dochargeddotproduct==1){
        sprintf(filename, "%s/chargedhisto.%i-%i", directory, mincycle, cycle);
        output_dotproduct(filename, dpbins, chargedhisto, chargedcount, dpbinwidth);
        sprintf(filename, "%s/chargedouterhisto.%i-%i", directory, mincycle, cycle);
        output_dotproduct(filename, dpbins, chargedouterhisto, chargedoutercount, dpbinwidth);
    }
	return 0;
}
	

