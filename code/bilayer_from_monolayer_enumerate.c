#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"
#include "mymath.h"
#include "peptoid_functions.h"
#include "search_functions.h"

int main(int argc, char *argv[]){
		
	int charlength=200, i, j, k, chemistrymax=20, sidechaintypesmax=10, chaintypes, *Nchainsoftype, *chainlengthoftype, **chainchemistry, Nchains, Nmonomers, Nnodes, *chainlength, Nnodetypes, leafsymmetry, terminuscode, biasaminotrans;
	double kelvins, saltmolarity, temperature, newenergy;
	char **chemistry_char;
	linkedlistset linkset;
	chemistry_char=xmalloc(sizeof(char *)*chemistrymax);
	for(i=0;i<chemistrymax;i++) chemistry_char[i]=xmalloc(sizeof(char)*charlength);
    double_triple box_dimension, shift;
	coord *coordarray;
	monomernodes **monomerid;
	declare_array(char, inputfile, charlength);	
	declare_array(char, outputfile, charlength);	
	bonded_params my_bonded_params;
	nonbonded_params my_nonbonded_params;
	energycomponents my_energycomponents, bestmy_energycomponents;
	monomertypes monomercount;
	monomercount.type=malloc(sizeof(double)*sidechaintypesmax+1);

    biasaminotrans=loadparam(argc, "-biasaminotrans", argv, "1").i;

	my_bonded_params.K10=13595.3;
	my_bonded_params.K11=-15881.3;
	my_bonded_params.K12=6946.51;
	my_bonded_params.K13=-1348.29;
	my_bonded_params.K14=97.9797;
	my_bonded_params.K20=2.36497;
	my_bonded_params.K21=-3.87853;
	my_bonded_params.K22=0.463066;
	my_bonded_params.K23=-3.33253;
	my_bonded_params.K24=5.34237;
	my_bonded_params.kl=28.2832;
	my_bonded_params.r0l=3.36096;
	my_bonded_params.sl=1.95114;
	my_bonded_params.kr=12.6891;
	my_bonded_params.r0r=1.51312;
	my_bonded_params.sr=-4.10588;
 	my_bonded_params.sidechain=xmalloc(sizeof(sidechain_params)*sidechaintypesmax+1);
	my_bonded_params.sidechain[1].k1=23.086;
	my_bonded_params.sidechain[1].rperp0=0.343266;
	my_bonded_params.sidechain[1].rpar0=1.60786;
	my_bonded_params.sidechain[1].r0=3.44926;
	my_bonded_params.sidechain[1].J10=13.942-1.88336+2.84809-0.853415;		//	+ J200 + J300 - numerical minimum
	my_bonded_params.sidechain[1].J11=-82.0459;
	my_bonded_params.sidechain[1].J12=186.411;
	my_bonded_params.sidechain[1].J13=-164.93;
	my_bonded_params.sidechain[1].J14=49.0435;
	my_bonded_params.sidechain[1].J201=-2.02912;
	my_bonded_params.sidechain[1].J202=12.0487;
	my_bonded_params.sidechain[1].J203=-10.0584;
	my_bonded_params.sidechain[1].J204=3.44353;
	my_bonded_params.sidechain[1].J205=-0.536649;
	my_bonded_params.sidechain[1].J206=0.0315484;
	my_bonded_params.sidechain[1].J220=1.61338-13.0286;			//	+ J302
	my_bonded_params.sidechain[1].J221=-4.38724;
	my_bonded_params.sidechain[1].J222=1.24471;
	my_bonded_params.sidechain[1].J240=0.663689+33.5248;		//	+ J304
	my_bonded_params.sidechain[1].J311=9.4368;
	my_bonded_params.sidechain[1].J313=-45.2463;
	my_bonded_params.sidechain[1].J320=-2.36074;
	my_bonded_params.sidechain[1].J322=24.2649;
	my_bonded_params.sidechain[1].J331=-6.30323;
	my_bonded_params.sidechain[1].J340=0.711079;
	my_bonded_params.sidechain[2].k1=74.463;
	my_bonded_params.sidechain[2].rperp0=-0.891764;
	my_bonded_params.sidechain[2].rpar0=0.999481;
	my_bonded_params.sidechain[2].r0=2.93893;
    if(biasaminotrans==1){
        my_bonded_params.sidechain[2].J10=424.643;
        my_bonded_params.sidechain[2].J11=-1844.99;
    }
    else{
        my_bonded_params.sidechain[2].J10=418.586;
        my_bonded_params.sidechain[2].J11=-1837.99;
    }
	my_bonded_params.sidechain[2].J12=2964.05;
	my_bonded_params.sidechain[2].J13=-2079.92;
	my_bonded_params.sidechain[2].J14=537.819;
	my_bonded_params.sidechain[2].k2=12.8785;
	my_bonded_params.sidechain[2].r20=-0.0587984;
	my_bonded_params.sidechain[2].r21=-1.024;
	my_bonded_params.sidechain[2].r22=0.476617;
	my_bonded_params.sidechain[2].k3=50.2083;
	my_bonded_params.sidechain[2].r30=1.3901;
	my_bonded_params.sidechain[2].r31=1.55449;
	my_bonded_params.sidechain[2].r32=-0.155606;
	my_bonded_params.sidechain[3].k1=51.4346;
	my_bonded_params.sidechain[3].rperp0=0.389628;
	my_bonded_params.sidechain[3].rpar0=1.63953;
	my_bonded_params.sidechain[3].r0=1.98906-0.0574095;         //  minus minimum due to npar<1 preventing reaching rpar, rperp minimum
	my_bonded_params.sidechain[3].J10=193.754;
	my_bonded_params.sidechain[3].J11=-937.961;
	my_bonded_params.sidechain[3].J12=1671.22;
	my_bonded_params.sidechain[3].J13=-1279.85;
	my_bonded_params.sidechain[3].J14=355.337;
	my_bonded_params.sidechain[3].k2=30.9392;
	my_bonded_params.sidechain[3].r20=-2.08987;
	my_bonded_params.sidechain[3].r21=1.19511;
	my_bonded_params.sidechain[3].r22=-0.081282;
	my_bonded_params.sidechain[3].k3=23.0679;
	my_bonded_params.sidechain[3].r30=1.94004;
	my_bonded_params.sidechain[3].r31=1.58382;
	my_bonded_params.sidechain[3].r32=-0.210233;
 	my_bonded_params.sidechainorientationtype=xmalloc(sizeof(int)*sidechaintypesmax+1);
	my_bonded_params.sidechainorientationtype[1]=1;		//	phenyl
	my_bonded_params.sidechainorientationtype[2]=0;		//	amino
	my_bonded_params.sidechainorientationtype[3]=0;		//	carboxyl
    
    my_nonbonded_params.rhard0=loadparam(argc, "-backboneradius", argv, "2.05").f;
	my_nonbonded_params.one0.solvationenergy=10.1;
	my_nonbonded_params.one0.z0=-0.5;
	my_nonbonded_params.one0.uinterface=loadparam(argc, "-backbonesurfaceenergy", argv, "1.9").f;
	my_nonbonded_params.one0.zinterface=-0.5;
	my_nonbonded_params.one0.sigmainterface=0.6;
	my_nonbonded_params.rhard1=0.5*0.54*5.6;							//	phenyl
	my_nonbonded_params.one1.solvationenergy=0.76;						//	toluene
	my_nonbonded_params.one1.z0=0;
	my_nonbonded_params.one1.uinterface=loadparam(argc, "-phenylsurfaceenergy", argv, "3.6").f;	//	toluene
    my_nonbonded_params.one1.zinterface=1.2;
	my_nonbonded_params.one1.sigmainterface=loadparam(argc, "-phenylinterfacewidth", argv, "1.6").f;						//	toluene
	my_nonbonded_params.p11vac.eps0=1.8;								//	phenyl-phenyl
    my_nonbonded_params.p11vac.sigma0=loadparam(argc, "-phenylsigma", argv, "6.2").f;
	my_nonbonded_params.p11vac.xi=0.5;
	my_nonbonded_params.p11vac.Q=0.6;
	my_nonbonded_params.p11vac.chi=(0.54*0.54-1.)/(0.54*0.54+1.);			//	(kappa^2 - 1)/(kappa^2 + 1)
	my_nonbonded_params.p11vac.chiprime=(1.-3.6)/(1.+3.6);					//	(1 - kappaprime^(1/mu))/(1 + kappaprime^(1/mu)), fixing mu=1, nu=-2
    my_nonbonded_params.phenylcutoff2=8.49668*8.49668;              //	numerical solution to U=-0.1 eV (cutoff_v2.nb)
    
    
	my_nonbonded_params.p11vac.code=loadparam(argc, "-phenylinteractioncode", argv, "2").i;		//	0: vacuum; 1: solvated; 2: mixture
	my_nonbonded_params.p11vac.p11solv.eps0=2;								//	phenyl-phenyl
	my_nonbonded_params.p11vac.p11solv.sigma0=my_nonbonded_params.p11vac.sigma0;
	my_nonbonded_params.p11vac.p11solv.xi=0.5;
	my_nonbonded_params.p11vac.p11solv.chi=(0.54*0.54-1.)/(0.54*0.54+1.);			//	(kappa^2 - 1)/(kappa^2 + 1)
	my_nonbonded_params.p11vac.p11solv.chiprime=(1.-10)/(1.+10);					//	(1 - kappaprime^(1/mu))/(1 + kappaprime^(1/mu)), fixing mu=1, nu=-2
	my_nonbonded_params.p11vac.p11solv.w=2.8;
	my_nonbonded_params.p11vac.p11solv.amprep=1;
	my_nonbonded_params.p11vac.p11solv.chirep=(1.-25)/(1.+25);						//	(1 - kappaprime^(1/mu))/(1 + kappaprime^(1/mu)), fixing mu=1, nu=-2
	my_nonbonded_params.p11vac.p11solv.ampattr=0.3;
	my_nonbonded_params.p11vac.p11solv.chiattr=(1.-4)/(1.+4);						//	(1 - kappaprime^(1/mu))/(1 + kappaprime^(1/mu)), fixing mu=1, nu=-2
	my_nonbonded_params.p11vac.p11solv.interpolatemiddle=my_nonbonded_params.one1.zinterface;
	my_nonbonded_params.p11vac.p11solv.interpolatewidth=my_nonbonded_params.one1.sigmainterface;			//	should scale with, but could have prefactor
    if(my_nonbonded_params.p11vac.code>0) my_nonbonded_params.phenylcutoff2=10.98*10.98;							//	numerical solution to U=-0.1 eV (toluene_toluene.nb)
    
    
	my_nonbonded_params.solvationparams.shortrangepermittivity=20;      //  not 10
	my_nonbonded_params.solvationparams.waterradius=1.5;               //   not 0.6
	my_nonbonded_params.solvationparams.firstshellfraction=0.01;        //  not 0.3
	my_nonbonded_params.solvationparams.saturationrange=4;      //   not 2.5
	my_nonbonded_params.solvationparams.interfacethickness=1;		//	fit logistic function to water profile from Shaytan et al, J. Comp. Chem. 31 204 (2010)
    my_nonbonded_params.p2.charge=1;								//	amino group
	my_nonbonded_params.p2.dipole=0.48;
	my_nonbonded_params.p2.rhard=1.82;
	my_nonbonded_params.p2.solvationenergy=71.3;
    my_nonbonded_params.p2.shift=loadparam(argc, "-aminoshift", argv, "0.696").f;
	my_nonbonded_params.one2.solvationenergy=71.3;
	my_nonbonded_params.one2.z0=-1.5;
	my_nonbonded_params.one2.uinterface=0;
	my_nonbonded_params.one2.zinterface=0;
	my_nonbonded_params.one2.sigmainterface=0;
	my_nonbonded_params.p3.charge=-1;								//	carboxyl group
	my_nonbonded_params.p3.dipole=-0.81;
	my_nonbonded_params.p3.rhard=2.06;
	my_nonbonded_params.p3.solvationenergy=79.9;
    my_nonbonded_params.p3.shift=loadparam(argc, "-carboxylshift", argv, "0.151").f;
	my_nonbonded_params.one3.solvationenergy=79.9;
	my_nonbonded_params.one3.z0=-3.5;
	my_nonbonded_params.one3.uinterface=loadparam(argc, "-carboxylsurfaceenergy", argv, "1.5").f;
	my_nonbonded_params.one3.zinterface=-1.5;						//	Jungwirth, based on toluene water profile
	my_nonbonded_params.one3.sigmainterface=0.6;					//	Jungwirth, based on toluene water profile

	leafsymmetry=loadparam(argc, "-leafsymmetry", argv, "1").i;
	strcpy(inputfile, loadparam(argc, "-inputfile", argv, "0").s);
	strcpy(outputfile, loadparam(argc, "-outputfile", argv, "0").s);
	double firstx=loadparam(argc, "-firstx", argv, "0").f;
	double firsty=loadparam(argc, "-firsty", argv, "0").f;
	double firstz=loadparam(argc, "-firstz", argv, "5").f;
	double lastx=loadparam(argc, "-lastx", argv, "2").f;
	double lasty=loadparam(argc, "-lasty", argv, "2").f;
	double lastz=loadparam(argc, "-lastz", argv, "20").f;
	double xspace=loadparam(argc, "-xspace", argv, "0.1").f;
	double yspace=loadparam(argc, "-yspace", argv, "0.1").f;
	double zspace=loadparam(argc, "-zspace", argv, "0.1").f;

	//	System parameters
	
	kelvins=loadparam(argc, "-temperature", argv, "300").f;			//	in degrees Kelvin; default room T
    my_bonded_params.factor=loadparam(argc, "-bondedfactor", argv, "1.").f;
    my_nonbonded_params.p11vac.factor=loadparam(argc, "-phenylfactor", argv, "1.5").f;
    my_nonbonded_params.p2.factor=sqrt(loadparam(argc, "-chargedfactor", argv, "1.").f);
    my_nonbonded_params.p3.factor=sqrt(loadparam(argc, "-chargedfactor", argv, "1.").f);
	saltmolarity=loadparam(argc, "-saltmolarity", argv, "0.2").f;
	my_nonbonded_params.cutoff2=pow(loadparam(argc, "-cutoff", argv, "13.27").f, 2);			//	numerical solution to lambda_Debye*W(Coulombconstant/(waterpermittivity*threshold*lambda_Debye)) (cutoff.nb)
	loadseries(argc, "-chemistry", argv, chemistry_char);			
	terminuscode=loadparam(argc, "-terminuscode", argv, "0").i;
	my_nonbonded_params.solvationparams.interface=loadparam(argc, "-interface", argv, "0").i;
	my_nonbonded_params.vacuumthickness=loadparam(argc, "-vacuumthickness", argv, "20").f;	
	
	//syntax e.g.: -chemistry chaintypes patterntype1(0:alternating) Ntype1 length1(total monomers) nonpolar1 firstpolar1 secondpolar1 ...
	//or: -chemistry chaintypes patterntype1(1:block) Ntype1 polarblocks1 nonpolar1 polar1block1type polar1block1length(dimers) ...

	//	Derived parameters
	
	temperature=onedegreeinkcals*kelvins;
	my_nonbonded_params.solvationparams.debyelength=sqrt(waterpermittivity*temperature/(4*M_PI*coulombconstant*2*saltmolarity*avogadrosnumberover10raised27));
	parse_chemistry(chemistry_char, &chaintypes, &Nchainsoftype, &chainlengthoftype, &chainchemistry, &Nmonomers, &Nchains, &monomercount, sidechaintypesmax, &Nnodetypes, terminuscode);
	declare_array(double, avheight, Nnodetypes);
	declare_array(double, bestavheight, Nnodetypes);
	my_energycomponents.solvation=malloc(sizeof(double)*Nnodetypes+1);
	my_energycomponents.interface=malloc(sizeof(double)*Nnodetypes+1);
	bestmy_energycomponents.solvation=malloc(sizeof(double)*Nnodetypes+1);
	bestmy_energycomponents.interface=malloc(sizeof(double)*Nnodetypes+1);

	//	Configuration parameters
	
	bilayer config, newconfig;
	monolayer my_monolayer;
			
	config.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
	config.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;		//	0 perp, 1 par, 2 alt (initially)
	config.offsetfraction=0.5;

	input_monolayer(inputfile, &my_monolayer, chainlengthoftype[0], Nnodetypes, 0);	
	shift=my_monolayer.backboneoffset[2][1];
	for(i=0;i<Nnodetypes;i++){
		for(j=0;j<Nnodetypes;j++){
			(my_monolayer.backboneoffset[i][j].x)-=(shift.x);
			(my_monolayer.backboneoffset[i][j].y)-=(shift.y);
			(my_monolayer.backboneoffset[i][j].z)-=(shift.z);
		}
	}
	
	config.terminusspacing=my_monolayer.terminusspacing;
	config.monomerspacing=my_monolayer.monomerspacing;
	config.interchainspacing=my_monolayer.interchainspacing;
	config.monomerscaleoffsetfraction=my_monolayer.monomerscaleoffsetfraction;

	double leafspacing=loadparam(argc, "-leafspacing", argv, "10").f;
	double leafxoffsetfrac=loadparam(argc, "-leafxoffsetfrac", argv, "1").f;
	double leafyoffsetfrac=loadparam(argc, "-leafyoffsetfrac", argv, "1").f;

	config.length=chainlengthoftype[0]*config.monomerspacing+config.terminusspacing;
	config.offset=config.length*config.offsetfraction+config.monomerscaleoffsetfraction*config.monomerspacing;
    
    allocate_3d_tensor_nozero(double_triple, config.backboneoffset, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, config.backbonenz, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, config.backbonephi, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor_nozero(double_triple, config.sidechainrelativepos, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, config.sidechainnpar, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, config.sidechainnphi, 2, Nnodetypes, Nnodetypes);
	allocate_3d_tensor_nozero(double_triple, newconfig.backboneoffset, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, newconfig.backbonenz, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, newconfig.backbonephi, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor_nozero(double_triple, newconfig.sidechainrelativepos, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, newconfig.sidechainnpar, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, newconfig.sidechainnphi, 2, Nnodetypes, Nnodetypes);
 
    for(i=0;i<2;i++){
        config.backboneoffset[i][0][2].z=(1-i*2)*my_monolayer.backboneoffset[0][2].z+(i-0.5)*leafspacing;
        config.backboneoffset[i][0][3].z=(1-i*2)*my_monolayer.backboneoffset[0][3].z+(i-0.5)*leafspacing;
        config.backboneoffset[i][2][1].z=(1-i*2)*my_monolayer.backboneoffset[2][1].z+(i-0.5)*leafspacing;
        config.backboneoffset[i][3][1].z=(1-i*2)*my_monolayer.backboneoffset[3][1].z+(i-0.5)*leafspacing;						//	mirror image
        config.sidechainrelativepos[i][2][1].z=my_monolayer.sidechainrelativepos[2][1].z;
		config.sidechainrelativepos[i][3][1].z=my_monolayer.sidechainrelativepos[3][1].z;
        config.sidechainrelativepos[i][0][2].z=my_monolayer.sidechainrelativepos[0][2].z;
        config.sidechainrelativepos[i][0][3].z=my_monolayer.sidechainrelativepos[0][3].z;
		if((leafsymmetry==-1)&&(i==1)){
			config.backboneoffset[i][0][2].x=-my_monolayer.backboneoffset[0][2].x+i*leafxoffsetfrac*config.monomerspacing;
			config.backboneoffset[i][0][3].x=-my_monolayer.backboneoffset[0][3].x+i*leafxoffsetfrac*config.monomerspacing;
			config.backboneoffset[i][2][1].x=-my_monolayer.backboneoffset[2][1].x+i*leafxoffsetfrac*config.monomerspacing;
			config.backboneoffset[i][3][1].x=-my_monolayer.backboneoffset[3][1].x+i*leafxoffsetfrac*config.monomerspacing;
		}
		else{
			config.backboneoffset[i][0][2].x=my_monolayer.backboneoffset[0][2].x+i*leafxoffsetfrac*config.monomerspacing;
			config.backboneoffset[i][0][3].x=my_monolayer.backboneoffset[0][3].x+i*leafxoffsetfrac*config.monomerspacing;
			config.backboneoffset[i][2][1].x=my_monolayer.backboneoffset[2][1].x+i*leafxoffsetfrac*config.monomerspacing;
			config.backboneoffset[i][3][1].x=my_monolayer.backboneoffset[3][1].x+i*leafxoffsetfrac*config.monomerspacing;
		}
        config.backboneoffset[i][0][2].y=my_monolayer.backboneoffset[0][2].y+i*leafyoffsetfrac*config.interchainspacing;
        config.backboneoffset[i][0][3].y=my_monolayer.backboneoffset[0][3].y+i*leafyoffsetfrac*config.interchainspacing;
        config.backboneoffset[i][2][1].y=my_monolayer.backboneoffset[2][1].y+i*leafyoffsetfrac*config.interchainspacing;
        config.backboneoffset[i][3][1].y=my_monolayer.backboneoffset[3][1].y+i*leafyoffsetfrac*config.interchainspacing;		//	0.5: alternating; 1: unlike charges face each other across bilayer
        config.backbonenz[i][0][2]=my_monolayer.backbonenz[0][2];
        config.backbonenz[i][0][3]=my_monolayer.backbonenz[0][3];
        config.backbonenz[i][2][1]=my_monolayer.backbonenz[2][1];
        config.backbonenz[i][3][1]=my_monolayer.backbonenz[3][1];                                      //  perpendicular, keeping phi at 0
		if((i==1)&&(leafsymmetry==-1)){
			config.backbonephi[i][0][2]=M_PI-my_monolayer.backbonephi[0][2];
			config.backbonephi[i][0][3]=M_PI-my_monolayer.backbonephi[0][3];
			config.backbonephi[i][2][1]=M_PI-my_monolayer.backbonephi[2][1];
			config.backbonephi[i][3][1]=M_PI-my_monolayer.backbonephi[3][1];
		}
		else{
			config.backbonephi[i][0][2]=my_monolayer.backbonephi[0][2];
			config.backbonephi[i][0][3]=my_monolayer.backbonephi[0][3];
			config.backbonephi[i][2][1]=my_monolayer.backbonephi[2][1];
			config.backbonephi[i][3][1]=my_monolayer.backbonephi[3][1];
		}			
		if((i==1)&&(leafsymmetry==-1)){
			config.sidechainrelativepos[i][2][1].x=-my_monolayer.sidechainrelativepos[2][1].x;
			config.sidechainrelativepos[i][3][1].x=-my_monolayer.sidechainrelativepos[3][1].x;
			config.sidechainrelativepos[i][0][2].x=-my_monolayer.sidechainrelativepos[0][2].x;
			config.sidechainrelativepos[i][0][3].x=-my_monolayer.sidechainrelativepos[0][3].x;
		}
		else{
			config.sidechainrelativepos[i][2][1].x=my_monolayer.sidechainrelativepos[2][1].x;
			config.sidechainrelativepos[i][3][1].x=my_monolayer.sidechainrelativepos[3][1].x;
			config.sidechainrelativepos[i][0][2].x=my_monolayer.sidechainrelativepos[0][2].x;
			config.sidechainrelativepos[i][0][3].x=my_monolayer.sidechainrelativepos[0][3].x;
		}
		config.sidechainrelativepos[i][2][1].y=my_monolayer.sidechainrelativepos[2][1].y;
		config.sidechainrelativepos[i][3][1].y=my_monolayer.sidechainrelativepos[3][1].y;
        config.sidechainrelativepos[i][0][2].y=my_monolayer.sidechainrelativepos[0][2].y;
        config.sidechainrelativepos[i][0][3].y=my_monolayer.sidechainrelativepos[0][3].y;
        config.sidechainnpar[i][2][1]=my_monolayer.sidechainnpar[2][1];
		config.sidechainnpar[i][3][1]=my_monolayer.sidechainnpar[3][1];
        config.sidechainnpar[i][0][2]=my_monolayer.sidechainnpar[0][2];
        config.sidechainnpar[i][0][3]=my_monolayer.sidechainnpar[0][3];
		if((i==1)&&(leafsymmetry==-1)){
			config.sidechainnphi[i][2][1]=M_PI-my_monolayer.sidechainnphi[2][1];
			config.sidechainnphi[i][3][1]=M_PI-my_monolayer.sidechainnphi[3][1];
			config.sidechainnphi[i][0][2]=M_PI-my_monolayer.sidechainnphi[0][2];
			config.sidechainnphi[i][0][3]=M_PI-my_monolayer.sidechainnphi[0][3];
		}
		else{
			config.sidechainnphi[i][2][1]=my_monolayer.sidechainnphi[2][1];
			config.sidechainnphi[i][3][1]=my_monolayer.sidechainnphi[3][1];
			config.sidechainnphi[i][0][2]=my_monolayer.sidechainnphi[0][2];
			config.sidechainnphi[i][0][3]=my_monolayer.sidechainnphi[0][3];
		}
    }
 
	//	Initialize
	
	associate_chains_and_nodes(Nmonomers, Nchains, chaintypes, Nchainsoftype, chainlengthoftype, chainchemistry, &Nnodes, &coordarray, &monomerid, &chainlength, my_bonded_params);

    //	Set up linked lists (new)

    double maxtranslate=0, smallest;
	linkset.phenylphenyl.mincellwidth=sqrt(my_nonbonded_params.phenylcutoff2)+maxtranslate;
	linkset.chargedcharged.mincellwidth=sqrt(my_nonbonded_params.cutoff2)+maxtranslate;
	linkset.shortrange.mincellwidth=2.*my_nonbonded_params.rhard0;
	if(2.*my_nonbonded_params.rhard1>linkset.shortrange.mincellwidth) linkset.shortrange.mincellwidth=2.*my_nonbonded_params.rhard1;
	if(2.*my_nonbonded_params.p2.rhard>linkset.shortrange.mincellwidth) linkset.shortrange.mincellwidth=2.*my_nonbonded_params.p2.rhard;
	if(2.*my_nonbonded_params.p3.rhard>linkset.shortrange.mincellwidth) linkset.shortrange.mincellwidth=2.*my_nonbonded_params.p3.rhard;
	linkset.shortrange.mincellwidth+=maxtranslate;
	linkset.phenylphenyl.maxneighbors=((int) 4*M_PI*linkset.phenylphenyl.mincellwidth/my_nonbonded_params.rhard1)+1;		//	kissing number with self
	smallest=my_nonbonded_params.p2.rhard;
	if(my_nonbonded_params.p3.rhard<smallest) smallest=my_nonbonded_params.p3.rhard;
	linkset.chargedcharged.maxneighbors=((int) 4*M_PI*linkset.chargedcharged.mincellwidth/smallest)+1;
	smallest=my_nonbonded_params.rhard0;
	if(my_nonbonded_params.rhard1<smallest) smallest=my_nonbonded_params.rhard1;
	if(my_nonbonded_params.p2.rhard<smallest) smallest=my_nonbonded_params.p2.rhard;
	if(my_nonbonded_params.p3.rhard<smallest) smallest=my_nonbonded_params.p3.rhard;
	linkset.shortrange.maxneighbors=((int) 4*M_PI*linkset.shortrange.mincellwidth/smallest)+1;

	FILE *outp;
	outp=fopen(outputfile, "w");
	copybilayer(config, &newconfig, Nnodetypes);
    int first=1;
 	for(leafspacing=firstz;leafspacing<=lastz;leafspacing+=zspace){
		for(leafxoffsetfrac=firstx;leafxoffsetfrac<=lastx;leafxoffsetfrac+=xspace){
			for(leafyoffsetfrac=firsty;leafyoffsetfrac<=lasty;leafyoffsetfrac+=yspace){

                for(i=0;i<2;i++){
                    newconfig.backboneoffset[i][0][2].z=(1-i*2)*my_monolayer.backboneoffset[0][2].z+(i-0.5)*leafspacing;
                    newconfig.backboneoffset[i][0][3].z=(1-i*2)*my_monolayer.backboneoffset[0][3].z+(i-0.5)*leafspacing;
                    newconfig.backboneoffset[i][2][1].z=(1-i*2)*my_monolayer.backboneoffset[2][1].z+(i-0.5)*leafspacing;
                    newconfig.backboneoffset[i][3][1].z=(1-i*2)*my_monolayer.backboneoffset[3][1].z+(i-0.5)*leafspacing;						//	mirror image
					if((leafsymmetry==-1)&&(i==1)){
						newconfig.backboneoffset[i][0][2].x=-my_monolayer.backboneoffset[0][2].x+i*leafxoffsetfrac*config.monomerspacing;
						newconfig.backboneoffset[i][0][3].x=-my_monolayer.backboneoffset[0][3].x+i*leafxoffsetfrac*config.monomerspacing;
						newconfig.backboneoffset[i][2][1].x=-my_monolayer.backboneoffset[2][1].x+i*leafxoffsetfrac*config.monomerspacing;
						newconfig.backboneoffset[i][3][1].x=-my_monolayer.backboneoffset[3][1].x+i*leafxoffsetfrac*config.monomerspacing;
					}
					else{
						newconfig.backboneoffset[i][0][2].x=my_monolayer.backboneoffset[0][2].x+i*leafxoffsetfrac*config.monomerspacing;
						newconfig.backboneoffset[i][0][3].x=my_monolayer.backboneoffset[0][3].x+i*leafxoffsetfrac*config.monomerspacing;
						newconfig.backboneoffset[i][2][1].x=my_monolayer.backboneoffset[2][1].x+i*leafxoffsetfrac*config.monomerspacing;
						newconfig.backboneoffset[i][3][1].x=my_monolayer.backboneoffset[3][1].x+i*leafxoffsetfrac*config.monomerspacing;
					}
                    newconfig.backboneoffset[i][0][2].y=my_monolayer.backboneoffset[0][2].y+i*leafyoffsetfrac*newconfig.interchainspacing;
                    newconfig.backboneoffset[i][0][3].y=my_monolayer.backboneoffset[0][3].y+i*leafyoffsetfrac*newconfig.interchainspacing;
                    newconfig.backboneoffset[i][2][1].y=my_monolayer.backboneoffset[2][1].y+i*leafyoffsetfrac*newconfig.interchainspacing;
                    newconfig.backboneoffset[i][3][1].y=my_monolayer.backboneoffset[3][1].y+i*leafyoffsetfrac*newconfig.interchainspacing;		//	0.5: alternating; 1: unlike charges face each other across bilayer
                }

                
                initialize_periodic_bilayer(Nmonomers, Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, newconfig, leafsymmetry);
				my_nonbonded_params.solvationparams.interfaceheight1=0.5*my_nonbonded_params.vacuumthickness;
				my_nonbonded_params.solvationparams.interfaceheight2=box_dimension.z-0.5*my_nonbonded_params.vacuumthickness;
				
                //	Update linked lists (new)
                
                if(first!=1){
                    free_3d_tensor((linkset.phenylphenyl.core.head), linkset.phenylphenyl.core.cellsperside.x, linkset.phenylphenyl.core.cellsperside.y);
                    free(linkset.phenylphenyl.core.list);
                    free(linkset.phenylphenyl.reverselist);
                    free(linkset.phenylphenyl.core.cell);
                    free(linkset.phenylphenyl.core.number_neighbors);
                    free_matrix(linkset.phenylphenyl.core.neighbor, Nnodes);
                    free_matrix(linkset.phenylphenyl.core.pair_energy, Nnodes);
                    free_3d_tensor((linkset.chargedcharged.core.head), linkset.chargedcharged.core.cellsperside.x, linkset.chargedcharged.core.cellsperside.y);
                    free(linkset.chargedcharged.core.list);
                    free(linkset.chargedcharged.reverselist);
                    free(linkset.chargedcharged.core.cell);
                    free(linkset.chargedcharged.core.number_neighbors);
                    free_matrix(linkset.chargedcharged.core.neighbor, Nnodes);
                    free_matrix(linkset.chargedcharged.core.pair_energy, Nnodes);
                    free_3d_tensor((linkset.shortrange.core.head), linkset.shortrange.core.cellsperside.x, linkset.shortrange.core.cellsperside.y);
                    free(linkset.shortrange.core.list);
                    free(linkset.shortrange.reverselist);
                    free(linkset.shortrange.core.cell);
                    free(linkset.shortrange.core.number_neighbors);
                    free_matrix(linkset.shortrange.core.neighbor, Nnodes);
                }
                first=0;
                
                configure_cells_struct(&(linkset.phenylphenyl), box_dimension);
                configure_cells_struct(&(linkset.chargedcharged), box_dimension);
                configure_cells_struct(&(linkset.shortrange), box_dimension);           //  recalc cellsperside, cellwidth
                allocate_linklist_pairenergies(&(linkset.phenylphenyl), Nnodes);
                allocate_linklist_pairenergies(&(linkset.chargedcharged), Nnodes);
                allocate_linklist(&(linkset.shortrange), Nnodes);
                cellconstructdoublylinked(coordarray, Nnodes, &(linkset.phenylphenyl), 1, 1);
                cellconstructdoublylinked(coordarray, Nnodes, &(linkset.chargedcharged), 2, 3);
                cellconstructdoublylinked(coordarray, Nnodes, &(linkset.shortrange), 0, 3);
                constructneighborlist_pairenergies_phenyl(&(linkset.phenylphenyl), Nnodes, coordarray, box_dimension, my_nonbonded_params);
                constructneighborlist_pairenergies_charged(&(linkset.chargedcharged), Nnodes, coordarray, box_dimension, my_nonbonded_params);
                constructneighborlist(&(linkset.shortrange), Nnodes, coordarray, box_dimension);

                newenergy=total_energy_cellstruct(coordarray, my_bonded_params, my_nonbonded_params, Nchains, chainlength, monomerid, box_dimension, linkset);
				fprintf(outp, "%f\t%f\t%f\t%f\n", leafxoffsetfrac, leafyoffsetfrac, leafspacing, newenergy);
			}
		}
	}
	fclose(outp);
				
	return 0;

}
