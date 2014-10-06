
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
		
	int charlength=200, i, j, chemistrymax=20, sidechaintypesmax=10, chaintypes, *Nchainsoftype, *chainlengthoftype, **chainchemistry, Nchains, Nmonomers, Nnodes, *chainlength, Nnodetypes, movieatoms, allatoms, cgatoms, cycles, t, movetype, moviecycle, calccycle, calcsperoutput, bestframe=0, bestallatomframe=0, chargetype, lastchargetype, selftype, fixedoffset, fixedinterchainspacing, last, self, xchains, terminuscode, abort, biasaminotrans;
	double kelvins, saltmolarity, temperature, oldenergy, newenergy, acceptanceprob, lowestenergy, randval, backbonezzigzag, backbonenz, backbonephi, minrelativeterminusspacing, maxrelativeterminusspacing, mininterchainspacing;
	char **chemistry_char;
	linkedlistset linkset;
	chemistry_char=xmalloc(sizeof(char *)*chemistrymax);
	for(i=0;i<chemistrymax;i++) chemistry_char[i]=xmalloc(sizeof(char)*charlength);
    double_triple box_dimension, bestavboxdimension;
	coord *coordarray;
	monomernodes **monomerid;
	declare_array(char, inputfile, charlength);	
	declare_array(char, coordfile, charlength);
	declare_array(char, base, charlength);
	bonded_params my_bonded_params;
	nonbonded_params my_nonbonded_params;
    cgparams mycgparams;
    allatomparams myallatomparams;
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
	
	//	initialize sidechains according to biased potential in trans configuration, then reset these parameters later
	
	my_bonded_params.sidechain[2].J10=424.643;
	my_bonded_params.sidechain[2].J11=-1844.99;

	my_bonded_params.sidechain[2].J12=2964.05;
	my_bonded_params.sidechain[2].J13=-2079.92;
	my_bonded_params.sidechain[2].J14=537.819;
	my_bonded_params.sidechain[2].k2=12.8785;
	my_bonded_params.sidechain[2].r20=-0.0587984;
	my_bonded_params.sidechain[2].r21=-1.024;
	my_bonded_params.sidechain[2].r22=0.476619;
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

    mycgparams.backrad=my_nonbonded_params.rhard0;
    mycgparams.phenrad=my_nonbonded_params.rhard1;
    mycgparams.aminrad=my_nonbonded_params.p2.rhard;
    mycgparams.carbrad=my_nonbonded_params.p3.rhard;
    
    myallatomparams.Cbetaspacing=1.43;
    myallatomparams.NterminusHdepth=-0.5*0.997;
    myallatomparams.NterminusHwidth=-0.5*sqrt(3.)*0.997;
    myallatomparams.Cnyldepth=-0.5*1.345;
    myallatomparams.Cnylwidth=-0.5*sqrt(3.)*1.345;
    myallatomparams.Onyldepth=-0.5*sqrt(3.)-1.23;                   //  assuming trans
    myallatomparams.Onylwidth=-0.5*sqrt(3.)*0.997;
    myallatomparams.CterminusHdepth=0.997*cos(117./180.*M_PI);
    myallatomparams.CterminusHwidth=0.997*sin(117./180.*M_PI);;
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

	strcpy(inputfile, loadparam(argc, "-inputfile", argv, "0").s);
	strcpy(base, loadparam(argc, "-base", argv, "../monolayer").s);

	//	System parameters
	
	kelvins=loadparam(argc, "-temperature", argv, "300").f;			//	in degrees Kelvin; default room T
    my_bonded_params.factor=loadparam(argc, "-bondedfactor", argv, "1.").f;
    my_nonbonded_params.p11vac.factor=loadparam(argc, "-phenylfactor", argv, "1.5").f;
    my_nonbonded_params.p2.factor=sqrt(loadparam(argc, "-chargedfactor", argv, "1.").f);
    my_nonbonded_params.p3.factor=sqrt(loadparam(argc, "-chargedfactor", argv, "1.").f);
	saltmolarity=loadparam(argc, "-saltmolarity", argv, "0.2").f;
	my_nonbonded_params.cutoff2=pow(loadparam(argc, "-cutoff", argv, "14.1151").f, 2);			//	numerical solution to U=-0.1 eV for salt 0.2 M (cutoff_v2.nb)
	loadseries(argc, "-chemistry", argv, chemistry_char);
	terminuscode=loadparam(argc, "-terminuscode", argv, "0").i;
	xchains=loadparam(argc, "-xchains", argv, "1").i;
	my_nonbonded_params.solvationparams.interface=loadparam(argc, "-interface", argv, "1").i;		
	my_nonbonded_params.vacuumthickness=loadparam(argc, "-vacuumthickness", argv, "20").f;	
	
	//syntax e.g.: -chemistry chaintypes patterntype1(0:alternating) Ntype1 length1(total monomers) nonpolar1 firstpolar1 secondpolar1 ...
	//or: -chemistry chaintypes patterntype1(1:block) Ntype1 polarblocks1 nonpolar1 polar1block1type polar1block1length(dimers) ...

	//	Schedule
	
	cycles=loadparam(argc, "-cycles", argv, "10000").i;
	moviecycle=loadparam(argc, "-moviecycle", argv, "100").i;
	calccycle=loadparam(argc, "-calccycle", argv, "10").i;
	calcsperoutput=loadparam(argc, "-calcsperoutput", argv, "10").i;
	
	//	Derived parameters
	
	temperature=onedegreeinkcals*kelvins;
	my_nonbonded_params.solvationparams.debyelength=sqrt(waterpermittivity*onedegreeinkcals*300/(4*M_PI*coulombconstant*2*saltmolarity*avogadrosnumberover10raised27));
	parse_chemistry(chemistry_char, &chaintypes, &Nchainsoftype, &chainlengthoftype, &chainchemistry, &Nmonomers, &Nchains, &monomercount, sidechaintypesmax, &Nnodetypes, terminuscode);
	declare_array(double, avheight, Nnodetypes);
	declare_array(double, bestavheight, Nnodetypes);
	my_energycomponents.solvation=malloc(sizeof(double)*Nnodetypes+1);
	my_energycomponents.interface=malloc(sizeof(double)*Nnodetypes+1);
	bestmy_energycomponents.solvation=malloc(sizeof(double)*Nnodetypes+1);
	bestmy_energycomponents.interface=malloc(sizeof(double)*Nnodetypes+1);

	//	Configuration parameters
	
	monolayer config, delta, newconfig;
	
    fixedoffset=loadparam(argc, "-fixedoffset", argv, "0").i;
    fixedinterchainspacing=loadparam(argc, "-fixedinterchainspacing", argv, "0").i;
	config.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
	config.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;		//	0 perp, 1 par, 2 alt (initially)
	minrelativeterminusspacing=loadparam(argc, "-minrelativeterminusspacing", argv, "-10").f;
	maxrelativeterminusspacing=loadparam(argc, "-maxrelativeterminusspacing", argv, "10").f;
	mininterchainspacing=loadparam(argc, "-mininterchainspacing", argv, "0").f;
	config.offsetfraction=0.5;
	if(*inputfile!='0'){
		input_monolayer(inputfile, &config, chainlengthoftype[0], Nnodetypes, 0);
		config.interchainspacing+=loadparam(argc, "-extrainterchainspacing", argv, "0").f;;
	}
	if(*inputfile=='0'){
		config.monomerspacing=loadparam(argc, "-monomerspacing", argv, "3.35341").f;
        if(terminuscode==1){
            config.terminusspacing=loadparam(argc, "-terminusspacing", argv, "3.35341").f;
            config.monomerscaleoffsetfraction=loadparam(argc, "-monomeroffsetfraction", argv, "0").f;
        }
        else{
            config.terminusspacing=loadparam(argc, "-terminusspacing", argv, "0").f;
            config.monomerscaleoffsetfraction=loadparam(argc, "-monomeroffsetfraction", argv, "1").f;
        }
		if(config.monomerspacing+config.terminusspacing<2*my_nonbonded_params.rhard0){
			printf("because of rhard0, increasing terminuspacing from %f to ", config.terminusspacing);
			config.terminusspacing+=2*config.monomerspacing;
			printf("%f\n", config.terminusspacing);
			config.monomerscaleoffsetfraction+=1;
		}
		config.interchainspacing=loadparam(argc, "-interchainspacing", argv, "4.5").f;
		backbonezzigzag=loadparam(argc, "-backbonezzigzag", argv, "1.54313").f;
		backbonenz=loadparam(argc, "-backbonenz", argv, "0.837134").f;
		backbonephi=loadparam(argc, "-backbonephi", argv, "1.94271").f;
	}
	  
	if(*inputfile=='0'){
		allocate_matrix_nozero(double_triple, config.backboneoffset, Nnodetypes, Nnodetypes);
		allocate_matrix(double, config.backbonenz, Nnodetypes, Nnodetypes);
		allocate_matrix(double, config.backbonephi, Nnodetypes, Nnodetypes);
		allocate_matrix_nozero(double_triple, config.sidechainrelativepos, Nnodetypes, Nnodetypes);
		allocate_matrix(double, config.sidechainnpar, Nnodetypes, Nnodetypes);
		allocate_matrix(double, config.sidechainnphi, Nnodetypes, Nnodetypes);

        for(i=0;i<2;i++){
            for(j=2;j<4;j++){
                if(i==1){
                    last=0;
                    self=j;
                }
                if(i==0){
                    last=j;
                    self=1;
                }
                config.backbonenz[last][self]=backbonenz;
                config.backbonephi[last][self]=backbonephi;
                config.sidechainrelativepos[last][self]=rsol(my_bonded_params.sidechain[self]);
                config.sidechainnpar[last][self]=nparhardsol(my_bonded_params.sidechainorientationtype[self], &(config.sidechainrelativepos[last][self]), my_bonded_params.sidechain[self]);
                config.sidechainnphi[last][self]=nphihardsol(my_bonded_params.sidechainorientationtype[self], config.sidechainnpar[last][self], config.sidechainrelativepos[last][self], my_bonded_params.sidechain[self]);
                
                config.backboneoffset[last][self].z=my_nonbonded_params.one1.zinterface-config.sidechainrelativepos[2][1].z;
                if(last==0) config.backboneoffset[last][self].z-=backbonezzigzag;
                config.backboneoffset[last][self].x=config.backboneoffset[last][self].y=0;
            }
        }
	}
	
	//	now set interactions to prescribed values, that may favor amino in cis or trans
	
	if(biasaminotrans==1){
        my_bonded_params.sidechain[2].J10=424.643;
        my_bonded_params.sidechain[2].J11=-1844.99;
    }
    else{
        my_bonded_params.sidechain[2].J10=418.586;
        my_bonded_params.sidechain[2].J11=-1837.99;
    }	

	config.length=chainlengthoftype[0]*config.monomerspacing+config.terminusspacing;
	config.offset=config.length*config.offsetfraction+config.monomerscaleoffsetfraction*config.monomerspacing;

	allocate_matrix_nozero(double_triple, newconfig.backboneoffset, Nnodetypes, Nnodetypes);
	allocate_matrix(double, newconfig.backbonenz, Nnodetypes, Nnodetypes);
	allocate_matrix(double, newconfig.backbonephi, Nnodetypes, Nnodetypes);
	allocate_matrix_nozero(double_triple, newconfig.sidechainrelativepos, Nnodetypes, Nnodetypes);
	allocate_matrix(double, newconfig.sidechainnpar, Nnodetypes, Nnodetypes);
	allocate_matrix(double, newconfig.sidechainnphi, Nnodetypes, Nnodetypes);
	allocate_matrix_nozero(double_triple, delta.backboneoffset, Nnodetypes, Nnodetypes);
	allocate_matrix(double, delta.backbonenz, Nnodetypes, Nnodetypes);
	allocate_matrix(double, delta.backbonephi, Nnodetypes, Nnodetypes);
	allocate_matrix_nozero(double_triple, delta.sidechainrelativepos, Nnodetypes, Nnodetypes);
	allocate_matrix(double, delta.sidechainnpar, Nnodetypes, Nnodetypes);
	allocate_matrix(double, delta.sidechainnphi, Nnodetypes, Nnodetypes);
	delta.monomerspacing=0.1;
	delta.monomerscaleoffsetfraction=0.05;
	delta.terminusspacing=0.1;
	delta.interchainspacing=0.1;
    for(i=0;i<Nnodetypes;i++){
        for(j=0;j<Nnodetypes;j++){
            delta.backboneoffset[i][j].x=delta.backboneoffset[i][j].y=delta.backboneoffset[i][j].z=0.1;
            delta.backbonenz[i][j]=0.05;
            delta.backbonephi[i][j]=0.1;
            delta.sidechainrelativepos[i][j].x=delta.sidechainrelativepos[i][j].y=delta.sidechainrelativepos[i][j].z=0.1;
            delta.sidechainnpar[i][j]=0.05;
            delta.sidechainnphi[i][j]=0.1;
        }
    }

	//	Initialize
	
	associate_chains_and_nodes(Nmonomers, Nchains, chaintypes, Nchainsoftype, chainlengthoftype, chainchemistry, &Nnodes, &coordarray, &monomerid, &chainlength, my_bonded_params);
	cgatoms=count_cg_atoms(Nchains, chainlength, monomerid, coordarray);
	allatoms=count_allatom_atoms(Nchains, chainlength, monomerid, coordarray);
		
	//	Set up bookkeeping
	
	initialize_periodic_monolayer_xchains(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, config, xchains);
 	my_nonbonded_params.solvationparams.interfaceheight1=0.5*my_nonbonded_params.vacuumthickness;
	my_nonbonded_params.solvationparams.interfaceheight2=box_dimension.z-0.5*my_nonbonded_params.vacuumthickness;

	declare_array(char, bestcgmoviefile, charlength);
	declare_array(char, bestcgsourcefile, charlength);
	declare_array(char, bestallatommoviefile, charlength);
	declare_array(char, bestallatomsourcefile, charlength);
	declare_array(char, besttotalenergyfile, charlength);
	declare_array(char, bestconfigfile, charlength);
	sprintf(bestcgmoviefile, "%s/c%i.best.cg.xyz", base, cycles);
	sprintf(bestcgsourcefile, "%s/c%i.best.cg.vmd.source", base, cycles);
	sprintf(bestallatommoviefile, "%s/c%i.best.allatom.xyz", base, cycles);
	sprintf(bestallatomsourcefile, "%s/c%i.best.allatom.vmd.source", base, cycles);
	sprintf(besttotalenergyfile, "%s/c%i.best.energy", base, cycles);
    initialize_energycomponents_nofiles(&my_energycomponents, Nnodetypes);
    initialize_energycomponents_nofiles(&bestmy_energycomponents, Nnodetypes);
    deletefile(besttotalenergyfile);
	sprintf(bestconfigfile, "%s/c%i.best.config", base, cycles);
    
	output_xyz_cg(Nchains, chainlength, monomerid, cgatoms, coordarray, box_dimension, bestcgmoviefile, bestcgsourcefile, mycgparams, &bestframe);
    output_xyz_allatom(Nchains, chainlength, monomerid, allatoms, coordarray, box_dimension, bestallatommoviefile, bestallatomsourcefile, myallatomparams, &bestallatomframe);
	
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
	configure_cells_struct(&(linkset.phenylphenyl), box_dimension);
	configure_cells_struct(&(linkset.chargedcharged), box_dimension);
	configure_cells_struct(&(linkset.shortrange), box_dimension);
	allocate_linklist_pairenergies(&(linkset.phenylphenyl), Nnodes);
	allocate_linklist_pairenergies(&(linkset.chargedcharged), Nnodes);
	allocate_linklist(&(linkset.shortrange), Nnodes);
	cellconstructdoublylinked(coordarray, Nnodes, &(linkset.phenylphenyl), 1, 1);
	cellconstructdoublylinked(coordarray, Nnodes, &(linkset.chargedcharged), 2, 3);
	cellconstructdoublylinked(coordarray, Nnodes, &(linkset.shortrange), 0, 3);
	constructneighborlist_pairenergies_phenyl(&(linkset.phenylphenyl), Nnodes, coordarray, box_dimension, my_nonbonded_params);
	constructneighborlist_pairenergies_charged(&(linkset.chargedcharged), Nnodes, coordarray, box_dimension, my_nonbonded_params);
	constructneighborlist(&(linkset.shortrange), Nnodes, coordarray, box_dimension);
	
	oldenergy=total_energy_cellstruct(coordarray, my_bonded_params, my_nonbonded_params, Nchains, chainlength, monomerid, box_dimension, linkset);
	
	calc_cellstruct_leafid(coordarray, my_bonded_params, my_nonbonded_params, Nchains, chainlength, monomerid, box_dimension, &bestmy_energycomponents, linkset, bestavheight, &bestavboxdimension);
	output_timeseries(besttotalenergyfile, &bestmy_energycomponents, 0, monomercount, 1., Nnodetypes, bestavheight, &bestavboxdimension, 0);

    lowestenergy=oldenergy;
    output_monolayer(bestconfigfile, config, oldenergy, 0);

	for(t=0;t<cycles;){
        if(fixedoffset==1){
			if(fixedinterchainspacing==1){
				movetype=rand_int(13);
				if(movetype>=2) movetype+=2;
				else if(movetype>=1) movetype++;
			}	
			else{
				movetype=rand_int(14);
				if(movetype>=1) movetype++;
			}
        }
		else if(fixedinterchainspacing==1){
            movetype=rand_int(14);
            if(movetype>=3) movetype++;
        }		
        else movetype=rand_int(15);
        chargetype=rand_int(2);
        lastchargetype=rand_int(2)+2;
        if(chargetype==0){
            selftype=lastchargetype;
            lastchargetype=0;
        }
        else{
            selftype=1;
        }
		copymonolayer(config, &newconfig, Nnodetypes);

		abort=0;
		if(movetype==0){
			newconfig.monomerspacing+=(2*rand_double-1)*delta.monomerspacing;
			if(newconfig.terminusspacing/newconfig.monomerspacing<minrelativeterminusspacing) abort=1;
			if(newconfig.terminusspacing/newconfig.monomerspacing>maxrelativeterminusspacing) abort=1;
		}
		if(movetype==1) newconfig.monomerscaleoffsetfraction+=(2*rand_double-1)*delta.monomerscaleoffsetfraction;
		if(movetype==2){
			newconfig.terminusspacing+=(2*rand_double-1)*delta.terminusspacing;
			if(newconfig.terminusspacing/newconfig.monomerspacing<minrelativeterminusspacing) abort=1;
			if(newconfig.terminusspacing/newconfig.monomerspacing>maxrelativeterminusspacing) abort=1;
		}
		if(movetype<3){
			newconfig.length=chainlengthoftype[0]*newconfig.monomerspacing+newconfig.terminusspacing;
			newconfig.offset=newconfig.length*newconfig.offsetfraction+newconfig.monomerscaleoffsetfraction*newconfig.monomerspacing;
		}
		if(movetype==3){
			newconfig.interchainspacing+=(2*rand_double-1)*delta.interchainspacing;
			if(newconfig.interchainspacing<mininterchainspacing) abort=1;
		}
		if(movetype==4) newconfig.backboneoffset[lastchargetype][selftype].x+=(2*rand_double-1)*delta.backboneoffset[lastchargetype][selftype].x;
		if(movetype==5) newconfig.backboneoffset[lastchargetype][selftype].y+=(2*rand_double-1)*delta.backboneoffset[lastchargetype][selftype].y;
		if(movetype==6) newconfig.backboneoffset[lastchargetype][selftype].z+=(2*rand_double-1)*delta.backboneoffset[lastchargetype][selftype].z;
        if(movetype==7){
            newconfig.backbonenz[lastchargetype][selftype]+=(2*rand_double-1)*delta.backbonenz[lastchargetype][selftype];
            reflect(&(newconfig.backbonenz[lastchargetype][selftype]));
        }
        if(movetype==8){
            newconfig.backbonephi[lastchargetype][selftype]+=(2*rand_double-1)*delta.backbonephi[lastchargetype][selftype];
            fmod((newconfig.backbonephi[lastchargetype][selftype]), (2*M_PI));
        }
		if(movetype==9) newconfig.sidechainrelativepos[lastchargetype][selftype].x+=(2*rand_double-1)*delta.sidechainrelativepos[lastchargetype][selftype].x;
		if(movetype==10) newconfig.sidechainrelativepos[lastchargetype][selftype].y+=(2*rand_double-1)*delta.sidechainrelativepos[lastchargetype][selftype].y;
		if(movetype==11) newconfig.sidechainrelativepos[lastchargetype][selftype].z+=(2*rand_double-1)*delta.sidechainrelativepos[lastchargetype][selftype].z;
        if(movetype==12){
            newconfig.sidechainnpar[lastchargetype][selftype]+=(2*rand_double-1)*delta.sidechainnpar[lastchargetype][selftype];
            reflect(&(newconfig.sidechainnpar[lastchargetype][selftype]));
        }
        if(movetype==13){
            newconfig.sidechainnphi[lastchargetype][selftype]+=(2*rand_double-1)*delta.sidechainnphi[lastchargetype][selftype];
            fmod((newconfig.sidechainnphi[lastchargetype][selftype]), (2*M_PI));
        }
		if(movetype==14){
			for(i=0;i<Nnodetypes;i++){
				for(j=0;j<Nnodetypes;j++){
					randval=rand_double;
					newconfig.backboneoffset[i][j].z+=(2*randval-1)*delta.backboneoffset[i][j].z;
				}
			}
		}
		initialize_periodic_monolayer(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, newconfig);
		my_nonbonded_params.solvationparams.interfaceheight1=0.5*my_nonbonded_params.vacuumthickness;
		my_nonbonded_params.solvationparams.interfaceheight2=box_dimension.z-0.5*my_nonbonded_params.vacuumthickness;
		
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

		if(abort==0){
			newenergy=total_energy_cellstruct(coordarray, my_bonded_params, my_nonbonded_params, Nchains, chainlength, monomerid, box_dimension, linkset);
			
			acceptanceprob=exp(-(newenergy-oldenergy)/temperature);
			if((acceptanceprob>1)||(rand_double<acceptanceprob)){
				copymonolayer(newconfig, &config, Nnodetypes);
				oldenergy=newenergy;
				if(oldenergy<lowestenergy){
					output_xyz_cg(Nchains, chainlength, monomerid, cgatoms, coordarray, box_dimension, bestcgmoviefile, bestcgsourcefile, mycgparams, &bestframe);
					output_xyz_allatom(Nchains, chainlength, monomerid, allatoms, coordarray, box_dimension, bestallatommoviefile, bestallatomsourcefile, myallatomparams, &bestallatomframe);    //  haven't written com version
					calc_cellstruct_leafid(coordarray, my_bonded_params, my_nonbonded_params, Nchains, chainlength, monomerid, box_dimension, &bestmy_energycomponents, linkset, bestavheight, &bestavboxdimension);
					output_timeseries(besttotalenergyfile, &bestmy_energycomponents, t+1, monomercount, 1., Nnodetypes, bestavheight, &bestavboxdimension, 0);
					output_monolayer(bestconfigfile, config, oldenergy, t+1);
					lowestenergy=oldenergy;
				}
			}
		}
		t++;
	}
				
	return 0;

}
