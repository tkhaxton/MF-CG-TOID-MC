#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"
#include "mymath.h"
#include "peptoid_functions.h"

int main(int argc, char *argv[]){
	long long int cycles, t;
	int charlength=400, i, j, maxseconds, chemistrymax=20, sidechaintypesmax=10, chaintypes, *Nchainsoftype, *chainlengthoftype, **chainchemistry, Nchains, Nmonomers, Nnodes, *chainlength, allatoms, cgatoms, moviecycle, frame, cgframe, allatomframe, calccycle, calcsperoutput, Nnodetypes, flag, ictype, qvals, leafsymmetry=1, mover, reset, resetframe, xchains, glycinebackbone, configNnodetypes, maxsurgery, pulling, pullchain, pullcalccycle, pullcalcsperoutput, max_neighbors, max_frustrated_links, *polymernumberneighbors, **polymerneighbor, *interior_member, *exterior_member, *in_cluster, *assigned, twoleavesflag=0, multiplesheetsflag=0, Nsheets=0, runningtime, terminuscode, inputflag, replicates, dovmmc, move3d=0, comcount=0, biasaminotrans;
	int_double *frustrated_link;
	double density, kelvins, saltmolarity, maxtranslate, maxrotate, maxrotatevmmc, translatefreq, rotatefreq, reptatefreq, surgeryfreq, temperature, myrand, maxlogaspect, qspace, aspectfreq, backbonezzigzag, backbonenz, backbonephi, smallest, molarity, wholepolymertranslatefreq, wholepolymerrotatefreq, chargedbackbonezoffset, leafspacing, leafxoffsetfrac, leafyoffsetfrac, residueshift=0, pullforce, wholebilayertranslatefreq, shiftbilayergapfreq, newheight, newwidth, newdepth, maxrforvmmc;
	char **chemistry_char;
	linkedlistset linkset;
	chemistry_char=xcalloc(chemistrymax, sizeof(char *));
	for(i=0;i<chemistrymax;i++) chemistry_char[i]=xcalloc(charlength, sizeof(char));
	double_triple box_dimension, avboxdimension, stack, avcom;
	coord *coordarray, *newcoordarray, *reversecoordarray;
	monomernodes **monomerid;
    pressureparam my_pressureparam;
	declare_array(char, inputfile, charlength);
	declare_array(char, coordfile, charlength);
	declare_array(char, base, charlength);
	declare_array(char, cgmoviefile, charlength);
	declare_array(char, cgmoviesourcefile, charlength);
	declare_array(char, allatommoviefile, charlength);
	declare_array(char, allatommoviesourcefile, charlength);
	declare_array(char, backbonesourcefile, charlength);
	declare_array(char, trajectoryfile, charlength);
	declare_array(char, comfile, charlength);
	declare_array(char, totalenergyfile, charlength);
	declare_array(char, pullfile, charlength);
	bonded_params my_bonded_params;
	nonbonded_params my_nonbonded_params;
    cgparams mycgparams;
    allatomparams myallatomparams;
	energycomponents my_energycomponents;
	monomertypes monomercount;
	monomercount.type=xcalloc(sidechaintypesmax, sizeof(double));
	polymer my_polymer;
 	monolayer my_monolayer;
    bilayer my_bilayer;
	time_t seconds, now;
    time(&now);
    srand((unsigned int) now);
	reptationparams leftreptation, rightreptation;
    declare_array_nozero(sidereptationparams, sidereptation, sidechaintypesmax);
    displacementparams displacement;
	   
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
 	my_bonded_params.sidechain=xcalloc(sidechaintypesmax, sizeof(sidechain_params));
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

 	my_bonded_params.sidechainorientationtype=xcalloc(sidechaintypesmax, sizeof(int));
	my_bonded_params.sidechainorientationtype[1]=1;		//	phenyl
	my_bonded_params.sidechainorientationtype[2]=0;		//	amino
	my_bonded_params.sidechainorientationtype[3]=0;		//	carboxyl
    
    leftreptation.rside0=3.;
    leftreptation.rforward0=0.;
    leftreptation.range=3.1;
    leftreptation.range2=leftreptation.range*leftreptation.range;
    leftreptation.maxconsecutivedirectordotproduct=0;
    rightreptation.rside0=3.;
    rightreptation.rforward0=0.;
    rightreptation.range=3.1;
    rightreptation.range2=rightreptation.range*rightreptation.range;
    rightreptation.maxconsecutivedirectordotproduct=0;
	
	glycinebackbone=loadparam(argc, "-glycinebackbone", argv, "0").i;   //  1 indicates glycine backbone to simulate phenylalanine monomer
	
    my_nonbonded_params.rhard0=loadparam(argc, "-backboneradius", argv, "2.05").f;
	my_nonbonded_params.one0.solvationenergy=10.1;
	my_nonbonded_params.one0.z0=-0.5;
	my_nonbonded_params.one0.uinterface=loadparam(argc, "-backbonesurfaceenergy", argv, "1.9").f;
	my_nonbonded_params.one0.zinterface=-0.5;
	my_nonbonded_params.one0.sigmainterface=0.6;
    if(glycinebackbone==1){
        my_nonbonded_params.one0.solvationenergy=150;
        my_nonbonded_params.one0.uinterface=0;
    }
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

    //  should be 9.85391*9.85391 for fine-tuned model (cutoff_v2.nb); U = -0.37 at wrong one
    
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
    
    //  should be 11.67 for fine-tuned model (toluene-toluene.nb); U=-0.42 at wrong one
    
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

	strcpy(base, loadparam(argc, "-base", argv, "../output").s);

	//	System parameters
	
 	kelvins=loadparam(argc, "-temperature", argv, "300").f;										//	in degrees Kelvin; default room T
    my_bonded_params.factor=loadparam(argc, "-bondedfactor", argv, "1.").f;
    my_nonbonded_params.p11vac.factor=loadparam(argc, "-phenylfactor", argv, "1.5").f;
    my_nonbonded_params.p2.factor=sqrt(loadparam(argc, "-chargedfactor", argv, "1.").f);
    my_nonbonded_params.p3.factor=sqrt(loadparam(argc, "-chargedfactor", argv, "1.").f);
    my_pressureparam.surfacepressure=loadparam(argc, "-surfacepressure", argv, "0").f;			//	kcal/mol/angstroms^2
	my_pressureparam.ypressure=loadparam(argc, "-ypressure", argv, "0").f;						//	kcal/mol/angstroms^2
	my_pressureparam.pressuretype=loadparam(argc, "-pressuretype", argv, "0").i;
	my_pressureparam.normalforceperunitarea=loadparam(argc, "-normalforceperunitarea", argv, "0").f;
    if((my_pressureparam.normalforceperunitarea!=0)&&(my_pressureparam.pressuretype!=7)&&(my_pressureparam.pressuretype!=0)) my_exit("normal force only set up for pressuretype 0 or 7");
	saltmolarity=loadparam(argc, "-saltmolarity", argv, "0.2").f;
	my_nonbonded_params.cutoff2=pow(loadparam(argc, "-cutoff", argv, "14.1151").f, 2);			//	numerical solution to U=-0.1 eV for salt 0.2 M (cutoff_v2.nb)
	loadseries(argc, "-chemistry", argv, chemistry_char);
	terminuscode=loadparam(argc, "-terminuscode", argv, "0").i;

	//syntax e.g.: -chemistry chaintypes patterntype1(0:alternating) Ntype1 length1(total monomers) nonpolar1 firstpolar1 secondpolar1 ...
	//or: -chemistry chaintypes patterntype1(1:block) Ntype1 polarblocks1 nonpolar1 polar1block1type polar1block1length(dimers) ...
    //or: -chemistry chaintypes patterntype1(2:not necessarily alternating phobic-philic) Ntype1 blocks1 type1 length1
	
 	my_nonbonded_params.solvationparams.interface=loadparam(argc, "-interface", argv, "0").i;
	my_nonbonded_params.vacuumthickness=loadparam(argc, "-vacuumthickness", argv, "40").f;	
	
	//	Box defined either by density or by monolayer parameters; density=0 indicates to use monolayer parameters
	
	strcpy(inputfile, loadparam(argc, "-inputfile", argv, "0").s);
    ictype=loadparam(argc, "-ictype", argv, "0").i;
    inputflag=loadparam(argc, "-inputflag", argv, "0").i;           //  1 for input from v0
	density=loadparam(argc, "-density", argv, "0.001").f;			//	monomers per angstroms^{-3}
	molarity=loadparam(argc, "-molarity", argv, "0.00002").f;
    
	//	Algorithm parameters
	
	maxrotate=loadparam(argc, "-maxrotate", argv, "0.8").f;
	maxtranslate=loadparam(argc, "-maxtranslate", argv, "0.4").f;
	maxlogaspect=loadparam(argc, "-maxlogaspect", argv, "0.002").f;
	translatefreq=loadparam(argc, "-translatefreq", argv, "0.5").f;
	wholepolymertranslatefreq=loadparam(argc, "-wholepolymertranslatefreq", argv, "0").f;     //  relative to one individual move per site on a polymer
	wholepolymerrotatefreq=loadparam(argc, "-wholepolymerrotatefreq", argv, "0").f;				//  relative to one individual move per site on a polymer
	maxrotatevmmc=loadparam(argc, "-maxrotatevmmc", argv, "0.8").f;
	maxrforvmmc=loadparam(argc, "-maxrforvmmc", argv, "8").f;								//	for self-consistent selfoverlap calc; (chainlength[0])^(1/3) times this
	wholebilayertranslatefreq=loadparam(argc, "-wholebilayertranslatefreq", argv, "0").f;     //  relative to one individual move per site on a bilayer
    shiftbilayergapfreq=loadparam(argc, "-shiftbilayergapfreq", argv, "0").f;
	aspectfreq=loadparam(argc, "-aspectfreq", argv, "0.15").f;                                //  relative to a whole cycle
    reptatefreq=loadparam(argc, "-reptatefreq", argv, "0").f;                                 //  relative to one individual move per site on a polymer
    maxsurgery=loadparam(argc, "-maxsurgery", argv, "1").i;                                   //  should set equal to xchains
    surgeryfreq=loadparam(argc, "-surgeryfreq", argv, "0").f;                               //  relative to one individual move per site on a row
    pulling=loadparam(argc, "-pulling", argv, "0").i;
    if(pulling==1){
        pullchain=0;
        pullforce=loadparam(argc, "-pullforce", argv, "0").f;
    }
	dovmmc=loadparam(argc, "-dovmmc", argv, "0").i;
	
	//	Schedule
	
	cycles=loadparam(argc, "-cycles", argv, "1000000").i;
	moviecycle=loadparam(argc, "-moviecycle", argv, "100000").i;
	calccycle=loadparam(argc, "-calccycle", argv, "100").i;
	calcsperoutput=loadparam(argc, "-calcsperoutput", argv, "10").i;
	qvals=loadparam(argc, "-qvals", argv, "25").i;
	qspace=loadparam(argc, "-qspace", argv, "0.02").f;
	maxseconds=(int) (3600*loadparam(argc, "-maxhours", argv, "0").f);
	reset=loadparam(argc, "-reset", argv, "0").i;
	resetframe=loadparam(argc, "-resetframe", argv, "0").i;
	pullcalccycle=loadparam(argc, "-calccycle", argv, "10").i;
	pullcalcsperoutput=loadparam(argc, "-calcsperoutput", argv, "100").i;

	//	Derived parameters
	
	temperature=onedegreeinkcals*kelvins;
	my_nonbonded_params.solvationparams.debyelength=sqrt(waterpermittivity*temperature/(4*M_PI*coulombconstant*2*saltmolarity*avogadrosnumberover10raised27));

	//	Using this aspect of linked list for initializing multiple bilayers

    linkset.shortrange.mincellwidth=2.*my_nonbonded_params.rhard0;
	if(2.*my_nonbonded_params.rhard1>linkset.shortrange.mincellwidth) linkset.shortrange.mincellwidth=2.*my_nonbonded_params.rhard1;
	if(2.*my_nonbonded_params.p2.rhard>linkset.shortrange.mincellwidth) linkset.shortrange.mincellwidth=2.*my_nonbonded_params.p2.rhard;
	if(2.*my_nonbonded_params.p3.rhard>linkset.shortrange.mincellwidth) linkset.shortrange.mincellwidth=2.*my_nonbonded_params.p3.rhard;

    if(ictype>0){
	
		//	Initialize
	
        if(!((my_nonbonded_params.solvationparams.interface==0)&&(ictype==4))&&!((ictype==7))){
            parse_chemistry(chemistry_char, &chaintypes, &Nchainsoftype, &chainlengthoftype, &chainchemistry, &Nmonomers, &Nchains, &monomercount, sidechaintypesmax, &Nnodetypes, terminuscode);
            associate_chains_and_nodes(Nmonomers, Nchains, chaintypes, Nchainsoftype, chainlengthoftype, chainchemistry, &Nnodes, &coordarray, &monomerid, &chainlength, my_bonded_params);
        }
		if(my_nonbonded_params.solvationparams.interface==0){
            if(ictype==1){                              //  dilute lattice (or single polymer) from backbone configuration
				aspectfreq=0;
				backbonezzigzag=loadparam(argc, "-backbonezzigzag", argv, "1.54313").f;
				backbonenz=loadparam(argc, "-backbonenz", argv, "0.837134").f;
				backbonephi=loadparam(argc, "-backbonephi", argv, "1.94271").f;
				initialize_polymer_config(&my_polymer, backbonezzigzag, backbonenz, backbonephi, Nnodetypes, my_bonded_params, 0);
				my_polymer.monomerspacing=loadparam(argc, "-monomerspacing", argv, "3.35341").f;
				initialize_2dlattice_fromconfig(Nmonomers, Nchains, density, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_polymer);
			}
			else if(ictype==3){                         //  bilayer from inputted configuration
                my_bilayer.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
                my_bilayer.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;
                my_bilayer.offsetfraction=0.5;
				residueshift=loadparam(argc, "-residueshift", argv, "0").f;
                input_bilayer(inputfile, &my_bilayer, chainlengthoftype[0], Nnodetypes, residueshift);
				xchains=loadparam(argc, "-xchains", argv, "1").i;
                initialize_periodic_bilayer_xchains(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_bilayer, leafsymmetry, xchains);
            }
			else if(ictype==4){							//	bulk solution, starting from copies of inputted isolated polymers ( CAN WORK FOR INTERFACE TOO )
				aspectfreq=0;
				Nchains=loadparam(argc, "-Nchains", argv, "1").i;
				input_and_copy_single_polymer(inputfile, &Nnodes, Nchains, &chainlength, &monomerid, &coordarray, &box_dimension, &(my_nonbonded_params.solvationparams.interface), &(my_nonbonded_params.vacuumthickness), &Nnodetypes, &monomercount, &t, &frame, molarity, my_nonbonded_params.cutoff2, &runningtime);
			}
			else if(ictype==5){							//	bulk solution, starting from free-floating bilayer fragment
				aspectfreq=0;
                my_bilayer.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;
                my_bilayer.offsetfraction=0.5;
                input_bilayer(inputfile, &my_bilayer, chainlengthoftype[0], Nnodetypes, residueshift);
				xchains=loadparam(argc, "-xchains", argv, "1").i;
                initialize_bilayer_fragment(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_bilayer, leafsymmetry, xchains, molarity);
			}
			else if(ictype==6){							//	bilayer starting from ordered backbone, default parameters from ordered_backbone.nb
				my_bilayer.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
                my_bilayer.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;
				my_bilayer.offsetfraction=0.5;
				my_bilayer.terminusspacing=0;
				my_bilayer.monomerspacing=loadparam(argc, "-monomerspacing", argv, "3.07739").f;
				my_bilayer.interchainspacing=loadparam(argc, "-interchainspacing", argv, "4.5").f;
				my_bilayer.monomerscaleoffsetfraction=loadparam(argc, "-monomeroffsetfraction", argv, "1").f;
				backbonezzigzag=loadparam(argc, "-backbonezzigzag", argv, "1.58157").f;
				backbonenz=loadparam(argc, "-backbonenz", argv, "0.87455").f;
				backbonephi=loadparam(argc, "-backbonephi", argv, "3.14159").f;
				chargedbackbonezoffset=loadparam(argc, "-chargedbackbonezoffset", argv, "-2").f;
				residueshift=loadparam(argc, "-residueshift", argv, "0").f;
				leafspacing=loadparam(argc, "-leafspacing", argv, "10").f;
				leafxoffsetfrac=loadparam(argc, "-leafxoffsetfrac", argv, "1").f;
				leafyoffsetfrac=loadparam(argc, "-leafyoffsetfrac", argv, "1").f;
				leafsymmetry=loadparam(argc, "-leafsymmetry", argv, "1").i;
				if(Nnodetypes<4) configNnodetypes=4;				//	create_ordered_bilayer assumes 4 node types
                else configNnodetypes=Nnodetypes;
                create_ordered_bilayer(&my_bilayer, chainlengthoftype[0], configNnodetypes, residueshift, chargedbackbonezoffset, backbonezzigzag, backbonenz, backbonephi, my_bonded_params, leafspacing, leafxoffsetfrac, leafyoffsetfrac, leafsymmetry);
				xchains=loadparam(argc, "-xchains", argv, "1").i;
                initialize_periodic_bilayer_xchains(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_bilayer, leafsymmetry, xchains);
			}
			else if(ictype==7){							// stacked bilayers
                replicates=loadparam(argc, "-replicates", argv, "2").i;
                stack.x=loadparam(argc, "-stackx", argv, "0").f;
                stack.y=loadparam(argc, "-stacky", argv, "0").f;
                stack.z=loadparam(argc, "-stackz", argv, "40").f;
                input_coords_and_replicate_stacks(inputfile, &Nnodes, &Nchains, &chainlength, &monomerid, &coordarray, &box_dimension, &(my_nonbonded_params.solvationparams.interface), &(my_nonbonded_params.vacuumthickness), &Nnodetypes, &monomercount, &t, &frame, &runningtime, replicates, stack, linkset.shortrange.mincellwidth);
            }
		}
		else{
            if(ictype==1) my_exit("not expecting interface 1, ictype 1");          
			else if(ictype==2){                                     //  monolayer from inputted configuration
                my_monolayer.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
                my_monolayer.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;
				if(chainlength[0]>2) my_monolayer.offsetfraction=0.5;
				else my_monolayer.offsetfraction=0;
				residueshift=loadparam(argc, "-residueshift", argv, "0").f;
                input_monolayer(inputfile, &my_monolayer, chainlengthoftype[0], Nnodetypes, residueshift);
 				xchains=loadparam(argc, "-xchains", argv, "1").i;
                initialize_periodic_monolayer_xchains(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_monolayer, xchains);
			}
			else if(ictype==3){                                     //  bilayer from inputted configuration
                my_bilayer.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
                my_bilayer.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;
                my_bilayer.offsetfraction=0.5;
				residueshift=loadparam(argc, "-residueshift", argv, "0").f;
                input_bilayer(inputfile, &my_bilayer, chainlengthoftype[0], Nnodetypes, residueshift);
				xchains=loadparam(argc, "-xchains", argv, "1").i;
                initialize_periodic_bilayer_xchains(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_bilayer, leafsymmetry, xchains);
            }
			else if(ictype==4){                                     //  monolayer from prescribed configuration
                my_monolayer.boxheight=loadparam(argc, "-boxheight", argv, "100").f;
                my_monolayer.phenylcode=loadparam(argc, "-phenylcode", argv, "1").f;
                if(chainlength[0]>1){
                    if(chainlength[0]>2) my_monolayer.offsetfraction=0.5;
                    else my_monolayer.offsetfraction=0;
					residueshift=loadparam(argc, "-residueshift", argv, "0").f;
                    input_monolayer(inputfile, &my_monolayer, chainlengthoftype[0], Nnodetypes, residueshift);
                }
                else{
					if(Nnodetypes<4) configNnodetypes=4;				//	even monomer uses indices [2][i] in my_monolayer
					else configNnodetypes=Nnodetypes;
					allocate_matrix_nozero(double_triple, my_monolayer.backboneoffset, configNnodetypes, configNnodetypes);
					allocate_matrix(double, my_monolayer.backbonenz, configNnodetypes, configNnodetypes);
					allocate_matrix(double, my_monolayer.backbonephi, configNnodetypes, configNnodetypes);
					allocate_matrix_nozero(double_triple, my_monolayer.sidechainrelativepos, configNnodetypes, configNnodetypes);
					allocate_matrix(double, my_monolayer.sidechainnpar, configNnodetypes, configNnodetypes);
					allocate_matrix(double, my_monolayer.sidechainnphi, configNnodetypes, configNnodetypes);
                    my_monolayer.interchainspacing=2*my_nonbonded_params.rhard0+0.1;
                    my_monolayer.length=my_nonbonded_params.p11vac.sigma0;
                    my_monolayer.terminusspacing=my_nonbonded_params.p11vac.sigma0-2*my_nonbonded_params.rhard0;
					my_monolayer.monomerspacing=my_monolayer.length-my_monolayer.terminusspacing;
                    my_monolayer.backboneoffset[2][1].z=my_nonbonded_params.solvationparams.interfacethickness;
                    my_monolayer.backboneoffset[2][1].x=my_monolayer.backboneoffset[2][1].y=0;
                    my_monolayer.backbonenz[2][1]=1;
                    my_monolayer.backbonephi[2][1]=1;
                    my_monolayer.sidechainrelativepos[2][1]=rsol(my_bonded_params.sidechain[1]);
                    my_monolayer.sidechainnpar[2][1]=nparhardsol(my_bonded_params.sidechainorientationtype[1], &(my_monolayer.sidechainrelativepos[2][1]), my_bonded_params.sidechain[1]);
                    my_monolayer.sidechainnphi[2][1]=nphihardsol(my_bonded_params.sidechainorientationtype[1], my_monolayer.sidechainnpar[2][1], my_monolayer.sidechainrelativepos[2][1], my_bonded_params.sidechain[1]);
                }
 				xchains=loadparam(argc, "-xchains", argv, "1").i;
                initialize_periodic_monolayer_xchains_bothinterfaces(Nchains, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_monolayer, xchains);
			}
			else if(ictype==5)  my_exit("not expecting interface 1, ictype 5");
            else if(ictype==6){                              //  dilute lattice (in xy instead of yz plane) (or single polymer) from backbone configuration
				backbonezzigzag=loadparam(argc, "-backbonezzigzag", argv, "1.5").f;
				backbonenz=loadparam(argc, "-backbonenz", argv, "1").f;
				backbonephi=loadparam(argc, "-backbonephi", argv, "0").f;
				initialize_polymer_config(&my_polymer, backbonezzigzag, backbonenz, backbonephi, Nnodetypes, my_bonded_params, my_nonbonded_params.one1.zinterface);
				my_polymer.monomerspacing=loadparam(argc, "-monomerspacing", argv, "3.5").f;
				initialize_2dxylattice_fromconfig(Nmonomers, Nchains, density, &box_dimension, chainlength, monomerid, coordarray, my_bonded_params, my_nonbonded_params, my_polymer);
			}
			else if(ictype==7){							//	exposing single polymer to interface
				aspectfreq=0;
				newheight=loadparam(argc, "-height", argv, "30").f;
				newwidth=loadparam(argc, "-width", argv, "30").f;
				newdepth=loadparam(argc, "-depth", argv, "15").f;
				expose_input_polymer_to_interface(inputfile, &Nnodes, &Nchains, &chainlength, &monomerid, &coordarray, &box_dimension, &(my_nonbonded_params.solvationparams.interface), &(my_nonbonded_params.vacuumthickness), &Nnodetypes, &monomercount, &t, &frame, &runningtime, newheight, newwidth, newdepth, my_nonbonded_params.vacuumthickness);
                move3d=1;
			}
		}
        t=0;frame=0;allatomframe=0;cgframe=0;
		newcoordarray=xcalloc(Nnodes, sizeof(coord));
        runningtime=0;
	}
	else{
		
		//	Input from file
		
		if(inputflag==1) input_coords_v0_bilayer(inputfile, &Nnodes, &Nchains, &chainlength, &monomerid, &coordarray, &box_dimension, &(my_nonbonded_params.solvationparams.interface), &(my_nonbonded_params.vacuumthickness), &Nnodetypes, &monomercount, &t, &frame, &runningtime);
		else input_coords(inputfile, &Nnodes, &Nchains, &chainlength, &monomerid, &coordarray, &box_dimension, &(my_nonbonded_params.solvationparams.interface), &(my_nonbonded_params.vacuumthickness), &Nnodetypes, &monomercount, &t, &frame, &runningtime);
		newcoordarray=xcalloc(Nnodes, sizeof(coord));
		if(reset==1){
			t=0;
			frame=0;
            cgframe=0;
            allatomframe=0;
            runningtime=0;
		}
		if(resetframe==1){
            frame=0;
            cgframe=0;
            allatomframe=0;
        }
        if(Nchains==1) aspectfreq=0;
	}
	cycles=t+cycles;
	allatomframe=frame;
	cgframe=frame;
    if(my_nonbonded_params.solvationparams.interface==0) move3d=1;      //  for vmmc
	
    //  Derived parameters
    
    for(i=0;i<Nnodes;i++){
        if(coordarray[i].leafid>0) twoleavesflag=1;
        if(coordarray[i].leafid>1) multiplesheetsflag=1;
        if(coordarray[i].leafid>Nsheets) Nsheets=coordarray[i].leafid;
    }
    Nsheets=(Nsheets)/2+1;
	my_nonbonded_params.solvationparams.interfaceheight1=0.5*my_nonbonded_params.vacuumthickness;
	my_nonbonded_params.solvationparams.interfaceheight2=box_dimension.z-0.5*my_nonbonded_params.vacuumthickness;
	cgatoms=count_cg_atoms(Nchains, chainlength, monomerid, coordarray);
	allatoms=count_allatom_atoms(Nchains, chainlength, monomerid, coordarray);
    aspectfreq=aspectfreq/(1.*Nnodes);
    wholepolymertranslatefreq=wholepolymertranslatefreq/(2.*chainlength[0]);
    wholepolymerrotatefreq=wholepolymerrotatefreq/(2.*chainlength[0]);
    wholebilayertranslatefreq=wholebilayertranslatefreq/(Nnodes/Nsheets);
    shiftbilayergapfreq=shiftbilayergapfreq/(Nnodes/Nsheets);
    reptatefreq=reptatefreq/(2.*chainlength[0]);
    surgeryfreq=surgeryfreq/(2.*chainlength[0]*maxsurgery);
    rotatefreq=1.-translatefreq-aspectfreq-reptatefreq-wholepolymertranslatefreq-wholepolymerrotatefreq-wholebilayertranslatefreq-surgeryfreq-shiftbilayergapfreq;
	maxrforvmmc*=pow(1.*chainlength[0], 1./3.);
	
	// Output files
	
    sprintf(cgmoviefile, "%s/vmd.cg.xyz", base);
	sprintf(cgmoviesourcefile, "%s/vmd.cg.source", base);
    sprintf(allatommoviefile, "%s/vmd.allatom.xyz", base);
	sprintf(allatommoviesourcefile, "%s/vmd.allatom.source", base);
    sprintf(comfile, "%s/com", base);
    sprintf(trajectoryfile, "%s/trajectory", base);
    sprintf(pullfile, "%s/pull", base);
	
	sprintf(totalenergyfile, "%s/timeseries", base);
	
	//	Output first frame early, so it is successful even if box is too small (coordinates for all-atom simulations)
	
	if(frame==0){
        output_trajectory(Nchains, chainlength, coordarray, monomerid, box_dimension, trajectoryfile, &frame, t);
        output_xyz_cg(Nchains, chainlength, monomerid, cgatoms, coordarray, box_dimension, cgmoviefile, cgmoviesourcefile, mycgparams, &cgframe);    //  haven't written com version
        output_xyz_allatom(Nchains, chainlength, monomerid, allatoms, coordarray, box_dimension, allatommoviefile, allatommoviesourcefile, myallatomparams, &allatomframe);    //  haven't written com version
    }	
    
	//	Set up linked lists
	
	linkset.phenylphenyl.mincellwidth=sqrt(my_nonbonded_params.phenylcutoff2)+maxtranslate;
	linkset.chargedcharged.mincellwidth=sqrt(my_nonbonded_params.cutoff2)+maxtranslate;
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

	// Set up vmmc
	
	if(dovmmc==1){
		max_neighbors=loadparam(argc, "-maxneighbors", argv, "50").i;
		if(max_neighbors>Nchains) max_neighbors=Nchains;
		max_frustrated_links=loadparam(argc, "-maxfrustratedlinks", argv, "20").i;
		if(max_frustrated_links<Nchains) max_frustrated_links=Nchains;
		allocate_array(int, polymernumberneighbors, Nchains);
		allocate_matrix(int, polymerneighbor, Nchains, max_neighbors);
		reversecoordarray=xcalloc(Nnodes, sizeof(coord));
		allocate_array(int, interior_member, Nchains);
		allocate_array(int, exterior_member, Nchains);
		allocate_array(int, in_cluster, Nchains);
		allocate_array(int_double, frustrated_link, Nchains);
		allocate_array(int, assigned, Nchains);
		update_neighbor_list(Nchains, chainlength, coordarray, my_nonbonded_params, linkset, monomerid, polymernumberneighbors, polymerneighbor, max_neighbors, assigned);
	}
	
	//	Calculation variables

	declare_array(double, avheight, Nnodetypes);
	my_energycomponents.solvation=xcalloc(Nnodetypes, sizeof(double));
	my_energycomponents.interface=xcalloc(Nnodetypes, sizeof(double));
	initialize_energycomponents_nofiles(&my_energycomponents, Nnodetypes);
    if(t==0){
        if(pulling==1) deletefile(pullfile);
        deletefile(totalenergyfile);
        deletefile(comfile);
    }
	avboxdimension.x=avboxdimension.y=avboxdimension.z=0;
    avcom.x=avcom.y=avcom.z=0;

	//	Initialize random number generator
	
	time(&seconds);

	//	Run simulations
	    
	for(;t<cycles;){
		for(i=0;i<Nnodes;i++){
			myrand=rand_double;
			mover=rand_int(Nnodes);
            while(coordarray[mover].nodetype<0) mover=rand_int(Nnodes);                    //  in case there is an empty node, e.g. for terminuscode 1
			if(myrand<translatefreq){
				if(coordarray[mover].nodetype==0){
                    if(pulling==1) mc_translate_single_hard_pull(mover, coordarray, maxtranslate, box_dimension, &(linkset.shortrange), my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength, pullchain, pullforce);
                    else mc_translate_single_hard(mover, coordarray, maxtranslate, box_dimension, &(linkset.shortrange), my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength);
                }
				else if(coordarray[mover].nodetype==1) mc_translate_single_phenyl(mover, coordarray, maxtranslate, box_dimension, &(linkset.shortrange), &(linkset.phenylphenyl), my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength);
				else mc_translate_single_charged(mover, coordarray, maxtranslate, box_dimension, &(linkset.shortrange), &(linkset.chargedcharged), my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength);
			}
            else{
                myrand-=translatefreq;
                if(myrand<rotatefreq){

                    //	shortrange (hard core) interactions are isotropic, so can ignore shortrange list

                    if(coordarray[mover].nodetype==0) mc_rotate_single_hard(mover, coordarray, maxrotate, box_dimension, linkset.shortrange.core, my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength);
                    else if(coordarray[mover].nodetype==1) mc_rotate_single_phenyl(mover, coordarray, maxrotate, box_dimension, linkset.phenylphenyl.core, my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength);
                    else if(coordarray[mover].nodetype>1) mc_rotate_single_charged(mover, coordarray, maxrotate, box_dimension, linkset.chargedcharged.core, my_bonded_params, my_nonbonded_params, temperature, monomerid, chainlength);
               }
                else{
                    myrand-=rotatefreq;
                    if(myrand<aspectfreq){
                        flag=mc_aspect_ratios_cellstruct(Nnodes, Nchains, coordarray, newcoordarray, &box_dimension, maxlogaspect, my_bonded_params, &my_nonbonded_params, chainlength, monomerid, linkset, my_pressureparam, temperature);
                        if(flag==1){
                            flag=change_cells_doublylinked(&(linkset.shortrange), box_dimension);
                            flag=change_cells_doublylinked(&(linkset.phenylphenyl), box_dimension);
                            flag=change_cells_doublylinked(&(linkset.chargedcharged), box_dimension);
                            cellconstructdoublylinked(coordarray, Nnodes, &(linkset.shortrange), 0, 3);
                            cellconstructdoublylinked(coordarray, Nnodes, &(linkset.phenylphenyl), 1, 1);
                            cellconstructdoublylinked(coordarray, Nnodes, &(linkset.chargedcharged), 2, 3);
                            constructneighborlist_pairenergies_phenyl(&(linkset.phenylphenyl), Nnodes, coordarray, box_dimension, my_nonbonded_params);
                            constructneighborlist_pairenergies_charged(&(linkset.chargedcharged), Nnodes, coordarray, box_dimension, my_nonbonded_params);
                            constructneighborlist(&(linkset.shortrange), Nnodes, coordarray, box_dimension);
                        }
                    }
                    else{
                        myrand-=aspectfreq;
                        if(myrand<wholepolymertranslatefreq){
							if(dovmmc==1) vmmc_translate_wholepolymer(Nchains, chainlength, newcoordarray, reversecoordarray, coordarray, maxtranslate, box_dimension, &linkset, my_bonded_params, my_nonbonded_params, temperature, monomerid, polymernumberneighbors, polymerneighbor, max_frustrated_links, linkset, max_neighbors, interior_member, exterior_member, in_cluster, frustrated_link, assigned, move3d);
                            else mc_translate_wholepolymer(Nchains, chainlength, newcoordarray, coordarray, maxtranslate, box_dimension, &linkset, my_bonded_params, my_nonbonded_params, temperature, monomerid);
                        }
                        else{
                            myrand-=wholepolymertranslatefreq;
							if(myrand<wholepolymerrotatefreq){
								if(dovmmc==1) vmmc_rotate_wholepolymer(Nchains, chainlength, newcoordarray, reversecoordarray, coordarray, maxrotatevmmc, box_dimension, &linkset, my_bonded_params, my_nonbonded_params, temperature, monomerid, polymernumberneighbors, polymerneighbor, max_frustrated_links, linkset, max_neighbors, interior_member, exterior_member, in_cluster, frustrated_link, assigned, move3d, maxrforvmmc);
							}
							else{
								myrand-=wholepolymerrotatefreq;
								if(myrand<wholebilayertranslatefreq){
									mc_translate_wholebilayer(Nnodes, Nsheets, chainlength, newcoordarray, coordarray, maxtranslate, box_dimension, &linkset, my_bonded_params, my_nonbonded_params, temperature, monomerid);
                                }
                                else{
                                    myrand-=wholebilayertranslatefreq;
                                    if(myrand<shiftbilayergapfreq){
                                        flag=mc_shiftbilayergap(Nnodes, Nsheets, chainlength, newcoordarray, coordarray, maxtranslate, &box_dimension, &linkset, my_bonded_params, my_nonbonded_params, temperature, monomerid, my_pressureparam.normalforceperunitarea, Nnodes/Nsheets);
                                        if(flag==1){
                                            flag=change_cells_doublylinked(&(linkset.shortrange), box_dimension);
                                            flag=change_cells_doublylinked(&(linkset.phenylphenyl), box_dimension);
                                            flag=change_cells_doublylinked(&(linkset.chargedcharged), box_dimension);
                                            cellconstructdoublylinked(coordarray, Nnodes, &(linkset.shortrange), 0, 3);
                                            cellconstructdoublylinked(coordarray, Nnodes, &(linkset.phenylphenyl), 1, 1);
                                            cellconstructdoublylinked(coordarray, Nnodes, &(linkset.chargedcharged), 2, 3);
                                            constructneighborlist_pairenergies_phenyl(&(linkset.phenylphenyl), Nnodes, coordarray, box_dimension, my_nonbonded_params);
                                            constructneighborlist_pairenergies_charged(&(linkset.chargedcharged), Nnodes, coordarray, box_dimension, my_nonbonded_params);
                                            constructneighborlist(&(linkset.shortrange), Nnodes, coordarray, box_dimension);
                                        }
                                    }
								}
                            }
                        }
                    }
                }
			}
        }
		t++;
		if(t%calccycle==0){
            calc_cellstruct_leafid(coordarray, my_bonded_params, my_nonbonded_params, Nchains, chainlength, monomerid, box_dimension, &my_energycomponents, linkset, avheight, &avboxdimension);
            if(Nchains==1) calc_com_trajectory(Nnodes, coordarray, box_dimension, &avcom, &comcount);

			if((t/calccycle)%calcsperoutput==0){
                time(&now);
				if(multiplesheetsflag==1) output_timeseries_sheets(totalenergyfile, &my_energycomponents, t, monomercount, 1./calcsperoutput, Nnodetypes, avheight, &avboxdimension, runningtime+now-seconds);
				else if(twoleavesflag==1) output_timeseries_leaves(totalenergyfile, &my_energycomponents, t, monomercount, 1./calcsperoutput, Nnodetypes, avheight, &avboxdimension, runningtime+now-seconds);
				else output_timeseries(totalenergyfile, &my_energycomponents, t, monomercount, 1./calcsperoutput, Nnodetypes, avheight, &avboxdimension, runningtime+now-seconds);
                if(t%moviecycle==0){
                    srand((unsigned int) now);
                    output_trajectory(Nchains, chainlength, coordarray, monomerid, box_dimension, trajectoryfile, &frame, t);
                    output_xyz_cg(Nchains, chainlength, monomerid, cgatoms, coordarray, box_dimension, cgmoviefile, cgmoviesourcefile, mycgparams, &cgframe);    //  haven't written com version
                    output_xyz_allatom(Nchains, chainlength, monomerid, allatoms, coordarray, box_dimension, allatommoviefile, allatommoviesourcefile, myallatomparams, &allatomframe);    //  haven't written com version
                    sprintf(coordfile, "%s/coord", base);
                    output_coords(coordfile, Nnodes, Nchains, chainlength, monomerid, coordarray, box_dimension, my_nonbonded_params.solvationparams.interface, my_nonbonded_params.vacuumthickness, Nnodetypes, monomercount, t, frame, runningtime+now-seconds);
                    if(Nchains==1) output_com(comfile, &avcom, &comcount, t);
                    if(maxseconds>0){
                        time(&now);
                        if((now-seconds)>maxseconds){
                            printf("aborting (time)\n");
                            return 0;
                        }
					}
                }
			}
		}
		if(pulling==1){
			if(t%pullcalccycle==0){
				calc_displacement(coordarray[monomerid[pullchain][0].backbone], coordarray[monomerid[pullchain][chainlength[pullchain]-1].backbone], &displacement, box_dimension.x);
				if((t/pullcalccycle)%pullcalcsperoutput==0) output_displacement_timeseries(pullfile, t, &displacement);
			}
		}

	}
	sprintf(coordfile, "%s/coord", base);
	output_coords(coordfile, Nnodes, Nchains, chainlength, monomerid, coordarray, box_dimension, my_nonbonded_params.solvationparams.interface, my_nonbonded_params.vacuumthickness, Nnodetypes, monomercount, t, frame, runningtime+now-seconds);
	return 0;
}
