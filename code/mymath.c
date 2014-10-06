#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mymath.h"

double_triple add_double_triple(double_triple a, double_triple b){
	double_triple result;
	result.x=(a).x+(b).x;
	result.y=(a).y+(b).y;
	result.z=(a).z+(b).z;
	return result;
}

double_triple scalar_multiply_double_triple(double_triple a, double s){
	double_triple result;
	result.x=s*a.x;
	result.y=s*a.y;
	result.z=s*a.z;
	return result;
}

double dot_product(double_triple a, double_triple b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

double norm(double_triple a){
	return sqrt(dot_product(a, a));
}

double_triple rand_unit_cube(){
	double_triple result;
	result.x=2*rand_double-1;
	result.y=2*rand_double-1;
	result.z=2*rand_double-1;
	return result;
}

double_triple rand_unit_ball(){
	double_triple result;
	double r2=1;
	while(r2>=1){
		result.x=2*rand_double-1;
		result.y=2*rand_double-1;
		result.z=2*rand_double-1;
		r2=dot_product(result, result);
	}
	return result;
}

double_triple rand_unit_sphere(){
	double_triple result;
	double r2=1;
	while(r2>=1){
		result.x=2*rand_double-1;
		result.y=2*rand_double-1;
		result.z=2*rand_double-1;
		r2=dot_product(result, result);
	}
	r2=1/(sqrt(r2));
	result.x*=r2;
	result.y*=r2;
	result.z*=r2;
	return result;
}

double_triple subtract_double_triple(double_triple a, double_triple b){
	double_triple result;
	result.x=a.x-b.x;
	result.y=a.y-b.y;
	result.z=a.z-b.z;
	return result;
}

void recenter_double_triple(double_triple *pa, double_triple b){
	recenter(((*pa).x), ((b).x));
	recenter(((*pa).y), ((b).y));
	recenter(((*pa).z), ((b).z));
}

void fmod_double_triple(double_triple *pa, double_triple b){
	fmod(((*pa).x), ((b).x));
	fmod(((*pa).y), ((b).y));
	fmod(((*pa).z), ((b).z));
}

void normalize(double_triple *pvec){
	double mynorm=norm(*pvec);
	(*pvec).x/=mynorm;
	(*pvec).y/=mynorm;
	(*pvec).z/=mynorm;
}

double_triple normed(double_triple vec){
	double_triple result=(vec);
    double oneovernorm=1./norm(vec);
    result.x*=oneovernorm;
    result.y*=oneovernorm;
    result.z*=oneovernorm;
    return result;
}

double_triple cross_product(double_triple a, double_triple b){
	double_triple result;
	result.x=a.y*b.z-a.z*b.y;
	result.y=a.z*b.x-a.x*b.z;
	result.z=a.x*b.y-a.y*b.x;
	return result;
}

double largest_cubic_root(double a3, double a2, double a1, double a0){
    
    // from http://www-old.me.gatech.edu/energy/andy_phd/appA.htm
    // except he was missing a 1/3 in definition of phi
    // corrected by looking at http://home.pipeline.com/~hbaker1/cubic3realroots.htm
    
    double p=a2/a3;
    double q=a1/a3;
    double r=a0/a3;
    double A=(3*q-p*p)/3;
    double B=(2*p*p*p-9*p*q+27*r)/27;
    double D=A*A*A/27+B*B/4;
    double M, N, y1, y2, y3, phi, x1, x2, x3, biggest, M3, N3, arg;
    if(D>0){                        //  one real root
        M3=-B/2+sqrt(D);
        if(M3<0){
            M=-pow(-M3, 1./3.);
        }
        else{
            M=pow(M3, 1./3.);
        }
        N3=-B/2-sqrt(D);
        if(N3<0){
            N=-pow(-N3, 1./3.);
        }
        else{
            N=pow(N3, 1./3.);
        }
        y1=M+N;
        x1=y1-p/3;
        biggest=x1;
    }
    else{                           //  three real roots
        if(B>0){
            arg=B*B/4/(-A*A*A/27);
            phi=acos(-sqrt(arg))/3;
        }
        else{
            arg=B*B/4/(-A*A*A/27);
            phi=acos(sqrt(arg))/3;
        }
        y1=2*sqrt(-A/3)*cos(phi);
        y2=2*sqrt(-A/3)*cos(phi+2./3.*M_PI);
        y3=2*sqrt(-A/3)*cos(phi+4./3.*M_PI);
        x1=y1-p/3;
        x2=y2-p/3;
        x3=y3-p/3;
        if(x1>x2){
            if(x1>x3) biggest=x1;
            else biggest=x3;
        }
        else{
            if(x2>x3) biggest=x2;
            else biggest=x3;
        }
    }
    return biggest;
}

double quartic_global_minimum(double p4, double p3, double p2, double p1){
    double a3=4*p4;
    double a2=3*p3;
    double a1=2*p2;
    double a0=p1;

    // solve cubic
    
    // from http://www-old.me.gatech.edu/energy/andy_phd/appA.htm
    // except he was missing a 1/3 in definition of phi
    // corrected by looking at http://home.pipeline.com/~hbaker1/cubic3realroots.htm
    
    double p=a2/a3;
    double q=a1/a3;
    double r=a0/a3;
    double A=(3*q-p*p)/3;
    double B=(2*p*p*p-9*p*q+27*r)/27;
    double D=A*A*A/27+B*B/4;
    double M, N, y1, y2, y3, phi, x1, x2, x3, lowest, M3, N3, arg;
    if(D>0){                        //  one real root
        M3=-B/2+sqrt(D);
        if(M3<0){
            M=-pow(-M3, 1./3.);
        }
        else{
            M=pow(M3, 1./3.);
        }
        N3=-B/2-sqrt(D);
        if(N3<0){
            N=-pow(-N3, 1./3.);
        }
        else{
            N=pow(N3, 1./3.);
        }
        y1=M+N;
        x1=y1-p/3;
        lowest=x1;
    }
    else{                           //  three real roots
        if(B>0){
            arg=B*B/4/(-A*A*A/27);
            phi=acos(-sqrt(arg))/3;
        }
        else{
            arg=B*B/4/(-A*A*A/27);
            phi=acos(sqrt(arg))/3;
        }
        y1=2*sqrt(-A/3)*cos(phi);
        y2=2*sqrt(-A/3)*cos(phi+2./3.*M_PI);
        y3=2*sqrt(-A/3)*cos(phi+4./3.*M_PI);
        x1=y1-p/3;
        x2=y2-p/3;
        x3=y3-p/3;
        if(p4*pow(x1, 4)+p3*pow(x1, 3)+p2*pow(x1, 2)+p1*x1<p4*pow(x2, 4)+p3*pow(x2, 3)+p2*pow(x2, 2)+p1*x2){
            if(p4*pow(x1, 4)+p3*pow(x1, 3)+p2*pow(x1, 2)+p1*x1<p4*pow(x3, 4)+p3*pow(x3, 3)+p2*pow(x3, 2)+p1*x3){
                lowest=x1;
            }
            else lowest=x3;
        }
        else{
            if(p4*pow(x2, 4)+p3*pow(x2, 3)+p2*pow(x2, 2)+p1*x2<p4*pow(x3, 4)+p3*pow(x3, 3)+p2*pow(x3, 2)+p1*x3){
                lowest=x2;
            }
            else lowest=x3;
        }
    }
    return lowest;
}

void matrix_multiply(double **a, double **b, double **c){
    c[0][0]=a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0];
    c[0][1]=a[0][0]*b[0][1]+a[0][1]*b[1][1]+a[0][2]*b[2][1];
    c[0][2]=a[0][0]*b[0][2]+a[0][1]*b[1][2]+a[0][2]*b[2][2];
    c[1][0]=a[1][0]*b[0][0]+a[1][1]*b[1][0]+a[1][2]*b[2][0];
    c[1][1]=a[1][0]*b[0][1]+a[1][1]*b[1][1]+a[1][2]*b[2][1];
    c[1][2]=a[1][0]*b[0][2]+a[1][1]*b[1][2]+a[1][2]*b[2][2];
    c[2][0]=a[2][0]*b[0][0]+a[2][1]*b[1][0]+a[2][2]*b[2][0];
    c[2][1]=a[2][0]*b[0][1]+a[2][1]*b[1][1]+a[2][2]*b[2][1];
    c[2][2]=a[2][0]*b[0][2]+a[2][1]*b[1][2]+a[2][2]*b[2][2];
}

void forward_matrix(double angle, double_triple axis_vector, double **forward){
	double costheta=cos(angle);
	double oneminuscos=1.0-cos(angle);
	double sintheta=sin(angle);
	forward[0][0]=costheta+axis_vector.x*axis_vector.x*oneminuscos;
	forward[0][1]=axis_vector.x*axis_vector.y*oneminuscos-axis_vector.z*sintheta;
	forward[0][2]=axis_vector.x*axis_vector.z*oneminuscos+axis_vector.y*sintheta;
	forward[1][0]=axis_vector.x*axis_vector.y*oneminuscos+axis_vector.z*sintheta;
	forward[1][1]=costheta+axis_vector.y*axis_vector.y*oneminuscos;
	forward[1][2]=axis_vector.y*axis_vector.z*oneminuscos-axis_vector.x*sintheta;
	forward[2][0]=axis_vector.x*axis_vector.z*oneminuscos-axis_vector.y*sintheta;
	forward[2][1]=axis_vector.y*axis_vector.z*oneminuscos+axis_vector.x*sintheta;
	forward[2][2]=costheta+axis_vector.z*axis_vector.z*oneminuscos;
}

void forward_and_backward_matrix(double angle, double_triple axis_vector, double **forward, double **backward){
	double costheta=cos(angle);
	double oneminuscos=1.0-cos(angle);
	double sintheta=sin(angle);
	forward[0][0]=costheta+axis_vector.x*axis_vector.x*oneminuscos;
	forward[0][1]=axis_vector.x*axis_vector.y*oneminuscos-axis_vector.z*sintheta;
	forward[0][2]=axis_vector.x*axis_vector.z*oneminuscos+axis_vector.y*sintheta;
	forward[1][0]=axis_vector.x*axis_vector.y*oneminuscos+axis_vector.z*sintheta;
	forward[1][1]=costheta+axis_vector.y*axis_vector.y*oneminuscos;
	forward[1][2]=axis_vector.y*axis_vector.z*oneminuscos-axis_vector.x*sintheta;
	forward[2][0]=axis_vector.x*axis_vector.z*oneminuscos-axis_vector.y*sintheta;
	forward[2][1]=axis_vector.y*axis_vector.z*oneminuscos+axis_vector.x*sintheta;
	forward[2][2]=costheta+axis_vector.z*axis_vector.z*oneminuscos;
	backward[0][0]=costheta+axis_vector.x*axis_vector.x*oneminuscos;
	backward[0][1]=axis_vector.x*axis_vector.y*oneminuscos+axis_vector.z*sintheta;
	backward[0][2]=axis_vector.x*axis_vector.z*oneminuscos-axis_vector.y*sintheta;
	backward[1][0]=axis_vector.x*axis_vector.y*oneminuscos-axis_vector.z*sintheta;
	backward[1][1]=costheta+axis_vector.y*axis_vector.y*oneminuscos;
	backward[1][2]=axis_vector.y*axis_vector.z*oneminuscos+axis_vector.x*sintheta;
	backward[2][0]=axis_vector.x*axis_vector.z*oneminuscos+axis_vector.y*sintheta;
	backward[2][1]=axis_vector.y*axis_vector.z*oneminuscos-axis_vector.x*sintheta;
	backward[2][2]=costheta+axis_vector.z*axis_vector.z*oneminuscos;
}

double_triple rotate_by_matrix(double_triple start, double **matrix){
	double_triple result;
	result.x=start.x*matrix[0][0]+start.y*matrix[0][1]+start.z*matrix[0][2];
	result.y=start.x*matrix[1][0]+start.y*matrix[1][1]+start.z*matrix[1][2];
	result.z=start.x*matrix[2][0]+start.y*matrix[2][1]+start.z*matrix[2][2];
	return result;
}



