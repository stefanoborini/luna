#include <math.h>
#include <stdio.h>
#include "luna.h"

main(int argc,char *argv[])
{
char filename[80],line[80],luogo[40];
FILE *fpun;

int anno,mese,giorno,FINESTRA;

if (argc<3) {
	printf("mancano parametri indispensabili\n");
	exit(0);
	}


if (!(fpun=fopen(argv[1],"r"))) {
	printf("\nFile %s non trovato\n\n",filename);
	exit(1);
	}

FINESTRA = atoi(argv[2]);

fscanf(fpun,"%d %d %d %s",&anno,&mese,&giorno,luogo);

while (!feof(fpun)) {
	if (anno<=0 || mese<=0 || giorno<=0) {
		fscanf(fpun,"%d %d %d %s",&anno,&mese,&giorno,luogo);
		continue;
		}

	converti_giuliano(giorno,mese,anno);
	printf("\nCalcolo per il %d/%d/%d (%.10g), luogo %s....\n",giorno,mese,anno,jd,luogo);
	jd_attuale=jd-FINESTRA;
	while (jd_attuale<=jd+FINESTRA) {
		printf("  JD = %.10g\n",jd_attuale);
		T = (jd_attuale-2415020.0)/36525;
		calcola_parametri_sole();
		calcola_parametri_luna();
		calcola_fase_lunare();
		calcola_separazione_angolare();
		jd_attuale+=0.1;
		}
	fscanf(fpun,"%d %d %d %s",&anno,&mese,&giorno,luogo);
	}
fclose(fpun);
}


converti_giuliano(int giorno,int mese,int anno)
{
int A,B;

if (!(mese>2)) {
	anno--;
	mese+=12;
	}

A=anno/100;
B=2-A+A/4;

jd= (int)(365.25*(anno+4716))+(int)(30.6*(mese+1))+giorno+B-1524.5;
}

calcola_parametri_sole()
{
double L,M,A,B,C,D,E,F,H,e;
double T2=pow(T,2);
double T3=pow(T,3);


M= 358.47583 + 35999.04975*T - 0.00015*T2 - 0.0000033*T3;
L= 279.69668 + 36000.76892*T + 0.0003025*T2;
e= 0.01675104-0.0000418*T-0.000000126*T2;
F=(1.919460 - 0.004789*T - 0.000014*T2) * dsin(M)
  + (0.020094 - 0.0001*T) * dsin(2*M)
  + 0.000293*dsin(3*M);

LV=L+F;
v=M+F;

/*calcolo delle correzioni*/

A=153.23+22518.7541*T;
B=216.57+45037.5082*T;
C=312.69+32964.3577*T;
D=350.74+445267.1142*T-0.00144*T2;
E=231.19+20.20*T;
H=353.40+65928.7155*T;

LV+=0.00134*dcos(A);
LV+=0.00154*dcos(B);
LV+=0.002*dcos(C);
LV+=0.00179*dsin(D);
LV+=0.00178*dsin(E);

/**************************************************/

/********** calcolo del raggio vettore  ***********/

R=(1.0000002*(1-pow(e,2)))/(1+e*dcos(v));

/****** correzioni al raggio vettore   *****/

R+=0.00000543*dsin(A);
R+=0.00001575*dsin(B);
R+=0.00001627*dsin(C);
R+=0.00003076*dcos(D);
R+=0.00000927*dsin(H);

}


calcola_parametri_luna()
{
double L1,M,M1,D,F,omega,e,e2,B,omega1,omega2;

double T2=pow(T,2);
double T3=pow(T,3);

L1 = 270.434164 + 481267.8831*T - 0.001133*T2 + 0.0000019*T3;
M= 358.475833 + 35999.0498*T - 0.00015*T2 - 0.0000033*T3;
M1 = 296.104608 + 477198.8491*T + 0.009192*T2 + 0.0000144*T3;
D  = 350.737486 + 445267.1142*T - 0.001436*T2 + 0.0000019*T3;
F  = 11.250889  + 483202.0251*T - 0.003211*T2 - 0.0000003*T3;
omega = 259.183275 - 1934.1420*T + 0.002078*T2 + 0.0000022*T3;

/* additivi correttivi */

L1+=0.000233*dsin(51.2+20.2*T);
M-=0.001778*dsin(51.2+20.2*T);
M1+=0.000817*dsin(51.2+20.2*T);
D+=0.002011*dsin(51.2+20.2*T);
L1+=0.003964*dsin(346.560+132.870*T-0.0091731*T2);
M1+=0.003964*dsin(346.560+132.870*T-0.0091731*T2);
D+=0.003964*dsin(346.560+132.870*T-0.0091731*T2);
F+=0.003964*dsin(346.560+132.870*T-0.0091731*T2);
L1+=0.001964*dsin(omega);
M1+=0.002541*dsin(omega);
D+=0.001964*dsin(omega);
F-=0.024691*dsin(omega);
F-=0.004328*dsin(omega+275.05-2.3*T);

/* fine additivi correttivi */

e = 1 - 0.002495*T - 0.00000752*T2;
e2 = pow(e,2);

lambda = L1 + 6.288750*dsin(M1)
         + 1.274018*dsin(2*D-M1)
         + 0.658309*dsin(2*D)
         + 0.213616*dsin(2*M1)
         - e*0.185596*dsin(M)
         - 0.114336*dsin(2*F)
         + 0.058793*dsin(2*D-2*M1)
         + e*0.057212*dsin(2*D-M-M1)
         + 0.053320*dsin(2*D+M1)
         + e*0.045874*dsin(2*D-M)
         + e*0.041024*dsin(M1-M)
         - 0.034718*dsin(D)
         - e*0.030465*dsin(M+M1)
         + 0.015326*dsin(2*D-2*F)
         - 0.012528*dsin(2*F+M1)
         - 0.010980*dsin(2*F-M1)
         + 0.010674*dsin(4*D-M1)
         + 0.010034*dsin(3*M1)
         + 0.008548*dsin(4*D-2*M1)
         - e*0.007910*dsin(M-M1+2*D)
         - e*0.006783*dsin(2*D+M)
         + 0.005162*dsin(M1-D)
         + e*0.005*dsin(M+D)
         + e*0.004049*dsin(M1-M+2*D)
         + 0.003996*dsin(2*M1+2*D)
         + 0.003862*dsin(4*D)
         + 0.003665*dsin(2*D-3*M1)
         + e*0.002695*dsin(2*M1-M)
         + 0.002602*dsin(M1-2*F-2*D)
         + e*0.002396*dsin(2*D-M-2*M1)
         - 0.002349*dsin(M1+D)
         + e2*0.002249*dsin(2*D-2*M)
         - e*0.002125*dsin(2*M1+M)
         - e2*0.002079*dsin(2*M)
         + e2*0.002059*dsin(2*D-M1-2*M)
         - 0.001773*dsin(M1+2*D-2*F)
         - 0.001595*dsin(2*F+2*D)
         + e*0.001220*dsin(4*D-M-M1)
         - 0.001110*dsin(2*M1+2*F)
         + 0.000892*dsin(M1-3*D)
         - e*0.000811*dsin(M+M1+2*D)
         + e*0.000761*dsin(4*D-M-2*M1)
         + e2*0.000717*dsin(M1-2*M)
         + e2*0.000704*dsin(M1-2*M-2*D)
         + e*0.000693*dsin(M-2*M1+2*D)
         + e*0.000598*dsin(2*D-M-2*F)
         + 0.000550*dsin(M1+4*D)
         + 0.000538*dsin(4*M1)
         + e*0.000521*dsin(4*D-M)
         + 0.000486*dsin(2*M1-D);

B = + 5.128189*dsin(F)
    + 0.280606*dsin(M1+F)
    + 0.277693*dsin(M1-F)
    + 0.173238*dsin(2*D-F)
    + 0.055413*dsin(2*D+F-M1)
    + 0.046272*dsin(2*D-F-M1)
    + 0.032573*dsin(2*D+F)
    + 0.017198*dsin(2*M1+F)
    + 0.009267*dsin(2*D+M1-F)
    + 0.008823*dsin(2*M1-F)
    + e*0.008247*dsin(2*D-M-F)
    + 0.004323*dsin(2*D-F-2*M1)
    + 0.0042*dsin(2*D+F+M1)
    + e*0.003372*dsin(F-M-2*D)
    + e*0.002472*dsin(2*D+F-M-M1)
    + e*0.002222*dsin(2*D+F-M)
    + e*0.002072*dsin(2*D-F-M-M1)
    + e*0.001877*dsin(F-M+M1)
    + 0.001828*dsin(4*D-F-M1)
    - e*0.001803*dsin(F+M)
    - 0.00175*dsin(3*F)
    + e*0.001570*dsin(M1-M-F)
    - 0.001487*dsin(F+D)
    - e*0.001481*dsin(F+M+M1)
    + e*0.001417*dsin(F-M-M1)
    + e*0.001350*dsin(F-M)
    + 0.001330*dsin(F-D)
    + 0.001106*dsin(F+3*M1)
    + 0.001020*dsin(4*D-F)
    + 0.000833*dsin(F+4*D-M1)
    + 0.000781*dsin(M1-3*F)
    + 0.000670*dsin(F+4*D-2*M1)
    + 0.000606*dsin(2*D-3*F)
    + 0.000597*dsin(2*D+2*M1-F)
    + e*0.000492*dsin(2*D+M1-M-F)
    + 0.000450*dsin(2*M1-F-2*D)
    + 0.000439*dsin(3*M1-F)
    + 0.000423*dsin(F+2*D+2*M1)
    + 0.000422*dsin(2*D-F-3*M1)
    - e*0.000367*dsin(M+F+2*D-M1)
    - e*0.000353*dsin(M+F+2*D)
    + 0.000331*dsin(F+4*D)
    + e*0.000317*dsin(2*D+F-M+M1)
    + e2*0.000306*dsin(2*D-2*M-F)
    - 0.000283*dsin(M1+3*F);


omega1=0.0004664*dcos(omega);
omega2=0.0000754*dcos(omega+275.05-2.30*T);
beta=B*(1-omega1-omega2);

pi_greco= + 0.950724
          + 0.051818*dcos(M1)
          + 0.009531*dcos(2*D-M1)
          + 0.007843*dcos(2*D)
          + 0.002824*dcos(2*M1)
          + 0.000857*dcos(2*D+M1)
          + e*0.000533*dcos(2*D-M)
          + e*0.000401*dcos(2*D-M-M1)
          + e*0.000320*dcos(M1-M)
          - 0.000271*dcos(D)
          - e*0.000264*dcos(M+M1)
          - 0.000198*dcos(2*F-M1)
          + 0.000173*dcos(3*M1)
          + 0.000167*dcos(4*D-M1)
          - e*0.000111*dcos(M)
          + 0.000103*dcos(4*D-2*M1)
          - 0.000084*dcos(2*M1-2*D)
          - e*0.000083*dcos(2*D+M)
          + 0.000079*dcos(2*D+2*M1)
          + 0.000072*dcos(4*D)
          + e*0.000064*dcos(2*D-M+M1)
          - e*0.000063*dcos(2*D+M-M1)
          + e*0.000041*dcos(M+D)
          + e*0.000035*dcos(2*M1-M)
          - 0.000033*dcos(3*M1-2*D)
          - 0.000030*dcos(M1+D)
          - 0.000029*dcos(2*F-2*D)
          - e*0.000029*dcos(2*M1+M)
          + e2*0.000026*dcos(2*D-2*M)
          - 0.000023*dcos(2*F-2*D+M1)
          + e*0.000019*dcos(4*D-M-M1);

distanza_luna=6378.14/dsin(pi_greco);


}

calcola_fase_lunare()
{
double cosd,d,M,M1;
float i,k;

double T2= pow(T,2);
double T3= pow(T,3);

M= 358.475833 + 35999.0498*T - 0.00015*T2 - 0.0000033*T3;
M1 = 296.104608 + 477198.8491*T + 0.009192*T2 + 0.0000144*T3;


cosd=dcos(lambda-LV)*dcos(beta);

d=arccos(cosd);

i=(float)(180-d-0.1468*((1-0.0549*dsin(M1))/(1-0.0167*dsin(M)))*dsin(d));

k=(1+dcos(i))/2;


printf("Frazione illuminata = %f\n",k);
}

calcola_separazione_angolare()
{
double T2=pow(T,2);
double T3=pow(T,3);
double epsilon=23.452294-0.0130125*T-0.00000164*T2+0.000000503*T3;
double ARsole,ARluna,decsole,decluna;
double separazione;

ARsole=arctan((dcos(epsilon)*dsin(LV))/dcos(LV));
decsole=arcsin(dsin(epsilon)*dsin(LV));

ARluna=arctan((dsin(lambda)*dcos(epsilon)-dtan(beta)*dsin(epsilon))/dcos(lambda));

decluna=arcsin(dsin(beta)*dcos(epsilon)+dcos(beta)*dsin(epsilon)*dsin(lambda));


separazione=arccos(dsin(decsole)*dsin(decluna)+dcos(decsole)*dcos(decluna)*dcos(ARsole-ARluna));

printf(" Separazione angolare = %.10g\n",separazione);
printf(" Distanza sole = %.10g UA\n",R);
printf(" Distanza luna = %.10g Km\n",distanza_luna);
}



/************************************************************/
/************************************************************/
/************************************************************/






double dsin(double angle)
{
double angolo,angrad;
angolo=riduci(angle);
angrad=(angolo*PI)/180;
return sin(angrad);
}


double dcos(double angle)
{
double angolo,angrad;
angolo=riduci(angle);
angrad=(angolo*PI)/180;
return cos(angrad);
}

double dtan(double angle)
{
double angolo,angrad;
angolo=riduci(angle);
angrad=(angolo*PI)/180;
return tan(angrad);
}

double riduci(double angle)
{
int detrai=0;
double angolo;
if (angle>0) {
	while (detrai<angle) detrai+=360;
	detrai-=360;
	angolo=angle-detrai;
	}
else {
	double assoluto=abs(angle);
	while(detrai<assoluto) detrai+=360;
	angolo=angle+detrai;
	}
return angolo;
}

double arccos(double val)
{
double angrad,angolo;
angrad=acos(val);
angolo=180*(angrad/PI);
return angolo;
}

double arcsin(double val)
{
double angrad,angolo;
angrad=asin(val);
angolo=180*(angrad/PI);
return angolo;
}

double arctan(double val)
{
double angrad,angolo;
angrad=atan(val);
angolo=180*(angrad/PI);
return angolo;
}

