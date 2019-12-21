#include "pch.h"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#define PRES double
#define NXB 15
#define NYB 12
#define NX NXB*3
#define NY NYB*3
#define NYK2 NYB*2
#define REP 3000
#define DEL 100
#define AMAT 1.1f
#define TEM1 5.0f
#define TEM2 15.0f
#define HX 0.2f
#define HY 0.3f

using namespace std;

int main(int argc, char **argv)
{
	int i, j, k;
	int idt = 0;
	int ndt = 0;
	PRES T1 = TEM1, 
		T2 = TEM2, 
		h = HX, 
		r = HY, 
		a = AMAT, t0;
	PRES T[NY+1][NX+1], 
		TT[NY+1][NX+1];
	PRES rr;
	if (h < r)
		rr = h;
	else
		rr = r;
	PRES tau = 0.25f*rr*rr / a;
	PRES alf_1 = -h / r;
	PRES alf_2 = -r / h;
	PRES alf_3 = 0.5f * alf_2;
	PRES alf_4 = 0.5f * alf_1;
	PRES bet_1 = a * tau / (h*r);
	PRES bet_2 = 2.0f*bet_1;
	PRES bet_3 = 4.0f*bet_1;
	PRES bet_4 = 4.0f*bet_1/3;
	PRES gam_1 = -2.f*(alf_1 + alf_2);
	PRES gam_2 = -1.5*(alf_1 + alf_2);
	PRES gam_3 = -(alf_1 + alf_2);
	PRES gam_4 = -(alf_3 + alf_4);

	char filename[128];

	for (j = 0; j <= NY; j++) 
		for (i = 0; i <= NX; i++) 
			T[j][i] = 0.0f;

	for (i = 0; i <= NXB; i++)
	{
		T[NY][i] = T1; 
		TT[NY][i] = T1;
	}
	for (i = NXB * 2; i <= NX; i++)
	{
		T[0][i] = T2; 
		TT[0][i] = T2;
	}


	ofstream fout("c:/work/T1.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	for (j = 0; j <= NY; j++) 
		for (i = 0; i <= NX; i++) 
		{
			PRES w = T[j][i];
			fout.write((char*)&w, sizeof w);
		}
	fout.close();

	for (k = 0; k < REP; k++) 
	{
		for (j = 0; j <= NY; j++) 
		{
			for (i = 0; i <= NX; i++) 
			{
				t0 = T[j][i];
				//D
				if (i == NX && j == NY)
					TT[j][i] = t0 - bet_3 * (alf_4*T[j - 1][i] + alf_3 * T[j][i - 1] + gam_4 * t0);
				//E
				else if (i == 0 && j == NYB * 2)
					TT[j][i] = t0 - bet_3 * (alf_4*T[j + 1][i] + alf_3 * T[j][i + 1] + gam_4 * t0);
				//F
				else if (i == NXB && j == NYB * 2)
					TT[j][i] = t0 - bet_4 * (alf_4*T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_2 * t0);
				//G
				else if (i == NXB && j == NYB)
					TT[j][i] = t0 - bet_3 * (alf_4*T[j + 1][i] + alf_3 * T[j][i + 1] + gam_4 * t0);
				//H
				else if (i == NXB * 2 && j == NYB)
					TT[j][i] = t0 - bet_4 * (alf_4*T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_2 * t0);

				//горизонтали
				//BD
				else if (NXB < i && i < NX && j == NY)
					TT[j][i] = t0 - bet_2 * (alf_3*T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j - 1][i] + gam_3 * t0);

				//EF
				else if (0 < i && i < NXB && j == NYB * 2)
					TT[j][i] = t0 - bet_2 * (alf_3*T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_3 * t0);

				//GH
				else if (NXB < i && i < NXB * 2 && j == NYB)
					TT[j][i] = t0 - bet_2 * (alf_3*T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_3 * t0);

				//вертикали
				//AE
				else if (i == 0 && NYB * 2 < j && j < NY)
					TT[j][i] = t0 - bet_2 * (alf_4*T[j - 1][i] + alf_4 * T[j + 1][i] + alf_2 * T[j][i + 1] + gam_3 * t0);

				//FG
				else if (i == NXB && NYB < j && j < NYB * 2)
					TT[j][i] = t0 - bet_2 * (alf_4*T[j - 1][i] + alf_4 * T[j + 1][i] + alf_2 * T[j][i + 1] + gam_3 * t0);

				//HI
				else if (i == NXB * 2 && 0 < j && j < NYB)
					TT[j][i] = t0 - bet_2 * (alf_4*T[j - 1][i] + alf_4 * T[j + 1][i] + alf_2 * T[j][i + 1] + gam_3 * t0);

				//DJ
				else if (i == NX && 0 < j && j < NY)
					TT[j][i] = t0 - bet_2 * (alf_4*T[j - 1][i] + alf_4 * T[j + 1][i] + alf_2 * T[j][i - 1] + gam_3 * t0);

				//области
				//ABEF
				else if (0 < i && i <= NXB && NYB * 2 < j && j < NY)
					TT[j][i] = t0 - bet_1 * (alf_1*T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_1 * t0);

				//BCHG
				else if (NXB < i && i <= NXB * 2 && NYB < j && j < NY)
					TT[j][i] = t0 - bet_1 * (alf_1*T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_1 * t0);

				//CDJI
				else if (NXB * 2 < i && i < NX && 0 < j && j < NY)
					TT[j][i] = t0 - bet_1 * (alf_1*T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] + gam_1 * t0);

			}
		}

		for (j=0; j<=NY; j++) 
		{
			for (i=0; i<=NX; i++)
			{
				T[j][i] = TT[j][i];
			} 
		}

		idt++;
		if (idt == DEL) 
		{ 
			idt=0; 
			ndt++;
			sprintf_s(filename, sizeof(filename), "c:/work/T%d.dat", ndt+1);
			ofstream fout(filename,ios_base::out | ios_base::trunc | ios_base::binary);
			for (j = 0; j <= NY; j++) 
			{  
				for (i = 0; i <= NX; i++) 
				{
					PRES w = T[j][i];
					fout.write((char*)&w, sizeof w); 
				} 
			}
			fout.close(); 
		}
	}
	int n_x = NX + 1; 
	int n_y = NY + 1; 
	int n_k = ndt;
	ofstream fou("c:/work/Param.dat",ios_base::out | ios_base::trunc | ios_base::binary);
	fou.write((char*)&n_x, sizeof n_x);
	fou.write((char*)&n_y, sizeof n_y);
	fou.write((char*)&n_k, sizeof n_y);
	fou.close();
	for (j = 0; j <= 5; j++)
	{
		for (i = 0; i <= 4; i++)
		{
			cout <<setw (20)<< T[j][i];
		}
		cout << endl;
	}
}