/*****************************************************************************
 *                         Autocorrelation Function                          *
 *                              Pedro H Mendes                               *
 ****************************************************************************/

/*****************************************************************************
 *                               INCLUDES                                    *
 ****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*****************************************************************************
 *                               DEFINITIONS                                 *
 ****************************************************************************/
#define NDATA	1000000		//1e6
#define KMAX	400
#define NHEADER	2

/*****************************************************************************
 *                             MAIN PROGRAM                                  *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	int i, k, t;
	double col_data, col_extra, data[NDATA];
	double Oi, Oi2, var_Oi;
	double sum_Oi, sum_Oik, sum_OiOik;
	double Ak;
	FILE *data_file;

	Oi = 0.0;
	Oi2 = 0.0;

	data_file = fopen(argv[1], "r");

	for(i=0; i<NDATA; i++)
	{
		if(i < NHEADER)
		{
			fscanf(data_file, "%*[^\n]");
		}
			fscanf(data_file, "%d %lf %lf", &t, &col_extra, &col_data);
			col_data = fabs(col_data);
			data[i] = col_data;
			Oi += col_data;
			Oi2 += (col_data*col_data);
	}

	fclose(data_file);

	Oi = Oi/(1.*NDATA);
	Oi2 = Oi2/(1.*NDATA);
	var_Oi = Oi2 - (Oi*Oi);

	printf("%d\t%lf\n", 0, 1.);

	for(k=0; k<KMAX; k++)
	{
		sum_Oi = 0.;
		sum_Oik = 0.;
		sum_OiOik = 0.;

		for(i=0; i<(NDATA-k); i++)
		{
			sum_Oi += data[i];
			sum_Oik += data[i+k];
			sum_OiOik += (data[i]*data[i+k]);
		}

		sum_Oi = sum_Oi/(1.*(NDATA-k));
		sum_Oik = sum_Oik/(1.*(NDATA-k));
		sum_OiOik = sum_OiOik/(1.*(NDATA-k));

		Ak = (sum_OiOik - (sum_Oi*sum_Oik))/var_Oi;

		printf("%d\t%lf\n", k, Ak);
	}

	return 0;
}
