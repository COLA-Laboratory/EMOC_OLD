#include "wfg.h"

static void trimLine(char line[])
{
	int i = 0;

	while(line[i] != '\0')
	{
		if (line[i] == '\r' || line[i] == '\n')
		{
			line[i] = '\0';
			break;
		}
		i++;
	}
}

void printContents(FILECONTENTS *f)
{
	int i , j, k;
	for (i = 0; i < f->nFronts; i++)
	{
		printf("Front %d:\n", i+1);
		for (j = 0; j < f->fronts[i].nPoints; j++)
		{
			printf("\t");
			for ( k = 0; k < f->fronts[i].n; k++)
			{
				printf("%f ", f->fronts[i].points[j].objectives[k]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

FILECONTENTS *read_data(void *ptr)
{
	population_real *pop = (population_real *) ptr;

	int front = 0, point = 0, objective = 0;
	// init the struct
	FILECONTENTS *fc = malloc(sizeof(FILECONTENTS));
	fc->nFronts = 0;
	fc->fronts = NULL;
	front = fc->nFronts;
	fc->nFronts++;
	fc->fronts = realloc(fc->fronts, sizeof(FRONT) * fc->nFronts);
	fc->fronts[front].nPoints = 0;
	fc->fronts[front].points = NULL;
	// read the data

	int i,j;

	for ( i = 0; i < popsize; i++)
	{
		FRONT *f = &fc->fronts[front];
		point = f->nPoints;
		f->nPoints++;
		f->points = realloc(f->points, sizeof(POINT) * f->nPoints);
		f->n = 0;
		f->points[point].objectives = NULL;

		for ( j = 0; j < number_objective; j++)
		{
			POINT *p = &f->points[point];
			objective = f->n;
			f->n++;
			p->objectives = realloc(p->objectives, sizeof(OBJECTIVE) * f->n);
			p->objectives[objective] = pop->ind[i].obj[j];
		}
	}

/*
	while(fgets(line, sizeof line, fp) != NULL)
	{
		trimLine(line);
		if (strcmp(line, "#") == 0)
		{
			front = fc->nFronts;
			fc->nFronts++;
			fc->fronts = realloc(fc->fronts, sizeof(FRONT) * fc->nFronts);
			fc->fronts[front].nPoints = 0;
			fc->fronts[front].points = NULL;
		}
		else
		{
			FRONT *f = &fc->fronts[front];
			point = f->nPoints;
			f->nPoints++;
			f->points = realloc(f->points, sizeof(POINT) * f->nPoints);
			f->n = 0;
			f->points[point].objectives = NULL;
			char *tok = strtok(line, " \t\n");
			do
			{
				POINT *p = &f->points[point];
				objective = f->n;
				f->n++;
				p->objectives = realloc(p->objectives, sizeof(OBJECTIVE) * f->n);
				p->objectives[objective] = atof(tok);
			} while ((tok = strtok(NULL, " \t\n")) != NULL);
		}
	}

*/
	//fc->nFronts--;
	// for (int i = 0; i < fc->nFronts; i++) fc->fronts[i].n = fc->fronts[i].points[0].nObjectives;

	//printf("Read %d fronts\n", fc->nFronts);
	//printContents(fc);
	return fc;
}
