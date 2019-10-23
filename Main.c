#define _CRT_SECURE_NO_WARNINGS

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Main.h"

#define ROW 1
#define COL 0
#define ASCENDING 0
#define DESCENDING 1
#define DATA_TAG 0
#define MASTER 0
#define POINTS_NUM_OF_ATTRIBUTES 3
#define DIMENSIONS 2

int main(int argc, char *argv[])
{
	int myId, numOfProcess, n;
	int coords[DIMENSIONS];
	Point* pointsArr = NULL;
	Point myPoint;
	const char* FILE_NAME = "C:\\Users\\morso\\source\\repos\\HW2_mor_sarel\\HW2_mor_sarel\\PointsFile.txt";

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);
	MPI_Comm cart_comm;
	MPI_Datatype PointType;

	if (myId == MASTER)
	{
		int numOfPoints;
		pointsArr = readPointsFromFile(FILE_NAME, &numOfPoints);
		n = (int)(sqrt(numOfPoints));

		printf("The points matrix before sort is :\n ");
		printMat(pointsArr, n);
		printf("The distances matrix before sort is :\n ");
		printDistances(pointsArr, n);
	}
		MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD); //Send to all slaves the number n
		
		if (numOfProcess != n * n)//if num of points does not match num of process, Abort and exit program
		{
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

	cart_comm = createCartComm(myId, n, coords);
	createPointType(&PointType);

	MPI_Scatter(pointsArr, 1, PointType, &myPoint, 1, PointType, MASTER, cart_comm);
	shearSort(cart_comm, myId, &myPoint, PointType, n, coords);
	MPI_Gather(&myPoint, 1, PointType, pointsArr, 1, PointType, MASTER, cart_comm);

	if (myId == MASTER)
	{
		//Print Result
		printf("The points matrix before sort is :\n ");
		printMat(pointsArr, n);
		printf("The distances matrix before sort is :\n ");
		printDistances(pointsArr, n);

	}

	MPI_Finalize();
	return 0;
}

Point* readPointsFromFile(const char* fileName, int* numOfPoints)
{
	Point* pointsArr;
	FILE* fp;
	int i;

	fp = fopen(fileName, "r");
	if (!fp)
	{
		printf("\nim here\n");
		return NULL;
	}
	fscanf(fp, "%d", numOfPoints);
	if (*numOfPoints == 0)
	{
		fclose(fp);
		return NULL;
	}

	pointsArr = (Point*)malloc((*numOfPoints) * sizeof(Point));
	if (!pointsArr)
	{
		fclose(fp);
		return NULL;
	}

	for (i = 0; i < *numOfPoints; i++)
	{
		fscanf(fp, "%lf", &(pointsArr[i].x));
		fscanf(fp, "%lf", &(pointsArr[i].y));
		fscanf(fp, "%lf", &(pointsArr[i].z));
	}
	fclose(fp);
	return pointsArr;
}

void printMat(Point* mat, int n)
{
	int i, j;
	printf("\n");
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			printPoint(&mat[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

MPI_Comm createCartComm(int myId, int n, int *coords)
{
	MPI_Comm cartComm;
	int nDims = DIMENSIONS;
	int dim[DIMENSIONS], period[DIMENSIONS], reorder;
	dim[0] = n;     // num of columns
	dim[1] = n;     // num of rows
	period[0] = 0;  // cols are cyclic
	period[1] = 0; // rows are cyclic
	reorder = 1; // allows changing the order of processes ids

	MPI_Cart_create(MPI_COMM_WORLD, nDims, dim, period, reorder, &cartComm);
	MPI_Cart_coords(cartComm, myId, nDims, coords);

	return cartComm;
}

void shearSort(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int coords[])
{
	int i;

	for (i = 0; i < 2 * log(n) + 1; i++)
	{
		//Rows First
		if (coords[0] % 2 == 0)
		{
			oddEvenSort(comm_cart, myId, ROW, DESCENDING, myPoint, PointType, n, coords);
		}
		else
		{
			oddEvenSort(comm_cart, myId, ROW, ASCENDING, myPoint, PointType, n, coords);
		}

		//Then Cols...
		oddEvenSort(comm_cart, myId, COL, DESCENDING, myPoint, PointType, n, coords);
	}
}

void oddEvenSort(MPI_Comm comm_cart, int myId, int rowOrCol, int sortOrder, Point* myPoint, MPI_Datatype PointType, int n, int coords[])
{
	if (rowOrCol == ROW)
	{
		rowSort(comm_cart, myId, sortOrder, myPoint, PointType, n, coords);
	}
	else //Col were chosen
	{
		colSort(comm_cart, myId, sortOrder, myPoint, PointType, n, coords);
	}
}

void rowSort(MPI_Comm comm_cart, int myId, int sortOrder, Point* myPoint, MPI_Datatype PointType, int n, int coords[])
{
	int step;

	for (step = 0; step < n; step++)
	{
		rowScenario(comm_cart, myId, myPoint, PointType, n, step, sortOrder, coords);
	}
}

void rowScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int step, int sortOrder, int coords[])
{
	int leftId, rightId;
	int disp = 1;
	MPI_Cart_shift(comm_cart, ROW, disp, &leftId, &rightId);

	if (step % 2 == 0)//if the step is even
	{
		rowEvenScenario(comm_cart, myId, myPoint, PointType, n, rightId, leftId, sortOrder, coords);
	}
	else // step is odd
	{
		rowOddScenario(comm_cart, myId, myPoint, PointType, n, rightId, leftId, sortOrder, coords);
	}
}

void rowEvenScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int rightId, int leftId, int sortOrder, int coords[])
{
	if (coords[1] % 2 == 0) //process id is even in row
	{
		sendReciveAndCompare(comm_cart, myPoint, PointType, myId, rightId, sortOrder);
	}
	else //process id is odd in row
	{
		sendReciveAndCompare(comm_cart, myPoint, PointType, myId, leftId, sortOrder);
	}
}

void rowOddScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int rightId, int leftId, int sortOrder, int coords[])
{
	MPI_Status status;
	Point hisPoint;

	if (coords[1] % 2 == 0) //process id is even in row
	{
		if (coords[1] != 0) //check if it is not the first or last process - edge case
		{
			sendReciveAndCompare(comm_cart, myPoint, PointType, myId, leftId, sortOrder);
		}
	}
	else //process id is odd in row
	{
		if (coords[1] != n - 1) //check if it is not the first or last process - edge case
		{
			sendReciveAndCompare(comm_cart, myPoint, PointType, myId, rightId, sortOrder);
		}

	}
}

void colSort(MPI_Comm comm_cart, int myId, int sortOrder, Point* myPoint, MPI_Datatype PointType, int n, int coords[])
{
	int step;

	for (step = 0; step < n; step++)
	{
		colScenario(comm_cart, myId, myPoint, PointType, n, step, sortOrder, coords);
	}
}

void colScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int step, int sortOrder, int coords[])
{
	int upId, downId;
	int disp = 1;
	MPI_Cart_shift(comm_cart, COL, disp, &upId, &downId);

	if (step % 2 == 0)//step is even
	{
		colEvenScenario(comm_cart, myId, myPoint, PointType, n, downId, upId, sortOrder, coords);
	}
	else //step is odd
	{
		colOddScenario(comm_cart, myId, myPoint, PointType, n, downId, upId, sortOrder, coords);
	}
}
void colEvenScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int downId, int upId, int sortOrder, int coords[])
{
	if (coords[0] % 2 == 0)//process id is even in col
	{
		sendReciveAndCompare(comm_cart, myPoint, PointType, myId, downId, sortOrder);
	}
	else //process id is odd in col
	{
		sendReciveAndCompare(comm_cart, myPoint, PointType, myId, upId, sortOrder);
	}
}

void colOddScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int downId, int upId, int sortOrder, int coords[])
{
	MPI_Status status;
	Point hisPoint;

	if (coords[0] % 2 == 0) //process id is even in col
	{
		if (coords[0] != 0) //check if it is not the first or last process - edge case
		{
			sendReciveAndCompare(comm_cart, myPoint, PointType, myId, upId, sortOrder);
		}
	}

	else //process id is odd in col
	{
		if (coords[0] != n - 1) //check if it is not the first or last process - edge case
		{
			sendReciveAndCompare(comm_cart, myPoint, PointType, myId, downId, sortOrder);
		}
	}
}

void sendReciveAndCompare(MPI_Comm comm_cart, Point* myPoint, MPI_Datatype PointType, int myId, int hisId, int sortOrder)
{
	Point hisPoint;
	MPI_Status status;
	MPI_Send(myPoint, 1, PointType, hisId, DATA_TAG, comm_cart);
	MPI_Recv(&hisPoint, 1, PointType, hisId, DATA_TAG, comm_cart, &status);

	if (sortOrder == ASCENDING) //Ascending order sort
	{
		if (compareDistances(myPoint, &hisPoint) > 0 && myId < hisId)
		{
			*myPoint = hisPoint;
		}
		if (compareDistances(myPoint, &hisPoint) < 0 && myId > hisId)
		{
			*myPoint = hisPoint;
		}
	}
	else //decending order sort
	{
		if (compareDistances(myPoint, &hisPoint) < 0 && myId < hisId)
		{
			*myPoint = hisPoint;
		}
		if (compareDistances(myPoint, &hisPoint) > 0 && myId > hisId)
		{
			*myPoint = hisPoint;
		}
	}
}

void createPointType(MPI_Datatype* PointType)
{
	int blockLengths[POINTS_NUM_OF_ATTRIBUTES] = { 1,1,1 };
	MPI_Aint disp[POINTS_NUM_OF_ATTRIBUTES];
	MPI_Datatype types[POINTS_NUM_OF_ATTRIBUTES] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };

	disp[0] = offsetof(Point, x);
	disp[1] = offsetof(Point, y);
	disp[2] = offsetof(Point, z);



	MPI_Type_create_struct(POINTS_NUM_OF_ATTRIBUTES, blockLengths, disp, types, PointType);
	MPI_Type_commit(PointType);
}

double calculateDistanceFromZero(Point* point)
{
	return sqrt((point->x)*(point->x) + (point->y)*(point->y) + (point->z)*(point->z));
}

int compareDistances(Point* point1, Point* point2) // Return 1 if point1 distance is bigger then point2 distance, else return -1
{
	if (calculateDistanceFromZero(point1) > calculateDistanceFromZero(point2))
		return 1;
	else
		return -1;

}

void printPoint(Point* point)
{
	printf("(%lf,\t%lf,\t%lf)\t", point->x, point->y, point->z);
}

void printDistances(Point* pointArr, int n)
{
	int i, j;
	printf("\n");
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			double distance = calculateDistanceFromZero(&(pointArr[i * n + j]));
			printf("%lf\t", distance);
		}
		printf("\n");
	}
	printf("\n");

}