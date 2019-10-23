#pragma once


typedef struct
{
	double x;
	double y;
	double z;
}Point;

Point* readPointsFromFile(const char* fileName, int* numOfPoints);
void printMat(Point* mat, int n);
MPI_Comm createCartComm(int myId, int n, int *coords);
void shearSort(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int coords[]);
void oddEvenSort(MPI_Comm comm_cart, int myId, int rowOrCol, int sortOrder, Point* myPoint, MPI_Datatype PointType, int n, /*int index,*/ int coords[]);
void createPointType(MPI_Datatype* PointType);
double calculateDistanceFromZero(Point* point);
int compareDistances(Point* point1, Point* point2);
void rowSort(MPI_Comm comm_cart, int myId, int sortOrder, Point* myPoint, MPI_Datatype PointType, int n, int coords[]);
void rowScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int step, int sortOrder, int coords[]);
void rowEvenScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int rightId, int leftId, int sortOrder, int coords[]);
void rowOddScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int rightId, int leftId, int sortOrder, int coords[]);
void colSort(MPI_Comm comm_cart, int myId, int sortOrder, Point* myPoint, MPI_Datatype PointType, int n, int coords[]);
void colScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int step, int sortOrder, int coords[]);
void colEvenScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int downId, int upId, int sortOrder, int coords[]);
void colOddScenario(MPI_Comm comm_cart, int myId, Point* myPoint, MPI_Datatype PointType, int n, int downId, int upId, int sortOrder, int coords[]);
void sendReciveAndCompare(MPI_Comm comm_cart, Point* myPoint, MPI_Datatype PointType, int myId, int hisId, int sortOrder);
void printPoint(Point* point);

void printDistances(Point* pointArr, int n);
