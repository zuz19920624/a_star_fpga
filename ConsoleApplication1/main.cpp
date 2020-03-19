// ConsoleApplication1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <stack>
#include <set>
#include <ctime>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <opencv.hpp>
// A C++ Program to implement A* Search Algorithm 
using namespace std;

#define ROW 256 
#define COL 256
#define ROW_LOG2 8
#define COL_LOG2 8
#define CELL_NUM_IN_BLOCK 8 //the number of cell in block
#define BLOCK_NUM_IN_CACHE 15 // the number of blocks in cache
#define CELl_NUM_LOG2      3
#define NUM_OF_BLOCK_IN_X ( (ROW + CELL_NUM_IN_BLOCK - 1) / CELL_NUM_IN_BLOCK)
#define NUM_OF_BLOCK_IN_Y ( (COL + CELL_NUM_IN_BLOCK - 1) / CELL_NUM_IN_BLOCK)
#define SIZE_OF_CACHE (CELL_NUM_IN_BLOCK * CELL_NUM_IN_BLOCK * 4 * BLOCK_NUM_IN_CACHE) //counting in cells
#define SIZE_OF_CELL  10
#define PENALTY_OF_CACHE_MISS 10
#define TOTAL_MAPS 10000
#define TOTAL_PROBLEM 10
#define RANDOM    10
#define MAX_VALUE 0xffffffff

long long numOfOperations = 0;
long long current_numOfOps = 0;
long long maxiumOpenList = 0;
typedef long clock_t;
long long numOfCacheMiss = 0;
long long current_cacheMiss = 0;
double path_length = 0;

// Creating a shortcut for int, int pair type 
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type 
typedef pair<int, pair<int, int>> pPair;

// A structure to hold the neccesary parameters 
struct cell
{
	// Row and Column index of its parent 
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1 
	int parent_i, parent_j;
	// f = g + h 
	int f, g, h;
}; //16B

struct block
{
	cell c[CELL_NUM_IN_BLOCK][CELL_NUM_IN_BLOCK];
};

struct blockCache
{
	block b[2][BLOCK_NUM_IN_CACHE];   //2 way
	int tag[2][BLOCK_NUM_IN_CACHE];   //{block_j, block_i}, indicating the block inside
	int init[2][BLOCK_NUM_IN_CACHE];  //{block_j, block_i}, indicating the block inside
	int dirty[2][BLOCK_NUM_IN_CACHE]; //indicating the cell is dirty or not
	int lru[BLOCK_NUM_IN_CACHE];
};

blockCache bcache[4]; //4 banks

//functions related to cache

//transfer the block addr into cache addr
int transferBankAddr(int block_i, int block_j)
{
	return ((block_i % 2) + (block_j % 2) * 2 );
}

int transferBlockAddr(int block_i, int block_j)
{
	return ( ((block_i / 2) % BLOCK_NUM_IN_CACHE) + ((block_j / 2)  << (ROW_LOG2 - 1 - CELl_NUM_LOG2) ) ) % BLOCK_NUM_IN_CACHE;
}

//write the data of block into memory
void writeBlockIntoMem(int cache_addr, int way)
{

}
//read the data of block from memory into block cache
void copycell(cell* c1, cell* c2)
{
	c1->f = c2->f;
	c1->g = c2->g;
	c1->h = c2->h;
	c1->parent_i = c2->parent_i;
	c1->parent_j = c2->parent_j;
}

void writeCache(int cache_addr, int way, int block_i, int block_j, cell cellDetails[][COL], int bank_addr)
{
	for (int i = 0; i < CELL_NUM_IN_BLOCK; i++)
		for (int j = 0; j < CELL_NUM_IN_BLOCK; j++)
		{
			int cell_x = block_i * CELL_NUM_IN_BLOCK + i;
			int cell_y = block_j * CELL_NUM_IN_BLOCK + j;
			copycell(&bcache[bank_addr].b[way][cache_addr].c[i][j], &cellDetails[cell_x][cell_y]);
			
		}
	bcache[bank_addr].init[way][cache_addr] = 1;
	bcache[bank_addr].tag[way][cache_addr] = block_j * 128 + block_i;
}

void readBlockFromMem(int block_i, int block_j, cell cellDetails[][COL])
{
	numOfCacheMiss++;
	//read all cells from memory
	int bank_addr = transferBankAddr(block_i, block_j);
	int cache_addr = transferBlockAddr(block_i, block_j);
	//find the position in cache
	if (bcache[bank_addr].init[0][cache_addr] == 1 && bcache[bank_addr].dirty[0][cache_addr] == 0xffff)
	{
		writeBlockIntoMem(cache_addr, 0);
		//write data into current position
		writeCache(cache_addr, 0, block_i, block_j, cellDetails, bank_addr);
		bcache[bank_addr].lru[cache_addr] = 0;
	}
	else if (bcache[bank_addr].init[1][cache_addr] == 1 && bcache[bank_addr].dirty[1][cache_addr] == 0xffff)
	{
		writeBlockIntoMem(cache_addr, 1);
		//write data into current position
		writeCache(cache_addr, 1, block_i, block_j, cellDetails, bank_addr);
		bcache[bank_addr].lru[cache_addr] = 1;
	}
	else if (bcache[bank_addr].init[0][cache_addr] == 0)
	{
		//write data into current position
		writeCache(cache_addr, 0, block_i, block_j, cellDetails, bank_addr);
		bcache[bank_addr].lru[cache_addr] = 0;
	}
	else if (bcache[bank_addr].init[1][cache_addr] == 0)
	{
		//write data into current position
		writeCache(cache_addr, 1, block_i, block_j, cellDetails, bank_addr);
		bcache[bank_addr].lru[cache_addr] = 1;
	}
	else if(bcache[bank_addr].lru[cache_addr] == 0)
	{
		writeBlockIntoMem(cache_addr, 1);
		//write data into current position
		writeCache(cache_addr, 1, block_i, block_j, cellDetails, bank_addr);
		bcache[bank_addr].lru[cache_addr] = 1;
	}
	else if (bcache[bank_addr].lru[cache_addr] == 1)
	{
		writeBlockIntoMem(cache_addr, 0);
		//write data into current position
		writeCache(cache_addr, 0, block_i, block_j, cellDetails, bank_addr);
		bcache[bank_addr].lru[cache_addr] = 0;
	}
}

bool blockInCache(int block_i, int block_j, int* way)
{
	int bank_addr = transferBankAddr(block_i, block_j);
	int cache_addr = transferBlockAddr(block_i, block_j);
	int tag = block_j * 128 + block_i;
	if (tag == bcache[bank_addr].tag[0][cache_addr] && bcache[bank_addr].init[0][cache_addr])
	{
		*way = 0;
		return true;
	}
	else if (tag == bcache[bank_addr].tag[1][cache_addr] && bcache[bank_addr].init[1][cache_addr])
	{
		*way = 1;
		return true;
	}
	else
	{
		return false;
	}
}

//read the info of cell
void readCellFromBlockCache(int cell_x, int cell_y, cell cellDetails[][COL])
{
	int block_x = cell_x / CELL_NUM_IN_BLOCK;
	int block_y = cell_y / CELL_NUM_IN_BLOCK;
	int way;
	//judge the whether the block is in the cache
	if (blockInCache(block_x, block_y, &way) == true) //in cache
	{
		return;
	}
	else
	{
		readBlockFromMem(block_x, block_y, cellDetails);
		return;
	}
}

// A Utility Function to check whether given cell (row, col) 
// is a valid cell or not. 
bool isValid(int row, int col)
{
	// Returns true if row number and column number 
	// is in range 
	return (row >= 0) && (row < ROW) &&
		(col >= 0) && (col < COL);
}

// A Utility Function to check whether the given cell is 
// blocked or not 
bool isUnBlocked(int** grid, int row, int col)
{
	// Returns true if the cell is not blocked else false 
	if (grid[row][col] == 1)
		return (true);
	else
		return (false);
}

// A Utility Function to check whether destination cell has 
// been reached or not 
bool isDestination(int row, int col, Pair dest)
{
	if (row == dest.first && col == dest.second)
		return (true);
	else
		return (false);
}

// A Utility Function to calculate the 'h' heuristics. 
double calculateHValue(int row, int col, Pair dest)
{
	// Return using the distance formula 
	//return ((double)sqrt((row - dest.first) * (row - dest.first)
	//	+ (col - dest.second) * (col - dest.second)));

	//return the Chebyshev  distance
	int x_abs = abs(row - dest.first);
	int y_abs = abs(col - dest.second);
	return (x_abs > y_abs ? x_abs * 10 : y_abs * 10);
}

void writingmaps(int** cellDetails, char* fileName)
{
	ofstream ofile;
	ofile.open(fileName, 'w');
	if (ofile.is_open())
	{
		for(int i = 0; i < ROW; i++)
			for (int j = 0; j < COL; j++)
			{
				ofile << cellDetails[i][j];
			}
		ofile.close();
	}
}

// A Utility Function to trace the path from the source 
// to destination 
void tracePath(cell cellDetails[][COL], Pair dest)
{
	//printf("\nThe Path is ");
	int row = dest.first;
	int col = dest.second;

	stack<Pair> Path;

	while (!(cellDetails[row][col].parent_i == row
		&& cellDetails[row][col].parent_j == col))
	{
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row][col].parent_i;
		int temp_col = cellDetails[row][col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	while (!Path.empty())
	{
		pair<int, int> p = Path.top();
		Path.pop();
		//printf("-> (%d,%d) ", p.first, p.second);
	}

	path_length += cellDetails[dest.first][dest.second].g;

	return;
}

// A Function to find the shortest path between 
// a given source cell to a destination cell according 
// to A* Search Algorithm 
int aStarSearch(int** grid, Pair src, Pair dest)
{
	// If the source is out of range 
	if (isValid(src.first, src.second) == false)
	{
		//printf("Source is invalid\n");
		return 0;
	}

	// If the destination is out of range 
	if (isValid(dest.first, dest.second) == false)
	{
		//printf("Destination is invalid\n");
		return 0;
	}

	// Either the source or the destination is blocked 
	if (isUnBlocked(grid, src.first, src.second) == false ||
		isUnBlocked(grid, dest.first, dest.second) == false)
	{
		//printf("Source or the destination is blocked\n");
		return 0;
	}

	// If the destination cell is the same as source cell 
	if (isDestination(src.first, src.second, dest) == true)
	{
		//printf("We are already at the destination\n");
		return 0;
	}

	// Create a closed list and initialise it to false which means 
	// that no cell has been included yet 
	// This closed list is implemented as a boolean 2D array 
	bool closedList[ROW][COL];
	memset(closedList, false, sizeof(closedList));

	// Declare a 2D array of structure to hold the details 
	//of that cell 
	cell cellDetails[ROW][COL];

	int i, j;

	for (i = 0; i < ROW; i++)
	{
		for (j = 0; j < COL; j++)
		{
			cellDetails[i][j].f = MAX_VALUE;
			cellDetails[i][j].g = MAX_VALUE;
			cellDetails[i][j].h = MAX_VALUE;
			cellDetails[i][j].parent_i = -1;
			cellDetails[i][j].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node 
	i = src.first, j = src.second;
	cellDetails[i][j].f = 0;
	cellDetails[i][j].g = 0;
	cellDetails[i][j].h = 0;
	cellDetails[i][j].parent_i = i;
	cellDetails[i][j].parent_j = j;

	/*
	Create an open list having information as-
	<f, <i, j>>
	where f = g + h,
	and i, j are the row and column index of that cell
	Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	This open list is implenented as a set of pair of pair.*/
	set<pPair> openList;

	// Put the starting cell on the open list and set its 
	// 'f' as 0 
	openList.insert(make_pair(0, make_pair(i, j)));

	// We set this boolean value as false as initially 
	// the destination is not reached. 
	bool foundDest = false;
	
	int temp;

	while (!openList.empty())
	{
		numOfOperations++;
		pPair p = *openList.begin();
		if (maxiumOpenList < openList.size())
			maxiumOpenList = openList.size();


		// Remove this vertex from the open list 
		openList.erase(openList.begin());

		// Add this vertex to the closed list 
		i = p.second.first;
		j = p.second.second;
		closedList[i][j] = true;

		/*
			Generating all the 8 successor of this cell

				N.W N N.E
				\ | /
				\ | /
				W----Cell----E
					/ | \
				/ | \
				S.W S S.E

			Cell-->Popped Cell (i, j)
			N --> North	 (i-1, j)
			S --> South	 (i+1, j)
			E --> East	 (i, j+1)
			W --> West		 (i, j-1)
			N.E--> North-East (i-1, j+1)
			N.W--> North-West (i-1, j-1)
			S.E--> South-East (i+1, j+1)
			S.W--> South-West (i+1, j-1)*/

			// To store the 'g', 'h' and 'f' of the 8 successors 
		int gNew, hNew, fNew;

		//----------- 1st Successor (North) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j) == true)
		{
			readCellFromBlockCache(i-1, j, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j].parent_i = i;
				cellDetails[i - 1][j].parent_j = j;
				cellDetails[i - 1][j].g = cellDetails[i][j].g + 10;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}
			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j] == false &&
				isUnBlocked(grid, i - 1, j) == true)
			{
				gNew = cellDetails[i][j].g + 10;
				hNew = calculateHValue(i - 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j].f == MAX_VALUE ||
					cellDetails[i - 1][j].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i - 1, j)));

					// Update the details of this cell 
					cellDetails[i - 1][j].f = fNew;
					cellDetails[i - 1][j].g = gNew;
					cellDetails[i - 1][j].h = hNew;
					cellDetails[i - 1][j].parent_i = i;
					cellDetails[i - 1][j].parent_j = j;
				}
			}
		}

		//----------- 2nd Successor (South) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j) == true)
		{
			readCellFromBlockCache(i + 1, j, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j].parent_i = i;
				cellDetails[i + 1][j].parent_j = j;
				cellDetails[i + 1][j].g = cellDetails[i][j].g + 10;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}
			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j] == false &&
				isUnBlocked(grid, i + 1, j) == true)
			{
				gNew = cellDetails[i][j].g + 10;
				hNew = calculateHValue(i + 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j].f == MAX_VALUE ||
					cellDetails[i + 1][j].f > fNew)
				{
					openList.insert(make_pair(fNew, make_pair(i + 1, j)));
					// Update the details of this cell 
					cellDetails[i + 1][j].f = fNew;
					cellDetails[i + 1][j].g = gNew;
					cellDetails[i + 1][j].h = hNew;
					cellDetails[i + 1][j].parent_i = i;
					cellDetails[i + 1][j].parent_j = j;
				}
			}
		}

		//----------- 3rd Successor (East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i, j + 1) == true)
		{
			readCellFromBlockCache(i , j + 1, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i][j + 1].parent_i = i;
				cellDetails[i][j + 1].parent_j = j;
				cellDetails[i][j + 1].g = cellDetails[i][j].g + 10;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i][j + 1] == false &&
				isUnBlocked(grid, i, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 10;
				hNew = calculateHValue(i, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i][j + 1].f == MAX_VALUE ||
					cellDetails[i][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i, j + 1)));

					// Update the details of this cell 
					cellDetails[i][j + 1].f = fNew;
					cellDetails[i][j + 1].g = gNew;
					cellDetails[i][j + 1].h = hNew;
					cellDetails[i][j + 1].parent_i = i;
					cellDetails[i][j + 1].parent_j = j;
				}
			}
		}

		//----------- 4th Successor (West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i, j - 1) == true)
		{
			readCellFromBlockCache(i, j - 1, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i][j - 1].parent_i = i;
				cellDetails[i][j - 1].parent_j = j;
				cellDetails[i][j - 1].g = cellDetails[i][j].g + 10;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i][j - 1] == false &&
				isUnBlocked(grid, i, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 10;
				hNew = calculateHValue(i, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i][j - 1].f == MAX_VALUE ||
					cellDetails[i][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i, j - 1)));

					// Update the details of this cell 
					cellDetails[i][j - 1].f = fNew;
					cellDetails[i][j - 1].g = gNew;
					cellDetails[i][j - 1].h = hNew;
					cellDetails[i][j - 1].parent_i = i;
					cellDetails[i][j - 1].parent_j = j;
				}
			}
		}

		//----------- 5th Successor (North-East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j + 1) == true)
		{
			readCellFromBlockCache(i - 1, j + 1, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j + 1].parent_i = i;
				cellDetails[i - 1][j + 1].parent_j = j;
				cellDetails[i - 1][j + 1].g = cellDetails[i][j].g + 14;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j + 1] == false &&
				isUnBlocked(grid, i - 1, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 14;
				hNew = calculateHValue(i - 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j + 1].f == MAX_VALUE ||
					cellDetails[i - 1][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i - 1, j + 1)));

					// Update the details of this cell 
					cellDetails[i - 1][j + 1].f = fNew;
					cellDetails[i - 1][j + 1].g = gNew;
					cellDetails[i - 1][j + 1].h = hNew;
					cellDetails[i - 1][j + 1].parent_i = i;
					cellDetails[i - 1][j + 1].parent_j = j;
				}
			}
		}

		//----------- 6th Successor (North-West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j - 1) == true)
		{
			readCellFromBlockCache(i - 1, j - 1, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j - 1].parent_i = i;
				cellDetails[i - 1][j - 1].parent_j = j;
				cellDetails[i - 1][j - 1].g = cellDetails[i][j].g + 14;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j - 1] == false &&
				isUnBlocked(grid, i - 1, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 14;
				hNew = calculateHValue(i - 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j - 1].f == MAX_VALUE ||
					cellDetails[i - 1][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew, make_pair(i - 1, j - 1)));
					// Update the details of this cell 
					cellDetails[i - 1][j - 1].f = fNew;
					cellDetails[i - 1][j - 1].g = gNew;
					cellDetails[i - 1][j - 1].h = hNew;
					cellDetails[i - 1][j - 1].parent_i = i;
					cellDetails[i - 1][j - 1].parent_j = j;
				}
			}
		}

		//----------- 7th Successor (South-East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j + 1) == true)
		{
			readCellFromBlockCache(i + 1, j + 1, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j + 1].parent_i = i;
				cellDetails[i + 1][j + 1].parent_j = j;
				cellDetails[i + 1][j + 1].g = cellDetails[i][j].g + 14;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j + 1] == false &&
				isUnBlocked(grid, i + 1, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 14;
				hNew = calculateHValue(i + 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j + 1].f == MAX_VALUE ||
					cellDetails[i + 1][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i + 1, j + 1)));

					// Update the details of this cell 
					cellDetails[i + 1][j + 1].f = fNew;
					cellDetails[i + 1][j + 1].g = gNew;
					cellDetails[i + 1][j + 1].h = hNew;
					cellDetails[i + 1][j + 1].parent_i = i;
					cellDetails[i + 1][j + 1].parent_j = j;
				}
			}
		}

		//----------- 8th Successor (South-West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j - 1) == true)
		{
			readCellFromBlockCache(i + 1, j - 1, cellDetails);
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j - 1].parent_i = i;
				cellDetails[i + 1][j - 1].parent_j = j;
				cellDetails[i + 1][j - 1].g = cellDetails[i][j].g + 14;
				//printf("The destination cell is found\n");
				tracePath(cellDetails, dest);
				foundDest = true;
				return 1;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j - 1] == false &&
				isUnBlocked(grid, i + 1, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 14;
				hNew = calculateHValue(i + 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j - 1].f == MAX_VALUE ||
					cellDetails[i + 1][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i + 1, j - 1)));

					// Update the details of this cell 
					cellDetails[i + 1][j - 1].f = fNew;
					cellDetails[i + 1][j - 1].g = gNew;
					cellDetails[i + 1][j - 1].h = hNew;
					cellDetails[i + 1][j - 1].parent_i = i;
					cellDetails[i + 1][j - 1].parent_j = j;
				}
			}
		}
	}

	// When the destination cell is not found and the open 
	// list is empty, then we conclude that we failed to 
	// reach the destiantion cell. This may happen when the 
	// there is no way to destination cell (due to blockages) 
	if (foundDest == false);
		//printf("Failed to find the Destination Cell\n");

	return foundDest;
}


// Driver program to test above function 
int main()
{
	/* Description of the Grid-
	1--> The cell is not blocked
	0--> The cell is blocked */
	clock_t startTime, endTime;
	int total_time = 0;
	
	//int grid[ROW][COL];
	int** grid;
	//initialize grid
	grid = new int* [ROW];
	for (int i = 0; i < ROW; i++)
		grid[i] = new int[COL];

	//fixed start/stop cases
	for (int random = 0; random <= 50; random += 10)
	{
		for (int count = 0; count < TOTAL_MAPS; count++)
		{
			for (int i = 0; i < ROW; i++)
				for (int j = 0; j < COL; j++)
				{
					int k = rand() % 100;
					if (k >= random)
						grid[i][j] = 1;
					else
						grid[i][j] = 0;
				}

			int start_x = rand() % 255;
			int start_y = rand() % 255;
			int stop_x = rand() % 255;
			int stop_y = rand() % 255;

			// Source is the left-most bottom-most corner 
			Pair src = make_pair(0, 0);

			// Destination is the left-most top-most corner 
			Pair dest = make_pair(255, 255);

			current_numOfOps = numOfOperations;
			current_cacheMiss = numOfCacheMiss;
			startTime = clock();
			int findPath = aStarSearch(grid, src, dest);
			endTime = clock();




			if (!findPath)
			{
				count--;
				numOfOperations = current_numOfOps;
				numOfCacheMiss = current_cacheMiss;
			}
			else
			{
				char grid_file_name[200];
				sprintf_s(grid_file_name, "grid_file_worst_%d_%d", random, count);
				writingmaps(grid, grid_file_name);
				total_time += endTime - startTime;
			}
			printf("%d", count);
		}
		int cache_miss_pen;
		cache_miss_pen = numOfCacheMiss * (CELL_NUM_IN_BLOCK * CELL_NUM_IN_BLOCK + PENALTY_OF_CACHE_MISS) / TOTAL_MAPS;

		std::cout << "\ntotal operations : " << numOfOperations / TOTAL_MAPS << std::endl;
		std::cout << "total maps : " << TOTAL_MAPS << std::endl;
		std::cout << "total time : " << total_time / TOTAL_MAPS << std::endl;
		std::cout << "maxium open list size" << maxiumOpenList << std::endl;
		std::cout << "cache miss " << numOfCacheMiss / TOTAL_MAPS << std::endl;
		std::cout << "penalty " << cache_miss_pen << std::endl;
		std::cout << "total number of cells in block " << SIZE_OF_CACHE << std::endl;
		std::cout << "total cache size " << SIZE_OF_CACHE * SIZE_OF_CELL << std::endl;
		std::cout << "path length : " << path_length / 10.0 / TOTAL_MAPS << std::endl;

		

		std::string fileName = "result_worse_";
		char sRandom[10];
		ofstream ofile;
		sprintf_s(sRandom, "%d", random);
		fileName += sRandom;
		ofile.open(fileName, 'w');
		if (!ofile.is_open())
		{
			std::cout << "cannot open file " << fileName << std::endl;
			continue;
		}
		ofile << path_length / 10.0 / TOTAL_MAPS << ' '
			<< numOfOperations / TOTAL_MAPS << ' '
			<< total_time / TOTAL_MAPS << ' '
			<< (numOfOperations / TOTAL_MAPS + cache_miss_pen / 12) / 100.0 / 1000.0 << std::endl;
		ofile.close();

		numOfOperations = 0;
		numOfCacheMiss = 0;
		numOfOperations = 0;
		path_length = 0;
		total_time = 0;

	}

	//random start/stop cases
	for (int random = 0; random <= 50; random += 10)
	{
		for (int count = 0; count < TOTAL_MAPS; count++)
		{
			for (int i = 0; i < ROW; i++)
				for (int j = 0; j < COL; j++)
				{
					int k = rand() % 100;
					if (k >= random)
						grid[i][j] = 1;
					else
						grid[i][j] = 0;
				}

			int start_x = rand() % 255;
			int start_y = rand() % 255;
			int stop_x = rand() % 255;
			int stop_y = rand() % 255;

			// Source is the left-most bottom-most corner 
			Pair src = make_pair(start_x, start_y);

			// Destination is the left-most top-most corner 
			Pair dest = make_pair(stop_x, stop_y);

			current_numOfOps = numOfOperations;
			current_cacheMiss = numOfCacheMiss;
			startTime = clock();
			int findPath = aStarSearch(grid, src, dest);
			endTime = clock();




			if (!findPath)
			{
				count--;
				numOfOperations = current_numOfOps;
				numOfCacheMiss = current_cacheMiss;
			}
			else
			{
				char grid_file_name[200];
				sprintf_s(grid_file_name, "grid_file_random_%d_%d", random, count);
				writingmaps(grid, grid_file_name);
				total_time += endTime - startTime;
			}
			printf("%d", count);
		}
		int cache_miss_pen;
		cache_miss_pen = numOfCacheMiss * (CELL_NUM_IN_BLOCK * CELL_NUM_IN_BLOCK + PENALTY_OF_CACHE_MISS) / TOTAL_MAPS;

		std::cout << "\ntotal operations : " << numOfOperations / TOTAL_MAPS << std::endl;
		std::cout << "total maps : " << TOTAL_MAPS << std::endl;
		std::cout << "total time : " << total_time / TOTAL_MAPS << std::endl;
		std::cout << "maxium open list size" << maxiumOpenList << std::endl;
		std::cout << "cache miss " << numOfCacheMiss / TOTAL_MAPS << std::endl;
		std::cout << "penalty " << cache_miss_pen << std::endl;
		std::cout << "total number of cells in block " << SIZE_OF_CACHE << std::endl;
		std::cout << "total cache size " << SIZE_OF_CACHE * SIZE_OF_CELL << std::endl;
		std::cout << "path length : " << path_length / 10.0 / TOTAL_MAPS << std::endl;



		std::string fileName = "result_random";
		char sRandom[10];
		ofstream ofile;
		sprintf_s(sRandom, "%d", random);
		fileName += sRandom;
		ofile.open(fileName, 'w');
		if (!ofile.is_open())
		{
			std::cout << "cannot open file " << fileName << std::endl;
			continue;
		}
		ofile << path_length / 10.0 / TOTAL_MAPS << ' '
			<< numOfOperations / TOTAL_MAPS << ' '
			<< total_time / TOTAL_MAPS << ' '
			<< (numOfOperations / TOTAL_MAPS + cache_miss_pen / 12) / 100.0 / 1000.0 << std::endl;
		ofile.close();

		numOfOperations = 0;
		numOfCacheMiss = 0;
		numOfOperations = 0;
		path_length = 0;
		total_time = 0;

	}

	//loading maps
	std::string imageName("testmap.png");
	cv::Mat src;
	src = cv::imread(imageName);
	for (int random = 0; random <= 50; random += 10)
	{
		for (int count = 0; count < 1; count++)
		{
			for (int i = 0; i < ROW; i++)
				for (int j = 0; j < COL; j++)
				{
					int k = src.at<cv::Vec3b>(i, j)[0];
					if (k >= 128)
						grid[i][j] = 1;
					else
						grid[i][j] = 0;
				}

			int start_x = rand() % 255;
			int start_y = rand() % 255;
			int stop_x = rand() % 255;
			int stop_y = rand() % 255;

			// Source is the left-most bottom-most corner 
			Pair src = make_pair(0, 0);

			// Destination is the left-most top-most corner 
			Pair dest = make_pair(255, 255);

			current_numOfOps = numOfOperations;
			current_cacheMiss = numOfCacheMiss;
			startTime = clock();
			int findPath = aStarSearch(grid, src, dest);
			endTime = clock();




			if (!findPath)
			{
				count--;
				numOfOperations = current_numOfOps;
				numOfCacheMiss = current_cacheMiss;
			}
			else
			{
				char grid_file_name[200];
				sprintf_s(grid_file_name, "grid_file_worst_%d_%d", random, count);
				writingmaps(grid, grid_file_name);
				total_time += endTime - startTime;
			}
			printf("%d", count);
		}
		long long  cache_miss_pen;
		cache_miss_pen = numOfCacheMiss * (CELL_NUM_IN_BLOCK * CELL_NUM_IN_BLOCK + PENALTY_OF_CACHE_MISS) / TOTAL_MAPS;

		std::cout << "\ntotal operations : " << numOfOperations / TOTAL_MAPS << std::endl;
		std::cout << "total maps : " << TOTAL_MAPS << std::endl;
		std::cout << "total time : " << total_time / TOTAL_MAPS << std::endl;
		std::cout << "maxium open list size" << maxiumOpenList << std::endl;
		std::cout << "cache miss " << numOfCacheMiss / TOTAL_MAPS << std::endl;
		std::cout << "penalty " << cache_miss_pen << std::endl;
		std::cout << "total number of cells in block " << SIZE_OF_CACHE << std::endl;
		std::cout << "total cache size " << SIZE_OF_CACHE * SIZE_OF_CELL << std::endl;
		std::cout << "path length : " << path_length / 10.0 / TOTAL_MAPS << std::endl;



		std::string fileName = "result_rts_";
		char sRandom[10];
		ofstream ofile;
		sprintf_s(sRandom, "%d", random);
		fileName += sRandom;
		ofile.open(fileName, 'w');
		if (!ofile.is_open())
		{
			std::cout << "cannot open file " << fileName << std::endl;
			continue;
		}
		ofile << path_length / 10.0 / TOTAL_MAPS << ' '
			<< numOfOperations / TOTAL_MAPS << ' '
			<< total_time / TOTAL_MAPS << ' '
			<< (numOfOperations / TOTAL_MAPS + cache_miss_pen / 12) / 100.0 / 1000.0 << std::endl;
		ofile.close();

		numOfOperations = 0;
		numOfCacheMiss = 0;
		numOfOperations = 0;
		path_length = 0;
		total_time = 0;

	}



	return(0);
}


// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
