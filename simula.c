#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define HEAD 0
#define TAIL -1
#define VSIT -2
#define EPTY -3
#define AAPATH -4
#define CIRCLE -5
#define ABPATH -6
#define GENERATION 10000
#define GENENUM 1000
#define runno 30

void encodeAdjacency(int* genome, int genen, int* adjacency);
int calculateDCJdistance(int* adjacency_a, int* adjacency_b, int length);
int getStartValue(int* adjacency, int length);
int getRturnValueOfSecond(int* adjacency, int length, int st);
int getRturnValueOfFirst(int* adjacency, int length, int st);
void searchCutNode(int *target_node, int* adjacency_a, int* adjacency_b, int length);
int getIndex(int val, int* adjacency_a, int length);
void simulated_annealing(int** passgenomes, int* median_genome, int length, int* bestMScore);
void createInitialTrueGenome(int** original_list, int length);
void createCandiateMedianGenome(int* previousMedianGenome, int** original_list, int length);
int calculateMedianScore(int* median_genome, int** original_list, int length);
double acceptanceProbability(int energy, int newEnergy, double temperature);
void decodeAdjacency(int* median_genome,  int* medianGenome_adj, int circularnum, int* indices, int length);
void DCJOperation(int* adjacency_a, int* adjacency_b, int length);

int checkTwoGenomeSame(int* genome_a, int* genome_b, int length);
void printAdjacency(int* adjacency, int length);
void adjacencyDecode(int* adjacency, int* median_genome, int length);
int checkNodeInCircle(int *cutNode, int* adjacency_a, int* adjacency_b, int length);
int getReturnFromAdj(int value, int* adjacency, int length);
void splitBigCircleIntoSmall(int* adjacency_a, int* adjacency_b, int length, int *cutNode);
void produceSmallCircle(int* target_node, int* adjacency_a, int* adjacency_b, int length);
void medinaOptimizaiton(int** original_list, int* median_genome, int length, int* bestMScore);
void createOptimizeCandiateMedianGenome(int* previousMedianGenome, int** original_list, int length);

int* copy_adj_a;
int* copy_adj_b;
int* temp_initial_genome;
int** candidateInitialTrue;
int* initialTrueGenome;
int* updateTrueGenome;
int** candidateMedianGenme;
int* currentTrueGenome;
int* bestGenome;
int* best_scorelist;
int* current_scorelist;
int* adjacency_a_index;
int* adjacency_b_index;
int* adjacency_index;
int* median_genome;
int* temp_genome;
int* converge;
int* final_medianGenome;
/*change starts*/
int* circleadjacnecy;
/*change ends*/
int** original_list;
double avg_optindex=0;

void encodeAdjacency(int* genome, int genen, int* adjacency){
	int i;
	
	//head
	adjacency[0] = HEAD;

	for (i=0; i<genen; i++)
	{
		if (genome[i]>0)
		{
			adjacency[i*2+1] = genome[i] * 2 - 1;
			adjacency[i*2+2] = genome[i] * 2;
		}
		else
		{
			adjacency[i*2+1] = genome[i] * 2 * (-1);
			adjacency[i*2+2] = genome[i] * 2 * (-1) - 1;
		}
	}

	//Tail
	adjacency[genen*2 + 1] = TAIL;
}

int getStartValue(int* adjacency, int length)
{
	int i;
	for(i=0; i<length; i++)
	{
		if (adjacency[i] != VSIT && adjacency[i] != HEAD && adjacency[i] != TAIL)
		{
			int val = adjacency[i];
			adjacency[i] = VSIT;
			return val;
		}
	}
	return EPTY;
}

int getRturnValueOfSecond(int* adjacency, int length, int st)
{
	int i;

	if (st == TAIL)
	{
		i = adjacency_b_index[length-1];
	}
	else
	{
		i = adjacency_b_index[st];
	}


	if (adjacency[i] == st)
	{
		int val;
		if (i%2 == 0)
		{
			val = adjacency[i+1];
			adjacency[i] = VSIT;
			if (val != TAIL) adjacency[i+1] = VSIT;
		}
		else
		{
			val = adjacency[i-1];
			adjacency[i] = VSIT;
			if (val != HEAD) adjacency[i-1] = VSIT;
			}
			if (val == HEAD || val == TAIL)
			{
				return ABPATH;
			}
			else
			{
				return val;
			}
	}
	

	return EPTY;
}

int getRturnValueOfFirst(int* adjacency, int length, int st)
{
	int i;

	/*debug

	printf("************\n");
	for (i=0; i<length; i++)
	{
		printf("---adjacency[i]---%d\n", adjacency[i]);
	}
	printf("^^^^^^^^^^^^^\n");
	*/

	if (st == TAIL)
	{
		i = adjacency_a_index[length-1];
	}
	else
	{
		i = adjacency_a_index[st];
	}


	if (adjacency[i] == st)
	{
		int val; 
		if (i%2 == 0)
		{
			val = adjacency[i+1];
			adjacency[i] = VSIT;
			// printf("val: %d\n", val);
			if (val != TAIL) adjacency[i+1] = VSIT;	
		}
		else
		{
			val = adjacency[i-1];
			adjacency[i] = VSIT;
			if (val != HEAD) adjacency[i-1] = VSIT;
		}

		if(val == TAIL || val == HEAD)
		{
			return AAPATH;
		}
		else if(val == VSIT)
		{
			return CIRCLE;
		}
		return val;
	}
	return EPTY;
}
//By using DCJ sorting to generate teh neighbour
void createCandiateMedianGenome(int* previousMedianGenome, int** original_list, int length)
{
	int i, k;
	for (i=0; i<3; i++)
	{
		//copy previousMedianGenome to a temperary genome 
		for (k=0; k<length; k++)
		{
			temp_initial_genome[k] = previousMedianGenome[k];
		}

		int distance = calculateDCJdistance(temp_initial_genome, original_list[i], length);
		//int distance1 = calculateDCJdistance(temp_initial_genome, original_list[i+1], length);
		//int distance2 = calculateDCJdistance(temp_initial_genome, original_list[i+2], length);

		if (distance > 0)
		{
			DCJOperation(temp_initial_genome, original_list[i], length);
		}

		/*
		int boolean = checkTwoGenomeSame(temp_initial_genome, original_list[i], length);
		if (boolean == 0)
		{
			DCJOperation(temp_initial_genome, original_list[i], length);
		}
		*/
		//store the new genome to candidatemedianGenome
		for (k=0; k<length; k++)
		{
			candidateMedianGenme[i][k] = temp_initial_genome[k];
		}
	}
}

void createInitialTrueGenome(int** original_list, int length)
{
	int i,j,k,l;
	l=0;
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			if (i!=j)
			{
				//copy original one genome to temp_initial_genome, includes circularnum and indices
				for (k=0; k<length; k++)
				{
					temp_initial_genome[k] = original_list[i][k];
				}
				int distance =  calculateDCJdistance(temp_initial_genome, original_list[j], length);
				distance = distance / 10;
				while(distance>0)
				{
					DCJOperation(temp_initial_genome, original_list[j], length);
					distance--;
				}
				//copy the temp_initial_genome to 
				for (k=0; k<length; k++)
				{
					candidateInitialTrue[l][k] = temp_initial_genome[k];
				}
				l++;
			}
		}
	}
}

void simulated_annealing(int** original_list, int* median_genome, int length, int* bestMScore)
{
	//step1:to initialize the population
	int i,r, temperature, formerState, currentState, bestMedianScore, M, index;
	int optindex[runno], N;
	double coolingRate;
	int sum_optindex=0;


	createInitialTrueGenome(original_list, length);
	//generate random value from 0-6;
	r = rand() % 6;	

	//copy to initialTrueGenome includes the circularnum and indices
	for (i=0; i<length; i++)
	{
		initialTrueGenome[i] = candidateInitialTrue[r][i];
	}

	//simulate annealing
	temperature = 10;

	coolingRate = 0.9;

	//store the initialTrueGenome to updateTrueGenome
	for (i=0; i<length; i++)
	{
		updateTrueGenome[i] = initialTrueGenome[i];
	}

	formerState = 0;
	currentState = 3*length;

	formerState = calculateMedianScore(initialTrueGenome, original_list, length);
	//printf("initialValue: %d", formerState);
	bestMedianScore = 3*length;

//step2:update the median genome and fitness
	M = GENERATION; //change the generation at two positions
	N = M;

	for(int run=0;run<runno;run++)
	{
		while(M>0)
		{
			index = -1;
			//(1)update the median genome crroding to theupdateTrueGenome
			createCandiateMedianGenome(updateTrueGenome, original_list, length);
			index = rand() % 3;	
			//calculate the fitness for each canditate median genome
			currentState = calculateMedianScore(candidateMedianGenme[index], original_list, length);
			current_scorelist[N-M] = currentState;
			//copy to currentTrueGenome

			for (i=0; i<length; i++)
			{
				currentTrueGenome[i] = candidateMedianGenme[index][i];
			}
			// decide if the current solution is better than the former state
			if (acceptanceProbability(formerState, currentState, temperature)> ((double) rand() / (RAND_MAX)))
			{
				formerState = currentState;
				if (currentState < bestMedianScore)
				{
					//add bestTrueGenome here
					for(i=0; i<length; i++)
					{
						bestGenome[i] = currentTrueGenome[i];
					}
					bestMedianScore = currentState;
				}
				//pass value from currentTrueGenome to updateTrueGenome
				for (i=0; i<length; i++)
				{
					updateTrueGenome[i] = currentTrueGenome[i];
				}
			}
			converge[GENERATION-M] = currentState;
			//converge[GENERATION-M] = bestMedianScore;
			temperature = temperature*coolingRate;
			M--;
		}
		//step3:save the gbest and fitness
	for (i=0; i<N; i++)
	{
		if (converge[i] == bestMedianScore)
		{
			optindex[run] = i;
			break;
		}
	}
	}
	

	for (int i = 0; i < runno; i++)
	{
		/* code */
		sum_optindex=sum_optindex+optindex[i];
	}
	avg_optindex=sum_optindex/runno;
	// printf("Index: %d\n", optindex);
	// printf("Best: %d\n", bestMedianScore);
	// printf("BEST NUM: %d\n", bestGenNum);
	//decodeAdjacency(median_genome, bestGenome, bestGenNum, bestGenomeIndices, length);
	*bestMScore = bestMedianScore;

	// for (i=0; i<bestGenNum; i++)
	// {
	// 	printf("IN DICE: %d %d\n", bestGenomeIndices[2*i], bestGenomeIndices[2*i+1]);
	// }

	//printf("BestScore: %d, optindex: %d\n", bestMedianScore, optindex);
	// for (i=0; i<bestGenNum; i++)
	// {
	// 	printf("indices0: %d, indices1: %d\n", bestGenomeIndices[i*2], bestGenomeIndices[i*2+1]);
	// }
	// // //print out the best genome
	// for(i=0; i<length; i++)
	// {
	// 	printf("%d ", bestGenome[i]);
	// }
	// printf("\n");
}

int calculateDCJdistance(int* adjacency_a, int* adjacency_b, int length)
{
	int i;
	//record the index of adjacency
	for (i=0; i<length; i++)
	{
		if (adjacency_a[i] == TAIL)
		{
			adjacency_a_index[length-1] = i;
		}
		else
		{
			adjacency_a_index[adjacency_a[i]] = i;
		}
	}

	for (i=0; i<length; i++)
	{
		if (adjacency_b[i] == TAIL)
		{
			adjacency_b_index[length-1] = i;
		}
		else
		{
			adjacency_b_index[adjacency_b[i]] = i;
		}
	}
	//deep copy

	for(i=0; i<length; i++)
	{
		copy_adj_a[i] = adjacency_a[i];
		copy_adj_b[i] = adjacency_b[i];
	}

	int ab_path = 0;
	int circle = 0;
	while(1){
		int st = getStartValue(copy_adj_a, length);
		if (st == EPTY)
		{
			break;
		}
		while (1)
		{
			//printf("bbbbb\n");
			int rv = getRturnValueOfSecond(copy_adj_b, length, st);
			//printf("rv: %d\n", rv);
			if (rv == ABPATH)
			{
				ab_path++;
				break;
			}
			st = getRturnValueOfFirst(copy_adj_a, length, rv);
			//printf("in st: %d\n", st);
			if (st == AAPATH)
			{
				break;
			}
			else if (st == CIRCLE)
			{
				circle++;
				break;
			}
			else if (st == EPTY)
			{
				break;
			}
		}
		//printf(".....\n");
	}
	// printf("circle: %d\n", circle);
	int gNum = length / 2 - 1;
	int distance = gNum - circle - ab_path/2;
	return distance;
}


void searchCutNode(int *target_node, int* adjacency_a, int* adjacency_b, int length)
{
	int r,i; 
	int temp_node[2];
	while(1)
	{
		r = rand() % (length);

		/* find a feasible target node*/
		if (r%2 == 0)
		{
			temp_node[0] = adjacency_a[r];
			temp_node[1] = adjacency_a[r+1];
		}
		else
		{
			temp_node[0] = adjacency_a[r-1];
			temp_node[1] = adjacency_a[r];
		}
		//check if the values are equal to adjacency b
		for (i=0; i<length; i++)
		{
			if (temp_node[0] == adjacency_b[i])
			{
				if (i%2 == 0)
				{
					if(temp_node[1] != adjacency_b[i+1])
					{
						target_node[0] = adjacency_b[i];
						target_node[1] = adjacency_b[i+1];
						return;
					}
				}
				else
				{
					if(temp_node[1] != adjacency_b[i-1])
					{
						target_node[0] = adjacency_b[i-1];
						target_node[1] = adjacency_b[i];
						return;
					}
				}
			}
		}
	}
}

/*new change start*/

int checkNodeInCircle(int *cutNode, int* adjacency_a, int* adjacency_b, int length)
{
	int first = cutNode[0];
	int second = cutNode[1];
	int returned_b, returned_a, start_a, start_b, step;
	step = 0;
	start_b = first;
	while(1)
	{
		returned_b = getReturnFromAdj(start_b, adjacency_b, length);
		step++;
		if (returned_b == HEAD || returned_b == TAIL) return -1;
		start_a = returned_b;
		returned_a = getReturnFromAdj(start_a, adjacency_a, length);
		step++;
		if (returned_b == HEAD || returned_b == TAIL) return -1;
		if (returned_a == first) return step;
		start_b = returned_a;
	}
	printf("error in checkNodeInCircle function!\n");
	return -1;
}


int getReturnFromAdj(int value, int* adjacency, int length)
{
	int i;
	for (i=0; i<length; i++)
	{
		if (adjacency[i] == value)
		{
			if (i%2==0)
			{
				return adjacency[i+1];
			}
			else
			{
				return adjacency[i-1];
			}
		}
	}
	printf("getReturnFromAdj function error!\n");
	return -1;
}

/*new change ends*/

int getIndex(int val, int* adjacency_a, int length)
{
	int i;
	for (i=0; i<length; i++)
	{
		if (val == adjacency_a[i])
		{
			return i;
		}
	}
	return -1;
}

// int findAdjVal(int val, int* adjacency_b, int length)
// {
// 	int i;
// 	for(i=0; i<length; i++)
// 	{
// 		if (adjacency_b[i] == val)
// 		{
// 			if (i%2 == 0)
// 			{
// 				return adjacency_b[i+1];
// 			}
// 			else
// 			{
// 				return adjacency_b[i-1];
// 			}
// 		}
// 	}
// 	return -1;
// }


int calculateMedianScore(int* median_genome, int** original_list, int length)
{
	int i,j;
	int distance = 0;
	for (i=0; i<3; i++)
	{
		// printf("582\n");
		distance += calculateDCJdistance(median_genome, original_list[i], length);
		// printf("584\n");
	}
	return distance;
}

void decodeAdjacency(int* median_genome, int* medianGenome_adj, int circularnum, int* indices, int length)
{
	int i, j, t, tail_index;
	if (circularnum == 0)
	{
		for (i=1; i<length-1; i=i+2)
		{
			if (medianGenome_adj[i] % 2 == 0)
			{
				median_genome[i/2] = (medianGenome_adj[i] / 2)*(-1);
			}
			else
			{
				median_genome[i/2] = medianGenome_adj[i+1] / 2;
			}
		}
	}
	else
	{
		//first decode the head tail list
		tail_index = indices[0] - 1;
		for (i=1; i<tail_index; i=i+2)
		{
			if (medianGenome_adj[i] % 2 == 0)
			{
				median_genome[i/2] = (medianGenome_adj[i] / 2)*(-1);
			}
			else
			{
				median_genome[i/2] = medianGenome_adj[i+1] / 2;
			}
		}
		for (i=0; i<circularnum; i++)
		{
			//move the last ele to head, and every other ele move to right one cell
			t = medianGenome_adj[indices[i*2+1]];
			for (j= indices[i*2]; j < indices[i*2+1]; j++)
			{
				medianGenome_adj[indices[i*2+1]+indices[i*2]-j]= medianGenome_adj[indices[i*2+1]+indices[i*2]-j-1];
			}
			medianGenome_adj[indices[i*2]] = t;
			for (j= indices[i*2]; j< indices[i*2+1]; j=j+2)
			{
				if (medianGenome_adj[j] % 2 == 0)
				{
					median_genome[j/2-1] = (medianGenome_adj[j] / 2)*(-1);
				}
				else
				{
					median_genome[j/2-1] = (medianGenome_adj[j+1] / 2);
				}
			}
		}
	}
}



double acceptanceProbability(int energy, int newEnergy, double temperature)
{
	if (newEnergy < energy)
	{
		return 1;
	}
	//double k3=0.0;
	//k3=exp((energy - newEnergy)/temperature);
	//return k3;
	return exp((energy - newEnergy)/temperature);
}

void readDataFromFile(int *arr, char *fileName)
{
    char buffer[12000];
    char *record;
    char *line;

    //open read in stream
    // char *path = "./data/";
    // char str[20];
    // strcpy(str, path);
    // strcat(str, fileName);
    FILE *fstream = fopen(fileName, "r");

    if(fstream == NULL)
    {
        printf("File opening failed!\n");
        exit(0);
    }

    char *ptr = "C:";
    int index = 0;
	while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
	{
		record = strtok(line," ");
		if (!strcmp(ptr, record))
		{
			while(record != NULL)
			{
				record = strtok(NULL, " ");
				if (record != NULL && *record != '\n')
				{
					*(arr+index) = atoi(record);
					index++;
				}
			}
		}
	}
   	fclose(fstream);
}

int checkTwoGenomeSame(int* genome_a, int* genome_b, int length)
{
	int i, j;
	for (i=0; i<length; i++)
	{
		if (genome_a[i] != genome_b[i])
		{
			return 0;
		}
	}
	return 1;
}


/*new change starts*/
void splitBigCircleIntoSmall(int* adjacency_a, int* adjacency_b, int length, int *cutNode)
{
	int start = cutNode[0];
	int index = 2;
	circleadjacnecy[0] = cutNode[0];
	circleadjacnecy[1] = cutNode[1];
	int start_a, start_b, returned_a, returned_b, first, second, first_adj, second_adj, first_index, second_index, i;
	start_b = start;
	while(1)
	{
		returned_b = getReturnFromAdj(start_b, adjacency_b, length);
		start_a = returned_b;
		returned_a = getReturnFromAdj(start_a, adjacency_a, length);
		if (returned_a == start) break;
		circleadjacnecy[index++] = start_a;
		circleadjacnecy[index++] = returned_a;
		start_b = returned_a;
	}
	// printf("index: %d\n", index);
	first = rand() % (index/2);
	while(1)
	{
		second = rand() % (index/2);
		if (first != second) break;
	}
	// printf("first: %d\n", first);
	// printf("second: %d\n", second);
	first_adj = circleadjacnecy[first*2];
	second_adj = circleadjacnecy[second*2];
	// find the corresponding two nodes
	for (i=0; i<length;i ++)
	{
		if (adjacency_a[i] == first_adj)
		{
			first_index = i;
		}
		if (adjacency_a[i] == second_adj)
		{
			second_index = i;
		}
	}
	// printf("first_index: %d\n", first_index);
	// printf("second_index: %d\n", second_index);
	int distance, new_distance;
	distance = calculateDCJdistance(adjacency_a, adjacency_b, length);
	int temp, delta;
	if (first_index%2==0 && second_index%2==0)
	{
		temp = adjacency_a[first_index+1];
		adjacency_a[first_index+1] = adjacency_a[second_index];
		adjacency_a[second_index] = temp;
		new_distance = calculateDCJdistance(adjacency_a, adjacency_b, length);
		delta = distance - new_distance;
		// printf("758: %d\n", delta);
		if (delta  == 1)
		{
			return;
		}
		else if (delta == 0)
		{
			temp = adjacency_a[first_index+1];
			adjacency_a[first_index+1] = adjacency_a[second_index+1];
			adjacency_a[second_index+1] = temp;
		}
	}
	else if (first_index%2==0 && second_index%2==1)
	{
		temp = adjacency_a[first_index+1];
		adjacency_a[first_index+1] = adjacency_a[second_index];
		adjacency_a[second_index] = temp;
		new_distance = calculateDCJdistance(adjacency_a, adjacency_b, length);
		delta = distance - new_distance;
		// printf("777: %d\n", delta);
		if (delta  == 1)
		{
			return;
		}
		else if (delta == 0)
		{
			temp = adjacency_a[first_index+1];
			adjacency_a[first_index+1] = adjacency_a[second_index-1];
			adjacency_a[second_index-1] = temp;
		}
	}
	else if (first_index%2==1 && second_index%2==0)
	{
		temp = adjacency_a[first_index-1];
		adjacency_a[first_index-1] = adjacency_a[second_index];
		adjacency_a[second_index] = temp;	
		new_distance = calculateDCJdistance(adjacency_a, adjacency_b, length);
		delta = distance - new_distance;
		// printf("795: %d\n", delta);
		if (delta  == 1)
		{
			return;
		}
		else if (delta == 0)
		{
			temp = adjacency_a[first_index-1];
			adjacency_a[first_index-1] = adjacency_a[second_index+1];
			adjacency_a[second_index+1] = temp;
		}
	}
	else
	{
		temp = adjacency_a[first_index-1];
		adjacency_a[first_index-1] = adjacency_a[second_index];
		adjacency_a[second_index] = temp;	
		new_distance = calculateDCJdistance(adjacency_a, adjacency_b, length);
		delta = distance - new_distance;
		// printf("815: %d\n", delta);
		if (delta  == 1)
		{
			return;
		}
		else if (delta == 0)
		{
			temp = adjacency_a[first_index-1];
			adjacency_a[first_index-1] = adjacency_a[second_index-1];
			adjacency_a[second_index-1] = temp;
		}
	}
}

void produceSmallCircle(int* target_node, int* adjacency_a, int* adjacency_b, int length)
{
	int i, fn_index, sn_index, temp;
	for (i=0; i<length; i++)
	{
		if (adjacency_a[i] == target_node[0])
		{
			fn_index = i;
		}
		if (adjacency_a[i] == target_node[1])
		{
			sn_index = i;
		}
	}
	if (fn_index % 2 == 0)
	{
		temp = adjacency_a[fn_index+1];
		adjacency_a[fn_index+1] = adjacency_a[sn_index];
		adjacency_a[sn_index] = temp;
	}
	else
	{
		temp = adjacency_a[fn_index-1];
		adjacency_a[fn_index-1] = adjacency_a[sn_index];
		adjacency_a[sn_index] = temp;
	}
}

void DCJOperation(int* adjacency_a, int* adjacency_b, int length)
{
	int i, fn_index, sn_index, temp;
	int target_node[2];
	searchCutNode(target_node, adjacency_a, adjacency_b, length);
	if (checkNodeInCircle(target_node, adjacency_a, adjacency_b, length) <= 4)
	{
		// printf("856\n");
		produceSmallCircle(target_node, adjacency_a, adjacency_b, length);
	}
	else 
	{
		// printf("861\n");
		splitBigCircleIntoSmall(adjacency_a, adjacency_b, length, target_node);
	}
	
}

/*new change ends*/

void printAdjacency(int* adjacency, int length)
{ 	int i;
	for (i=0; i<length; i++)
	{
		printf("%d ", adjacency[i]);
	}
	printf("\n");
}


void adjacencyDecode(int* adjacency, int* median_genome, int length)
{
	int i, index, start, start_index, val, flag;
	for (i=0; i<length; i++)
	{
		if (adjacency[i] == TAIL)
		{
			adjacency_index[length-1] = i;
		}
		else
		{
			adjacency_index[adjacency[i]] = i;
		}
	}
	//copy every adjacency into temp_genome
	for (i=0; i<length; i++)
	{
		temp_genome[i] = adjacency[i];
	}

	index = 0;

	start = HEAD;
	start_index = adjacency_index[start];
	temp_genome[start_index] = EPTY;
	while(1)
	{
		if (start_index % 2 == 0)
		{
			val = adjacency[start_index+1];
			temp_genome[start_index+1] = EPTY;
		}
		else
		{
			val = adjacency[start_index-1];
			temp_genome[start_index-1] = EPTY;
		}
		if (val == TAIL)
		{
			flag = 0;
			for (i=0; i<length; i++)
			{
				if (temp_genome[i] != EPTY)
				{
					flag = 1;
					start = temp_genome[i];
					temp_genome[i] = EPTY;
					start_index = adjacency_index[start];
					break;
				}
			}
			if (flag == 0) break;
		}
		else
		{
			if (val % 2 == 0)
			{
				median_genome[index++] = (val/2)*(-1);
				start = val - 1;
			}
			else
			{
				median_genome[index++] = (val+1) / 2;
				start = val + 1;
			}
			start_index = adjacency_index[start];
			if (temp_genome[start_index] == EPTY)
			{
				flag = 0;
				for (i=0; i<length; i++)
				{
					if (temp_genome[i] != EPTY)
					{
						flag = 1;
						start = temp_genome[i];
						temp_genome[i] = EPTY;
						start_index = adjacency_index[start];
						break;
					}
				}
				if (flag == 0) break;
			}
			else
			{
				temp_genome[start_index] = EPTY;
			}
		}
	}

}


/*new changes starts*/

void medinaOptimizaiton(int** original_list, int* median_genome, int length, int* bestMScore)
{
	int repeat, i, score, current_score, index, initialScore;
	repeat = 100;
	index = -1;
	// printf("989\n");
	initialScore = calculateMedianScore(median_genome, original_list, length);
	while(repeat>0)
	{
		score = length*3;
		// printf("993\n");
		createOptimizeCandiateMedianGenome(median_genome, original_list, length);
		for(i=0; i<3; i++)
		{
			// printf("998\n");
			current_score = calculateMedianScore(candidateMedianGenme[i], original_list, length);
			// printf("1000\n");
			if (current_score < score)
			{
				score = current_score;
				index = i;
			}
		}
		if (score < initialScore)
		{
			initialScore = score;
			for (i=0; i<length; i++)
			{
				bestGenome[i] = candidateMedianGenme[index][i];
				median_genome[i] = candidateMedianGenme[index][i];
			}
			*bestMScore = score;
		}
		repeat--;
	}
}

void createOptimizeCandiateMedianGenome(int* previousMedianGenome, int** original_list, int length)
{
	int i, k, step;
	for (i=0; i<3; i++)
	{
		//copy previousMedianGenome to a temperary genome 
		for (k=0; k<length; k++)
		{
			temp_initial_genome[k] = previousMedianGenome[k];
		}

		int distance = calculateDCJdistance(temp_initial_genome, original_list[i], length);

		if (distance > 0)
		{
			step = distance - 1;
			while (step>0)
			{
				DCJOperation(temp_initial_genome, original_list[i], length);
				step--;
			}
		}

		/*
		int boolean = checkTwoGenomeSame(temp_initial_genome, original_list[i], length);
		if (boolean == 0)
		{
			DCJOperation(temp_initial_genome, original_list[i], length);
		}
		*/
		//store the new genome to candidatemedianGenome
		for (k=0; k<length; k++)
		{
			candidateMedianGenme[i][k] = temp_initial_genome[k];
		}
	}	
}

/*new change ends*/

int main(int argc, char **argv)
{
	clock_t begin, end;
	double time_spent;

	begin = clock();

	srand(time(NULL));
	int i, bestMScore, optMScore, j, bestGNum, optMGenomeNum,tail_index,t;
	int genen = GENENUM; /* gene number*/
	int* genome = (int *)malloc(genen*sizeof(int));

	/*read data from files*/
	int *readingenes = (int *)malloc(genen*3*sizeof(int));
	/*	
	char* p;
	char* readinFileName;
	p = strtok(argv[1], "/");
	while(p)
	{
		p = strtok(NULL, "/");
		if(p)
		readinFileName = p;	
	}
	
	//printf("%s\n", readinFileName);
	*/
	readDataFromFile(readingenes, argv[1]);
	//print out read in data
	// for (i=0; i<genen*3; i++)
	// {
	// 	printf("%d ", readingenes[i]);
	// }

	int *g = (int *)malloc(genen*sizeof(int));
	int *g2 = (int *)malloc(genen*sizeof(int));
	int *g3 = (int *)malloc(genen*sizeof(int));

	for (i=0; i<genen; i++)
	{
		g[i] = readingenes[i];
	}

	for (i=0; i<genen; i++)
	{
		g2[i] = readingenes[i+genen];
	}

	for (i=0; i<genen; i++)
	{
		g3[i] = readingenes[i+genen*2];
	}

	free(readingenes);


	//int g[200] = {1,2,3,4,5,6,7,8,-180,60,-155,-154,-153,-152,190,191,192,193,-110,-109,-108,-99,-98,-97,-96,85,86,87,88,89,90,91,92,93,94,95,-84,-131,-130,-129,-128,-157,-156,61,62,63,-53,-52,32,33,34,177,178,-137,-136,-135,-134,-133,-132,-83,-82,-81,-80,-79,-78,-102,-164,-163,-40,103,104,-199,-198,-197,-196,-195,-194,111,112,179,-183,-45,-118,-51,-50,-49,-48,-47,-46,101,-175,-174,165,166,167,168,169,170,171,172,173,176,35,-127,-126,-125,-124,-123,-122,-121,-120,-119,-44,-43,-42,25,26,27,140,141,142,143,144,145,146,147,148,-21,-20,-19,-139,-138,113,114,115,116,117,-31,-30,158,159,160,161,162,41,-24,-23,-22,149,150,151,-189,-188,75,76,77,-39,-38,-37,-36,-29,-28,-18,-17,-16,-15,-14,-13,-12,-11,-10,185,186,187,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,54,55,56,57,58,59,184,-9,181,182,-100,-107,-106,-105,200};
	//int g2[200] = {1,2,3,4,88,-161,-160,-159,-158,-157,-156,-155,-154,-153,-152,190,178,113,114,87,5,6,7,8,9,10,11,12,13,14,15,16,17,168,169,170,171,79,80,81,82,83,84,85,86,115,116,-90,-89,162,163,164,-173,-172,-78,-77,109,110,-193,-192,-191,-177,53,54,55,56,-24,-23,-22,149,150,151,-189,-188,75,76,-108,-107,-106,21,-148,-147,-146,-179,-112,-111,194,195,196,197,198,199,-39,-38,-37,-36,-35,-176,102,103,104,-121,28,29,-18,-167,-166,-165,174,175,-101,46,47,48,49,120,30,31,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,-183,-182,-181,-180,184,-25,42,43,44,45,-100,-99,-98,-97,95,96,-94,-93,-92,-91,117,118,119,50,51,52,-34,-33,-32,-122,105,-20,-19,-27,-26,185,186,187,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-41,-40,200};
	//int g3[200] = {1,2,-23,-22,149,150,-21,-20,-19,-67,-66,-65,-64,199,-39,-38,-37,-36,-35,-176,102,103,104,-121,-29,-110,-109,77,78,148,-187,-186,-185,26,27,68,69,70,71,72,73,74,151,12,13,14,15,177,178,113,114,115,116,-90,-89,162,163,164,193,-28,-18,91,92,93,40,41,-24,3,4,88,-161,-160,-159,-158,-157,-156,-155,-154,-153,-152,190,191,192,-173,-172,9,10,-179,-112,-111,194,195,196,32,33,34,87,5,6,7,-129,-128,-127,-126,-125,-124,-123,-31,-30,-120,-49,-48,-60,-59,-58,-57,-56,-55,-54,-53,136,-141,-140,-139,-138,-137,61,166,167,117,118,119,50,51,52,-86,-85,-84,181,182,183,-17,-16,-135,63,-198,-197,-122,105,106,107,108,-76,-75,188,189,-11,146,147,-42,25,-184,180,-83,-82,-81,-80,-79,142,143,144,145,168,169,170,171,-8,130,131,132,-46,101,-175,-174,165,62,-134,133,47,43,44,45,-100,-99,-98,-97,95,96,-94,200};


	int length = genen*2+2;//memory size
	int* adjacency_a = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(g, genen, adjacency_a);


	int* adjacency_b = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(g2, genen, adjacency_b);

	int* adjacency_c = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(g3, genen, adjacency_c);

	/*chnage starts*/
	circleadjacnecy = (int *)malloc((length)*sizeof(int));
	/*change ends*/

	//print out the prev list
	// for (i=0; i<length; i++)
	// {
	// 	printf("%d ", adjacency_a[i]);
	// }

	copy_adj_a = (int *)malloc(length*sizeof(int));
	copy_adj_b = (int *)malloc(length*sizeof(int));
	temp_initial_genome = (int *)malloc(length*sizeof(int));
	candidateInitialTrue = (int **)malloc(6 * sizeof(int *));
	initialTrueGenome = (int *)malloc(length*sizeof(int));
	updateTrueGenome = (int *)malloc(length * sizeof(int));
	currentTrueGenome = (int *)malloc(length* sizeof(int));
	bestGenome = (int *)malloc(length* sizeof(int));
	best_scorelist = (int *)malloc(GENERATION* sizeof(int));
	current_scorelist = (int *)malloc(GENERATION* sizeof(int));
	adjacency_a_index = (int *)malloc(length* sizeof(int));
	adjacency_b_index = (int *)malloc(length* sizeof(int));
	adjacency_index = (int *)malloc(length* sizeof(int));
	median_genome = (int *)malloc(genen*sizeof(int));
	temp_genome = (int *)malloc(length*sizeof(int));
	converge = (int *)malloc(GENERATION*sizeof(int));
	final_medianGenome = (int *)malloc(GENENUM*sizeof(int));

	original_list = (int **)malloc(3*sizeof(int*));


	/*
	int dist = calculateDCJdistance(adjacency_a, adjacency_b, length);


	while(dist>0)
	{
		DCJOperation(adjacency_a, adjacency_b, length);
		dist = calculateDCJdistance(adjacency_a, adjacency_b, length);
		printf("dist: %d\n", dist);
	}

	printf("print result\n");
	for (i=0; i<length; i++)
	{
		printf("%d ", adjacency_a[i]);
	}
	printf("\n");
	adjacencyDecode(adjacency_a, median_genome, length);
	for (i=0; i<GENENUM; i++)
	{
		printf("%d ", median_genome[i]);
	}
	*/

	for (i=0; i<6; i++)
    {
    	candidateInitialTrue[i] = (int *)malloc(length * sizeof(int));
    }
    candidateMedianGenme = (int **)malloc( 3* sizeof(int *));
    for (i=0; i<3; i++)
    {
    	original_list[i] = (int *)malloc(length * sizeof(int));
    	candidateMedianGenme[i] = (int *)malloc(length * sizeof(int));
    }
	for (i=0; i<length; i++)
	{
		original_list[0][i] = adjacency_a[i];
		original_list[1][i] = adjacency_b[i];
		original_list[2][i] = adjacency_c[i];
	}

	int *bestMedian = (int *)malloc(length*sizeof(int));
	int *optIndices = (int *)malloc(genen*sizeof(int));
	int* dupBestMedian = (int *)malloc(length*sizeof(int));
	int *trueGenome = (int *)malloc(genen*sizeof(int));
	int* trueGenomeAdj = (int *)malloc((length)*sizeof(int));
	int *optMedian = (int *)malloc((length)*sizeof(int));
	optMScore = genen*3;
	bestMScore = 0;
	bestGNum = 0;


	// int abdistance = calculateDCJdistance(original_list[0], original_list[1], length);
	// printf("abdistance: %d\n", abdistance);
	// int numbers = 10;
	// while(abdistance>0)
	// {
	// 	DCJOperation(original_list[0], original_list[1], length);
	// 	abdistance = calculateDCJdistance(original_list[0], original_list[1], length);
	// 	printf("distance: %d\n", abdistance);
	// 	// numbers--;
	// }

	// exit(0);

	//simulated_annealing(original_list, median_genome, length, &bestMScore);
	/*
	printf("best score: %d\n", bestMScore);
	for (i=0; i<length; i++)
	{
		printf("%d ", bestGenome[i]);
	}
	printf("\n");
	*/
	
	// for (i=0; i<GENENUM; i++)
	// {
	// 	printf ("%d ", final_medianGenome[i]);
	// }
	// printf("\n");
	
	// printf("best score: %d\n", bestMScore);	
	// end = clock();
	// time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	// printf("time: %.2f\n", time_spent);
	//optMGenomeNum = 0;
	

	//run ten times select the best one
	
	for (i=0; i<1; i++)
	{
		simulated_annealing(original_list, median_genome, length, &bestMScore);
		if (bestMScore < optMScore)
		{
			optMScore = bestMScore;
			for (j=0; j<length; j++)
			{
				bestMedian[j] = bestGenome[j];
				//optMedian[j] = bestGenome[j];
			}
			for (j=0; j<GENERATION; j++)
			{
				best_scorelist[j] = current_scorelist[j];
			}
		}
	}

	for(i=0; i<length; i++)
	{
		optMedian[i] = bestMedian[i];
	}
	medinaOptimizaiton(original_list, optMedian, length, &bestMScore);
	for(i=0; i<length; i++)
	{
		bestMedian[i] = bestGenome[i];
	}
	
	
	//copy best median genome
	for (i=0; i<length; i++)
	{
		dupBestMedian[i] = bestMedian[i];
	}

	// printf("OUT %d\n", optMGenomeNum);  

 //    printf("optMScore: %d\n", optMScore);

    // for (i=0; i<optMGenomeNum; i++)
    // {
    // 	printf("optIndices: %d %d", optIndices[i*2], optIndices[i*2+1]);
    // }

    // for (i=0; i<length; i++)
    // {
    // 	printf("%d ", bestMedian[i]);
    // }

    //calculate best median to real distance
    
    for (i=0; i<genen; i++)
    {
    	trueGenome[i] = i+1;
    }

    
	encodeAdjacency(trueGenome, genen, trueGenomeAdj);
	

    
    //write the result to out file
    char *out = ".out";
    char str[40];
    strcpy(str, argv[1]);
    strcat(str, out);

    FILE *f = fopen(str, "w");
    if (f == NULL)
    {
       printf("Error opening file!\n");
       exit(1);
    }

   // print the Test betsGenome
    int TestbestGenome[length];
	fprintf(f,"->TestbestGenome:\n");
  	for(i=0; i<length; i++)
	{
		TestbestGenome[i] = bestMedian[i];
		fprintf(f, "%d ", TestbestGenome[i]);
	}
	fprintf(f,"\n");

    fprintf(f, "Optimal Median Score: %d \n", optMScore);

    fprintf(f, "Medain Genome:\n");

    adjacencyDecode(bestMedian, final_medianGenome, length);

	fprintf(f,"->C:\n");
	for (i=0; i<GENENUM; i++)
	{
		fprintf(f, "%d ", final_medianGenome[i]);
	}
	fprintf(f,"\n");

    fprintf(f, "To true distance: %d \n", calculateDCJdistance(dupBestMedian,trueGenomeAdj, length));
    fprintf(f, "Convergence:\n");
    for (i=0; i<GENERATION; i++)
    {
		fprintf(f, "%d ", best_scorelist[i]);
    }
    fprintf(f, "\n");
    fprintf(f, "average_coveragence:%.2f\n",avg_optindex);

    end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//printf("time: %.2f\n", time_spent);
	fprintf(f, "running time: %.2f\n", time_spent);
	fclose(f);
	

	/* multiple chromsome to genome    

    int* medianAdj = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(median_genome, genen, medianAdj);

	printf("merge Chromesome: %d\n", calculateMedianScore(medianAdj, original_list, length));
	
	free(medianAdj);
	*/ 
	
    /*
	int dist = calculateDCJdistance(adjacency_a, adjacency_b, length);

	printf("distance: %d\n", dist);
	int circularNum = 0;
	int *indices = (int*)calloc(2*genen, sizeof(int));

	//int current_distance = 0;
	while(dist>0)
	{	
		printf("start while\n");
		for (i=0; i<length; i++)
		{
			printf("%d ", adjacency_a[i]);
		}
		printf("\n");
		printf("before DCJOperation\n");
		DCJOperation(adjacency_a, adjacency_b, length, &circularNum, indices);
		printf("after DCJOperation:\n");
		printf("circularNum: %d\n", circularNum);
		for (i=0; i<circularNum; i++)
		{
			printf("indice0: %d, indice1: %d\n", indices[i*2], indices[i*2+1]);
		}
		for (i=0; i<length; i++)
		{
			printf("%d ", adjacency_a[i]);
		}
		printf("\n");
		dist = calculateDCJdistance(adjacency_a, adjacency_b, length);
		printf("current: %d\n", dist);
		// printf("after calculateDCJdistance:\n");
		// for (i=0; i<length; i++)
		// {
		// 	printf("%d ", adjacency_a[i]);
		// }
		// printf("end while\n");
	}
	*/
	// for (i=0; i<length; i++)
	// {
	// 	printf("%d ", adjacency_a[i]);
	// }

	//DCJOperation(adjacency_a, adjacency_b, length, &circularNum, indices);

	//print out the new list
	// for (i=0; i<length; i++)
	// {
	// 	printf("%d ", adjacency_a[i]);
	// }

	free(dupBestMedian);
	free(trueGenomeAdj);
	free(trueGenome);
	free(optIndices);
	free(median_genome);
	free(genome);
	free(adjacency_a);
	free(adjacency_b);
	free(adjacency_c);
    free(copy_adj_a);
	free(copy_adj_b);
	free(temp_initial_genome);
	free(initialTrueGenome);
	free(updateTrueGenome);
	free(currentTrueGenome);
	for (i=0; i<6; i++)
    {
    	free(candidateInitialTrue[i]); 
    }
    for (i=0; i<3; i++)
    {
    	free(original_list[i]);
    	free(candidateMedianGenme[i]);
    }
    free(original_list);
    free(candidateMedianGenme);
	free(candidateInitialTrue);
	free(bestGenome);
	free(g);
	free(g2);
	free(g3);
	free(bestMedian);
	free(best_scorelist); 
	free(current_scorelist);
	free(adjacency_a_index);
	free(adjacency_b_index);
	free(adjacency_index);
	free(temp_genome);
	free(converge);
	free(final_medianGenome);
	/*add new starts*/
	free(circleadjacnecy);
	free(optMedian);
	/*add new ends*/

	return 0;
}
