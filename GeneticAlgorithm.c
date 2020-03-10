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
#define GENERATION 100
#define GENENUM 1000
#define GENOMELENGTH 2002

void encodeAdjacency(int* genome, int genen, int* adjacency);
int calculateDCJdistance(int* adjacency_a, int* adjacency_b, int length);
int getStartValue(int* adjacency, int length);
int getRturnValueOfSecond(int* adjacency, int length, int st);
int getRturnValueOfFirst(int* adjacency, int length, int st);
void DCJOperation(int* adjacency_a, int* adjacency_b, int length, int *circularNum, int *indices);
void searchCutNode(int *target_node, int* adjacency_a, int* adjacency_b, int length);
int getIndex(int val, int* adjacency_a, int length);
void simulated_annealing(int** passgenomes, int* median_genome, int length, int* bestMScore, int* bestGenomeNum);
void createInitialTrueGenome(int** original_list, int length);
void createCandiateMedianGenome(int* previousMedianGenome, int circular_num, int* indices, int** original_list, int length);
int calculateMedianScore(int* median_genome, int** original_list, int length);
double acceptanceProbability(int energy, int newEnergy, double temperature);
void decodeAdjacency(int* median_genome,  int* medianGenome_adj, int circularnum, int* indices, int length);

void generateInitialPool(int** original_list, int length);
void geneticAlgorithm(int** original_list, int length);
void quickSort(int low, int high);
void shuffle(int *array, size_t n);
void geneticAlgorithmForParents(int parent1_index, int parent2_index, int** original_list, int length);

int* copy_adj_a;
int* copy_adj_b;
int* temp_list;
int* reverse_list;
int* temp_list2;
int* internal_templist;
int* temp_htlist;
int* temp_initial_genome;
int* temp_indices;
int** candidateInitialTrue;
int* initilcandid_circular_num;
int** initilcandid_indices;
int* initialTrueGenome;
int* initialTrueGenomeIndices;
int* updateTrueGenome;
int* updateTrueGenomeIndices;
int** candidateMedianGenme;
int** candidMedian_indices;
int* candidMedian_circular_num;
int* currentTrueGenome;
int* currentTrueGenomeIndices;
int* bestGenome;
int* bestGenomeIndices;
int** genome_pool;
int** g_pool_indices;
int* g_pool_cnum;
// int *global_circularNum;
int *global_indices;
int *original_copy_a;
int *original_copy_b;
int *fitnessScore;
int *temporary_genome;
int *temporary_indices;
int *best_genome;
int *best_genomeIndices;
int *index_array;
int *parent1;
int *parent2;
int *parent1Indices;
int *parent2Indices;
int *parent1Copy;
int *parent2Copy;
int *parent1CopyIndices;
int *parent2CopyIndices;
int **candidateGenomes;
int **candidateIndices;

int* converge;

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

	for (i=0; i<length; i++)
	{
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
	for (i=0; i<length; i++)
	{
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
	}
	return EPTY;
}

int calculateDCJdistance(int* adjacency_a, int* adjacency_b, int length)
{
	int i;
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
		// printf("st %d\n", st);
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

	// printf("ab_path: %d\n", ab_path);
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

void DCJOperation(int* adjacency_a, int* adjacency_b, int length, int *circularNum, int *indices)
{
	//get the target adj node to cut
	int i;
	int target_node[2];
	searchCutNode(target_node, adjacency_a, adjacency_b, length);
	//get the indices
	int findex = getIndex(target_node[0], adjacency_a, length);
	int sindex = getIndex(target_node[1], adjacency_a, length);
	//printf("findex %d, sindex %d\n", findex, sindex);
	// printf("val: %d , val: %d\n", target_node[0], target_node[1]);
	//if findex> sindex swap
	if (findex > sindex)
	{
		int t = findex;
		findex = sindex;
		sindex = t;
		int t_node = target_node[0];
		target_node[0] = target_node[1];
		target_node[1] = t_node;
	}
	int tail_index = getIndex(TAIL, adjacency_a, length);
	// printf("tail: %d\n", tail_index);
	//case 1 between head and tail
	if (findex<=tail_index && sindex<=tail_index)
	{
		if(findex % 2 == 0 && sindex % 2 == 0)
		{
			//printf("318\n");
			for (i=0; i<sindex-findex; i++)
			{
				temp_list[i] = adjacency_a[i+findex+1];
			}
			//reverse
			for (i=0; i<sindex-findex; i++)
			{
				reverse_list[i] = temp_list[sindex-findex-1-i];
			}
			for (i=0; i<sindex-findex; i++)
			{
				adjacency_a[i+findex+1] = reverse_list[i];
			}
		}
		else if(findex % 2 == 0 && sindex % 2 == 1)
		{
			//printf("339\n");
			for (i=0; i<sindex-findex-1; i++)
			{
				temp_list[i] = adjacency_a[i+findex+1];
			}
			for (i=0; i<tail_index - sindex + 1; i++)
			{
				temp_list2[i] = adjacency_a[i+sindex];
			}	
			for (i=0; i<tail_index - sindex + 1; i++)
			{
				adjacency_a[findex+i+1] = temp_list2[i];
			}
			for (i=0; i<sindex-findex-2; i++)
			{
				adjacency_a[tail_index+findex-sindex+2+i] = temp_list[i+1];
			}
			adjacency_a[tail_index]  = temp_list[0];


			//move the right side two cells

			int current_circular = (*circularNum)-1;
			for (i=current_circular; i>=0; i--)
			{
				//move to rigthside two cells
				indices[i*2+2] = indices[i*2];
				indices[i*2+3] = indices[i*2+1];
			}

			/*
			indices[2*wchCircular_a] = f_circle_stt_index;
			indices[2*wchCircular_a+1] = s_circle_end_index - sindex + findex + 1;
			indices[2*wchCircular_a+2] = s_circle_end_index - sindex + findex + 2;
			indices[2*wchCircular_a+3] = s_circle_end_index;

			
			for (i=current_circular*2-1; i>2*wchCircular_a; i=i-2)
			{
				//move to rightside two cell
				indices[i+2] = indices[i];
				indices[i+1] = indices[i-1];
			}
			indices[2*wchCircular_a] = s_circle_end_index+findex-sindex+2;
			indices[2*wchCircular_a+1] = s_circle_end_index;
		
			(*circularNum) ++;

			int current_circular;
			current_circular = (*circularNum);
			for (i=current_circular*2-1; i>0; i=i-2)
			{
				//move to rightside two cell
				indices[i+2] = indices[i];
				indices[i+1] = indices[i-1];
			}
			*/
			indices[0] = tail_index+findex-sindex+2;
			indices[1] = tail_index;
			(*circularNum) ++;
			// printf("++++++++++++++++++\n");
			// for (i=0; i<length; i++)
			// {
			// 	printf("%d ", adjacency_a[i]);
			// }
			// printf("------------------\n");
			// int cn = (*circularNum);
			// indices[cn*2-2] = tail_index+findex-sindex+2;
			// indices[cn*2-1] = tail_index;
		}
		else if(findex % 2 == 1 && sindex % 2 == 1)
		{
			//printf("384\n");
			for (i=0; i<sindex-findex; i++)
			{
				temp_list[i] = adjacency_a[i+findex];
			}
			//reverse
			for (i=0; i<sindex-findex; i++)
			{
				reverse_list[i] = temp_list[sindex-findex-1-i];
			}
			for (i=0; i<sindex-findex; i++)
			{
				adjacency_a[i+findex] = reverse_list[i];
			}
		}
		else if(findex % 2 == 1 && sindex % 2 == 0)
		{
			//printf("405\n");
			for (i=0; i<sindex-findex+1; i++)
			{
				temp_list[i] = adjacency_a[i+findex];
			}
			for (i=0; i<tail_index - sindex; i++)
			{
				temp_list2[i] = adjacency_a[i+sindex+1];
			}	
			for (i=0; i<tail_index - sindex; i++)
			{
				adjacency_a[findex+i] = temp_list2[i];
			}
			for (i=0; i<sindex-findex; i++)
			{
				adjacency_a[tail_index+findex-sindex+i] = temp_list[i+1];
			}
			adjacency_a[tail_index] = temp_list[0];

			//printf("428\n");
			/*
			int current_circular;
			printf("circularNum : %d\n", *circularNum);
			current_circular = (*circularNum);
			for (i=current_circular*2-1; i>0; i=i-2)
			{
				printf("433!!!\n");
				//move to rightside two cell
				indices[i+2] = indices[i];
				indices[i+1] = indices[i-1];
			}*/

			int current_circular = (*circularNum)-1;
			for (i=current_circular; i>=0; i--)
			{
				//move to rigthside two cells
				indices[i*2+2] = indices[i*2];
				indices[i*2+3] = indices[i*2+1];
			}

			indices[0] = tail_index+findex-sindex;
			indices[1] = tail_index;
			(*circularNum) ++;

			// printf("++++++++++++++++++\n");
			// for (i=0; i<length; i++)
			// {
			// 	printf("%d ", adjacency_a[i]);
			// }
			// printf("------------------\n");

			// (*circularNum) ++;
			// int cn = (*circularNum);
			// indices[cn*2-2] = tail_index+findex-sindex;
			// indices[cn*2-1] = tail_index;
		}
	}
	else if(findex<=tail_index && sindex > tail_index)
	{
		int circle_stt_index, circle_end_index;
		int whichCircular;
		//first find which interval does the sindex belongs to?
		for (i=0; i<(*circularNum); i++)
		{
			if (sindex>=indices[i*2] && sindex<=indices[i*2+1])
			{
				whichCircular = i;
				circle_stt_index = indices[i*2];
				circle_end_index = indices[i*2+1];
				break;
			}
		}
		if (findex % 2 == 0 && sindex % 2 == 0)
		{
			//printf("458\n");
			//get the first half
			for (i=0; i<sindex - circle_stt_index + 1; i++)
			{
				temp_list[i] = adjacency_a[sindex-i];
			}
			//get the second half
			for (i=0; i<circle_end_index-sindex; i++)
			{
				temp_list[sindex-circle_stt_index+1+i] = adjacency_a[circle_end_index-i];
			}
			//get tail of ht list
			for (i=0; i<tail_index-findex; i++)
			{
				temp_htlist[i] = adjacency_a[findex+1+i];
			}
			//store the internal values between the tail and circle_stt
			int l = circle_stt_index - tail_index-1;
			if (l>0) 
			{
				for (i=0; i<l; i++)
				{
					internal_templist[i] = adjacency_a[tail_index+i+1];
				}
			}

			//reassign to original list
			for (i=0; i<circle_end_index- circle_stt_index +1; i++)
			{
				adjacency_a[findex+1+i] = temp_list[i];
			}
			for (i=0; i<tail_index-findex; i++)
			{
				adjacency_a[circle_end_index - circle_stt_index + findex + 2 + i] = temp_htlist[i];
			}

			if (l>0)
			{
				for (i=0; i<l; i++)
				{
					adjacency_a[tail_index + circle_end_index - circle_stt_index + 2 + i] = internal_templist[i];
				}
				//update indices of circular before current one

				for (i=0; i<whichCircular; i++)
				{
					indices[i*2] = indices[i*2] + circle_end_index - circle_stt_index + 1;
					indices[i*2+1] = indices[i*2+1] + circle_end_index - circle_stt_index + 1;
				}

				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}

				(*circularNum)--;
			}
			else
			{	
				//for circular at right, move left two cells
				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}
				(*circularNum)--;
			}

			// int cn = (*circularNum);
			// while((indices[cn*2] != 0) && (indices[cn*2+1] != 0))
			// {
			// 	indices[cn*2 - 2] = indices[cn*2];
			// 	indices[cn*2 - 1] = indices[cn*2+1];
			// 	cn++;
			// }
			// indices[cn*2-2] = 0;
			// indices[cn*2-1] = 0;
			// (*circularNum) --;

		}
		else if (findex % 2 == 0 && sindex % 2 == 1)
		{
			//printf("541\n");
			//get the first half
			for (i=0; i<circle_end_index - sindex + 1; i++)
			{
				temp_list[i] = adjacency_a[sindex+i];
			}
			//get the second half
			for (i=0; i<sindex - circle_stt_index; i++)
			{
				temp_list[circle_end_index - sindex+1+i] = adjacency_a[circle_stt_index+i];
			}
			//get tail of ht list
			for (i=0; i<tail_index-findex; i++)
			{
				temp_htlist[i] = adjacency_a[findex+1+i];
			}
			//store the internal values between the tail and circle_stt
			int l = circle_stt_index - tail_index-1;
			if (l>0) 
			{
				for (i=0; i<l; i++)
				{
					internal_templist[i] = adjacency_a[tail_index+i+1];
				}
			}

			//reassign to original list
			for (i=0; i<circle_end_index- circle_stt_index +1; i++)
			{
				adjacency_a[findex+1+i] = temp_list[i];
			}
			for (i=0; i<tail_index-findex; i++)
			{
				adjacency_a[circle_end_index - circle_stt_index + findex + 2 + i] = temp_htlist[i];
			}

			if (l>0)
			{
				for (i=0; i<l; i++)
				{
					adjacency_a[tail_index + circle_end_index - circle_stt_index + 2 + i] = internal_templist[i];
				}
				//update indices of circular before current one
				for (i=0; i<whichCircular; i++)
				{
					indices[i*2] = indices[i*2] + circle_end_index - circle_stt_index + 1;
					indices[i*2+1] = indices[i*2+1] + circle_end_index - circle_stt_index + 1;
				}

				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}

				(*circularNum)--;
			}
			else
			{	
				//for circular at right, move left two cells
				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}
				(*circularNum)--;
			}

			// int cn = (*circularNum);
			// while((indices[cn*2] != 0) && (indices[cn*2+1] != 0))
			// {
			// 	indices[cn*2 - 2] = indices[cn*2];
			// 	indices[cn*2 - 1] = indices[cn*2+1];
			// 	cn++;
			// }
			// indices[cn*2-2] = 0;
			// indices[cn*2-1] = 0;
			// (*circularNum) --;
		}
		else if(findex % 2 == 1 && sindex % 2 == 0)
		{
			//printf("622\n");
			//get the first half
			for (i=0; i<circle_end_index - sindex; i++)
			{
				temp_list[i] = adjacency_a[sindex+i+1];
			}
			//get the second half
			for (i=0; i<sindex - circle_stt_index+1; i++)
			{
				temp_list[circle_end_index - sindex+i] = adjacency_a[circle_stt_index+i];
			}
			//get tail of ht list
			for (i=0; i<tail_index-findex+1; i++)
			{
				temp_htlist[i] = adjacency_a[findex+i];
			}
			//store the internal values between the tail and circle_stt
			int l = circle_stt_index - tail_index-1;
			if (l>0) 
			{
				for (i=0; i<l; i++)
				{
					internal_templist[i] = adjacency_a[tail_index+i+1];
				}
			}
			//reassign to original list
			for (i=0; i<circle_end_index- circle_stt_index +1; i++)
			{
				adjacency_a[findex+i] = temp_list[i];
			}
			for (i=0; i<tail_index-findex+1; i++)
			{
				adjacency_a[circle_end_index - circle_stt_index + findex + 1 + i] = temp_htlist[i];
			}

			if (l>0)
			{
				for (i=0; i<l; i++)
				{
					adjacency_a[tail_index + circle_end_index - circle_stt_index + 2 + i] = internal_templist[i];
				}
				//update indices of circular before current one
				for (i=0; i<whichCircular; i++)
				{
					indices[i*2] = indices[i*2] + circle_end_index - circle_stt_index + 1;
					indices[i*2+1] = indices[i*2+1] + circle_end_index - circle_stt_index + 1;
				}

				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}

				(*circularNum)--;
			}
			else
			{	
				//for circular at right, move left two cells
				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}
				(*circularNum)--;
			}

			// int cn = (*circularNum);
			// while((indices[cn*2] != 0) && (indices[cn*2+1] != 0))
			// {
			// 	indices[cn*2 - 2] = indices[cn*2];
			// 	indices[cn*2 - 1] = indices[cn*2+1];
			// 	cn++;
			// }
			// indices[cn*2-2] = 0;
			// indices[cn*2-1] = 0;
			// (*circularNum) --;
		}
		else if (findex % 2 == 1 && sindex % 2 == 1)
		{
			//printf("703\n");
			//get the first half
			for (i=0; i<sindex - circle_stt_index; i++)
			{
				temp_list[i] = adjacency_a[sindex-i-1];
			}
			//get the second half
			for (i=0; i<circle_end_index-sindex+1; i++)
			{
				temp_list[sindex - circle_stt_index+i] = adjacency_a[circle_end_index-i];
			}
			//get tail of ht list
			for (i=0; i<tail_index-findex+1; i++)
			{
				temp_htlist[i] = adjacency_a[findex+i];
			}
			//store the internal values between the tail and circle_stt
			int l = circle_stt_index - tail_index-1;
			if (l>0)
			{
				for (i=0; i<l; i++)
				{
					internal_templist[i] = adjacency_a[tail_index+i+1];
				}
			}
			//reassign to original list
			for (i=0; i<circle_end_index- circle_stt_index +1; i++)
			{
				adjacency_a[findex+i] = temp_list[i];
			}
			for (i=0; i<tail_index-findex+1; i++)
			{
				adjacency_a[circle_end_index - circle_stt_index + findex + 1 + i] = temp_htlist[i];
			}
			if (l>0)
			{
				for (i=0; i<l; i++)
				{
					adjacency_a[tail_index + circle_end_index - circle_stt_index + 2 + i] = internal_templist[i];
				}
				//update indices of circular before current one
				for (i=0; i<whichCircular; i++)
				{
					indices[i*2] = indices[i*2] + circle_end_index - circle_stt_index + 1;
					indices[i*2+1] = indices[i*2+1] + circle_end_index - circle_stt_index + 1;
				}

				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}

				(*circularNum)--;
			}
			else
			{	
				//for circular at right, move left two cells
				for (i=whichCircular+1; i<(*circularNum); i++)
				{
					indices[i*2-2] = indices[i*2];
					indices[i*2-1] = indices[i*2+1];
				}
				(*circularNum)--;
			}

			// int cn = (*circularNum);
			// while((indices[cn*2] != 0) && (indices[cn*2+1] != 0))
			// {
			// 	indices[cn*2 - 2] = indices[cn*2];
			// 	indices[cn*2 - 1] = indices[cn*2+1];
			// 	cn++;
			// }
			// indices[cn*2-2] = 0;
			// indices[cn*2-1] = 0;
			// (*circularNum)--;
		}
	}
	else if(findex > tail_index && sindex > tail_index)
	{
		int f_circle_stt_index, f_circle_end_index;
		int s_circle_stt_index, s_circle_end_index;
		int wchCircular_a, wchCircular_b;
		//int ii, jj;
		//check if these two node in the same circular
		for (i=0; i<(*circularNum); i++)
		{
			if (findex>=indices[i*2] && findex<=indices[i*2+1])
			{
				f_circle_stt_index = indices[i*2];
				f_circle_end_index = indices[i*2+1];
				wchCircular_a = i;
				break;
			}
		}

		for (i=0; i<(*circularNum); i++)
		{
			if (sindex>=indices[i*2] && sindex<=indices[i*2+1])
			{
				s_circle_stt_index = indices[i*2];
				s_circle_end_index = indices[i*2+1];
				wchCircular_b = i;
				break;
			}
		}

		if (wchCircular_a == wchCircular_b)
		{
			//two nodes at the same circle
			if(findex % 2 == 0 && sindex % 2 == 0)
			{
				//printf("815\n");
				for (i=0; i<sindex-findex; i++)
				{
					temp_list[i] = adjacency_a[i+findex+1];
				}
				//reverse
				for (i=0; i<sindex-findex; i++)
				{
					reverse_list[i] = temp_list[sindex-findex-1-i];
				}
				for (i=0; i<sindex-findex; i++)
				{
					adjacency_a[i+findex+1] = reverse_list[i];
				}
			}
			else if(findex % 2 == 0 && sindex % 2 == 1)
			{
				//printf("836\n");
				for (i=0; i<sindex-findex-1; i++)
				{
					temp_list[i] = adjacency_a[i+findex+1];
				}
				for (i=0; i<s_circle_end_index - sindex + 1; i++)
				{
					temp_list2[i] = adjacency_a[i+sindex];
				}	
				for (i=0; i<s_circle_end_index - sindex + 1; i++)
				{
					adjacency_a[findex+i+1] = temp_list2[i];
				}
				for (i=0; i<sindex-findex-2; i++)
				{
					adjacency_a[s_circle_end_index+findex-sindex+2+i] = temp_list[i+1];
				}
				adjacency_a[s_circle_end_index] = temp_list[0]; 
				//update the original

				//move the right side two cells

				int current_circular = (*circularNum)-1;
				for (i=current_circular; i>wchCircular_a; i--)
				{
					//move to rigthside two cells
					indices[i*2+2] = indices[i*2];
					indices[i*2+3] = indices[i*2+1];
				}

				indices[2*wchCircular_a] = f_circle_stt_index;
				indices[2*wchCircular_a+1] = s_circle_end_index - sindex + findex + 1;
				indices[2*wchCircular_a+2] = s_circle_end_index - sindex + findex + 2;
				indices[2*wchCircular_a+3] = s_circle_end_index;

				/*
				for (i=current_circular*2-1; i>2*wchCircular_a; i=i-2)
				{
					//move to rightside two cell
					indices[i+2] = indices[i];
					indices[i+1] = indices[i-1];
				}
				indices[2*wchCircular_a] = s_circle_end_index+findex-sindex+2;
				indices[2*wchCircular_a+1] = s_circle_end_index;
				*/
				(*circularNum) ++;

				// indices[ii] = s_circle_end_index - sindex + findex + 1;

				// (*circularNum) ++;
				// int cn = (*circularNum);
				// indices[cn*2-2] = s_circle_end_index+findex-sindex+2;
				// indices[cn*2-1] = s_circle_end_index;
				
			}
			else if(findex % 2 == 1 && sindex % 2 == 0)
			{
				//printf("882\n");
				for (i=0; i<sindex-findex+1; i++)
				{
					temp_list[i] = adjacency_a[i+findex];
				}
				for (i=0; i<s_circle_end_index - sindex; i++)
				{
					temp_list2[i] = adjacency_a[i+sindex+1];
				}	
				for (i=0; i<s_circle_end_index - sindex; i++)
				{
					adjacency_a[findex+i] = temp_list2[i];
				}
				for (i=0; i<sindex-findex; i++)
				{
					adjacency_a[s_circle_end_index+findex-sindex+i] = temp_list[i+1];
				}
				adjacency_a[s_circle_end_index] = temp_list[0];


				int current_circular = (*circularNum)-1;
				for (i=current_circular; i>wchCircular_a; i--)
				{
					//move to rigthside two cells
					indices[i*2+2] = indices[i*2];
					indices[i*2+3] = indices[i*2+1];
				}

				indices[2*wchCircular_a] = f_circle_stt_index;
				indices[2*wchCircular_a+1] = s_circle_end_index - sindex + findex - 1;
				indices[2*wchCircular_a+2] = s_circle_end_index - sindex + findex;
				indices[2*wchCircular_a+3] = s_circle_end_index;

				/*
				//move the right side two cells
				int current_circular = (*circularNum);
				for (i=current_circular*2-1; i>2*wchCircular_a; i=i-2)
				{
					//move to rightside two cell
					indices[i+2] = indices[i];
					indices[i+1] = indices[i-1];
				}
				indices[2*wchCircular_a] = s_circle_end_index+findex-sindex;
				indices[2*wchCircular_a+1] = s_circle_end_index;
				*/
				(*circularNum) ++;

				// indices[ii] = s_circle_end_index - sindex - 1 +findex;

				// (*circularNum) ++;
				// int cn = (*circularNum);
				// indices[cn*2-2] = s_circle_end_index+findex-sindex;
				// indices[cn*2-1] = s_circle_end_index;
			}
			else if(findex % 2 == 1 && sindex % 2 == 1)
			{
				//printf("926\n");
				for (i=0; i<sindex-findex; i++)
				{
					temp_list[i] = adjacency_a[i+findex];
				}
				//reverse
				for (i=0; i<sindex-findex; i++)
				{
					reverse_list[i] = temp_list[sindex-findex-1-i];
				}
				for (i=0; i<sindex-findex; i++)
				{
					adjacency_a[i+findex] = reverse_list[i];
				}
			}
		}
		else 
		{
			if(findex % 2 == 0 && sindex % 2 == 0)
			{
				//printf("950\n");
				//get the first half
				for (i=0; i<sindex - s_circle_stt_index + 1; i++)
				{
					temp_list[i] = adjacency_a[sindex-i];
				}
				//get the second half
				for (i=0; i<s_circle_end_index-sindex; i++)
				{
					temp_list[sindex-s_circle_stt_index+1+i] = adjacency_a[s_circle_end_index-i];
				}
				//get tail of ht list
				for (i=0; i<f_circle_end_index-findex; i++)
				{
					temp_htlist[i] = adjacency_a[findex+i+1];
				}
				//store the internal values between the tail and circle_stt
				int l = s_circle_stt_index - f_circle_end_index-1;
				if (l>0) 
				{
					for (i=0; i<l; i++)
					{
						internal_templist[i] = adjacency_a[f_circle_end_index+i+1];
					}
				}
				//assign to original list
				for (i=0; i<s_circle_end_index - s_circle_stt_index + 1; i++)
				{
					adjacency_a[findex+1+i] = temp_list[i];
				}
				for (i=0; i<f_circle_end_index-findex; i++)
				{
					adjacency_a[s_circle_end_index - s_circle_stt_index+findex+2+i] = temp_htlist[i];
				}

				if (l>0)
				{
					for (i=0; i<l; i++)
					{
						adjacency_a[s_circle_end_index - s_circle_stt_index+2 +f_circle_end_index + i] = internal_templist[i];
					}
					//update indices of circular before current one
					// i = 0;
					// while(i<wchCircular_a)
					// {
					// 	indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// 	indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// 	i++;
					// }
					// //for wchCircular_a
					// indices[wchCircular_a*2] = indices[wchCircular_a*2];
					// indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index+findex + 1;
					// //for circular at right, move left two cells
					// for (i=wchCircular_a*2+2; i<(*circularNum)*2; i=i+2)
					// {
					// 	indices[i-2] = indices[i];
					// 	indices[i-1] = indices[i+1];
					// }

					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index + 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					for (i=wchCircular_a+1; i<wchCircular_b; i++)
					{
						indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index + 1;
						indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index + 1;
					}

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}

					(*circularNum)--;
				}
				else
				{	
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index+ 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					// for (i=wchCircular_a+1; i<wchcircular_b; i++)
					// {
					// 	indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// 	indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// }

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}

				//indices[ii] = s_circle_end_index - s_circle_stt_index+1 + f_circle_end_index;

				//refresh the indices
			}
			else if (findex % 2 == 0 && sindex % 2 == 1)
			{
				//printf("1030\n");
				//get the first half
				for (i=0; i<s_circle_end_index - sindex + 1; i++)
				{
					temp_list[i] = adjacency_a[sindex+i];
				}
				//get the second half
				for (i=0; i<sindex -s_circle_stt_index ; i++)
				{
					temp_list[s_circle_end_index - sindex+1+i] = adjacency_a[s_circle_stt_index+i];
				}

				//get tail of ht list
				for (i=0; i<f_circle_end_index-findex; i++)
				{
					temp_htlist[i] = adjacency_a[findex+i+1];
				}
				//store the internal values between the tail and circle_stt
				int l = s_circle_stt_index - f_circle_end_index-1;
				if (l>0)
				{
					for (i=0; i<l; i++)
					{
						internal_templist[i] = adjacency_a[f_circle_end_index+i+1];
						
					}
				}
				//assign to original list
				for (i=0; i<s_circle_end_index - s_circle_stt_index + 1; i++)
				{
					adjacency_a[findex+1+i] = temp_list[i];
				}
				for (i=0; i < f_circle_end_index-findex; i++)
				{
					adjacency_a[s_circle_end_index - s_circle_stt_index + findex +2 + i] = temp_htlist[i]; 
				}

				//refresh the indices
				if (l>0)
				{
					for (i=0; i<l; i++)
					{
						adjacency_a[s_circle_end_index - s_circle_stt_index +2 +f_circle_end_index + i] = internal_templist[i];
					}
					//update indices of circular before current one
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index+ 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					for (i=wchCircular_a+1; i<wchCircular_b; i++)
					{
						indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index+ 1;
						indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+ 1;
					}

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}
				else
				{	
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index + 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					// for (i=wchCircular_a+1; i<wchcircular_b; i++)
					// {
					// 	indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// 	indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// }

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}
			}
			else if (findex % 2 == 1 && sindex % 2 == 0)
			{
				//printf("1108\n");
				//get the first half
				for (i=0; i<s_circle_end_index - sindex ; i++)
				{
					temp_list[i] = adjacency_a[sindex+i+1];
				}
				//get the second half
				for (i=0; i<sindex -s_circle_stt_index+1 ; i++)
				{
					temp_list[s_circle_end_index - sindex+i] = adjacency_a[s_circle_stt_index+i];
				}
				//get tail of ht list
				for (i=0; i<f_circle_end_index-findex+1; i++)
				{
					temp_htlist[i] = adjacency_a[findex+i];
				}
				//store the internal values between the tail and circle_stt
				int l = s_circle_stt_index - f_circle_end_index-1;
				if (l>0) 
				{	
					for (i=0; i<l; i++)
					{
						internal_templist[i] = adjacency_a[f_circle_end_index+i+1];
						
					}
				}
				//reassign to original list
				for (i=0; i<s_circle_end_index - s_circle_stt_index + 1; i++)
				{
					adjacency_a[findex+i] = temp_list[i];
				}
				for (i=0; i < f_circle_end_index-findex+1; i++)
				{
					adjacency_a[s_circle_end_index - s_circle_stt_index+findex+1+i] = temp_htlist[i];
				}

				//refresh the indices
				if (l>0)
				{
					for (i=0; i<l; i++)
					{
						adjacency_a[s_circle_end_index - s_circle_stt_index+ f_circle_end_index+2 + i] = internal_templist[i];
					}
					//update indices of circular before current one
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index + 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					for (i=wchCircular_a+1; i<wchCircular_b; i++)
					{
						indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index + 1;
						indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index + 1;
					}

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}
				else
				{	
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index + 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					// for (i=wchCircular_a+1; i<wchcircular_b; i++)
					// {
					// 	indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// 	indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// }

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}
			}
			else if (findex % 2 == 1 && sindex % 2 == 1)
			{
				//printf("1185\n");
				//get the first half
				for (i=0; i<sindex - s_circle_stt_index; i++)
				{
					temp_list[i] = adjacency_a[sindex-i-1];
				}
				//get the second half
				for (i=0; i<s_circle_end_index-sindex+1; i++)
				{
					temp_list[sindex-s_circle_stt_index+i] = adjacency_a[s_circle_end_index-i];
				}
				//get tail of ht list
				for (i=0; i<f_circle_end_index-findex+1; i++)
				{
					temp_htlist[i] = adjacency_a[findex+i];
				}
				//store the internal values between the tail and circle_stt
				int l = s_circle_stt_index - f_circle_end_index-1;
				if (l>0) 
				{
					for (i=0; i<l; i++)
					{
						internal_templist[i] = adjacency_a[f_circle_end_index+i+1];
						
					}
				}
				//assign to original list
				for (i=0; i<s_circle_end_index - s_circle_stt_index + 1; i++)
				{
					adjacency_a[findex+i] = temp_list[i];
				}
				for (i=0; i<f_circle_end_index-findex+1; i++)
				{
					adjacency_a[s_circle_end_index - s_circle_stt_index+findex+1+i] = temp_htlist[i];
				}

				if (l>0)
				{
					for (i=0; i<l; i++)
					{
						adjacency_a[s_circle_end_index - s_circle_stt_index+ f_circle_end_index+2 + i] = internal_templist[i];
					}
					//update indices of circular before current one
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index + 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					for (i=wchCircular_a+1; i<wchCircular_b; i++)
					{
						indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index + 1;
						indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+ 1;
					}

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}
				else
				{	
					//for wchCircular_a
					indices[wchCircular_a*2] = indices[wchCircular_a*2];
					indices[wchCircular_a*2+1] = indices[wchCircular_a*2+1] +  s_circle_end_index - s_circle_stt_index + 1;
					//keep all indices left wchCircular_a, and increase same index value between circular_a and circular_b
					//for all right wchcircular_b, move left two cells

					// for (i=wchCircular_a+1; i<wchcircular_b; i++)
					// {
					// 	indices[i*2] = indices[i*2] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// 	indices[i*2+1] = indices[i*2+1] + s_circle_end_index - s_circle_stt_index+findex + 1;
					// }

					for (i=wchCircular_b+1; i<(*circularNum); i++)
					{
						indices[i*2-2] = indices[i*2];
						indices[i*2-1] = indices[i*2+1];
					}
					(*circularNum)--;
				}
			}
		}
	}
}

/*genetic algorithm*/

void generateInitialPool(int** original_list, int length)
{
	int i,j,k,l, distance, steps, n, counter, global_circularNum;
	counter = 0;
	//sorting a to b
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			if (i != j)
			{
				int* temp_original_a = original_list[i];
				int* temp_original_b = original_list[j];
				//calculate the DCJ distance between a and b
				distance = calculateDCJdistance(temp_original_a, temp_original_b, length);
				//generate 1/10 2/10 3/10 4/10 5/10 6/10 distance
				int d[6] = {0};
				for(k=0; k<6; k++)
				{
					d[k] = (int) (distance * (k+1) * 0.1);
				}
				for (k=0; k<6; k++)
				{
					steps = d[k];
					if (steps > 0)
					{
						n = 50;
						while (n > 0)
						{
							steps = d[k];
							global_circularNum = 0;
							//deep copy original a & original b
							for (l=0; l<length; l++)
							{
								original_copy_a[l] = original_list[i][l];
								original_copy_b[l] = original_list[j][l];
							}
							while(steps > 0)
							{
								DCJOperation(original_copy_a, original_copy_b, length, &global_circularNum, global_indices);
								steps--;
							}
							//copy into the initial pool
							for (l=0; l<length; l++)
							{
								genome_pool[counter][l] = original_copy_a[l];
								g_pool_indices[counter][l] = global_indices[l];
							}
							g_pool_cnum[counter] = global_circularNum;
							n--;
							counter++;
						}
					}
				}
			}
		}
	}
}
// Accroding to the plosone
void geneticAlgorithm(int** original_list, int length)
{
	int i,best_score,score,low, high,val;
	best_score = calculateDCJdistance(original_list[0], original_list[1], length)
	+ calculateDCJdistance(original_list[0], original_list[2], length) 
	+ calculateDCJdistance(original_list[1], original_list[2], length);
	best_score = best_score / 2;
	//compute the fitness score
	for (i=0; i<1800; i++)
	{
		int* temp_genome = genome_pool[i];
		score = 0;
		score = calculateDCJdistance(temp_genome, original_list[0], length)
		+ calculateDCJdistance(temp_genome, original_list[1], length)
		+ calculateDCJdistance(temp_genome, original_list[2], length);
		score = GENENUM - (score - best_score);
		fitnessScore[i] = score;
	}

	low = 0;
	high = 1799; // the index of last genome in initial pool.
	quickSort(low, high);
	// for (i=0; i<1800; i++)
	// {
	// 	val = 
	// 	calculateDCJdistance(genome_pool[i], original_list[0], length)+
	// 	calculateDCJdistance(genome_pool[i], original_list[1], length)+
	// 	calculateDCJdistance(genome_pool[i], original_list[2], length);
		
	// 	printf ("val: %d\n", GENENUM - val + best_score);
	// }

	// for (i=0; i<1800; i++)
	// {
	// 	printf ("fitnessScore: %d\n", fitnessScore[i]);
	// }
	for (i=0; i<1620; i++)
	{
		index_array[i] = 180 + i;
	}

	//(1)selection:randomly pick two genomes as parents
	shuffle(index_array, 1620);

	for (i=0; i<1620; i=i+2)
	{
		geneticAlgorithmForParents(index_array[i], index_array[i+1], original_list, length);
	}
}


void geneticAlgorithmForParents(int parent1_index, int parent2_index, int** original_list, int length)
{
	int i, j, k, t_val, p1_cn, p2_cn, p1_fScore, p2_fScore, dist, steps, p1_cnCopy, p2_cnCopy;
	int d_atom, d_p1toa, d_p2toa, d_btom, d_p1tob, d_p2tob, d_ctom, d_p1toc, d_p2toc;
	int largest_index, difference, largest_j, differ, best_score, score, temp, best_first, best_second;
	//copy value from initialpool to parents
	p1_cn = g_pool_cnum[parent1_index];
	p2_cn = g_pool_cnum[parent2_index];
	for (i=0; i<length; i++)
	{
		parent1[i] = genome_pool[parent1_index][i];
		parent2[i] = genome_pool[parent2_index][i];
		parent1Indices[i] = g_pool_indices[parent1_index][i];
		parent2Indices[i] = g_pool_indices[parent2_index][i]; 
	}

	//(2)cross over:according to the better one:from worse to the better
	//DCJOperation(temp_mediangenome, given_genomes)
	p1_fScore = fitnessScore[parent1_index];
	p2_fScore = fitnessScore[parent2_index];

	dist = calculateDCJdistance(parent1, parent2, length);
	if (p1_fScore > p2_fScore)
	{
		// generate C1
		steps = 0;
		if (dist > 0) { steps = rand() % dist + 1;}
		while (steps > 0)
		{
			DCJOperation(parent2, parent1, length, &p2_cn, parent2Indices);
			steps--;
		}
	}
	else 
	{
		// generate C1
		steps = 0;
		if (dist > 0) { steps = rand() % dist + 1;}
		while (steps > 0)
		{
			DCJOperation(parent1, parent2, length, &p1_cn, parent1Indices);
			steps--;	
		}
	}

	//mutation
	p1_cnCopy = g_pool_cnum[parent1_index];
	p2_cnCopy = g_pool_cnum[parent2_index];
	for (i=0; i<length; i++)
	{
		parent1Copy[i] = genome_pool[parent1_index][i];
		parent2Copy[i] = genome_pool[parent2_index][i];
		parent1CopyIndices[i] = g_pool_indices[parent1_index][i];
		parent2CopyIndices[i] = g_pool_indices[parent2_index][i];
	}	
	int p1toLeave[3] = {0};
	int p2toLeave[3] = {0};

	//lower bound d1m = (d12+d13-d23)/2
	d_atom = (calculateDCJdistance(original_list[0], original_list[1], length)+calculateDCJdistance(original_list[0], original_list[2], length)-calculateDCJdistance(original_list[1], original_list[2], length))/2;
	d_p1toa = calculateDCJdistance(parent1Copy, original_list[0], length);
	d_p2toa = calculateDCJdistance(parent2Copy, original_list[0], length);
	p1toLeave[0] = abs(d_p1toa - d_atom);
	p2toLeave[0] = abs(d_p2toa - d_atom);

	//lower bound d2m = (d12+d23-d13)/2
	d_btom = (calculateDCJdistance(original_list[0], original_list[1], length)+calculateDCJdistance(original_list[1], original_list[2], length)-calculateDCJdistance(original_list[0], original_list[2], length))/2;
	d_p1tob = calculateDCJdistance(parent1Copy, original_list[1], length);
	d_p2tob = calculateDCJdistance(parent2Copy, original_list[1], length);
	p1toLeave[1] = abs(d_p1tob - d_btom);
	p2toLeave[1] = abs(d_p2tob - d_btom);

	//lower bound d3m = (d13+d23-d12)/2
	d_ctom = (calculateDCJdistance(original_list[0], original_list[2], length)+calculateDCJdistance(original_list[1], original_list[2], length)-calculateDCJdistance(original_list[0], original_list[1], length))/2;
	d_p1toc = calculateDCJdistance(parent1Copy, original_list[2], length);
	d_p2toc = calculateDCJdistance(parent2Copy, original_list[2], length);
	p1toLeave[2] = abs(d_p1toc - d_ctom);
	p2toLeave[2] = abs(d_p2toc - d_ctom);

	largest_index = 0;
	difference = 0;
	for (i=0; i<3; i++)
	{
		if (p1toLeave[i] > difference)
		{
			difference = p1toLeave[i];
			largest_index = i;
		}
	}

	largest_j = 0;
	differ = 0;
	for (i=0; i<3; i++)
	{
		if (p2toLeave[i] > differ)
		{
			differ = p2toLeave[i];
			largest_j = i;
		}
	}
    // DCJ sortinh m steps away from parent1Copy to original_list[largest_index]:最大差异
	dist = calculateDCJdistance(parent1Copy, original_list[largest_index], length);
	steps = 0;
	if (dist > 0) {steps = rand() % dist + 1;}
	while (steps > 0)
	{
		DCJOperation(parent1Copy, original_list[largest_index], length, &p1_cnCopy, parent1CopyIndices);
		steps--;
	}

	dist = calculateDCJdistance(parent2Copy, original_list[largest_j], length);
	steps = 0;
	if (dist > 0) {steps = rand() % dist + 1;}
	while (steps > 0)
	{
		DCJOperation(parent2Copy, original_list[largest_j], length, &p2_cnCopy, parent2CopyIndices);
		steps--;
	}

	//add all four candidates into a array
	int cnarray[4] = {0};

	for (i=0; i<length; i++)
	{
		candidateGenomes[0][i] = parent1[i];
		candidateGenomes[1][i] = parent2[i];
		candidateGenomes[2][i] = parent1Copy[i];
		candidateGenomes[3][i] = parent2Copy[i];	

		
		candidateIndices[0][i] = parent1Indices[i];
		candidateIndices[1][i] = parent2Indices[i];
		candidateIndices[2][i] = parent1CopyIndices[i];
		candidateIndices[3][i] = parent2CopyIndices[i];

	}
	


	cnarray[0] = p1_cn;
	cnarray[1] = p2_cn;
	cnarray[2] = p1_cnCopy;
	cnarray[3] = p2_cnCopy;

	best_score = calculateDCJdistance(original_list[0], original_list[1], length) + calculateDCJdistance(original_list[0], original_list[2], length) + calculateDCJdistance(original_list[1], original_list[2], length);
	best_score /= 2;

	int candidatesFitScore[4] = {0};
	score = 0;
	for (i=0; i<4; i++)
	{
		score = calculateDCJdistance(candidateGenomes[i], original_list[0], length) + calculateDCJdistance(candidateGenomes[i], original_list[1], length) + calculateDCJdistance(candidateGenomes[i], original_list[2], length);
		score = GENENUM - (score - best_score);
		candidatesFitScore[i] = score; 
	}

	// int canFitScore[4] = {0};
	// //store the original 
	// for (i=0; i<4; i++)
	// {
	// 	canFitScore[i] = candidatesFitScore[i];
	// }

	//sort the fitness score
	for (i=0; i<4; i++)
	{
		for (j=i+1; j<4; j++)
		{
			if (candidatesFitScore[j] > candidatesFitScore[i])
			{
				//swap
				temp = candidatesFitScore[i];
				candidatesFitScore[i] = candidatesFitScore[j];
				candidatesFitScore[j] = temp;

				for (k=0; k<length; k++)
				{
					temporary_genome[k] = candidateGenomes[i][k];
					temporary_indices[k] = candidateIndices[i][k];
				}

				t_val = cnarray[i];
				cnarray[i] = cnarray[j];
				cnarray[j] = t_val;

				for (k=0; k<length; k++)
				{
					candidateGenomes[i][k] = candidateGenomes[j][k];
					candidateIndices[i][k] = candidateIndices[j][k];
				}

				for (k=0; k<length; k++)
				{
					candidateGenomes[j][k] = temporary_genome[k];
					candidateIndices[j][k] = temporary_indices[k];
				}
			}
		}
	}


	//copy these best two into the pool
	g_pool_cnum[parent1_index] = cnarray[0];
	g_pool_cnum[parent2_index] = cnarray[1];
	for(i=0; i<length; i++)
	{
		genome_pool[parent1_index][i] = candidateGenomes[0][i];
		g_pool_indices[parent1_index][i] = candidateIndices[0][i];
		genome_pool[parent2_index][i] = candidateGenomes[1][i];
		g_pool_indices[parent2_index][i] = candidateIndices[1][i];
	}
}

void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

void quickSort(int low, int high)
{
	int middle, pivot, i, j, fitscore_temp,k;
	if (low >= high)
		return;

	//pick the pivot
	middle = low + (high - low) / 2;
	pivot = fitnessScore[middle];
	i = low, j = high;
	while (i <= j)
	{
		while (fitnessScore[i] > pivot)
		{
			i++;
		}

		while (fitnessScore[j] < pivot)
		{
			j--;
		}

		if (i <= j)
		{
			//swap 
			fitscore_temp = fitnessScore[i];
			fitnessScore[i] = fitnessScore[j];
			fitnessScore[j] = fitscore_temp;

			for (k=0; k<GENOMELENGTH; k++)
			{
				temporary_genome[k] = genome_pool[i][k];
				temporary_indices[k] = g_pool_indices[i][k];
			}
			
			int tval = g_pool_cnum[i];

			for (k=0; k<GENOMELENGTH; k++)
			{
				genome_pool[i][k] = genome_pool[j][k];
				g_pool_indices[i][k] = g_pool_indices[j][k];
			}

			g_pool_cnum[i] = g_pool_cnum[j];

			for (k=0; k<GENOMELENGTH; k++)
			{
				genome_pool[j][k] = temporary_genome[k];
				g_pool_indices[j][k] = temporary_indices[k];
			}

			g_pool_cnum[j] = tval;

			i++;
			j--;
		}
	}

	if (low < j)
	{
		quickSort(low, j);
	}

	if (high > i)
	{
		quickSort(i, high);
	}
}

int calculateMedianScore(int* median_genome, int** original_list, int length)
{
	int i;
	int distance = 0;
	for (i=0; i<3; i++)
	{
		distance += calculateDCJdistance(median_genome, original_list[i], length);
	}
	return distance;
}
/*
void simulated_annealing(int** original_list, int* median_genome, int length, int* bestMScore, int* bestGenomeNum)
{
	int i,r,iniTrueGeNum, temperature, updateGeNum, formerState, currentState, bestMedianScore, M, index, currTrueGeNum;
	int optindex, N, bestGenNum;
	double coolingRate;
	createInitialTrueGenome(original_list, length);
	//generate random value from 0-6;
	r = rand() % 6;	

	//copy to initialTrueGenome includes the circularnum and indices
	for (i=0; i<length; i++)
	{
		initialTrueGenome[i] = candidateInitialTrue[r][i];
	}
	iniTrueGeNum = initilcandid_circular_num[r];
	for (i=0; i< iniTrueGeNum; i++)
	{
		initialTrueGenomeIndices[i*2] = initilcandid_indices[r][i*2];
		initialTrueGenomeIndices[i*2+1] = initilcandid_indices[r][i*2+1];
	}

	//simulate annealing
	temperature = 10;

	coolingRate = 0.9;

	//store the initialTrueGenome to updateTrueGenome
	for (i=0; i<length; i++)
	{
		updateTrueGenome[i] = initialTrueGenome[i];
	}
	updateGeNum = iniTrueGeNum;
	for (i=0; i<updateGeNum; i++)
	{
		updateTrueGenomeIndices[i*2] = initialTrueGenomeIndices[i*2];
		updateTrueGenomeIndices[i*2+1] = initialTrueGenomeIndices[i*2+1];
	}

	formerState = 0;
	currentState = 3*length;

	formerState = calculateMedianScore(initialTrueGenome, original_list, length);
	//printf("initialValue: %d", formerState);
	bestMedianScore = 3*length;

	M = GENERATION; //change the generation at two positions
	N = M;
	int *converge = (int *)malloc(M*sizeof(int));
	while(M>1)
	{
		index = -1;
		createCandiateMedianGenome(updateTrueGenome, updateGeNum, updateTrueGenomeIndices, original_list, length);
		index = rand() % 3;	
		currentState = calculateMedianScore(candidateMedianGenme[index], original_list, length);
		//copy to currentTrueGenome
		for (i=0; i<length; i++)
		{
			currentTrueGenome[i] = candidateMedianGenme[index][i];
		}
		currTrueGeNum = candidMedian_circular_num[index];
		for(i=0; i<currTrueGeNum; i++)
		{
			currentTrueGenomeIndices[i*2] = candidMedian_indices[index][i*2];
			currentTrueGenomeIndices[i*2+1] = candidMedian_indices[index][i*2+1];
		}
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
				bestGenNum = currTrueGeNum;
				for (i=0; i<bestGenNum; i++)
				{
					bestGenomeIndices[i*2] = currentTrueGenomeIndices[i*2];
					bestGenomeIndices[i*2+1] = currentTrueGenomeIndices[i*2+1];
				}
				bestMedianScore = currentState;
			}
			//pass value from currentTrueGenome to updateTrueGenome
			for (i=0; i<length; i++)
			{
				updateTrueGenome[i] = currentTrueGenome[i];
			}
			updateGeNum = currTrueGeNum;
			for (i=0; i<updateGeNum; i++)
			{
				updateTrueGenomeIndices[i*2] = currentTrueGenomeIndices[i*2];
				updateTrueGenomeIndices[i*2+1] = currentTrueGenomeIndices[i*2+1];
			}
		}
		converge[GENERATION-M] = currentState;
		temperature = temperature*coolingRate;
		M--;
	}

	for (i=0; i<N; i++)
	{
		if (converge[i] == bestMedianScore)
		{
			optindex = i;
			break;
		}
	}
	*bestGenomeNum = bestGenNum;
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
	free(converge);
}
*/

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

void createCandiateMedianGenome(int* previousMedianGenome, int circular_num, int* indices, int** original_list, int length)
{
	int i, k, circularNum;
	for (i=0; i<3; i++)
	{
		//copy previousMedianGenome to a temporary genome includes circular num and indices
		for (k=0; k<length; k++)
		{
			temp_initial_genome[k] = previousMedianGenome[k];
		}
		circularNum = circular_num;
		for (k=0; k<circularNum; k++)
		{
			temp_indices[k*2] = indices[k*2];
			temp_indices[k*2+1] = indices[k*2+1];
		}
		int distance = calculateDCJdistance(temp_initial_genome, original_list[i], length);
		if (distance > 0)
		{
			DCJOperation(temp_initial_genome, original_list[i], length, &circularNum, temp_indices);
		}
		//store the new genome to candidatemedianGenome
		for (k=0; k<length; k++)
		{
			candidateMedianGenme[i][k] = temp_initial_genome[k];
		}
		candidMedian_circular_num[i] = circularNum;
		for (k=0; k<circularNum; k++)
		{
			candidMedian_indices[i][k*2] = temp_indices[k*2];
			candidMedian_indices[i][k*2+1] = temp_indices[k*2+1];
		}
	}
}

void createInitialTrueGenome(int** original_list, int length)
{
	int i,j,k,l,circularNum;
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
				circularNum = 0;
				//reset the temp_indices to zero
				for (k=0; k<length; k++)
				{
					temp_indices[k] = 0;
				} 
				int distance =  calculateDCJdistance(temp_initial_genome, original_list[j], length);
				distance = distance / 10;
				while(distance>0)
				{
					DCJOperation(temp_initial_genome, original_list[j], length, &circularNum, temp_indices);
					distance--;
				}
				//copy the temp_initial_genome to 
				for (k=0; k<length; k++)
				{
					candidateInitialTrue[l][k] = temp_initial_genome[k];
				}
				initilcandid_circular_num[l] = circularNum;
				for (k=0; k<circularNum; k++)
				{
					initilcandid_indices[l][k*2] = temp_indices[k*2] ;
					initilcandid_indices[l][k*2+1] = temp_indices[k*2+1] ;
				}
				l++;
			}
		}
	}
}
/*
double acceptanceProbability(int energy, int newEnergy, double temperature)
{
	if (newEnergy < energy)
	{
		return 1;
	}
	return exp((energy - newEnergy)/temperature);
}
*/

void readDataFromFile(int *arr, char *fileName)
{
	char buffer[6000];
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

int main(int argc, char *argv[])
{
	clock_t begin, end;
	double time_spent;

	begin = clock();

	srand(time(NULL));
	int i, bestMScore, optMScore, j, bestGNum, optMGenomeNum,tail_index,t, generation, length,tempMScore;
	int k=0,temp_optindex=0;
	int genen = GENENUM; /* gene number*/
	int* genome = (int *)malloc(genen*sizeof(int));

	/*read data from files*/
	int *readingenes = (int *)malloc(genen*3*sizeof(int));
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

	length = genen*2+2;
	int* adjacency_a = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(g, genen, adjacency_a);

	int* adjacency_b = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(g2, genen, adjacency_b);

	int* adjacency_c = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(g3, genen, adjacency_c);

	//printf("dist: %d\n", calculateDCJdistance(adjacency_a, adjacency_b, genen*2+2));

	//print out the prev list
	// for (i=0; i<length; i++)
	// {
	// 	printf("%d ", adjacency_a[i]);
	// }

	copy_adj_a = (int *)malloc(length*sizeof(int));
	copy_adj_b = (int *)malloc(length*sizeof(int));
	temp_list = (int*) malloc(length*sizeof(int));
	temp_list2 = (int*) malloc(length*sizeof(int));
	reverse_list = (int*) malloc(length*sizeof(int));
	internal_templist = (int *)malloc(length*sizeof(int));
	temp_htlist = (int *)malloc(length*sizeof(int));
	temp_initial_genome = (int *)malloc(length*sizeof(int));
	temp_indices = (int*)calloc(length, sizeof(int));
	initilcandid_indices = (int **)malloc(6 * sizeof(int *));
	candidateInitialTrue = (int **)malloc(6 * sizeof(int *));
	initialTrueGenome = (int *)malloc(length*sizeof(int));
	initilcandid_circular_num = (int *)malloc(6*sizeof(int));
	initialTrueGenomeIndices = (int *)malloc(length * sizeof(int));
	updateTrueGenome = (int *)malloc(length * sizeof(int));
	candidMedian_circular_num = (int *)malloc(3*sizeof(int));
	candidMedian_indices = (int **)malloc(3*sizeof(int*));
	updateTrueGenomeIndices = (int *)malloc(length * sizeof(int));
	currentTrueGenome = (int *)malloc(length* sizeof(int));
	currentTrueGenomeIndices = (int *)malloc(length*sizeof(int));
	bestGenome = (int *)malloc(length* sizeof(int));
	bestGenomeIndices =(int *)malloc(length* sizeof(int));

	int** original_list = (int **)malloc(3*sizeof(int*));
	int*  median_genome = (int *)malloc(genen*sizeof(int));

	for (i=0; i<6; i++)
    {
    	initilcandid_indices[i] = (int *)malloc(length * sizeof(int));
    	candidateInitialTrue[i] = (int *)malloc(length * sizeof(int));
    }
    candidateMedianGenme = (int **)malloc( 3* sizeof(int *));
    for (i=0; i<3; i++)
    {
    	original_list[i] = (int *)malloc(length * sizeof(int));
    	candidMedian_indices[i] = (int *)malloc(length * sizeof(int));
    	candidateMedianGenme[i] = (int *)malloc(length * sizeof(int));
    }
	for (i=0; i<length; i++)
	{
		original_list[0][i] = adjacency_a[i];
		original_list[1][i] = adjacency_b[i];
		original_list[2][i] = adjacency_c[i];
	}

	int *bestMedian = (int *)malloc(length*sizeof(int));
	int *optIndices = (int *)malloc(length*sizeof(int));


	/*genetic algorithm*/
	genome_pool = (int **)malloc(1800*sizeof(int*));
	for (i=0; i<1800; i++)
	{
		genome_pool[i] = (int *)malloc(length * sizeof(int));
	}

	g_pool_indices = (int **)malloc(1800*sizeof(int *));
	for (i=0; i<1800; i++)
	{
		g_pool_indices[i] = (int *)malloc(length * sizeof(int));
	}
	g_pool_cnum = (int *)malloc(1800 * sizeof(int));

	global_indices = (int *)malloc(length * sizeof(int));
	// global_circularNum = (int *)malloc(sizeof(int));

	original_copy_a = (int *) malloc(length * sizeof(int));
	original_copy_b = (int *) malloc(length * sizeof(int));

	fitnessScore = (int *) malloc(1800 * sizeof(int));

	temporary_genome = (int *) malloc(length * sizeof(int));
	temporary_indices = (int *) malloc(length * sizeof(int));

	best_genome = (int *) malloc(length * sizeof(int));
	best_genomeIndices = (int *)malloc(length * sizeof(int));

	index_array = (int *)malloc(1620*sizeof(int));

	parent1 = (int *)malloc(length*sizeof(int));
	parent2 = (int *)malloc(length*sizeof(int));
	parent1Indices = (int *)malloc(length*sizeof(int));
	parent2Indices = (int *)malloc(length*sizeof(int));

	parent1Copy = (int *)malloc(length*sizeof(int));
	parent2Copy = (int *)malloc(length*sizeof(int));
	parent1CopyIndices = (int *)malloc(length*sizeof(int));
	parent2CopyIndices = (int *)malloc(length*sizeof(int));

	candidateGenomes = (int **)malloc(4*sizeof(int *));
	candidateIndices = (int **)malloc(4*sizeof(int *));

	converge = (int *)malloc(GENERATION*sizeof(int));

	for (i=0; i<4; i++)
	{
		candidateGenomes[i] = (int *)malloc(length*sizeof(int));
		candidateIndices[i] = (int *)malloc(length*sizeof(int));
	}

	// printf("program start!");
	// int * indicesss = (int *)malloc(length*sizeof(int));
	// int a = 0;
	// int ddd = 0;
	// ddd = calculateDCJdistance(original_list[0], original_list[1], length);
	// printf("%d \n", ddd);
	// while (ddd > 0)
	// {
	// 	DCJOperation(original_list[0], original_list[1], length, &a, indicesss);
	// 	ddd = calculateDCJdistance(original_list[0], original_list[1], length);
	// 	printf("%d \n", ddd);
	// }

	// exit(0);
	// printf("d: %d\n", calculateDCJdistance(original_list[0], original_list[1], length));
	// exit(0);

	generateInitialPool(original_list, length);

	generation = GENERATION;
	
	while (generation>0)
	{
		geneticAlgorithm(original_list, length);
		printf("bestscore: %d\n", fitnessScore[0]);
		printf("MedianScore: %d\n", calculateDCJdistance(genome_pool[0],original_list[0], length)+calculateDCJdistance(genome_pool[0],original_list[1], length) +calculateDCJdistance(genome_pool[0],original_list[2], length));
		tempMScore = calculateDCJdistance(genome_pool[0],original_list[0], length)+calculateDCJdistance(genome_pool[0],original_list[1], length) +calculateDCJdistance(genome_pool[0],original_list[2], length);
		converge[GENERATION-generation] = tempMScore;		
		generation--;
	}

	optMScore =  calculateDCJdistance(genome_pool[0],original_list[0], length)+calculateDCJdistance(genome_pool[0],original_list[1], length) +calculateDCJdistance(genome_pool[0],original_list[2], length);
	
	for(k=0; k<GENERATION; k++)
	{
		if (converge[k] == optMScore)
		{
			temp_optindex = k+1;
			break;
		}
	}
	/*
	optMScore = genen*3;
	bestMScore = 0;
	bestGNum = 0;
	//optMGenomeNum = 0;
	//run ten times select the best one
	for (i=0; i<10; i++)
	{
		simulated_annealing(original_list, median_genome, length, &bestMScore, &bestGNum);
		if (bestMScore < optMScore)
		{
			optMScore = bestMScore;
			optMGenomeNum = bestGNum;
			for (j=0; j<length; j++)
			{
				bestMedian[j] = bestGenome[j];
			}
			for (j=0; j<optMGenomeNum; j++)
			{
				optIndices[j*2] = bestGenomeIndices[j*2];
				optIndices[j*2+1] = bestGenomeIndices[j*2+1];	
			}
		}
	}
	*/
	int* dupBestMedian = (int *)malloc(length*sizeof(int));
	//copy best median genome
	for (i=0; i<length; i++)
	{
		dupBestMedian[i] = genome_pool[0][i];
		bestMedian[i] = genome_pool[0][i];
		optIndices[i] = g_pool_indices[0][i];
	}
	optMGenomeNum = g_pool_cnum[0];

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
    int *trueGenome = (int *)malloc(genen*sizeof(int));
    for (i=0; i<genen; i++)
    {
    	trueGenome[i] = i+1;
    }

    int* trueGenomeAdj = (int *)malloc((length)*sizeof(int));
	encodeAdjacency(trueGenome, genen, trueGenomeAdj);

    /*write the result to out file*/
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

    fprintf(f, "Optimal Median Score: %d \n", optMScore);

    fprintf(f, "Medain Genome:\n");
    if (optMGenomeNum == 0)
    {
    	fprintf(f,"->C:\n");
		for (i=1; i<length-1; i=i+2)
		{
			if (bestMedian[i] % 2 == 0)
			{
				fprintf(f, "%d ", (bestMedian[i] / 2)*(-1));
			}
			else
			{
				fprintf(f, "%d ", bestMedian[i+1] / 2);
			}
		}
		fprintf(f,"\n");
    }
    else
    {
    	//first decode the head tail list
		tail_index = optIndices[0] - 1;
		fprintf(f, "->C:\n");
		for (i=1; i<tail_index; i=i+2)
		{
			if (bestMedian[i] % 2 == 0)
			{
				fprintf(f, "%d ", (bestMedian[i] / 2)*(-1));
			}
			else
			{
				fprintf(f, "%d ", bestMedian[i+1] / 2);
			}
		}
		fprintf(f, "\n");
		for (i=0; i<optMGenomeNum; i++)
		{
			fprintf(f, "->C:\n");
			//move the last ele to head, and every other ele move to right one cell
			//t = 0;
			t = bestMedian[optIndices[i*2+1]];
			for (j= optIndices[i*2]; j < optIndices[i*2+1]; j++)
			{
				bestMedian[optIndices[i*2+1]+optIndices[i*2]-j]= bestMedian[optIndices[i*2+1]+optIndices[i*2]-j-1];
			}
			bestMedian[optIndices[i*2]] = t;
			for (j= optIndices[i*2]; j< optIndices[i*2+1]; j=j+2)
			{
				if (bestMedian[j] % 2 == 0)
				{
					fprintf(f, "%d ", (bestMedian[j] / 2)*(-1));
				}
				else
				{
					fprintf(f, "%d ", bestMedian[j+1] / 2);
				}
			}
			fprintf(f, "\n");
		}	
    }

    fprintf(f, "To true distance: %d \n", calculateDCJdistance(dupBestMedian,trueGenomeAdj, length));
    fprintf(f, "Convergence:\n");
    for (i=0; i<GENERATION; i++)
    {
		fprintf(f, "%d ", converge[i]);
    }
    fprintf(f, "\n");
  
	
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
	free(temp_list);
	free(temp_list2);
	free(internal_templist);
	free(temp_htlist);
	free(reverse_list);
	free(temp_initial_genome);
	free(temp_indices);
	free(initialTrueGenome);
	free(initilcandid_circular_num);
	free(initialTrueGenomeIndices);
	free(updateTrueGenome);
	free(candidMedian_circular_num);
	free(updateTrueGenomeIndices);
	free(currentTrueGenome);
	free(currentTrueGenomeIndices);
	for (i=0; i<6; i++)
    {
    	free(initilcandid_indices[i]);
    	free(candidateInitialTrue[i]); 
    }
    for (i=0; i<3; i++)
    {
    	free(original_list[i]);
    	free(candidMedian_indices[i]);
    	free(candidateMedianGenme[i]);
    }
    free(original_list);
    free(candidMedian_indices);
    free(candidateMedianGenme);
    free(initilcandid_indices);
	free(candidateInitialTrue);
	free(bestGenomeIndices);
	free(bestGenome);
	free(g);
	free(g2);
	free(g3);
	free(bestMedian);

	/*genetic algo*/
	for (i=0; i<1800; i++)
	{
		free(genome_pool[i]);
		free(g_pool_indices[i]);
	}
	free(genome_pool);
	free(g_pool_indices);
	free(g_pool_cnum);
	free(global_indices);
	// free(global_circularNum); 
	free(original_copy_a);
	free(original_copy_b);
	free(temporary_genome);
	free(temporary_indices);
	free(best_genome);
	free(best_genomeIndices);
	free(index_array);
	free(parent1); 
	free(parent2); 
	free(parent1Indices); 
	free(parent2Indices);
	free(parent1Copy);
	free(parent2Copy);
	free(parent1CopyIndices); 
	free(parent2CopyIndices);
	free(fitnessScore);

	free(converge);

	for (i=0; i<4; i++)
	{
		free(candidateGenomes[i]);
		free(candidateIndices[i]);
	}
	free(candidateGenomes);
	free(candidateIndices);

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//printf("time: %.2f\n", time_spent);
	fprintf(f, "average_coveragence:%d\n",temp_optindex);
	fprintf(f, "running time: %.2f\n", time_spent);
	fclose(f);

	return 0;
}
