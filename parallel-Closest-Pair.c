#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <stddef.h>
#include <math.h>
#include <time.h>
#define  MASTER		0


typedef struct point_struct {
    int x, y;
}point;

///////////////////////////////// sort by X_coordinate //////////////////////
void b_s_x(int n, struct point_struct p_x2[]){
    int t1=0,t2=0,i,j;

	///////////////////////sorting x bye first-coordinate
	for(i=0 ; i<n-1 ; i++){
		for(j=0 ; j< n-i-1 ; j++){
			if(p_x2[j].x > p_x2[j+1].x){
				t1 = p_x2[j].x;
				t2 = p_x2[j].y;

				p_x2[j].x = p_x2[j+1].x;
				p_x2[j].y = p_x2[j+1].y;

				p_x2[j+1].x = t1;
				p_x2[j+1].y = t2;
			}
		}
	}

}
//////////////////////////////////sort by y-coordinate /////////////////////
void b_s_y(int n, struct point_struct p_y2[]){
    int t1=0,t2=0,i,j;

	///////////////////////// sorting y by 2rd-coordinate
	for(i=0 ; i<n-1 ; i++){
		for(j=0 ; j< n-i-1 ; j++){
			if(p_y2[j].y > p_y2[j+1].y){
				t1 = p_y2[j].x;
				t2 = p_y2[j].y;

				p_y2[j].x = p_y2[j+1].x;
				p_y2[j].y = p_y2[j+1].y;

				p_y2[j+1].x = t1;
				p_y2[j+1].y = t2;
			}
		}
	}

}

////////////////////////////////// Distance ////////////////////////////////
double distance(int x1, int x2, int y1, int y2){
    return  sqrt( pow((y1 - x1),2) + pow( (y2 - x2),2)  );
}
////////////////////////////////// brute-forth /////////////////////////////
double brute_forth(int x1,int y1,int x2,int y2,int x3,int y3){
	double d1_2=distance( x1, y1, x2, y2),
	       d1_3=distance( x1, y1, x3, y3),
	       d2_3=distance( x2, y2, x3, y3),
	       d = d1_2;
	       if(d < d1_3)
		  d = d1_3;
	       else if(d<d2_3)
		  d = d2_3;
	return d; 
}

///////////////////////////////// boundary_check ///////////////////////////
double boundary_check(int f, int n, struct point_struct p_y4[], int mid_x, double min_dist){
    struct point_struct y_strip[n];
    int y_strip_len, i, j=0, k=0;
    double dist=min_dist;

    for(j=0, i=f; i<n ; i++){
        if( abs(p_y4[i].x - mid_x) < min_dist){
            y_strip[j].x = p_y4[i].x;
            y_strip[j].y = p_y4[i].y;
            j++;
        }
    }

    y_strip_len = j;

    double d;
    for(j=0 ; j< y_strip_len-1 ; j++){
        k = j+1;
        while( (k<y_strip_len-1) && abs(y_strip[k].y - y_strip[j].y)<min_dist){
            d = distance(y_strip[k].x, y_strip[k].y, y_strip[j].x, y_strip[j].y);
            if(dist >= d)
                dist = d;
            k++;
        }
    }

    return dist;

}


///////////////////////////////// Closest_pair /////////////////////////////
double Closest_Pair(int taskid,int f, int l, int n, struct point_struct p_x3[]){
    int i;
    if(n == 1){
        double d=9999999;
        return d;
    }
    if(n == 2){
        double d=distance(p_x3[f].x, p_x3[f].y, p_x3[l].x, p_x3[l].y);
        return d;
    }
    if(n == 3){
        double d=brute_forth(p_x3[f].x, p_x3[f].y, p_x3[f+1].x, p_x3[f+1].y, p_x3[l].x, p_x3[l].y);
        return d;
    }

    int mid;
    mid = n/2;

    double dist_l,dist_r;

    dist_l = Closest_Pair(taskid, f, f+mid-1, mid, p_x3);
    dist_r = Closest_Pair(taskid, f+mid, l, n-mid, p_x3);

    double min_dist;
    if(dist_l <= dist_r)
        min_dist = dist_l;
    else
        min_dist = dist_r;
    return min_dist;
}



int main(int argc, char *argv[])
{
    int numtasks, taskid; 
    int n = atoi(argv[1]);
    const int num_item =2;
    int blocklengths[2] = {1,1};
    int seed = clock();
    int i,j;
    clock_t endt,start;
    srand(seed);

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    double inicio = MPI_Wtime( );

    if(n % numtasks != 0){
	printf("points number not dividable by number of processors\n");
	MPI_Finalize();
	return -1;
    }
    MPI_Status status;

    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype mpi_point_type;
    MPI_Aint offsets[2];

    offsets[0] = offsetof(point, x);
    offsets[1] = offsetof(point, y);

    MPI_Type_create_struct(num_item, blocklengths, offsets, types, &mpi_point_type);
    MPI_Type_commit(&mpi_point_type);
    
    point s[n];
    struct point_struct *p_x = (struct point_struct*)malloc(n*sizeof(struct point_struct));
    
    double *dist_closest_pair = (double*)malloc(numtasks*sizeof(double));
    int offset[numtasks], share_len = (n/numtasks);
    offset[taskid] = taskid * share_len;
    
    if (taskid == MASTER){
	for(i = 0; i < n ; i++) {
		s[i].x =  rand()%1000;
		s[i].y =  rand()%1000;
	}

	for(i=0 ; i<n ; i++){
		p_x[i].x = s[i].x;
		p_x[i].y = s[i].y;
	}
	start  = clock();
  b_s_x(n,p_x);

    }//MASTER

    MPI_Scatter(&p_x[0], share_len, mpi_point_type, &p_x[offset[taskid]], share_len, mpi_point_type, MASTER , MPI_COMM_WORLD);

    
    for(i=0 ;i < numtasks; i++){
	if(taskid == i ){
	    dist_closest_pair[taskid] = Closest_Pair(taskid,offset[taskid], offset[taskid]+share_len-1, share_len, p_x);
	}
    }

    MPI_Gather(&dist_closest_pair[taskid] , 1, MPI_DOUBLE, &dist_closest_pair[taskid] , 1 , MPI_DOUBLE , MASTER , MPI_COMM_WORLD);

    if(taskid == MASTER){
	point p_y[2*share_len];
	int x[numtasks-1];
	for(i=0 ;i< numtasks-1 ; i++){
	    x[i]= (i*share_len)+share_len;
	} 
	double d_boundary[numtasks-1], d_min_proc=dist_closest_pair[0];
	for(i=1 ; i<numtasks ; i++){
	    if(d_min_proc > dist_closest_pair[i])
		d_min_proc=dist_closest_pair[i];
	}

	for(i=0 ; i<numtasks-1 ; i++){
	    for(j=x[i]-share_len ; j<x[i]+share_len ; j++){
		p_y[j] = p_x[j];
	    } 
	    b_s_y(2*share_len,p_y);
	    d_boundary[i] = boundary_check(x[i]-share_len, 2*share_len, p_y, x[i], d_min_proc );
	}

	double D_min = d_min_proc;	
	for(i=0 ; i<numtasks-1 ; i++){
	    if(d_boundary[i] < D_min )
		D_min = d_boundary[i];
	}

    double final = MPI_Wtime( );
    double time = final - inicio ;
	
	printf("\n minimum distanse is : %f.\n Tempo de execução de: %f segundos\n",D_min, time);
   }

    MPI_Finalize();
    return 0;
}











