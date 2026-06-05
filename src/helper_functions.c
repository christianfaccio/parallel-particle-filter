#include "helper_functions.h"
#include <stdio.h>
#include <stdlib.h>

/* Grow a dynamic array to hold at least count+1 elements.
 * Returns the (possibly moved) buffer, or NULL on allocation failure
 * (in which case the original buffer is untouched and must be freed by the caller). */
static void *grow_if_needed(void *arr, int count, int *capacity, size_t elem_size) {
	if (count < *capacity) return arr;
        int new_cap = (*capacity == 0) ? 64 : *capacity * 2;
        void *p = realloc(arr, (size_t)new_cap * elem_size);
        if (!p) return NULL;
        *capacity = new_cap;
        return p;
}

bool read_map_data(const char *filename, Map *map) {
	FILE *f = fopen(filename, "r");
      	if (!f) return false;

      	map->landmark_list = NULL;
      	map->count = 0;
      	map->capacity = 0;

      	float x, y;
      	int id;
      	while (fscanf(f, "%f %f %d", &x, &y, &id) == 3) {
	  	void *tmp = grow_if_needed(map->landmark_list, map->count, &map->capacity, sizeof *map->landmark_list);
	  	if (!tmp) { free(map->landmark_list); fclose(f); return false; }
	  	map->landmark_list = tmp;

	  	map->landmark_list[map->count].x_f  = x;
	  	map->landmark_list[map->count].y_f  = y;
	  	map->landmark_list[map->count].id_i = id;
	  	map->count++;
      	}
      	fclose(f);
      	return true;
}

bool read_control_data(const char *filename, control_s **position_meas, int *count) {
	FILE *f = fopen(filename, "r");
        if (!f) return false;

      	control_s *arr = NULL;
      	int n = 0, cap = 0;

      	double velocity, yawrate;
      	while (fscanf(f, "%lf %lf", &velocity, &yawrate) == 2) {
        	void *tmp = grow_if_needed(arr, n, &cap, sizeof *arr);
          	if (!tmp) { free(arr); fclose(f); return false; }
          	arr = tmp;

          	arr[n].velocity = velocity;
          	arr[n].yawrate  = yawrate;
          	n++;
      	}
      	fclose(f);
      	*position_meas = arr;
      	*count = n;
      	return true;
}

bool read_gt_data(const char *filename, ground_truth **gt, int *count) {
        FILE *f = fopen(filename, "r");
        if (!f) return false;

      	ground_truth *arr = NULL;
      	int n = 0, cap = 0;

      	double x, y, azimuth;
      	while (fscanf(f, "%lf %lf %lf", &x, &y, &azimuth) == 3) {
        	void *tmp = grow_if_needed(arr, n, &cap, sizeof *arr);
          	if (!tmp) { free(arr); fclose(f); return false; }
          	arr = tmp;

          	arr[n].x     = x;
          	arr[n].y     = y;
          	arr[n].theta = azimuth;
          	n++;
	}
      	fclose(f);
      	*gt = arr;
      	*count = n;
      	return true;
}

bool read_landmark_data(const char *filename, LandmarkObs **observations, int *count) {
	FILE *f = fopen(filename, "r");
      	if (!f) return false;

      	LandmarkObs *arr = NULL;
      	int n = 0, cap = 0;

      	double local_x, local_y;
      	while (fscanf(f, "%lf %lf", &local_x, &local_y) == 2) {
        	void *tmp = grow_if_needed(arr, n, &cap, sizeof *arr);
        	if (!tmp) { free(arr); fclose(f); return false; }
          	arr = tmp;

          	arr[n].id = 0;       /* set during data association, as in the original */
          	arr[n].x  = local_x;
          	arr[n].y  = local_y;
          	n++;
      	}
      	fclose(f);
      	*observations = arr;
      	*count = n;
	return true;
}
