/*
 *  map.h 
 */

#ifndef MAP_H_
#define MAP_H_

typedef struct {
    int   id_i;  			/* Landmark ID */
    float x_f;   			/* Landmark x-position in map (global coords) */
    float y_f;   			/* Landmark y-position in map (global coords) */
} single_landmark_s;

typedef struct {
    single_landmark_s *landmark_list; 	/* dynamically grown array */
    int count;                       	/* number of landmarks in use */
    int capacity;                     	/* allocated slots */
} Map;

#endif /* MAP_H_ */
