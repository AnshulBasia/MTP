#include "Graph.h"

int main(int argc, char *argv[]){

	
	char *graphFile = argv[1];
	char *featFile = argv[2];
	char *queryFile = argv[3];
	int dir_control = atoi(argv[4]);
	

	srand(time(NULL));
	Graph *tg = new Graph(graphFile, dir_control);
	
}
